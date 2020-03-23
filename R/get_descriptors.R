get_descriptors <- function(vcf_filename, read_length, fasta_filename) {
  # Check required arguments
  defined <- ls()
  passed <- names(as.list(match.call())[-1])
  check_arguments(defined, passed)
  # Get mutations ids and samplename
  mutation_ids <- get_mut_id(vcf_filename)
  samplename <- get_samplename(vcf_filename)
  # Descriptors from the VCF file
  descriptors_vcf <- extract_vcf( mutation_ids = mutation_ids, vcf_filename = vcf_filename, read_length = read_length, samplename = samplename)
  # Descriptors from FASTA file
  descriptors_fasta <- extract_fa(mutation_ids = mutation_ids, fasta_filename = fasta_filename, k = 2)
}

#' Title
#'
#' @param vcf_filename
#' @param mutation_ids
#' @param read_length
#' @param samplename
#'
#' @return
#' @export
#' @import tibble
#' @import dplyr
#'
#' @examples
extract_vcf <- function(vcf_filename, mutation_ids, read_length, samplename) {
  ## Relate each descriptor with its argument
  args_to_vcf_funs <- list(
    allele.freq = "get_AF",
    alt.bases = "get_alt_bases",
    ref.bases = "get_ref_bases",
    ref.allele = "get_ref_allele",
    alt.allele = "get_alt_allele"
    #base.qual = "get_base_qual",
    #base.qual.frac = "get_base_qual_fraction",
    #frag.length = "get_fraglength",
    #pos.from.end = "get_pos_from_end",
    #map.qual = "get_map_qual",
    #FdeamC = "get_FdeamC",
    # SB = "get_SB"
  )
  depth <- get_depth(vcf_filename = vcf_filename, samplename = samplename)
  raw_vcf_descriptors <- sapply(names(args_to_vcf_funs), function(x) {
    fun <- get(args_to_vcf_funs[[x]])
    fun(vcf_filename, depth = depth, read_length = read_length, samplename = samplename)
  })
  raw_vcf_descriptors <- lapply(rapply(raw_vcf_descriptors, enquote, how = "unlist"), eval)
  vcf_descriptors <- bind_rows(raw_vcf_descriptors) %>%
    add_column(id = mutation_ids)
}

#' Extract descriptors from fasta file
#'
#' \code{extract_fa} extracts the descriptors of the data that are derived from a reference genome fasta file
#'
#' @param mutation_ids Tibble containing in each row the mutation id of the variants to be subjected to filtering. Mutation id format is chr:pos.
#' @param fasta_filename Reference genome data was aligned to.
#' @param k Number of bases to the right and to the left of the variant locus to get the genomic sequence from.
#'
#' @return Tibble containing the base composition of k positions before and k positions after the variant locus. Columns are mutation_id, current, bases, before_2, before_1, after_1, after_2, before and after.
#' @export
#' @import dplyr
#' @import tidyr
#' @examples
extract_fa <- function(mutation_ids, fasta_filename, k = 2) {
  message(sprintf("Variant context of +- %s positions will be extracted.", k))
  # Get descriptors related to surrounding bases
  surr_bases_raw <- get_surr_bases(fasta_filename = fasta_filename, mutation_ids = mutation_ids, k = k)
  ## Tidy and format surr_bases_raw
  odd_idx <- seq(from = 1, to = length(surr_bases_raw), by = 2)
  even_idx <- seq(from = 2, to = length(surr_bases_raw), by = 2)
  base_mat <- data.frame(region = surr_bases_raw[odd_idx],
                             bases = surr_bases_raw[even_idx]) %>%
    as_tibble()
  ## Separate bases in before and after & obtain mutation id_from region
  ## The following lines are specific to k = 2. They won't work for another choice of k
  base_mat %>%
    separate(col = bases, into = c("before", "after"), sep = 3, remove = FALSE) %>%
    separate(col = before, into = c("before", "current"), sep = -1, remove = TRUE) %>%
    separate(col = before, into = c("before_2", "before_1"), sep = 1, remove = FALSE) %>%
    separate(col = after, into = c("after_1", "after_2"), sep = 1, remove = FALSE) %>%
    mutate(chr = sub(":.*", "", region) %>%
             sub(">", "", .)) %>%
    mutate(pos = sub(".*:", "", region) %>%
             sub("-.*", "", .) %>%
             as.numeric() %>%
             magrittr::add(k)) %>%
    mutate(id = paste(chr, pos, sep = ":")) %>%
    select(id, current, bases, before_2, before_1, after_1, after_2, before, after) -> surr_bases
  ## Convert to uppercase
  cols_toupper <- setdiff(colnames(surr_bases), c("id", "bases"))
  surr_bases %>%
    mutate_at(cols_toupper, toupper) -> surr_bases
  return(surr_bases)
}


#' Get surrounding bases
#'
#' Returns the base composition of the fragment that surrounds each of the provided genomic variants. The region goes from the locus - \strong{k} to locus + \strong{k}, including the locus itself.
#'
#' @param fasta_filename Reference genome data was aligned to.
#' @param mutation_ids Tibble containing in each row the mutation id of the variants to be subjected to filtering. Mutation id format is samplename:chr:pos.
#' @param k Number of bases to the right and to the left of the variant locus to get the genomic sequence from.
#'
#' @return Genomic bases \strong{k} positions before and after each variant in \strong{mutation_ids}, including the variant locus itself.
#' @export
#' @import dplyr
#' @examples
get_surr_bases <- function(fasta_filename, mutation_ids, k) {
  tmpdir <- "/tmp"
  tmp_filename <- file.path(tmpdir, "mut_regions.temp")
  write_region_file(mutation_ids = mutation_ids, k = k, tmp_filename = tmp_filename)
  cmd_run <- sprintf("samtools faidx %s -r %s", fasta_filename, tmp_filename)
  bases <- system(cmd_run, intern = TRUE)
  cmd_rm <- sprintf("rm %s", tmp_filename)
  system(cmd_rm)
  return(bases)
}

#' Write temporary locus region file
#'
#' Creates the file needed by samtools faidx to get the surrounding bases of a locus
#'
#' @param mutation_ids Vector containing the mutation id of the variants to be subjected to filtering. Mutation id format is samplename:chr:pos.
#' @param k Number of bases to the right and to the left of the variant locus to get the genomic sequence from.
#' @param tmp_filename Path of the temporary file.
#'
#' @return None
#' @export
#' @import dplyr
#' @importFrom utils write.table
#' @examples
write_region_file <- function(mutation_ids, k, tmp_filename) {
  strsplit(mutation_ids, ":") %>%
    data.frame() %>%
    t() %>%
    data.frame() %>%
    magrittr::set_colnames(c("chr", "pos")) %>%
    mutate_at(1, as.character) %>%
    mutate_at(2, as.character) %>%
    mutate_at(2, as.numeric) %>%
    mutate(start = pos - k, end = pos + k) %>%
    mutate(region = sprintf("%s:%i-%i", chr, start, end)) %>%
    select(region) %>%
    write.table(., file = tmp_filename, sep = "\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
}

#' Get sample name from vcf file
#'
#' It extracts from the VCF filename the name given to the sample in which variants were called
#'
#' It corresponds to the tumor field in the vcf filename. In tumor-only variant calling, that is the only sample. In tumor-normal variant calling, it corresponds to the tumor sample
#'
#' @param vcf_filename VCF filename to extract the sample name from
#'
#' @return Name given to the sample the vcf filename contains variants from
#' @export
#'
#' @examples
get_samplename <- function(vcf_filename) {
  cmd <- sprintf("grep '##tumor_sample=' %s", vcf_filename)
  samplename_line <- system(cmd, intern = TRUE)
  samplename <- gsub("##tumor_sample=", "", samplename_line)
  return(samplename)
}

#' Get mutation id
#'
#' Parses a vcf file and creates a vector of the mutation ids of the variants to be subjected to filtering
#'
#' @param vcf_filename Variant call vcf filename
#'
#' @return Vector containing the mutation id of the variants to be subjected to filtering. Mutation id format is samplename:chr:pos.
#' @export
#'
#' @examples
get_mut_id <- function(vcf_filename) {
  cmd <- sprintf("bcftools query -i 'FILTER=\"PASS\"' -f '%%CHROM:%%POS\n' %s", vcf_filename)
  mutation_ids <- system(cmd, intern = TRUE)
  return(mutation_ids)
}

#' Get depth
#'
#' Gets read depth at variant locus
#'
#' @param vcf_filename Variant call vcf filename
#' @param samplename Name given to the sample the vcf filename contains variants from, ordered in order of appearance in the vcf file
#' @param ...
#'
#' @return Vector containing the read depth at each variant locus
#' @export
#' @import dplyr
#'
#' @examples
get_depth <- function(vcf_filename, samplename, ...) {
  cmd <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%AD]\n'", samplename, vcf_filename) # FORMAT fields need the brackets, INFO dont
  depth <- system(cmd, intern = TRUE) %>%
    as_tibble() %>%
    separate(value, into = c("ref", "alt"), sep = ",") %>% # this will produce NAs when we have tri/cuatriallelic sites. However it's not a problem because we later discard those sites and only keep C:G>T:A
    mutate_all(funs(as.numeric)) %>%
    mutate(depth = ref + alt) %>%
    pull(depth)
  return(depth)
}

#' Get allele frequency
#'
#' @inheritParams get_depth
#'
#' @return Vector containing the allele frequency (AF) of each variant, ordered in order of appearance in the vcf file
#' @export
#'
#' @examples
get_AF <- function(vcf_filename, samplename, ...) {
  cmd <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%AF]\n'", samplename, vcf_filename)
  AF <- system(cmd, intern = TRUE)
  AF <- as.numeric(AF)
  return(AF)
}

#' Get reference allele
#'
#' @param vcf_filename Variant call vcf filename
#' @param ...
#'
#' @return Vector containing the reference allele of each variant, ordered in order of appearance in the vcf file
#' @export
#'
#' @examples
get_ref_allele <- function(vcf_filename, ...) {
  cmd <- sprintf("bcftools query -i 'FILTER=\"PASS\"' -f '%%REF\n' %s", vcf_filename)
  ref_alleles <- system(cmd, intern = TRUE)
  return(ref_alleles)
}

#' Get alternate allele
#'
#' @inheritParams get_ref_allele
#'
#' @return Vector containing the alternate allele of each variant, ordered in order of appearance in the vcf file
#' @export
#'
#' @examples
get_alt_allele <- function(vcf_filename, ...) {
  cmd <- sprintf("bcftools query -i 'FILTER=\"PASS\"' -f '%%ALT\n' %s", vcf_filename)
  alt_alleles <- system(cmd, intern = TRUE)
  return(alt_alleles)
}

#' Get number of alternate alleles
#'
#' Get number of reads that support the alternate allele
#'
#' @inheritParams get_depth
#' @param depth Vector containing the read depth at each variant locus
#' @param ...
#'
#' @return Vector containing the number of reads that support the alternate allele
#' @export
#'
#' @examples
get_alt_bases <- function(vcf_filename, depth, samplename, ...) {
  cmd <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%AD{1}]\n'", samplename, vcf_filename)
  alt_bases <- system(cmd, intern = TRUE)
  alt_bases <- as.numeric(alt_bases)
  alt_base_ratio <- alt_bases/depth
  return(list(alt_bases, alt_base_ratio))
}

#' Get number of reference alleles
#'
#' @inheritParams get_alt_bases
#'
#' @return Vector containing the number of reads that support the reference allele
#' @export
#'
#' @examples
get_ref_bases <- function(vcf_filename, depth, samplename, ...) {
  cmd <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%AD{0}]\n'", samplename, vcf_filename)
  ref_bases <- system(cmd, intern = TRUE)
  ref_bases <- as.numeric(ref_bases)
  ref_base_ratio <- ref_bases/depth
  return(list(ref_bases, ref_base_ratio))
}
