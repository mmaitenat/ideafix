#' Get descriptors
#'
#' Get all the C:G > T:A variant descriptors needed for running the filtering.
#'
#' Get all the C:G > T:A variant descriptors needed for running the filtering. The descriptors are obtained from the variant vcf file and reference genome fasta file. The function can only run on a single vcf file.
#'
#' @param vcf_filename Variant call vcf filename
#' @param fasta_filename Reference genome data was aligned to.
#'
#' @return A tibble containing the descriptors of the C:G > T:A variants needed for running the filtering.
#' @export
#'
#' @examples
get_descriptors <- function(vcf_filename, fasta_filename) {
  # Check required arguments
  defined <- ls()
  passed <- names(as.list(match.call())[-1])
  check_arguments(defined = defined, passed = passed)
  # Get mutations ids and samplename
  mutation_ids <- get_mut_id(vcf_filename)
  samplename <- get_samplename(vcf_filename)
  # Descriptors from the VCF file
  message("Extracting variants from vcf file...")
  descriptors_vcf <- extract_vcf( mutation_ids = mutation_ids, vcf_filename = vcf_filename, samplename = samplename)
  message("Finished extracting variants from vcf file.")
  # Descriptors from FASTA file
  message("Extracting variants from fasta file...")
  descriptors_fasta <- extract_fa(mutation_ids = mutation_ids, fasta_filename = fasta_filename, k = 2)
  message("Finished extracting variants from fasta file.")
  # Join descriptors from vcf and FASTA
  all_descriptors <- full_join(descriptors_vcf, descriptors_fasta, by = "id")
  # Filter out those that are not C:G>T:A
  all_descriptors %>%
    mutate(mut = paste(ref.allele, alt.allele, sep = ":")) %>%
    filter(mut %in% c("C:T", "G:A")) %>%
    select(-mut) -> all_descriptors
  # Fix fragment length == 0 (actually this is a coding that indicates that the two reads map to different chromosomes)
  all_descriptors <- all_descriptors %>%
    mutate(norm.pos.from.end = replace(norm.pos.from.end, both.reads.aligned == 0, NA),
           frag.length = replace(frag.length, both.reads.aligned == 0, NA))
  # Rename columns
  all_descriptors <- rename(all_descriptors,
                            alt.bases = alt.bases1,
                            norm.alt.bases = alt.bases2,
                            ref.bases = ref.bases1,
                            norm.ref.bases = ref.bases2,
                            base.qual = base.qual1,
                            base.qual.frac = base.qual2,
                            FdeamC = ROB1,
                            SOB = ROB2,
                            SBGuo = SB1,
                            SBGATK = SB2)
  message("All descriptors have been extracted.")
  return(all_descriptors)
}

#' Extract descriptors from vcf file
#'
#' Extract model descriptors that are found in vcf files
#'
#' @param vcf_filename VCF filename to extract the sample name from
#' @param mutation_ids Vector that contains the mutation ids of the variants to be subjected to filtering. Mutation id format is chr:pos.
#' @param samplename Name given to the sample the vcf filename contains variants from, ordered in order of appearance in the vcf file
#'
#' @return A tibble containing the descriptors of the variants that are extracted from vcf files.
#' @export
#' @import tibble
#' @import dplyr
#'
#' @examples
extract_vcf <- function(vcf_filename, mutation_ids, samplename) {
  ## Relate each descriptor with its argument
  args_to_vcf_funs <- list(
    allele.freq = "get_AF",
    alt.bases = "get_alt_bases",
    ref.bases = "get_ref_bases",
    ref.allele = "get_ref_allele",
    alt.allele = "get_alt_allele",
    base.qual = "get_base_qual",
    frag.length = "get_fraglength",
    pos.from.end = "get_pos_from_end",
    map.qual = "get_map_qual",
    ROB = "get_RPOB",
    SB = "get_SB"
  )
  depth <- get_depth(vcf_filename = vcf_filename, samplename = samplename)
  raw_vcf_descriptors <- sapply(names(args_to_vcf_funs), function(x) {
    fun <- get(args_to_vcf_funs[[x]])
    fun(vcf_filename, depth = depth, samplename = samplename)
  })
  raw_vcf_descriptors <- lapply(rapply(raw_vcf_descriptors, enquote, how = "unlist"), eval)
  vcf_descriptors <- bind_rows(raw_vcf_descriptors) %>%
    add_column(id = mutation_ids)
  # Add some new features: both_reads_aligned, norm_pos_from_end
  vcf_descriptors %>%
    mutate(both.reads.aligned = ifelse(frag.length != 0, 1, 0),
           norm.pos.from.end = pos.from.end/frag.length) -> vcf_descriptors
  return(vcf_descriptors)
}

#' Extract descriptors from fasta file
#'
#' \code{extract_fa} extracts the descriptors of the data that are derived from a reference genome fasta file
#'
#' @param mutation_ids Vector that contains the mutation ids of the variants to be subjected to filtering. Mutation id format is chr:pos.
#' @param fasta_filename Reference genome data was aligned to.
#' @param k Number of bases to the right and to the left of the variant locus to get the genomic sequence from.
#'
#' @return Tibble containing the base composition of k positions before and k positions after the variant locus. Columns are mutation_id, current, bases, before.2, before.1, after.1, after.2, before and after.
#' @export
#' @import dplyr
#' @import tidyr
#' @import purrr
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
    separate(col = before, into = c("before.2", "before.1"), sep = 1, remove = FALSE) %>%
    separate(col = after, into = c("after.1", "after.2"), sep = 1, remove = FALSE) %>%
    mutate(chr = sub(":.*", "", region) %>%
             sub(">", "", .)) %>%
    mutate(pos = sub(".*:", "", region) %>%
             sub("-.*", "", .) %>%
             as.numeric() %>%
             magrittr::add(k)) %>%
    mutate(id = paste(chr, pos, sep = ":")) %>%
    select(id, current, bases, before.2, before.1, after.1, after.2, before, after) -> surr_bases
  # Obtain variable is.repeat.region. The variable will be TRUE if any of the position before, the position itself or the position after are lowercase (= repeat region)
  surr_bases %>%
    mutate(is.repeat.region = surr_bases %>%
             select(c("current", "before.1", "after.1")) %>%
             purrr::pmap_lgl(~any(check_lowercase(.)))) -> fa_descriptors
  ## Convert to uppercase
  cols_toupper <- setdiff(colnames(fa_descriptors), c("id", "bases", "is.repeat.region"))
  fa_descriptors %>%
    mutate_at(cols_toupper, toupper) -> fa_descriptors
  return(fa_descriptors)
}


#' Get surrounding bases
#'
#' Returns the base composition of the fragment that surrounds each of the provided genomic variants. The region goes from the locus - \strong{k} to locus + \strong{k}, including the locus itself.
#'
#' @param fasta_filename Reference genome data was aligned to.
#' @param mutation_ids Vector that contains the mutation ids of the variants to be subjected to filtering. Mutation id format is samplename:chr:pos.
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
#' @param ... Not to be checked
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
#' @param ... Not to be checked
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
#' @param ... Not to be checked
#'
#' @return List containing two vectors: the number of reads that support the alternate allele and ...
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
#' Get number of reads that support the reference allele
#'
#' @inheritParams get_alt_bases
#'
#' @return List containing two vectors: the number of reads that support the reference allele and ...
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

#' Get base quality
#'
#' Get median Phred quality of the bases that support the alternate allele.
#'
#' @inheritParams get_depth
#'
#' @return List containing two vectors: median base quality of the alternate allele and ratio between median base quality of the alternate allele and median base quality of the reference allele
#' @export
#'
#' @examples
get_base_qual <- function(vcf_filename, samplename, ...) {
  cmd_tumor <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%MBQ{1}]\n'", samplename, vcf_filename)
  cmd_normal <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%MBQ{0}]\n'", samplename, vcf_filename)
  base_qual_tumor <- system(cmd_tumor, intern = TRUE) %>%
    as.numeric()
  base_qual_normal <- system(cmd_normal, intern = TRUE) %>%
    as.numeric()
  base_qual_frac <- base_qual_tumor/base_qual_normal
  return(list(base_qual_tumor, base_qual_frac))
}

#' Get variant median fragment length
#'
#' Get median fragment length for reads that support the alternate allele
#'
#' @inheritParams get_depth
#'
#' @return Vector containing the median median fragment length of each variant, ordered in order of appearance in the vcf file
#' @export
#'
#' @examples
get_fraglength <- function(vcf_filename, samplename, ...) {
  cmd <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%MFRL{1}]\n'", samplename, vcf_filename)
  # Median fragment length for reads that support the alternate allele
  fraglen <- system(cmd, intern = TRUE) %>%
    as.numeric()
  return(fraglen)
}

#' Get distance from end of read
#'
#' Get median value of the position of the variant (starting from the end of the read) in reads that support the alternate allele.
#'
#' @inheritParams get_depth
#'
#' @return Vector containing the median distance from end of read of each variant, ordered in order of appearance in the vcf file
#' @export
#'
#' @examples
get_pos_from_end <- function(vcf_filename, samplename, ...) {
  cmd <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%MPOS]\n'", samplename, vcf_filename)
  pos_from_end <- system(cmd, intern = TRUE) %>%
    as.numeric() # NAs will appear for multiallelic variants
  return(pos_from_end)
}

#' Get mapping quality
#'
#' Get median Phred mapping quality (MAPQ) of the bases that support the alternate allele
#'
#' @inheritParams get_depth
#'
#' @return Vector containing the mapping quality of each variant, ordered in order of appearance in the vcf file
#' @export
#'
#' @examples
get_map_qual <- function(vcf_filename, samplename, ...) {
  cmd <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%MMQ]\n'", samplename, vcf_filename)
  map_qual <- system(cmd, intern = TRUE) %>%
    as.numeric() # NAs will appear for multiallelic variants
  return(map_qual)
}

#' Get read pair orientation bias
#'
#' Get two read pair orientation bias (RPOB) metrics: FdeamC and strand orientation bias
#'
#' @inheritParams get_depth
#'
#' @return List containing two vectors: for each variant, FdeamC and SOB
#' @export
#'
#' @examples
get_RPOB <- function(vcf_filename, samplename, ...) {
  ref_allele <-get_ref_allele(vcf_filename)
  alt_allele <- get_alt_allele(vcf_filename)
  F1R2_alt_cmd <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%F1R2{1}]\n'", samplename, vcf_filename)
  F1R2_alt <- system(F1R2_alt_cmd, intern = TRUE) %>%
    as.numeric()
  F2R1_alt_cmd <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%F2R1{1}]\n'", samplename, vcf_filename)
  F2R1_alt <- system(F2R1_alt_cmd, intern = TRUE) %>%
    as.numeric()
  denom <- F2R1_alt + F1R2_alt
  numer <- F1R2_alt # default numerator
  GA_idx <- ref_allele == "G" & alt_allele == "A"
  numer[GA_idx] <- F2R1_alt[GA_idx] # change if mutation is G>A
  FdeamC <- numer/denom
  # does not apply to not FFPE mutations
  notFFPE_idx <- !(paste(ref_allele, alt_allele, sep = ":") %in% c("C:T", "G:A"))
  FdeamC[notFFPE_idx] <- NA
  SOB <- (F1R2_alt - F2R1_alt)/denom
  return(list(FdeamC, SOB))
}

#' Get strand bias
#'
#' Get two strand bias (SB) metrics: SB defined by GATK and SB defined by Guo et al.
#'
#' @inheritParams get_depth
#'
#' @return List containing two vectors: for each variant, strand bias defined by GATK and strand bias defined by Guo et al.
#' @export
#'
#' @examples
get_SB <- function(vcf_filename, samplename, ...) {
  cmd_SB <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%SB]\n'", samplename, vcf_filename)
  SB_all <- system(cmd_SB, intern = TRUE) %>%
    as_tibble() %>%
    separate(value, into = c("ref_fw", "ref_rev", "alt_fw", "alt_rev"), sep = ",") %>%
    mutate_all(funs(as.numeric)) %>%
    mutate_all(funs(replace(., . == 0, 0.0001))) %>% # I am putting pseudozeros not to divide by 0
    mutate(SB = calc_SB_Guo(a = ref_fw, b = alt_fw, e = ref_rev, d = alt_rev),
           SBGATK = calc_SB_GATK(a = ref_fw, b = alt_fw, e = ref_rev, d = alt_rev))
  return(list(pull(SB_all, SB), pull(SB_all, SBGATK)))
}


#' Calculate Guo strand bias
#'
#' Calculate variant strand bias as defined by Guo et al (SB-Guo)
#'
#' This strand bias value was defined in \emph{The effect of strand bias in Illumina short-read sequencing data}. It can be accessed here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3532123/ Its value ranges from 0 to infinity.
#'
#' @param a Number of reads supporting the reference allele mapping to the forward strand
#' @param b Number of reads supporting the alternate allele mapping to the forward strand
#' @param e Number of reads supporting the reference allele mapping to the forward strand
#' @param d Number of reads supporting the alternate allele mapping to the reverse strand
#'
#' @return Vector containing SB-Guo of each variant, ordered in order of appearance in the vcf file
#' @export
#'
#' @examples
calc_SB_Guo <- function(a, b, e, d) {
  denom1 <- a + b
  denom2 <- e + d
  denom3 <- a + b + e + d
  numer <- abs((b / denom1) - (d / denom2))
  denom <- (b + d) / denom3
  return(numer / denom)
}

#' Calculate GATK strand bias
#'
#' Calculate variant strand bias as defined by GATK (SB-GATK)
#'
#' This strand bias value was defined by GATK and described in \emph{The effect of strand bias in Illumina short-read sequencing data}. It can be accessed here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3532123/ Its value ranges from 0 to infinity.
#'
#'
#' @inheritParams calc_SB_Guo
#'
#' @return Vector containing SB-GATK of each variant, ordered in order of appearance in the vcf file
#' @export
#'
#' @examples
calc_SB_GATK <- function (a, b, e, d) {
  numer1 <- (b * e) / ((a + b) * (e + d))
  numer2 <- (d * a) / ((e + d) * (a + b))
  denom <- (a + e) / (a + b + e + d)
  return(pmax(numer1/denom, numer2/denom))
}

#' Check case
#'
#' Check whether a character is lowercase
#'
#' @param char Character representing a genomic base
#'
#' @return Logical value indicating whether \strong{char} is lowercase.
#' @export
#'
#' @examples
check_lowercase <- function(char) {
  char == tolower(char)
}
