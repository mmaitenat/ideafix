#' Get variant descriptors
#'
#' \code{get_descriptors} obtains a set of descriptors of the C:G > T:A variants in \code{vcf_filename} that are relevant to their classification into deaminations or non-deaminations. These descriptors are extracted from \code{vcf_filename} and \code{fasta_filename}, or calculated using data retrieved from them.
#'
#'
#' @param vcf_filename character string naming the path to the input vcf, i.e. the vcf file containing the variants to classify. This file must have been generated with Mutect2, either in tumor only or tumor/normal mode with strand bias annotation enabled.
#' @param fasta_filename character string naming the path to the reference genome FASTA file the sequencing data was aligned to.
#'
#' @return Tibble containing the descriptors of the C:G > T:A variants in \code{vcf_filename} needed for their classification.
#' @export
#' @import dplyr
#' @details The returned tibble contains the values of each C:G > T:A variant for the following descriptors, divided in columns: VAF, number of alternate bases, normalized number of alternate bases, number of reference bases, normalized number of reference bases, reference allele, alternate allele, base quality, base quality fraction, fragment length, median position from read end, mapping quality, FDeamC, SOB, SB-GUO, SB-GATK, normalized median position from read end, base two positions before, base one position before, base two positions after, base one position after, dinucleotide before and dinucleotide after.
#'
#' For further detail in each of them, see each of the individual help files.
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
  descriptors_vcf <- extract_vcf(mutation_ids = mutation_ids, vcf_filename = vcf_filename, samplename = samplename)
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
  # Frequency cutoff
  AF_cutoff <- 0.3
  all_descriptors <- all_descriptors %>%
    filter(allele.freq <= AF_cutoff)
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

#' Retrieve descriptors from vcf file
#'
#'\code{extract_vcf} obtains from \code{vcf_filename}, a set of descriptors for the C:G > T:A variants defined in \code{mutation_ids}. These descriptors are relevant to their classification into deaminations or non-deaminations and may be directly extracted from \code{vcf_filename} or calculated using data retrieved from it.
#'
#' @param vcf_filename character string naming the path to the input vcf, i.e. the vcf file containing the variants in \code{mutation_ids}. This file must have been generated with Mutect2, either in tumor only or tumor/normal mode with strand bias annotation enabled.
#' @param mutation_ids character vector containing the ids of the loci to get the descriptors of. Id format is CHR:POS and can be obtained by calling the function \code{get_mut_id} on \code{vcf_filename}.
#' @param samplename character string naming the sample in \code{vcf_filename}. This must match the name given to the sample when running Mutect2. It can be obtained by calling \code{get_samplename} on \code{vcf_filename}.
#'
#' @return Tibble containing a set of descriptors for the C:G > T:A variants in \code{mutation_ids} and extracted or calculated using \code{vcf_filename}.
#' @import tibble
#' @import dplyr
#' @details The returned tibble contains the values of each C:G > T:A variant for the following descriptors, divided in columns: VAF, number of alternate bases, normalized number of alternate bases, number of reference bases, normalized number of reference bases, reference allele, alternate allele, base quality, base quality fraction, fragment length, median position from read end, normalized median position from read end, mapping quality, FDeamC, SOB, SB-GUO and SB-GATK.
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

#' Retrieve descriptors from fasta file
#'
#'\code{extract_fa} retrieves from \code{fasta_filename} a set of descriptors for the C:G > T:A variants defined in \code{mutation_ids}. These descriptors are relevant to their classification into deaminations or non-deaminations.
#'
#' @param mutation_ids character vector containing the ids of the loci to get the descriptors of. Id format is CHR:POS.
#' @param fasta_filename character string naming the path to the reference genome FASTA file the sequencing data was aligned to.
#' @param k integer with the number of bases to the right and to the left of the loci to get the genomic sequence from.
#'
#' @return Tibble containing a set of descriptors related to the genomic base sequence ranging from locus - \code{k} to locus + \code{k} for each locus in \code{mutation_ids}. These descriptors, divided in columns, are: base two positions before, base one position before, base two positions after, base one position after, dinucleotide before and dinucleotide after.
#'
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

#' Retrieve surrounding base sequence of locus
#'
#' \code{get_surr_bases} retrieves the base sequence of the \[locus - \code{k}, locus + \code{k}\] region of each locus in \code{mutation_ids} from \code{fasta_filename}.
#'
#' @param fasta_filename character string naming the path to the reference genome FASTA file the sequencing data was aligned to.
#' @param mutation_ids character vector containing the ids of the loci to get the surrounding bases of. Id format is CHR:POS and can be obtained by calling the function \code{get_mut_id} on the vcf file containing those loci.
#' @param k integer with the number of bases to the right and to the left of the loci to get the genomic sequence from.
#'
#' @return Character vector with the genomic base sequence ranging from locus - \code{k} to locus + \code{k} for each locus in \code{mutation_ids}.
#'
#' @details The sequence of the surrounding bases is retrieved from \code{fasta_filename} using samtools faidx tool. In the process, \code{get_surr_bases} creates a temporary file necessary for this tool to run.
#'
#' @examples
get_surr_bases <- function(fasta_filename, mutation_ids, k) {
  tmpdir <- "/tmp"
  tmp_filename <- file.path(tmpdir, "mut_regions.temp")
  write_region_file(ids = mutation_ids, k = k, tmp_filename = tmp_filename)
  cmd_run <- sprintf("samtools faidx %s -r %s", fasta_filename, tmp_filename)
  bases <- system(cmd_run, intern = TRUE)
  cmd_rm <- sprintf("rm %s", tmp_filename)
  system(cmd_rm)
  return(bases)
}

#' Write temporary loci region file
#'
#' \code{write_region_file} creates the file needed by samtools faidx to get the \code{k} surrounding bases of the loci defined in \code{mutation_ids}.
#'
#' @param ids character vector containing the ids of the loci to get the surrounding bases of. Id format is CHR:POS and can be obtained by calling the function \code{get_mut_id} on the vcf file containing those loci.
#' @param k integer with the number of bases to the right and to the left of the loci to get the genomic sequence from.
#' @param tmp_filename character string naming the path to the file to write.
#'
#' @return None
#' @import dplyr
#' @importFrom utils write.table
#' @examples
write_region_file <- function(ids, k, tmp_filename) {
  strsplit(ids, ":") %>%
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

#' Retrieve sample name
#'
#' \code{get_samplename} retrieves from \code{vcf_filename} the name given to the sample in which variants were called. The retrieved sample nane corresponds to the tumor field in \code{vcf_filename}. In tumor-only variant calling, that is the only sample. In tumor-normal variant calling, it corresponds to the tumor sample.
#'
#' @param vcf_filename character string naming the path to the vcf file to extract the sample name from.
#'
#' @return Character string with the name given to the sample when running Mutect2 on the data.
#'
#' @examples
get_samplename <- function(vcf_filename) {
  cmd <- sprintf("grep '##tumor_sample=' %s", vcf_filename)
  samplename_line <- system(cmd, intern = TRUE)
  samplename <- gsub("##tumor_sample=", "", samplename_line)
  return(samplename)
}

#' Retrieve mutation id
#'
#' This function parses \code{vcf_filename} and retrieves the mutation id (CHR:POS) of each C:G > T:A variant in the file.
#'
#' @param vcf_filename character string naming the path to the input vcf, i.e. the vcf file containing the variants to classify.
#'
#' @return Character vector with the mutation id, i.e. CHR:POS, of the C:G > T:A variants in \code{vcf_filename}.
#'
#' @examples
get_mut_id <- function(vcf_filename) {
  cmd <- sprintf("bcftools query -i 'FILTER=\"PASS\"' -f '%%CHROM:%%POS\n' %s", vcf_filename)
  mutation_ids <- system(cmd, intern = TRUE)
  return(mutation_ids)
}

#' Retrieve read depth
#'
#' This function retrieves the read depth at each C:G > T:A variant locus present in \code{vcf_filename}. This data is retrieved from \code{vcf_filename}.
#'
#' @param vcf_filename character string naming the path to the input vcf, i.e. the vcf file containing the variants to classify. This file must have been generated with Mutect2, either in tumor only or tumor/normal mode with strand bias annotation enabled.
#' @param samplename character string naming the sample in \code{vcf_filename}. This must match the name given to the sample when running Mutect2. It can be obtained by calling \code{get_samplename} on \code{vcf_filename}.
#' @param ... additional arguments not to be checked.
#'
#' @return Numeric vector with the read depth at each C:G > T:A variant locus present in \code{vcf_filename}, ordered by appearance in the file.
#' @import dplyr
#' @import tibble
#' @import tidyr
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

#' Retrieve allele frequency
#'
#' This function retrieves the allele frequency of each C:G > T:A variant present in \code{vcf_filename}. This data is retrieved from \code{vcf_filename}.
#'
#' @inheritParams get_depth
#'
#' @return Numeric vector with the allele frequency (AF) of each C:G > T:A variant present in \code{vcf_filename}, ordered by appearance in the file.
#'
#' @examples
get_AF <- function(vcf_filename, samplename, ...) {
  cmd <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%AF]\n'", samplename, vcf_filename)
  AF <- system(cmd, intern = TRUE)
  AF <- as.numeric(AF)
  return(AF)
}

#' Retrieve reference allele
#'
#' This function retrieves the reference allele of each C:G > T:A variant present in \code{vcf_filename}. This data is retrieved from \code{vcf_filename}.
#'
#' @param vcf_filename character string naming the path to the input vcf, i.e. the vcf file containing the variants to classify. This file must have been generated with Mutect2, either in tumor only or tumor/normal mode.
#' @param ... additional arguments not to be checked.
#'
#' @return Character vector with the reference allele of each C:G > T:A variant present in \code{vcf_filename}, ordered by appearance in the file.
#'
#' @examples
get_ref_allele <- function(vcf_filename, ...) {
  cmd <- sprintf("bcftools query -i 'FILTER=\"PASS\"' -f '%%REF\n' %s", vcf_filename)
  ref_alleles <- system(cmd, intern = TRUE)
  return(ref_alleles)
}

#' Retrieve alternate allele
#'
#' This function retrieves the alternate allele of each C:G > T:A variant present in \code{vcf_filename}. This data is retrieved from \code{vcf_filename}.
#'
#' @inheritParams get_ref_allele
#'
#' @return Character vector with the alternate allele of each C:G > T:A variant present in \code{vcf_filename}, ordered by appearance in the file.
#'
#' @examples
get_alt_allele <- function(vcf_filename, ...) {
  cmd <- sprintf("bcftools query -i 'FILTER=\"PASS\"' -f '%%ALT\n' %s", vcf_filename)
  alt_alleles <- system(cmd, intern = TRUE)
  return(alt_alleles)
}

#' Retrieve number of alternate alleles
#'
#' This function retrieves the number of reads that support the alternate allele of each C:G > T:A variant present in \code{vcf_filename}. This data is retrieved from \code{vcf_filename}.
#'
#' @inheritParams get_depth
#' @param depth numeric vector containing the read depth at each C:G > T:A variant locus
#' @param ... additional arguments not to be checked.
#'
#' @return List containing two numeric vectors: the number of reads that support the alternate allele and the number of reads that support the alternate allele divided by the total number of reads at that site. Each element in the vectors corresponds to a C:G > T:A variant in \code{vcf_filename} and the elements are ordered by appearance in \code{vcf_filename}.
#' @details \code{get_alt_bases} first retrieves the number of reads that support the alternate allele of each C:G > T:A variant present in \code{vcf_filename}. It then calculates the ratio between this value and the total number of reads at the site and reports both values, which is provided in \code{depth}.
#'
#' @examples
get_alt_bases <- function(vcf_filename, depth, samplename, ...) {
  cmd <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%AD{1}]\n'", samplename, vcf_filename)
  alt_bases <- system(cmd, intern = TRUE)
  alt_bases <- as.numeric(alt_bases)
  alt_base_ratio <- alt_bases/depth
  return(list(alt_bases, alt_base_ratio))
}

#' Retrieve number of reference alleles
#'
#' This function retrieves the number of reads that support the reference allele of each C:G > T:A variant present in \code{vcf_filename}. This data is retrieved from \code{vcf_filename}.
#'
#' @inheritParams get_alt_bases
#'
#' @return List containing two numeric vectors: the number of reads that support the reference allele and the number of reads that support the reference allele divided by the total number of reads at that site. Each element in the vectors corresponds to a C:G > T:A variant in \code{vcf_filename} and the elements are ordered by appearance in \code{vcf_filename}.
#' @details \code{get_ref_bases} first retrieves the number of reads that support the reference allele of each C:G > T:A variant present in \code{vcf_filename}. It then calculates the ratio between this value and the total number of reads at the site and reports both values, which is provided in \code{depth}.
#'
#' @examples
get_ref_bases <- function(vcf_filename, depth, samplename, ...) {
  cmd <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%AD{0}]\n'", samplename, vcf_filename)
  ref_bases <- system(cmd, intern = TRUE)
  ref_bases <- as.numeric(ref_bases)
  ref_base_ratio <- ref_bases/depth
  return(list(ref_bases, ref_base_ratio))
}

#' Retrieve base quality
#'
#' This function retrieves the median Phred quality of the bases that support the alternate allele of each C:G > T:A variant present in \code{vcf_filename}. This data is retrieved from \code{vcf_filename}.
#'
#' @inheritParams get_depth
#'
#' @return List containing two numeric vectors: the median Phred base quality of the alternate allele and the ratio between the median Phred base quality of the alternate allele and the median Phred base quality of the reference allele. Each element in the vectors corresponds to a C:G > T:A variant in \code{vcf_filename} and the elements are ordered by appearance in \code{vcf_filename}.
#' @import dplyr
#' @details \code{get_base_qual} first retrieves the median Phred quality of the bases that support the alternate allele of each C:G > T:A variant present in \code{vcf_filename}. It then retrieves this same value for the bases that support the reference allele of the same set of variants. It finally calculates the ratio between the two of them and reports two values: the median Phred quality of the bases that support the alternate allele and the calculated ratio.
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

#' Retrieve variant median fragment length
#'
#' This function retrieves the median fragment length of the reads that support the alternate allele of each C:G > T:A variant present in \code{vcf_filename}.This data is retrieved from \code{vcf_filename}.
#'
#' @inheritParams get_depth
#'
#' @return Numeric vector with the median median fragment length of each C:G > T:A variant present in \code{vcf_filename}, ordered by appearance in the file.
#'
#' @import dplyr
#'
#' @examples
get_fraglength <- function(vcf_filename, samplename, ...) {
  cmd <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%MFRL{1}]\n'", samplename, vcf_filename)
  # Median fragment length for reads that support the alternate allele
  fraglen <- system(cmd, intern = TRUE) %>%
    as.numeric()
  return(fraglen)
}

#' Retrieve variant distance from end of read
#'
#' \code{get_pos_from_end} retrieves, for each C:G > T:A variant present in \code{vcf_filename}, the median value of the position of the variant (starting from the end of the read) in those reads that support the alternate allele. This data is retrieved from \code{vcf_filename}.
#'
#' @inheritParams get_depth
#'
#' @return Numeric vector containing the median distance from end of read of each C:G > T:A variant in \code{vcf_filename}, ordered by appearance in the file.
#' @import dplyr
#'
#' @examples
get_pos_from_end <- function(vcf_filename, samplename, ...) {
  cmd <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%MPOS]\n'", samplename, vcf_filename)
  pos_from_end <- system(cmd, intern = TRUE) %>%
    as.numeric() # NAs will appear for multiallelic variants
  return(pos_from_end)
}

#' Retrieve mapping quality
#'
#' \code{get_map_qual} retrieves, for each C:G > T:A variant present in \code{vcf_filename}, the median Phred mapping quality (MAPQ) of the bases that support the alternate allele. This data is retrieved from \code{vcf_filename}.
#'
#' @inheritParams get_depth
#'
#' @return Numeric vector containing the MAPQ mapping quality of each C:G > T:A variant in \code{vcf_filename}, ordered by appearance in the file.
#' @import dplyr
#'
#' @examples
get_map_qual <- function(vcf_filename, samplename, ...) {
  cmd <- sprintf("bcftools view -s %s %s | bcftools query -i 'FILTER=\"PASS\"' -f '[%%MMQ]\n'", samplename, vcf_filename)
  map_qual <- system(cmd, intern = TRUE) %>%
    as.numeric() # NAs will appear for multiallelic variants
  return(map_qual)
}

#' Calculate read pair orientation bias
#'
#' \code{get_RPOB} calculates, for each C:G > T:A variant present in \code{vcf_filename}, two read pair orientation bias (RPOB) metrics: FdeamC and strand orientation bias (SOB). The calculation is based on data retrieved from \code{vcf_filename}.
#'
#' @inheritParams get_depth
#'
#' @return List containing two numeric vectors: FDeamC and SOB. Each element in the vectors corresponds to a C:G > T:A variant in \code{vcf_filename} and the elements are ordered by appearance in \code{vcf_filename}.
#'
#' @import dplyr
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

#' Calculate strand bias
#'
#' This function calculates, for each C:G > T:A variant present in \code{vcf_filename}, two strand bias (SB) metrics: SB defined by GATK (SB-GATK) and SB defined by Guo et. al. (SB-Guo). The calculation is based on data retrieved from \code{vcf_filename}.
#'
#' @inheritParams get_depth
#'
#' @return List of two numeric vectors. Each element in the vectors corresponds to a C:G > T:A variant in \code{vcf_filename}. The elements are ordered by appearance in \code{vcf_filename}. The first vector stores SB-GATK values and the second one SB-Guo values.
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @details See \code{\link{calc_SB_Guo}} and \code{\link{calc_SB_GATK}} for further detail on these SB values.
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
#' This function calculates variant strand bias as defined by Guo et. al. (SB-Guo).
#'
#' @param a Number of reads supporting the reference allele mapping to the forward strand.
#' @param b Number of reads supporting the alternate allele mapping to the forward strand.
#' @param e Number of reads supporting the reference allele mapping to the forward strand.
#' @param d Number of reads supporting the alternate allele mapping to the reverse strand.
#'
#' @return Numeric element with the SB-Guo value of the variant.
#' @details This strand bias value was defined in \emph{The effect of strand bias in Illumina short-read sequencing data}. It can be accessed here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3532123/. It ranges from 0 to infinity.
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
#' This function calculates variant strand bias as defined by GATK (SB-GATK).
#'
#' @inheritParams calc_SB_Guo
#'
#' @return Numeric element with the SB-GATK value of the variant.
#' @details This strand bias value was defined by GATK and described in \emph{The effect of strand bias in Illumina short-read sequencing data}. It can be accessed here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3532123/. It ranges from 0 to infinity.
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
#' \code{check_lowercase} checks whether a character is lowercase.
#'
#' @param char character string
#'
#' @return Logical value indicating whether \strong{char} is lowercase.
#'
#' @examples
check_lowercase <- function(char) {
  char == tolower(char)
}
