#' Classify variants
#'
#' \code{run_ideafix} classifies the C:G > T:A variants found in \code{vcf_filename} into deaminations or non-deaminations.
#'
#' @param vcf_filename character string naming the path to the input vcf, i.e. the vcf file containing the variants to classify. This file must have been generated with Mutect2, either in tumor only or tumor/normal mode with strand bias annotation enabled.
#' @param fasta_filename character string naming the path to the reference genome FASTA file the sequencing data was aligned to.
#' @param algorithm character string naming the algorithm to use to classify the variants. Can be "RF" or "XGBoost". Defaults to "RF".
#' @param outformat character string indicating the output file format. Can be "tsv" or "vcf". Defaults to "vcf".
#' @param outfolder character string naming the folder to write the output file to. Defaults to current working directory (\code{getwd}).
#' @param outname character string naming the output filename. Defaults to "ideafix_labels.txt" or "ideafix_labels.vcf", depending on \code{outformat} value .
#'
#' @return
#' @export
#' @details \code{run_ideafix} classifies the C:G > T:A variants found in \code{vcf_filename} into deaminations or non-deaminations using an RF or an XGBoost model. It first extracts a set of descriptors into a tibble either from \code{vcf_filename} or \code{fasta_filename}, or builds them using data retrieved from them using \code{get_descriptors}. These descriptors include VAF, number of alternate bases, normalized number of alternate bases, number of reference bases, normalized number of reference bases, reference allele, alternate allele, base quality, base quality fraction, fragment length, median position from read end, mapping quality, FDeamC, SOB, SB-GUO, SB-GATK, normalized median position from read end, base two positions before, base one position before, base two positions after, base one position after, dinucleotide before and dinucleotide after. Then, it classifies the variants based on those descriptors using the function \code{classify_variants}. These results are finally output to a text file, either to a new tsv file or appended to the input vcf file \code{vcf_filename}.
#'
#' @examples
run_ideafix <- function(vcf_filename, fasta_filename, algorithm = "RF", outformat = "vcf", outfolder = ".", outname = "ideafix_labels") {
  descriptors <- get_descriptors(vcf_filename = vcf_filename, fasta_filename = fasta_filename)
  predictions <- classify_variants(variant_descriptors = descriptors, algorithm = algorithm)
  annotate_deaminations(classification = predictions, format = outformat, outname = file.path(outfolder, outname), vcf_filename = vcf_filename)
  return(predictions)
}
