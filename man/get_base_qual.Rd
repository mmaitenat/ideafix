% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_descriptors.R
\name{get_base_qual}
\alias{get_base_qual}
\title{Retrieve base quality}
\usage{
get_base_qual(vcf_filename, samplename, ...)
}
\arguments{
\item{vcf_filename}{character string naming the path to the input vcf, i.e. the vcf file containing the variants to classify. This file must have been generated with Mutect2, either in tumor only or tumor/normal mode with strand bias annotation enabled.}

\item{samplename}{character string naming the sample in \code{vcf_filename}. This must match the name given to the sample when running Mutect2. It can be obtained by calling \code{get_samplename} on \code{vcf_filename}.}

\item{...}{additional arguments not to be checked.}
}
\value{
List containing two numeric vectors: the median Phred base quality of the alternate allele and the ratio between the median Phred base quality of the alternate allele and the median Phred base quality of the reference allele. Each element in the vectors corresponds to a C:G > T:A variant in \code{vcf_filename} and the elements are ordered by appearance in \code{vcf_filename}.
}
\description{
This function retrieves the median Phred quality of the bases that support the alternate allele of each C:G > T:A variant present in \code{vcf_filename}. This data is retrieved from \code{vcf_filename}.
}
\details{
\code{get_base_qual} first retrieves the median Phred quality of the bases that support the alternate allele of each C:G > T:A variant present in \code{vcf_filename}. It then retrieves this same value for the bases that support the reference allele of the same set of variants. It finally calculates the ratio between the two of them and reports two values: the median Phred quality of the bases that support the alternate allele and the calculated ratio.
}
