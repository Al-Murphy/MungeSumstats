% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_on_ref_genome.R
\name{check_on_ref_genome}
\alias{check_on_ref_genome}
\title{Ensure all SNPs are on the reference genome}
\usage{
check_on_ref_genome(
  sumstats_dt,
  path,
  ref_genome,
  on_ref_genome,
  indels = indels,
  rsids,
  imputation_ind,
  log_folder_ind,
  check_save_out,
  tabix_index,
  nThread,
  log_files,
  dbSNP,
  dbSNP_tarball
)
}
\arguments{
\item{path}{Filepath for the summary statistics file to be formatted. A
dataframe or datatable of the summary statistics file can also be passed
directly to MungeSumstats using the path parameter.}

\item{ref_genome}{name of the reference genome used for the GWAS ("GRCh37" or
"GRCh38"). Argument is case-insensitive. Default is NULL which infers the
reference genome from the data.}

\item{on_ref_genome}{Binary Should a check take place that all SNPs are on
the reference genome by SNP ID. Default is TRUE.}

\item{indels}{Binary does your Sumstats file contain Indels? These don't
exist in our reference file so they will be excluded from checks if this
value is TRUE. Default is TRUE.}

\item{imputation_ind}{Binary Should a column be added for each imputation
step to show what SNPs have imputed values for differing fields. This
includes a field denoting SNP allele flipping (flipped). On the flipped
value, this denoted whether the alelles where switched based on
MungeSumstats initial choice of A1, A2 from the input column headers and thus
may not align with what the creator intended.\strong{Note} these columns will be
in the formatted summary statistics returned. Default is FALSE.}

\item{log_folder_ind}{Binary Should log files be stored containing all
filtered out SNPs (separate file per filter). The data is outputted in the
same format specified for the resulting sumstats file. The only exception to
this rule is if output is vcf, then log file saved as .tsv.gz. Default is
FALSE.}

\item{tabix_index}{Index the formatted summary statistics with
\href{http://www.htslib.org/doc/tabix.html}{tabix} for fast querying.}

\item{nThread}{Number of threads to use for parallel processes.}

\item{log_files}{list of log file locations}

\item{dbSNP}{version of dbSNP to be used for imputation (144 or 155). See
\code{dbSNP_tarball} for different versions of dbSNP (including newer releases).}

\item{dbSNP_tarball}{Pass local versions of dbSNP in tarball format. Default
of NULL uses the dbSNP version passed in \code{dbSNP} parmeter. \code{dbSNP_tarball}
was enabled to help with dbSNP versions >=156, after the decision to no
longer provide dbSNP releases as bioconductor packages. dbSNP 156 tarball is
available here: http://149.165.171.124/SNPlocs/.}
}
\value{
A list containing two data tables:
\itemize{
\item \code{sumstats_dt}{
  : the modified summary statistics data table object}
\item \code{rsids}{
  : snpsById, filtered to SNPs of interest if loaded already. Or else NULL}
\item \code{log_files}{
  : log file list}
}
}
\description{
Ensure all SNPs are on the reference genome
}
\keyword{internal}
