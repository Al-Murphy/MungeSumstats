# Import full genome-wide GWAS summary statistics from Open GWAS

Requires internet access to run.

## Usage

``` r
import_sumstats(
  ids,
  vcf_dir = tempdir(),
  vcf_download = TRUE,
  save_dir = tempdir(),
  write_vcf = FALSE,
  download_method = "download.file",
  quiet = TRUE,
  force_new = FALSE,
  force_new_vcf = FALSE,
  nThread = 1,
  parallel_across_ids = FALSE,
  ...
)
```

## Arguments

- ids:

  List of Open GWAS study IDs (e.g. `c("prot-a-664", "ieu-b-4760")`).

- vcf_dir:

  Where to download the original VCF from Open GWAS. *WARNING:* This is
  set to [`tempdir()`](https://rdrr.io/r/base/tempfile.html) by default.
  This means the raw (pre-formatted) VCFs be deleted upon ending the R
  session. Change this to keep the raw VCF file on disk (e.g.
  `vcf_dir="./raw_vcf"`).

- vcf_download:

  Download the original VCF from Open GWAS.

- save_dir:

  Directory to save formatted summary statistics in.

- write_vcf:

  Whether to write as VCF (TRUE) or tabular file (FALSE).

- download_method:

  `"axel"` (multi-threaded) or `"download.file"` (single-threaded) .

- quiet:

  Run quietly.

- force_new:

  If a formatted file of the same names as `save_path` exists,
  formatting will be skipped and this file will be imported instead
  (default). Set `force_new=TRUE` to override this.

- force_new_vcf:

  Overwrite a previously downloaded VCF with the same path name.

- nThread:

  Number of threads to use for parallel processes.

- parallel_across_ids:

  If `parallel_across_ids=TRUE` and `nThread>1`, then each ID in `ids`
  will be processed in parallel.

- ...:

  Arguments passed on to
  [`format_sumstats`](https://al-murphy.github.io/MungeSumstats/reference/format_sumstats.md)

  `path`

  :   Filepath for the summary statistics file to be formatted. A
      dataframe or datatable of the summary statistics file can also be
      passed directly to MungeSumstats using the path parameter.

  `ref_genome`

  :   name of the reference genome used for the GWAS ("GRCh37" or
      "GRCh38"). Argument is case-insensitive. Default is NULL which
      infers the reference genome from the data.

  `convert_ref_genome`

  :   name of the reference genome to convert to ("GRCh37" or "GRCh38").
      This will only occur if the current genome build does not match.
      Default is not to convert the genome build (NULL).

  `chain_source`

  :   source of the chain file to use in liftover, if converting genome
      build ("ucsc" or "ensembl"). Note that the UCSC chain files
      require a license for commercial use. The Ensembl chain is used by
      default ("ensembl").

  `local_chain`

  :   Path to local chain file to use instead of downlaoding. Default of
      NULL i.e. no local file to use. NOTE if passing a local chain file
      make sure to specify the path to convert from and to the correct
      build like GRCh37 to GRCh38. We can not sense check this for local
      files. The chain file can be submitted as a gz file (as downloaed
      from source) or unzipped.

  `convert_small_p`

  :   Binary, should non-negative p-values \<= 5e-324 be converted to 0?
      Small p-values pass the R limit and can cause errors with
      LDSC/MAGMA and should be converted. Default is TRUE.

  `convert_large_p`

  :   Binary, should p-values \>1 be converted to 1? P-values \>1 should
      not be possible and can cause errors with LDSC/MAGMA and should be
      converted. Default is TRUE.

  `convert_neg_p`

  :   Binary, should p-values \<0 be converted to 0? Negative p-values
      should not be possible and can cause errors with LDSC/MAGMA and
      should be converted. Default is TRUE.

  `compute_z`

  :   Whether to compute Z-score column. Default is FALSE. This can be
      computed from Beta and SE with (Beta/SE) or P
      (Z:=sign(BETA)\*sqrt(stats::qchisq(P,1,lower=FALSE))). **Note**
      that imputing the Z-score from P for every SNP will not be
      perfectly correct and may result in a loss of power. This should
      only be done as a last resort. Use 'BETA' to impute by BETA/SE and
      'P' to impute by SNP p-value.

  `force_new_z`

  :   When a "Z" column already exists, it will be used by default. To
      override and compute a new Z-score column from P set
      `force_new_z=TRUE`.

  `compute_n`

  :   Whether to impute N. Default of 0 won't impute, any other integer
      will be imputed as the N (sample size) for every SNP in the
      dataset. **Note** that imputing the sample size for every SNP is
      not correct and should only be done as a last resort. N can also
      be inputted with "ldsc", "sum", "giant" or "metal" by passing one
      of these for this field or a vector of multiple. Sum and an
      integer value creates an N column in the output whereas giant,
      metal or ldsc create an Neff or effective sample size. If
      multiples are passed, the formula used to derive it will be
      indicated.

  `convert_n_int`

  :   Binary, if N (the number of samples) is not an integer, should
      this be rounded? Default is TRUE.

  `impute_beta`

  :   Binary, whether BETA should be imputed using other effect data if
      it isn't present in the sumstats. Note that this imputation is an
      approximation (for Z & SE approach) so could have an effect on
      downstream analysis. Use with caution. The different methods
      MungeSumstats will try and impute beta (in this order or priority)
      are:

      1.  log(OR) 2. Z x SE Default value is FALSE.

  `es_is_beta`

  :   Binary, whether to map ES to BETA. We take BETA to be any
      BETA-like value (including Effect Size). If this is not the case
      for your sumstats, change this to FALSE. Default is TRUE.

  `impute_se`

  :   Binary, whether the standard error should be imputed using other
      effect data if it isn't present in the sumstats. Note that this
      imputation is an approximation so could have an effect on
      downstream analysis. Use with caution. The different methods
      MungeSumstats will try and impute se (in this order or priority)
      are:

      1.  BETA / Z 2. abs(BETA/ qnorm(P/2)) Default is FALSE.

  `analysis_trait`

  :   If multiple traits were studied, name of the trait for analysis
      from the GWAS. Default is NULL.

  `ignore_multi_trait`

  :   If you have multiple traits (p-values) in the study but you want
      to ignorwe these and instead use a standard named p-value, set to
      TRUE. By default is FALSE which will check for multi-traits.

  `INFO_filter`

  :   numeric The minimum value permissible of the imputation
      information score (if present in sumstats file). Default 0.9.

  `FRQ_filter`

  :   numeric The minimum value permissible of the frequency(FRQ) of the
      SNP (i.e. Allele Frequency (AF)) (if present in sumstats file). By
      default no filtering is done, i.e. value of 0.

  `pos_se`

  :   Binary Should the standard Error (SE) column be checked to ensure
      it is greater than 0? Those that are, are removed (if present in
      sumstats file). Default TRUE.

  `effect_columns_nonzero`

  :   Binary should the effect columns in the data BETA,OR (odds
      ratio),LOG_ODDS,SIGNED_SUMSTAT be checked to ensure no SNP=0.
      Those that do are removed(if present in sumstats file). Default
      FALSE.

  `N_std`

  :   numeric The number of standard deviations above the mean a SNP's N
      is needed to be removed. Default is 5.

  `N_dropNA`

  :   Drop rows where N is missing.Default is TRUE.

  `chr_style`

  :   Chromosome naming style to use in the formatted summary statistics
      file ("NCBI", "UCSC", "dbSNP", or "Ensembl"). The NCBI and Ensembl
      styles both code chromosomes as `1-22, X, Y, MT`; the UCSC style
      is `chr1-chr22, chrX, chrY, chrM`; and the dbSNP style is
      `ch1-ch22, chX, chY, chMT`. Default is Ensembl.

  `rmv_chrPrefix`

  :   Is now deprecated, do. not use. Use chr_style instead - chr_style
      = 'Ensembl' will give the same result as rmv_chrPrefix=TRUE used
      to give.

  `rmv_chr`

  :   Chromosomes to exclude from the formatted summary statistics file.
      Use NULL if no filtering is necessary. Default is
      `c("X", "Y", "MT")` which removes all non-autosomal SNPs.

  `on_ref_genome`

  :   Binary Should a check take place that all SNPs are on the
      reference genome by SNP ID. Default is TRUE.

  `infer_eff_direction`

  :   Binary Should a check take place to ensure the alleles match the
      effect direction? Default is TRUE.

  `eff_on_minor_alleles`

  :   Binary Should MungeSumstats assume that the effects are
      majoritively measured on the minor alleles? Default is FALSE as
      this is an assumption that won't be appropriate in all cases.
      However, the benefit is that if we know the majority of SNPs have
      their effects based on the minor alleles, we can catch cases where
      the allele columns have been mislabelled.

  `strand_ambig_filter`

  :   Binary Should SNPs with strand-ambiguous alleles be removed.
      Default is FALSE.

  `allele_flip_check`

  :   Binary Should the allele columns be checked against reference
      genome to infer if flipping is necessary. Default is TRUE.

  `allele_flip_drop`

  :   Binary Should the SNPs for which neither their A1 or A2 base pair
      values match a reference genome be dropped. Default is TRUE.

  `allele_flip_z`

  :   Binary should the Z-score be flipped along with effect and FRQ
      columns like Beta? It is assumed to be calculated off the effect
      size not the P-value and so will be flipped i.e. default TRUE.

  `allele_flip_frq`

  :   Binary should the frequency (FRQ) column be flipped along with
      effect and z-score columns like Beta? Default TRUE.

  `bi_allelic_filter`

  :   Binary Should non-bi-allelic SNPs be removed. Default is TRUE.

  `flip_frq_as_biallelic`

  :   Binary Should non-bi-allelic SNPs frequency values be flipped as
      1-p despite there being other alternative alleles? Default is
      FALSE but if set to TRUE, this allows non-bi-allelic SNPs to be
      kept despite needing flipping.

  `snp_ids_are_rs_ids`

  :   Binary Should the supplied SNP ID's be assumed to be RSIDs. If
      not, imputation using the SNP ID for other columns like base-pair
      position or chromosome will not be possible. If set to FALSE, the
      SNP RS ID will be imputed from the reference genome if possible.
      Default is TRUE.

  `remove_multi_rs_snp`

  :   Binary Sometimes summary statistics can have multiple RSIDs on one
      row (i.e. related to one SNP), for example
      "rs5772025_rs397784053". This can cause an error so by default,
      the first RS ID will be kept and the rest removed e.g."rs5772025".
      If you want to just remove these SNPs entirely, set it to TRUE.
      Default is FALSE.

  `frq_is_maf`

  :   Conventionally the FRQ column is intended to show the minor/effect
      allele frequency (MAF) but sometimes the major allele frequency
      can be inferred as the FRQ column. This logical variable indicates
      that the FRQ column should be renamed to MAJOR_ALLELE_FRQ if the
      frequency values appear to relate to the major allele i.e. \>0.5.
      By default this mapping won't occur i.e. is TRUE.

  `indels`

  :   Binary does your Sumstats file contain Indels? These don't exist
      in our reference file so they will be excluded from checks if this
      value is TRUE. Default is TRUE.

  `drop_indels`

  :   Binary, should any indels found in the sumstats be dropped? These
      can not be checked against a reference dataset and will have the
      same RS ID and position as SNPs which can affect downstream
      analysis. Default is False.

  `drop_na_cols`

  :   A character vector of column names to be checked for missing
      values. Rows with missing values in any of these columns (if
      present in the dataset) will be dropped. If `NULL`, all columns
      will be checked for missing values. Default columns are SNP,
      chromosome, position, allele 1, allele2, effect columns
      (frequency, beta, Z-score, standard error, log odds, signed
      sumstats, odds ratio), p value and N columns.

  `dbSNP`

  :   version of dbSNP to be used for imputation (144 or 155). See
      `dbSNP_tarball` for different versions of dbSNP (including newer
      releases).

  `dbSNP_tarball`

  :   Pass local versions of dbSNP in tarball format. Default of NULL
      uses the dbSNP version passed in `dbSNP` parmeter. `dbSNP_tarball`
      was enabled to help with dbSNP versions \>=156, after the decision
      to no longer provide dbSNP releases as bioconductor packages.
      dbSNP 156 tarball is available here:
      http://149.165.171.124/SNPlocs/.

  `check_dups`

  :   whether to check for duplicates - if formatting QTL datasets this
      should be set to FALSE otherwise keep as TRUE. Default is TRUE.

  `sort_coordinates`

  :   Whether to sort by coordinates of resulting sumstats

  `save_path`

  :   File path to save formatted data. Defaults to
      `tempfile(fileext=".tsv.gz")`.

  `tabix_index`

  :   Index the formatted summary statistics with
      [tabix](http://www.htslib.org/doc/tabix.md) for fast querying.

  `return_data`

  :   Return `data.table`, `GRanges` or `VRanges` directly to user.
      Otherwise, return the path to the save data. Default is FALSE.

  `return_format`

  :   If return_data is TRUE. Object type to be returned
      ("data.table","vranges","granges").

  `ldsc_format`

  :   DEPRECATED, do not use. Use save_format="LDSC" instead.

  `save_format`

  :   Output format of sumstats. Options are NULL - standardised output
      format from MungeSumstats, LDSC - output format compatible with
      LDSC and openGWAS - output compatible with openGWAS VCFs. Default
      is NULL. **NOTE** - If LDSC format is used, the naming convention
      of A1 as the reference (genome build) allele and A2 as the effect
      allele will be reversed to match LDSC (A1 will now be the effect
      allele). See more info on this
      [here](https://groups.google.com/g/ldsc_users/c/S7FZK743w68). Note
      that any effect columns (e.g. Z) will be inrelation to A1 now
      instead of A2.

  `log_folder_ind`

  :   Binary Should log files be stored containing all filtered out SNPs
      (separate file per filter). The data is outputted in the same
      format specified for the resulting sumstats file. The only
      exception to this rule is if output is vcf, then log file saved as
      .tsv.gz. Default is FALSE.

  `log_mungesumstats_msgs`

  :   Binary Should a log be stored containing all messages and errors
      printed by MungeSumstats in a run. Default is FALSE

  `log_folder`

  :   Filepath to the directory for the log files and the log of
      MungeSumstats messages to be stored. Default is a temporary
      directory. Note the name of the log files (log messages and log
      outputs) are now the same as the name of the file specified in the
      save path parameter with the extension '\_log_msg.txt' and
      '\_log_output.txt' respectively.

  `imputation_ind`

  :   Binary Should a column be added for each imputation step to show
      what SNPs have imputed values for differing fields. This includes
      a field denoting SNP allele flipping (flipped). On the flipped
      value, this denoted whether the alelles where switched based on
      MungeSumstats initial choice of A1, A2 from the input column
      headers and thus may not align with what the creator
      intended.**Note** these columns will be in the formatted summary
      statistics returned. Default is FALSE.

  `mapping_file`

  :   MungeSumstats has a pre-defined column-name mapping file which
      should cover the most common column headers and their
      interpretations. However, if a column header that is in youf file
      is missing of the mapping we give is incorrect you can supply your
      own mapping file. Must be a 2 column dataframe with column names
      "Uncorrected" and "Corrected". See data(sumstatsColHeaders) for
      default mapping and necessary format.

## Value

Either a named list of data objects or paths, depending on the arguments
passed to `format_sumstats`.

## Examples

``` r
#only run the examples if user has internet access:
if(try(is.character(getURL("www.google.com")))==TRUE){
### Search by criteria
metagwas <- find_sumstats(
    traits = c("parkinson", "alzheimer"),
    min_sample_size = 5000
)
### Only use a subset for testing purposes
ids <- (dplyr::arrange(metagwas, nsnp))$id

### Default usage
## You can supply \code{import_sumstats()}
## with a list of as many OpenGWAS IDs as you want,
## but we'll just give one to save time.

## Call uses reference genome as default with more than 2GB of memory,
## which is more than what 32-bit Windows can handle so remove certain checks
## commented out down to runtime
# datasets <- import_sumstats(ids = ids[1])
}
#> Error in getURL("www.google.com") : could not find function "getURL"
```
