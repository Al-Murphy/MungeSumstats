# \`MungeSumstats\`: Getting started

## Citation

If you use the *MungeSumstats* package, please cite

[Murphy et al. MungeSumstats: A Bioconductor package for the
standardisation and quality control of many GWAS summary
statistics](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btab665/6380562).

## Overview

The *MungeSumstats* package is designed to facilitate the
standardisation of GWAS summary statistics as utilised in our Nature
Genetics paper¹.

The package is designed to handle the lack of standardisation of output
files by the GWAS community. There is a group who have now manually
standardised many GWAS: [R interface to the IEU GWAS database API •
ieugwasr](https://mrcieu.github.io/ieugwasr/) and
[gwasvcf](https://github.com/MRCIEU/gwasvcf) but because a lot of GWAS
remain closed access, these repositories are not all encompassing.

The [GWAS-Download
project](https://github.com/mikegloudemans/gwas-download) has collated
summary statistics from 200+ GWAS. This repository has been utilsed to
identify the most common formats, all of which can be standardised with
*MungeSumstats*.

Moreover, there is an emerging standard of VCF format for summary
statistics files with multiple, useful, associated R packages such as
*vcfR*. However, there is currently no method to convert VCF formats to
a standardised format that matches older approaches.

The *MungeSumstats* package standardises both VCF and the most common
summary statistic file formats to enable downstream integration and
analysis.

*MungeSumstats* also offers comprehensive Quality Control (QC) steps
which are important prerequisites for downstream analysis like Linkage
disequilibrium score regression (LDSC) and MAGMA.

Moreover, *MungeSumstats* is efficiently written resulting in all
reformatting and quality control checks completing in minutes for GWAS
summary statistics with 500k SNPs on a standard desktop machine. This
speed can be increased further by increasing the number of threads
(nThread) for `data.table` to use.

Currently *MungeSumstats* only works on data from humans, as it uses
human-based genome references.

## Aim

*MungeSumstats* will ensure that the all essential columns for analysis
are present and syntactically correct. Generally, summary statistic
files include (but are not limited to) the columns:

- SNP : SNP ID (rs IDs)
- CHR : Chromosome number
- BP : Base pair positions
- A1 : reference allele
- A2 : alternative allele
- Z : Z-score
- BETA : Effect size estimate relative to the alternative allele
- P : Unadjusted p-value for SNP
- SE : The standard error
- N : Sample size
- INFO: The imputation information score
- FRQ: The minor/effect allele frequency (MAF/EAF) of the SNP

*MungeSumstats* uses a mapping file to infer the inputted column names
(run `data("sumstatsColHeaders")` to view these). This mapping file is
far more comprehensive than any other publicly available munging tool
containing more than 200 unique mappings at the time of writing this
vignette. However, if your column headers are missing or if you want to
change the mapping, you can do so by passing your own mapping file (see
`format_sumstats(mapping_file)`).

*MungeSumstats* offers unmatched levels of quality control to ensure,
for example, consistency of allele assignment and direction of effects.
Tests run by *MungeSumstats* include:

- Check VCF format
- Check tab, space or comma delimited, zipped, csv or tsv file
- Check for header name synonyms
- Check for multiple models or traits in GWAS
- Check for uniformity in SNP ID - no mix of rs/missing rs/chr:bp
- Check for CHR:BP:A2:A1 all in one column
- Check for CHR:BP in one column
- Check for A1/A2 in one column
- Check if CHR and/or BP is missing (infer from reference genome)
- Check if SNP ID is missing (infer from reference genome)
- Check if A1 and/or A2 are missing (infer from reference genome)
- Check that vital columns are present (SNP,CHR,BP,P,A1,A2)
- Check for one signed/effect column (Z,OR,BETA,LOG_ODDS,SIGNED_SUMSTAT)
- Check for missing data
- Check for duplicated columns
- Check for small p-values (lower than 5e-324)
- Check N column is an integer
- Check for SNPs with N greater than 5 times standard dev. plus the mean
- Check SNPs are RS ID’s
- Check for uniformity of SNP ID format
- Check for duplicated rows, based on SNP ID
- Check for duplicated rows, based on base-pair position
- Check for SNPs on reference genome. Correct not found SNP IDs using
  CHR and BP (infer from reference genome)
- Check INFO score
- Check FRQ value
- Check FRQ is minor allele frequency (MAF)
- Check that the SNPs’ standard error (SE) is positive
- Check that SNPs’ effect columns (like BETA) aren’t equal to 0
- Check for strand-ambiguous SNPs
- Check for non-biallelic SNPs (infer from reference genome)
- Check for allele flipping
- Check for SNPs with nonstandard chromosome names
- Check for CHR format of X/Y chromosome - converting 23/24 to X/Y and M
  to MT
- Check for SNPs on excluded chromosomes (removes non-autosomal SNPs by
  default)
- Check for z-score (Z) and impute if missing
- Check for N and impute if missing
- Check output format is LDSC ready
- Check output format is IEU OpenGWAS ready
- Check and perform liftover to desired reference genome if necessary
- Check for indels in the sumstats and drop them if found (not run by
  default)

Users can specify which checks to run on their data. A **note** on the
allele flipping check: **MungeSumstats** infers the effect allele will
always be the A2 allele, this is the approach done for [IEU GWAS
VCF](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7805039/) and has such
also been adopted here. This inference is first from the inputted file’s
column headers however, the allele flipping check ensures this by
comparing A1, what should be the reference allele, to the reference
genome. If a SNP’s A1 DNA base doesn’t match the reference genome but
it’s A2 (what should be the alternative allele) does, the alleles will
be flipped along with the effect information (e.g. Beta, Odds Ratio,
signed summary statistics, FRQ, Z-score\*).

\*-by default the Z-score is assumed to be calculated off the effect
size not the P-value and so will be flipped. This can be changed by a
user.

If a test is failed, the user will be notified and if possible, the
input will be corrected. The QC steps from the checks above can also be
adjusted to suit the user’s analysis, see
[`MungeSumstats::format_sumstats`](https://al-murphy.github.io/MungeSumstats/reference/format_sumstats.md).

*MungeSumstats* can handle VCF, txt, tsv, csv file types or .gz/.bgz
versions of these file types. The package also gives the user the
flexibility to export the reformatted file as tab-delimited, VCF or R
native objects such as data.table, GRanges or VRanges objects. The
output can also be outputted in an **LDSC ready** format which means the
file can be fed directly into LDSC without the need for additional
munging. **NOTE** - If LDSC format is used, the naming convention of A1
as the reference (genome build) allele and A2 as the effect allele will
be reversed to match LDSC (A1 will now be the effect allele). See more
info on this
[here](https://groups.google.com/g/ldsc_users/c/S7FZK743w68). Note that
any effect columns (e.g. Z) will be inrelation to A1 now instead of A2.

Please read carefully through our [FAQ
Website](https://github.com/Al-Murphy/MungeSumstats/wiki/FAQ) to gain
insight on how best to run MungeSumstats on your data.

## Data

The *MungeSumstats* package contains small subsets of GWAS summary
statistics files. Firstly, on Educational Attainment by Okbay et al
2016: PMID: 27898078 PMCID: PMC5509058 DOI: 10.1038/ng1216-1587b.

Secondly, a VCF file (VCFv4.2) relating to the GWAS Amyotrophic lateral
sclerosis from ieu open GWAS project. Dataset: ebi-a-GCST005647:
<https://gwas.mrcieu.ac.uk/datasets/ebi-a-GCST005647/>

These datasets will be used to showcase *MungeSumstats* functionality.

## Running *MungeSumstats*

*MungeSumstats* is available on Bioconductor. To install the package on
Bioconductor run the following lines of code:

    if (!require("BiocManager")) install.packages("BiocManager")
    BiocManager::install("MungeSumstats")

Once installed, load the package:

``` r

library(MungeSumstats)
```

To standardise the summary statistics’ file format, simply call
[`format_sumstats()`](https://al-murphy.github.io/MungeSumstats/reference/format_sumstats.md)
passing in the path to your summary statistics file or directly pass the
summary statistics as a dataframe or datatable. You can specify which
genome build was used in the GWAS(GRCh37 or GRCh38) or, as default,
infer the genome build from the data.The reference genome is used for
multiple checks like deriving missing data such SNP/BP/CHR/A1/A2 and for
QC steps like removing non-biallelic SNPs, strand-ambiguous SNPs or
ensuring correct allele and direction of SNP effects. The path to the
reformatted summary statistics file can be returned by the function
call, the user can specify a location to save the file or the user can
return an R native object for the data: data.table, VRanges or GRanges
object.

Note that for a number of the checks implored by *MungeSumstats* a
reference genome is used. If your GWAS summary statistics file of
interest relates to *GRCh38*, you will need to install
`SNPlocs.Hsapiens.dbSNP155.GRCh38` and `BSgenome.Hsapiens.NCBI.GRCh38`
from Bioconductor as follows:

    #increase permissible time to download data, in case of slow internet access
    options(timeout=2000)
    BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38")
    BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")

If your GWAS summary statistics file of interest relates to *GRCh37*,
you will need to install `SNPlocs.Hsapiens.dbSNP155.GRCh37` and
`BSgenome.Hsapiens.1000genomes.hs37d5` from Bioconductor as follows:

    BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh37")
    BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")

These may take some time to install and are not included in the package
as some users may only need one of *GRCh37*/*GRCh38*.

The Educational Attainment by Okbay GWAS summary statistics file is
saved as a text document in the package’s external data folder so we can
just pass the file path to it straight to *MungeSumstats*.

**NOTE** - By default, Formatted results will be saved to
[`tempdir()`](https://rdrr.io/r/base/tempfile.html). This means all
formatted summary stats will be deleted upon ending the R session if not
copied to a local file path. Otherwise, to keep formatted summary stats,
change `save_path` ( e.g.`file.path('./formatted',basename(path))`), or
make sure to copy files elsewhere after processing (
e.g.`file.copy(save_path, './formatted/' )`.

``` r

eduAttainOkbayPth <- system.file("extdata","eduAttainOkbay.txt",
                                  package="MungeSumstats")
reformatted <- 
  MungeSumstats::format_sumstats(path=eduAttainOkbayPth,
                                 ref_genome="GRCh37")
```

    ## 
    ## 
    ## ******::NOTE::******
    ##  - Formatted results will be saved to `tempdir()` by default.
    ##  - This means all formatted summary stats will be deleted upon ending the R session.
    ##  - To keep formatted summary stats, change `save_path`  ( e.g. `save_path=file.path('./formatted',basename(path))` ),   or make sure to copy files elsewhere after processing  ( e.g. `file.copy(save_path, './formatted/' )`.
    ##  ********************

    ## Formatted summary statistics will be saved to ==>  /tmp/Rtmpt6oeYw/file393568a0e909.tsv.gz

    ## Warning: replacing previous import 'utils::findMatches' by
    ## 'S4Vectors::findMatches' when loading 'SNPlocs.Hsapiens.dbSNP155.GRCh37'

    ## Importing tabular file: /home/runner/work/_temp/Library/MungeSumstats/extdata/eduAttainOkbay.txt

    ## Checking for empty columns.

    ## Infer Effect Column

    ## First line of summary statistics file:

    ## MarkerName   CHR POS A1  A2  EAF Beta    SE  Pval    

    ## Allele columns are ambiguous, attempting to infer direction

    ## Can't infer allele columns from sumstats

    ## Standardising column headers.

    ## First line of summary statistics file:

    ## MarkerName   CHR POS A1  A2  EAF Beta    SE  Pval    

    ## Summary statistics report:
    ##    - 93 rows
    ##    - 93 unique variants
    ##    - 70 genome-wide significant variants (P<5e-8)
    ##    - 20 chromosomes

    ## Checking for multi-GWAS.

    ## Checking for multiple RSIDs on one row.

    ## Checking SNP RSIDs.

    ## Checking for merged allele column.

    ## Checking A1 is uppercase

    ## Checking A2 is uppercase

    ## Checking for incorrect base-pair positions

    ## Checking for missing data.

    ## Checking for duplicate columns.

    ## Checking for duplicated rows.

    ## INFO column not available. Skipping INFO score filtering step.

    ## Filtering SNPs, ensuring SE>0.

    ## Ensuring all SNPs have N<5 std dev above mean.

    ## Formatting CHR column (if necessary)- Removing chr prefix, converting 23/24 to X/Y and M to MT

    ## 47 SNPs (50.5%) have FRQ values > 0.5. Conventionally the FRQ column is intended to show the minor/effect allele frequency.
    ## The FRQ column was mapped from one of the following from the inputted  summary statistics file:
    ## FRQ, EAF, FREQUENCY, FRQ_U, F_U, MAF, FREQ, FREQ_TESTED_ALLELE, FRQ_TESTED_ALLELE, FREQ_EFFECT_ALLELE, FRQ_EFFECT_ALLELE, EFFECT_ALLELE_FREQUENCY, EFFECT_ALLELE_FREQ, EFFECT_ALLELE_FRQ, A2FREQ, A2FRQ, ALLELE_FREQUENCY, ALLELE_FREQ, ALLELE_FRQ, AF, MINOR_AF, EFFECT_AF, A2_AF, EFF_AF, ALT_AF, ALTERNATIVE_AF, INC_AF, A_2_AF, TESTED_AF, ALLELEFREQ, ALT_FREQ, EAF_HRC, EFFECTALLELEFREQ, FREQ.B, FREQ_EUROPEAN_1000GENOMES, FREQ_HAPMAP, FREQ_TESTED_ALLELE_IN_HRS, FRQ_U_113154, FRQ_U_31358, FRQ_U_344901, FRQ_U_43456, POOLED_ALT_AF, AF_ALT, AF.ALT, AF-ALT, ALT.AF, ALT-AF, A2.AF, A2-AF, AF.EFF, AF_EFF, ALL_AF

    ## As frq_is_maf=TRUE, the FRQ column will not be renamed. If the FRQ values were intended to represent major allele frequency,
    ## set frq_is_maf=FALSE to rename the column as MAJOR_ALLELE_FRQ and differentiate it from minor/effect allele frequency.

    ## Sorting coordinates with 'data.table'.

    ## Writing in tabular format ==> /tmp/Rtmpt6oeYw/file393568a0e909.tsv.gz

    ## Summary statistics report:
    ##    - 93 rows (100% of original 93 rows)
    ##    - 93 unique variants
    ##    - 70 genome-wide significant variants (P<5e-8)
    ##    - 20 chromosomes

    ## Done munging in 0.047 minutes.

    ## Successfully finished preparing sumstats file, preview:

    ## Reading header.

    ##           SNP   CHR       BP     A1     A2     FRQ   BETA    SE         P
    ##        <char> <int>    <int> <char> <char>   <num>  <num> <num>     <num>
    ## 1:   rs301800     1  8490603      T      C 0.17910  0.019 0.003 1.794e-08
    ## 2: rs11210860     1 43982527      A      G 0.36940  0.017 0.003 2.359e-10
    ## 3: rs34305371     1 72733610      A      G 0.08769  0.035 0.005 3.762e-14
    ## 4:  rs2568955     1 72762169      T      C 0.23690 -0.017 0.003 1.797e-08

    ## Returning path to saved data.

Here we know the summary statistics are based on the reference genome
GRCh37, GRCh38 can also be inputted. Moreover, if you are unsure of the
genome build, leave it as `NULL` and Mungesumstats will infer it from
the data.

Also note that the default dbSNP version used along with the reference
genome is the latest version available on Bioconductor (currently dbSNP
155) but older package versions are also available - Use the `dbSNP`
input parameter to control this.

Moreover, newer/different versions can be sued by passing a local
tarball package. Use the `dbSNP_tarball` parameter for this.
`dbSNP_tarball` was enabled to help with dbSNP versions \>=156, after
the decision to no longer provide dbSNP releases as bioconductor
packages. dbSNP 156 tarball is available
[here](http://149.165.171.124/SNPlocs/).

The arguments `format_sumstats` in that control the level of QC
conducted by *MungeSumstats* are:

- **convert_small_p** Binary, should `p-values < 5e-324` be converted to
  0? Small p-values pass the R limit and can cause errors with
  LDSC/MAGMA and should be converted. Default is TRUE.
- **convert_large_p** Binary, should p-values \>1 be converted to 1?
  P-values \>1 should not be possible and can cause errors with
  LDSC/MAGMA and should be converted. Default is TRUE.
- **convert_neg_p** Binary, should p-values \<0 be converted to 0?
  Negative p-values should not be possible and can cause errors with
  LDSC/MAGMA and should be converted. Default is TRUE.
- **compute_z** Whether to compute Z-score column from P. Default is
  FALSE. **Note** that imputing the Z-score for every SNP will not
  correct be perfectly correct and may result in a loss of power. This
  should only be done as a last resort.
- **force_new_z** When a “Z” column already exists, it will be used by
  default. To override and compute a new Z-score column from P set to
  TRUE.
- **compute_n** Whether to impute N. Default of 0 won’t impute, any
  other integer will be imputed as the N (sample size) for every SNP in
  the dataset. **Note** that imputing the sample size for every SNP is
  not correct and should only be done as a last resort. N can also be
  inputted with “ldsc”, “sum”, “giant” or “metal” by passing one of
  these for this field or a vector of multiple. Sum and an integer value
  creates an N column in the output whereas giant, metal or ldsc create
  an Neff or effective sample size. If multiples are passed, the formula
  used to derive it will be indicated.
- **convert_n_int** Binary, if N (the number of samples) is not an
  integer, should this be rounded? Default is TRUE. analysis_trait If
  multiple traits were studied, name of the trait for analysis from the
  GWAS. Default is NULL.
- **impute_beta** Binary, whether BETA should be imputed using other
  effect data if it isn’t present in the sumstats. Note that this
  imputation is an approximation so could have an effect on downstream
  analysis. Use with caution. The different methods MungeSumstats will
  try and impute beta (in this order or priority) are: 1. log(OR) 2. Z x
  SE. Default value is FALSE.
- **es_is_beta** Binary, whether to map ES to BETA. We take BETA to be
  any BETA-like value (including Effect Size). If this is not the case
  for your sumstats, change this to FALSE. Default is TRUE.
- **impute_se** Binary, whether the standard error should be imputed
  using other effect data if it isn’t present in the sumstats. Note that
  this imputation is an approximation so could have an effect on
  downstream analysis. Use with caution. The different methods
  MungeSumstats will try and impute se (in this order or priority)
  are: 1. BETA / Z
  2.  abs(BETA/ qnorm(P/2)). Default value is FALSE.
- **analysis_trait** If multiple traits were studied, name of the trait
  for analysis from the GWAS. Default is NULL.
- **INFO_filter** 0-1 The minimum value permissible of the imputation
  information score (if present in sumstatsfile). Default 0.9
- **FRQ_filter** 0-1 The minimum value permissible of the frequency(FRQ)
  of the SNP (i.e. Allele Frequency (AF)) (if present in sumstats file).
  By default no filtering is done, i.e. value of 0.  
- **pos_se** Binary Should the standard Error (SE) column be checked to
  ensure it is greater than 0? Those that are, are removed (if present
  in sumstats file). Default TRUE.
- **effect_columns_nonzero** Binary should the effect columns in the
  data BETA,OR (odds ratio),LOG_ODDS,SIGNED_SUMSTAT be checked to ensure
  no SNP=0. Those that do are removed(if present in sumstats file).
  Default TRUE.  
- **N_std** Numeric, the number of standard deviations above the mean a
  SNP’s N is needed to be removed. Default is 5. **N_dropNA** controls
  whether the SNPs with a missing N value are dropped or not (Default is
  TRUE).
- **N_dropNA** Drop rows where N is missing.Default is TRUE.
- **chr_style** Chromosome naming style to use in the formatted summary
  statistics file (“NCBI”, “UCSC”, “dbSNP”, or “Ensembl”). The NCBI and
  Ensembl styles both code chromosomes as `1-22, X, Y, MT`; the UCSC
  style is `chr1-chr22, chrX, chrY, chrM`; and the dbSNP style is
  `ch1-ch22, chX, chY, chMT`. Default is Ensembl.
- **rmv_chr** Chromosomes to exclude from the formatted summary
  statistics file. Use NULL if no filtering is necessary. Default is
  `c("X", "Y", "MT")` which removes all non-autosomal SNPs.
- **on_ref_genome** Binary, should a check take place that all SNPs are
  on the reference genome by SNP ID. Any SNPs not on the reference
  genome, will be corrected from the reference genome (if possible)
  using the chromosome and base pair position data. Default is TRUE
- **convert_ref_genome** name of the reference genome to convert to
  (“GRCh37” or “GRCh38”). This will only occur if the current genome
  build does not match. Default is not to convert the genome build
  (NULL).  
- **strand_ambig_filter** Binary, should SNPs with strand-ambiguous
  alleles be removed. Default is FALSE
- **allele_flip_check** Binary, should the allele columns be checked
  against reference genome to infer if flipping is necessary. Default is
  TRUE. **allele_flip_drop** controls whether the SNPs for which neither
  their A1 or A2 base pair values match a reference genome be dropped.
  Default is TRUE. **allele_flip_z** controls whether the Z-score value
  should be flipped along with effect and FRQ columns (e.g. Beta).
  Default is TRUE. **allele_flip_frq** controls whether the frequency
  (FRQ) value should be flipped along with effect and Z-score columns
  (e.g. Beta). Default is TRUE.
- **bi_allelic_filter** Binary, should non-biallelic SNPs be removed.
  Default is TRUE
- **flip_frq_as_biallelic** Binary, Should non-bi-allelic SNPs frequency
  values be flipped as 1-p despite there being other alternative
  alleles? Default is FALSE but if set to TRUE, this allows
  non-bi-allelic SNPs to be kept despite needing flipping.  
- **snp_ids_are_rs_ids** Binary, should the SNP IDs inputted be inferred
  as RS IDs or some arbitrary ID. Default is TRUE.  
- **remove_multi_rs_snp** Binary Sometimes summary statistics can have
  multiple RSIDs on one row (i.e. related to one SNP), for example
  “rs5772025_rs397784053”. This can cause an error so by default, the
  first RS ID will be kept and the rest removed e.g.”rs5772025”. If you
  want to just remove these SNPs entirely, set it to TRUE. Default is
  FALSE.
- **frq_is_maf** Binary, conventionally the FRQ column is intended to
  show the minor/effect allele frequency (MAF) but sometimes the major
  allele frequency can be inferred as the FRQ column. This logical
  variable indicates that the FRQ column should be renamed to
  MAJOR_ALLELE_FRQ if the frequency values appear to relate to the major
  allele i.e. \>0.5. By default mapping won’t occur i.e. is TRUE.
- **indels** Binary does your Sumstats file contain Indels? These don’t
  exist in our reference file so they will be excluded from checks if
  this value is TRUE. Further information -the reference dataset we use
  in MSS (dbSNP) does not include indels so any checks like is the SNP
  on the reference genome, attempts to impute any missing data for
  indels or check the direction of the effect columns can not be done
  for indels. Indels will be kept in the dataset if possible but certain
  situations (like if there is missing data) can cause an indel to be
  removed. See the printed information by MSS during your run to know if
  this affects you. Default is TRUE.
- **drop_indels** Binary should any indels found in the sumstats be
  dropped? These can not be checked against a reference dataset and will
  have the same RS ID and position as SNPs which can affect downstream
  analysis. Default is False.  
- **drop_na_cols** A character vector of column names to be checked for
  missing values. Rows with missing values in any of these columns (if
  present in the dataset) will be dropped. If `NULL`, all columns will
  be checked for missing values. Default columns are SNP, chromosome,
  position, allele 1, allele 2, effect columns (frequency, beta,
  Z-score, standard error, log odds, signed sumstats, odds ratio), p
  value and N columns.
- **dbSNP** The dbSNP version to use as a reference - defaults to the
  most recent version available (155). Note that with the 9x more SNPs
  in dbSNP 155 vs 144, run times will increase.
- **dbSNP_tarball** Pass local versions of dbSNP in tarball format.
  Default of NULL uses the dbSNP version passed in `dbSNP` parmeter.
  `dbSNP_tarball` was enabled to help with dbSNP versions \>=156, after
  the decision to no longer provide dbSNP releases as bioconductor
  packages. dbSNP 156 tarball is available here:
  <http://149.165.171.124/SNPlocs/>.
- **sort_coordinates** Whether to sort by coordinates of resulting
  sumstats.  
- **nThread** Number of threads to use for parallel processes.
- **write_vcf** Whether to write as VCF (TRUE) or tabular file (FALSE).
  While **tabix_index** is a binary input for whether to index the
  formatted summary statistics with
  [tabix](http://www.htslib.org/doc/tabix.md) for fast querying.
- **return_data** Return `data.table`, `GRanges` or `VRanges`directly to
  user. Otherwise, return the path to the save data. Default is FALSE.
- **return_format** If return_data is TRUE. Object type to be returned
  (“data.table”,“vranges”,“granges”).
- **save_format** Ensure that output format meets all requirements to be
  passed directly into LDSC (“ldsc”) without the need for additional
  munging or for IEU OpenGWAS format (“opengwas”) before saving as a
  VCF. **NOTE** - If LDSC format is used, the naming convention of A1 as
  the reference (genome build) allele and A2 as the effect allele will
  be reversed to match LDSC (A1 will now be the effect allele). See more
  info on this
  [here](https://groups.google.com/g/ldsc_users/c/S7FZK743w68). Note
  that any effect columns (e.g. Z) will be inrelation to A1 now instead
  of A2.
- **log_folder_ind** Should log files be stored containing all filtered
  out SNPs (separate file per filter). The data is outputted in the same
  format specified for the resulting sumstats file.
- **log_mungesumstats_msgs** Binary Should a log be stored containing
  all messages and errors printed by MungeSumstats in a run.
- **imputation_ind** Binary Should a column be added for each imputation
  step to show what SNPs have imputed values for differing fields. This
  includes a field denoting SNP allele flipping (flipped). On the
  flipped value, this denoted whether the alelles where switched based
  on MungeSumstats initial choice of A1, A2 from the input column
  headers and thus may not align with what the creator intended.**Note**
  these columns will be in the formatted summary statistics returned.
- **log_folder** File path to the directory for the log files and the
  log of MungeSumstats messages to be stored. Default is a temporary
  directory.
- **force_new** If a formatted file of the same names as exists,
  formatting will be skipped and this file will be imported instead
  (default). Set to override this.
- **mapping_file** MungeSumstats has a pre-defined column-name mapping
  file which should cover the most common column headers and their
  interpretations. However, if a column header that is in youf file is
  missing of the mapping we give is incorrect you can supply your own
  mapping file. Must be a 2 column dataframe with column names
  “Uncorrected” and “Corrected”. See `data(sumstatsColHeaders)` for
  default mapping and necessary format.

See `?MungeSumstats::format_sumstats()` for the full list of parameters
to control MungeSumstats QC and standardisation steps.

VCF files can also be standardised to the same format as other summary
statistic files. A subset of the Amyotrophic lateral sclerosis GWAS from
the ieu open GWAS project (a .vcf file) has been added to
*MungeSumstats* to demonstrate this functionality.Simply pass the path
to the file in the same manner you would for other summary statistic
files:

``` r

#save ALS GWAS from the ieu open GWAS project to a temp directory
ALSvcfPth <- system.file("extdata","ALSvcf.vcf", package="MungeSumstats")
```

``` r

reformatted_vcf <- 
  MungeSumstats::format_sumstats(path=ALSvcfPth, 
                                 ref_genome="GRCh37")
```

You can also get more information on the SNPs which have had data
imputed or have been filtered out by *MungeSumstats* by using the
`imputation_ind` and `log_folder_ind` parameters respectively. For
example:

``` r

#set
reformatted_vcf_2 <- 
  MungeSumstats::format_sumstats(path=ALSvcfPth,
                                 ref_genome="GRCh37",
                                 log_folder_ind=TRUE,
                                 imputation_ind=TRUE,
                                 log_mungesumstats_msgs=TRUE)
```

    ## Time difference of 0 secs

    ## Time difference of 0.2 secs

Check the file `snp_bi_allelic.tsv.gz` in the `log_folder` directory you
supply (by default a temp directory), for a list of SNPs removed as they
are non-bi-allelic. The text files containing the console output and
messages are also stored in the same directory.

Note you can also control the dbSNP version used as a reference dataset
by MungeSumstats using the `dbSNP` parameter. By default this will be
set to the most recent dbSNP version available (155).

Note that using `log_folder_ind` returns a list from `format_sumstats`
which includes the file locations of the differing classes of removed
SNPs. Using `log_mungesumstats_msgs` saves the messages to the console
to a file which is returned in the same list. Note that not all the
messages will also print to screen anymore when you set
`log_mungesumstats_msgs`:

``` r

names(reformatted_vcf_2)
```

    ## [1] "sumstats"  "log_files"

A user can load a file to view the excluded SNPs.

In this case, SNPs were filtered based on non-bi-allelic criterion:

``` r

print(reformatted_vcf_2$log_files$snp_bi_allelic)
```

    ## NULL

The different types of exclusion which lead to the names are explained
below:

- **snp_multi_rs_one_row** - Where the SNP (RS ID) contained more than
  one RS ID.
- **snp_missing_rs** - Where the SNP (RS ID) was missing the rs prefix.
  Note that these are only removed when other snps have an rs prefix.
- **snp_multi_colon** - Where the SNP ID has mutliple colons (“:”) in
  one SNP.
- **snp_not_found_from_bp_chr** - Where the RS ID was attempted to be
  imputed from the CHR and BP (Base-Pair) information, using the
  reference genome, but wasn’t successful.
- **chr_bp_not_found_from_snp** - Where the CHR and BP (Base-Pair) was
  attempted to be imputed from the SNP (RS ID), using the reference
  genome, but wasn’t successful.
- **alleles_not_found_from_snp** - Where the alleles (A1 and/or A2) was
  attempted to be imputed from the SNP (RS ID), using the reference
  genome, but wasn’t successful.
- **alleles_dont_match_ref_gen** - Where the alleles (A1 and/or A2)
  don’t match what’s on the reference genome.
- **missing_data** - Where there is data missing across the inputted
  columns.
- **dup_snp_id** - Where the SNP ID is duplicated in the input.
- **dup_base_pair_position** - Where the base-pair position is
  duplicated in the input.
- **info_filter** - SNP INFO value below the specified threshold.
- **se_neg** - SNPs SE (Standard Error) value is 0 or negative.
- **effect_col_zero** - SNPs effect column(s) value is zero e.g. BETA=0.
- **n_large** - SNPs N is N standard deviations greater than the mean.
- **n_null** - SNPs N is null.
- **chr_excl** - SNP has an unrecognized chromosome name or is on a
  chromosome that was specified to be excluded.
- **snp_strand_ambiguous** - SNP is strand ambiguous.
- **snp_bi_allelic** - SNP is not bi-allelic.
- **MungeSumstats_log_msg** - Text file of all messages to the console
  created during MungeSumstats run.
- **MungeSumstats_log_output** - Text file of all errors to the console
  created during MungeSumstats run.

Note to export to another type such as R native objects including
data.table, GRanges, VRanges or save as a VCF file, set
`return_data=TRUE` and choose your `return_format`:

``` r

#set
reformatted_vcf_2 <- 
  MungeSumstats::format_sumstats(path=ALSvcfPth,
                                 ref_genome="GRCh37", 
                                 log_folder_ind=TRUE,
                                 imputation_ind=TRUE,
                                 log_mungesumstats_msgs=TRUE,
                                 return_data=TRUE,
                                 return_format="GRanges")
```

Also you can now output a VCF compatible with [IEU
OpenGWAS](https://gwas.mrcieu.ac.uk/) format (Note that currently all
IEU OpenGWAS sumstats are GRCh37, MungeSumstats will throw a warning if
your data isn’t GRCh37 when saving):

``` r

#set
reformatted_vcf_2 <- 
  MungeSumstats::format_sumstats(path=ALSvcfPth,
                                 ref_genome="GRCh37", 
                                 write_vcf=TRUE,
                                 save_format ="openGWAS")
```

See our publication for further discussion of these checks and options:

[Murphy et al. MungeSumstats: A Bioconductor package for the
standardisation and quality control of many GWAS summary
statistics](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btab665/6380562).

## Extra Functionality

### Get genome builds

*MungeSumstats* also contains a function to quickly infer the genome
build of multiple summary statistic files. This can be called separately
to
[`format_sumstats()`](https://al-murphy.github.io/MungeSumstats/reference/format_sumstats.md)
which is useful if you want to just quickly check the genome build:

``` r

# Pass path to Educational Attainment Okbay sumstat file to a temp directory
eduAttainOkbayPth <- system.file("extdata", "eduAttainOkbay.txt",
                                  package = "MungeSumstats")
ALSvcfPth <- system.file("extdata","ALSvcf.vcf", package="MungeSumstats")
sumstats_list <- list(ss1 = eduAttainOkbayPth, ss2 = ALSvcfPth)

ref_genomes <- MungeSumstats::get_genome_builds(sumstats_list = sumstats_list)
```

### Liftover

*MungeSumstats* exposes the
[`liftover()`](https://al-murphy.github.io/MungeSumstats/reference/liftover.md)
function as a general utility for users.

Useful features include: - Automatic standardisation of genome build
names (i.e. “hg19”, “hg37”, and “GRCh37” will all be recognized as the
same genome build.) - Ability to specify `chrom_col` as well as both
`start_col` and `end_col` (for variants that span \>1bp). - Ability to
return in `data.table` or `GRanges` format. - Ability to specify which
chromosome format (e.g. “chr1” vs. 1) to return `GRanges` as.

``` r

sumstats_dt <- MungeSumstats::formatted_example()
```

    ## Standardising column headers.

    ## First line of summary statistics file:

    ## MarkerName   CHR POS A1  A2  EAF Beta    SE  Pval    

    ## Sorting coordinates with 'data.table'.

``` r

sumstats_dt_hg38 <- MungeSumstats::liftover(sumstats_dt = sumstats_dt, 
                                            ref_genome = "hg19",
                                            convert_ref_genome = "hg38")
```

    ## Performing data liftover from hg19 to hg38.

    ## First - Checking CHR format to avoid liftover issues with chain file.

    ## Formatting CHR column (if necessary)- Removing chr prefix, converting 23/24 to X/Y and M to MT

    ## Converting summary statistics to GenomicRanges.

    ## Downloading chain file...

    ## Downloading chain file from Ensembl.

    ## /tmp/Rtmpt6oeYw/GRCh37_to_GRCh38.chain.gz

    ## Reordering so first three column headers are SNP, CHR and BP in this order.

    ## Reordering so the fourth and fifth columns are A1 and A2.

``` r

knitr::kable(head(sumstats_dt_hg38))
```

| SNP        | CHR |       BP | A1  | A2  |     FRQ |   BETA |    SE |     P | IMPUTATION_gen_build |
|:-----------|:----|---------:|:----|:----|--------:|-------:|------:|------:|:---------------------|
| rs301800   | 1   |  8430543 | T   | C   | 0.17910 |  0.019 | 0.003 | 0e+00 | TRUE                 |
| rs11210860 | 1   | 43516856 | A   | G   | 0.36940 |  0.017 | 0.003 | 0e+00 | TRUE                 |
| rs34305371 | 1   | 72267927 | A   | G   | 0.08769 |  0.035 | 0.005 | 0e+00 | TRUE                 |
| rs2568955  | 1   | 72296486 | T   | C   | 0.23690 | -0.017 | 0.003 | 0e+00 | TRUE                 |
| rs1008078  | 1   | 90724174 | T   | C   | 0.37310 | -0.016 | 0.003 | 0e+00 | TRUE                 |
| rs61787263 | 1   | 98153158 | T   | C   | 0.76120 |  0.016 | 0.003 | 1e-07 | TRUE                 |

### Quick formatting

In some cases, users may not want to run the full munging pipeline
provided by  
[`MungeSumstats::format_sumstats`](https://al-murphy.github.io/MungeSumstats/reference/format_sumstats.md),
but still would like to take advantage of the file type conversion and
column header standardisation features. This will not be nearly as
robust as the full pipeline, but can still be helpful.

#### From disk

To do this, simply run the following:

``` r

eduAttainOkbayPth <- system.file("extdata", "eduAttainOkbay.txt",
                                  package = "MungeSumstats")
formatted_path <- tempfile(fileext = "_eduAttainOkbay_standardised.tsv.gz")


#### 1. Read in the data and standardise header names ####
dat <- MungeSumstats::read_sumstats(path = eduAttainOkbayPth, 
                                    standardise_headers = TRUE)
```

    ## Importing tabular file: /home/runner/work/_temp/Library/MungeSumstats/extdata/eduAttainOkbay.txt

    ## Checking for empty columns.

    ## Standardising column headers.

    ## First line of summary statistics file:

    ## MarkerName   CHR POS A1  A2  EAF Beta    SE  Pval    

``` r

knitr::kable(head(dat))
```

| SNP        | CHR |        BP | A1  | A2  |    FRQ |   BETA |    SE |     P |
|:-----------|:----|----------:|:----|:----|-------:|-------:|------:|------:|
| rs10061788 | 5   |  87934707 | A   | G   | 0.2164 |  0.021 | 0.004 | 0e+00 |
| rs1007883  | 16  |  51163406 | T   | C   | 0.3713 | -0.015 | 0.003 | 1e-07 |
| rs1008078  | 1   |  91189731 | T   | C   | 0.3731 | -0.016 | 0.003 | 0e+00 |
| rs1043209  | 14  |  23373986 | A   | G   | 0.6026 |  0.018 | 0.003 | 0e+00 |
| rs10496091 | 2   |  61482261 | A   | G   | 0.2705 | -0.018 | 0.003 | 0e+00 |
| rs10930008 | 2   | 161854736 | A   | G   | 0.7183 | -0.016 | 0.003 | 1e-07 |

``` r

#### 2. Write to disk as a compressed, tab-delimited, tabix-indexed file ####
formatted_path <- MungeSumstats::write_sumstats(sumstats_dt = dat,
                                                save_path = formatted_path,
                                                tabix_index = TRUE,
                                                write_vcf = FALSE,
                                                return_path = TRUE)   
```

    ## Sorting coordinates with 'data.table'.

    ## Writing in tabular format ==> /tmp/Rtmpt6oeYw/file393566488bf_eduAttainOkbay_standardised.tsv

    ## Writing uncompressed instead of gzipped to enable tabix indexing.

    ## Converting full summary stats file to tabix format for fast querying...

    ## Reading header.

    ## Ensuring file is bgzipped.

    ## Tabix-indexing file.

    ## Removing temporary .tsv file.

#### From `data.table`

If you already have your data imported as an `data.table`, you can also
standardise its headers like so:

``` r

#### Mess up some column names ####
dat_raw <- data.table::copy(dat)
data.table::setnames(dat_raw, c("SNP","CHR"), c("rsID","Seqnames"))
#### Add a non-standard column that I want to keep the casing for ####
dat_raw$Support <- runif(nrow(dat_raw))

dat2 <- MungeSumstats::standardise_header(sumstats_dt = dat_raw,
                                          uppercase_unmapped = FALSE, 
                                          return_list = FALSE )
```

    ## Standardising column headers.

    ## First line of summary statistics file:

    ## rsID Seqnames    BP  A1  A2  FRQ BETA    SE  P   Support 

    ## Returning unmapped column names without making them uppercase.

``` r

knitr::kable(head(dat2))
```

| SNP        | CHR |       BP | A1  | A2  |     FRQ |   BETA |    SE |     P |   Support |
|:-----------|:----|---------:|:----|:----|--------:|-------:|------:|------:|----------:|
| rs301800   | 1   |  8490603 | T   | C   | 0.17910 |  0.019 | 0.003 | 0e+00 | 0.9805397 |
| rs11210860 | 1   | 43982527 | A   | G   | 0.36940 |  0.017 | 0.003 | 0e+00 | 0.7415215 |
| rs34305371 | 1   | 72733610 | A   | G   | 0.08769 |  0.035 | 0.005 | 0e+00 | 0.0514463 |
| rs2568955  | 1   | 72762169 | T   | C   | 0.23690 | -0.017 | 0.003 | 0e+00 | 0.5302125 |
| rs1008078  | 1   | 91189731 | T   | C   | 0.37310 | -0.016 | 0.003 | 0e+00 | 0.6958239 |
| rs61787263 | 1   | 98618714 | T   | C   | 0.76120 |  0.016 | 0.003 | 1e-07 | 0.6885560 |

## Future Enhancements

The *MungeSumstats* package aims to be able to handle the most common
summary statistic file formats including VCF. If your file can not be
formatted by *MungeSumstats* feel free to report the bug on github:
<https://github.com/neurogenomics/MungeSumstats> along with your summary
statistic file header.

We also encourage people to edit the code to resolve their particular
issues too and are happy to incorporate these through pull requests on
github. If your summary statistic file headers are not recognised by
*MungeSumstats* but correspond to one of:

    SNP, BP, CHR, A1, A2, P, Z, OR, BETA, LOG_ODDS,
    SIGNED_SUMSTAT, N, N_CAS, N_CON, NSTUDY, INFO or FRQ 

feel free to update the
[`MungeSumstats::sumstatsColHeaders`](https://al-murphy.github.io/MungeSumstats/reference/sumstatsColHeaders.md)
following the approach in the data.R file and add your mapping. Then use
a pull request on github and we will incorporate this change into the
package.

A note on
[`MungeSumstats::sumstatsColHeaders`](https://al-murphy.github.io/MungeSumstats/reference/sumstatsColHeaders.md)
for summary statistic files with A0/A1. The mapping in
[`MungeSumstats::sumstatsColHeaders`](https://al-murphy.github.io/MungeSumstats/reference/sumstatsColHeaders.md)
converts A0 to A\*, this is a special case so that the code knows to map
A0/A1 to A1/A2 (ref/alt). The special case is needed since ordinarily A1
refers to the reference not the alternative allele.

A note on
[`MungeSumstats::sumstatsColHeaders`](https://al-murphy.github.io/MungeSumstats/reference/sumstatsColHeaders.md)
for summary statistic files with Effect Size (ES). By default, MSS takes
BETA to be any BETA-like value (including ES). This is coded into the
mapping file -
[`MungeSumstats::sumstatsColHeaders`](https://al-murphy.github.io/MungeSumstats/reference/sumstatsColHeaders.md).
If this isn’t the case for your sumstats, you can set the `es_is_beta`
parameter in
[`MungeSumstats::format_sumstats()`](https://al-murphy.github.io/MungeSumstats/reference/format_sumstats.md)
to FALSE to avoid this. Note this is done to try and capture most use
cases of MSS.

## Further functionality

See the [Open GWAS
vignette](https://neurogenomics.github.io/MungeSumstats/articles/OpenGWAS.html)
for how MungeSumstats can be used along with data from the MRC IEU Open
GWAS Project and also Mungesumstats’ functionality to handle lists of
summary statistics files.

## Session Information

    ## R version 4.6.1 (2026-06-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.5 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] MungeSumstats_1.21.2 BiocStyle_2.40.0    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.1                           
    ##  [2] dplyr_1.2.1                                
    ##  [3] blob_1.3.0                                 
    ##  [4] R.utils_2.13.0                             
    ##  [5] Biostrings_2.80.2                          
    ##  [6] bitops_1.1-0                               
    ##  [7] fastmap_1.2.0                              
    ##  [8] RCurl_1.98-1.20                            
    ##  [9] VariantAnnotation_1.58.0                   
    ## [10] GenomicAlignments_1.48.0                   
    ## [11] XML_3.99-0.25                              
    ## [12] digest_0.6.39                              
    ## [13] lifecycle_1.0.5                            
    ## [14] KEGGREST_1.52.2                            
    ## [15] RSQLite_3.53.3                             
    ## [16] magrittr_2.0.5                             
    ## [17] compiler_4.6.1                             
    ## [18] rlang_1.3.0                                
    ## [19] sass_0.4.10                                
    ## [20] tools_4.6.1                                
    ## [21] yaml_2.3.12                                
    ## [22] data.table_1.18.6.1                        
    ## [23] rtracklayer_1.72.0                         
    ## [24] knitr_1.52                                 
    ## [25] S4Arrays_1.12.1                            
    ## [26] bit_4.6.0                                  
    ## [27] curl_8.1.0                                 
    ## [28] DelayedArray_0.38.2                        
    ## [29] ieugwasr_1.2.0                             
    ## [30] abind_1.4-8                                
    ## [31] BiocParallel_1.46.0                        
    ## [32] BiocGenerics_0.58.1                        
    ## [33] desc_1.4.3                                 
    ## [34] R.oo_1.27.1                                
    ## [35] grid_4.6.1                                 
    ## [36] stats4_4.6.1                               
    ## [37] SummarizedExperiment_1.42.0                
    ## [38] cli_3.6.6                                  
    ## [39] rmarkdown_2.32                             
    ## [40] crayon_1.5.3                               
    ## [41] ragg_1.5.2                                 
    ## [42] generics_0.1.4                             
    ## [43] otel_0.2.0                                 
    ## [44] BSgenome.Hsapiens.1000genomes.hs37d5_0.99.1
    ## [45] httr_1.4.9                                 
    ## [46] rjson_0.2.23                               
    ## [47] DBI_1.3.0                                  
    ## [48] cachem_1.1.0                               
    ## [49] stringr_1.6.0                              
    ## [50] parallel_4.6.1                             
    ## [51] AnnotationDbi_1.74.0                       
    ## [52] BiocManager_1.30.27                        
    ## [53] XVector_0.52.0                             
    ## [54] restfulr_0.0.17                            
    ## [55] matrixStats_1.5.0                          
    ## [56] vctrs_0.7.3                                
    ## [57] Matrix_1.7-5                               
    ## [58] jsonlite_2.0.0                             
    ## [59] bookdown_0.48                              
    ## [60] IRanges_2.46.0                             
    ## [61] S4Vectors_0.50.3                           
    ## [62] bit64_4.8.6                                
    ## [63] GenomicFiles_1.48.0                        
    ## [64] systemfonts_1.3.2                          
    ## [65] GenomicFeatures_1.64.0                     
    ## [66] jquerylib_0.1.4                            
    ## [67] glue_1.8.1                                 
    ## [68] pkgdown_2.2.1                              
    ## [69] codetools_0.2-20                           
    ## [70] stringi_1.8.9                              
    ## [71] GenomeInfoDb_1.48.0                        
    ## [72] BiocIO_1.22.0                              
    ## [73] GenomicRanges_1.64.0                       
    ## [74] UCSC.utils_1.8.0                           
    ## [75] tibble_3.3.1                               
    ## [76] pillar_1.11.1                              
    ## [77] SNPlocs.Hsapiens.dbSNP155.GRCh37_0.99.24   
    ## [78] htmltools_0.5.9                            
    ## [79] Seqinfo_1.2.0                              
    ## [80] BSgenome_1.80.0                            
    ## [81] R6_2.6.1                                   
    ## [82] textshaping_1.0.5                          
    ## [83] evaluate_1.0.5                             
    ## [84] lattice_0.22-9                             
    ## [85] Biobase_2.72.0                             
    ## [86] R.methodsS3_1.8.2                          
    ## [87] png_0.1-9                                  
    ## [88] Rsamtools_2.28.0                           
    ## [89] cigarillo_1.2.1                            
    ## [90] memoise_2.0.1                              
    ## [91] bslib_0.12.0                               
    ## [92] SparseArray_1.12.3                         
    ## [93] xfun_0.62                                  
    ## [94] fs_2.1.0                                   
    ## [95] MatrixGenerics_1.24.0                      
    ## [96] pkgconfig_2.0.3

## References

1\.

Nathan G. Skene, T. E. B., Julien Bryois. Genetic identification of
brain cell types underlying schizophrenia. *Nature Genetics* (2018).
doi:[10.1038/s41588-018-0129-5](https://doi.org/10.1038/s41588-018-0129-5)
