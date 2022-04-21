## CHANGES IN VERSION 1.3.18

### New features
* Can now handle general remote sumstats not just IEU GWAS
* More column header mappings

## CHANGES IN VERSION 1.3.17

### New features
* Clean up of column header mapping file, including FREQUENCY given priority 
over MAF and addition of new CHR mappings.

## CHANGES IN VERSION 1.3.15

### Bug fixes

* Handle cases for multi-trait GWAS when P columns exists separate to the trait
specific P value so that when renaming occurs there isn't two P columns. 
Inputted P column will be renamed to 'P_input'
* Issue where 'check allele flip' wasn't running when the sumstats had all SNP 
IDs missing and incorrect direction of A1/A2 and effect columns has now been 
fixed.

## CHANGES IN VERSION 1.3.14

### New features

* `liftover` 
    - Now exported function.
    - Added args for more user flexibility.
    - Uses `GenomeInfoDb::mapGenomeBuilds` to standardise build names.
    - Warns users when mapped builds do not match one of the conversion options.
    - Choice to output as data.table or GRanges.
    - Added units tests for exported version.
* `standardise_sumstats_column_headers_crossplatform`
    - Exported as `standardise_header` while keeping the original function 
    name as an internal function (they call the same code).
    - Added unit tests for exported version.
* Added chunks to *Getting started` vignette
    - `liftover` tutorial
    - "Quick formatting" of headers and file formats.

### Bug fixes 

* `check_pos_se`: Remove extra `message()` call around string.
* `check_signed_col`: Remove extra `message()` call around string.
* `write_sumstats`
    - Added extra round of sorting when `tabix_index=TRUE` because this is 
    required for tabix.

## CHANGES IN VERSION 1.3.13

### New Features 

* Additional mappings for CHR
* Make A1, A2 upper-case

### Bug fixes 

* Bug fix for dealing with imputing SNP ID when there are indels

## CHANGES IN VERSION 1.3.11

### New Features 

* MungeSumstats can now handle Indels better. It will:
    - Not impute the RS ID of a SNP for an Indel
    - Not remove the Indel based on the RS ID not being present in the SNP ref
    dataset.
    - Not remove the Indel if it has the same base-pair location as a SNP in the
    sumstats.
* Can now handle vcfs with extensions .vcf.tsv, .vcf.tsv.gz and .vcf.tsv.bgz     

### Bug fixes 

* For non-bi-allelic SNP runs, no longer remove duplicated SNPs based on their 
base-pair position or their RS ID.

## CHANGES IN VERSION 1.3.9

### New Features 

* Exported functions. Added examples and unit tests:
    - `compute_nsize`
    - `standardise_sumstats_column_headers_crossplatform`
    - `formatted_example`
* New arguments:
    - `standardise_sumstats_column_headers_crossplatform`: 
    Added arg `uppercase_unmapped` to 
    to allow users to specify whether they want make the columns that could not be 
    mapped to a standard name uppercase (default=`TRUE` for backcompatibility).
    Added arg `return_list` to specify whether to return a named list 
    (default) or just the `data.table`.
    - `formatted_example`: 
    Added args `formatted` to specify whether the file should have its colnames standardised. 
    Added args `sorted` to specify whether the file should sort the data by coordinates. 
    Added arg `return_list` to specify whether to return a named list 
    (default) or just the `data.table`.
* Removed *codecode.yml* and *_pkgdown.yml* files (no longer necessary).
* Added Issues templates for Bugs and Feature requests. 
* Added `.datatable.aware=TRUE` to *.zzz* as [extra precaution](https://rdatatable.gitlab.io/data.table/articles/datatable-importing.html#optionally-import-data-table-suggests).  
* `vcf2df`: Documented arguments. 
* Made v2 of hex sticker: *inst/hex/hex.png*   

### Bug fixes 

* Regenerated the *gh-pages* branch after it accidentally got deleted. 
* Remove temporary *docs/* folder.  
* Updated GitHub Actions. 
* Updated *Dockerfile* so it doesn't run checks 
(this is now take care of by the GHA workflow). 
* Added Windows-specific folders to *.Rbuildignore*. 
* Made *to_GRanges.R* and *to_VRanges.R* file names lowercase
to be congruent with function names.


## CHANGES IN VERSION 1.3.7

### Bug fixes

* Bug in checking for bad characters in RSID fixed


## CHANGES IN VERSION 1.3.6

### New Features

* Columns Beta and Standard Error can now be imputed. However note that this 
imputation is an approximation so could have an effect on downstream 
analysis. Use with caution.


## CHANGES IN VERSION 1.3.5

### Bug fixes

* Flipping of Odds Ratio corrected (1/OR rather than -1*OR)

## CHANGES IN VERSION 1.3.4

### Bug fixes

* Issue downloading chain file resolved

## CHANGES IN VERSION 1.3.3

### New Features

* More mappings added to default mapping file.


## CHANGES IN VERSION 1.3.2

### Bug fixes

* Previously rsids with characters added (e.g. rs1234567w) would cause an error
when checking for the rsid on the reference genome. This has been fixed and 
the correct rsid will now be imputed from the reference genome for these cases.

## CHANGES IN VERSION 1.3.1

### New Features

* `import_sumstats`: Create individual folders for each GWAS dataset,
with a respective `logs` subfolder to avoid overwriting log files
when processing multiple GWAS.  
* `parse_logs`: **New function** to convert logs from one or more munged GWAS
into a `data.table`.  
* `list_sumstats`: **New function** to recursively search for local 
summary stats files previously munged with `MungeSumstats`.  
* Added new dataset `inst/extdata/MungeSumstats_log_msg.txt` 
to test logs files.  
* Added unit tests for `list_sumstats` and `parse_logs`. 
* Added new Docker vignette.  
* Updated GHA workflows using [r_workflows](https://github.com/neurogenomics/r_workflows).  
* Remove *docs/* folder as the website will now be pushed to 
the `gh-pages` branch automatically by new GHA workflow.  
* Made documentation in README more clear and concise.  
* Added checks for p-values >1 or <0 via args `convert_large_p` and 
`convert_neg_p`, respectively.
These are both handled by the new internal function `check_range_p_val`, 
which also reports the number of SNPs found meeting these criteria 
to the console/logs.  
* `check_small_p_val` records which SNPs were imputed in a  more robust way, 
by recording which SNPs met the criteria before making the changes (as opposed to inferred this info from which columns are 0 after making the changes). This 
function now only handles non-negative p-values, so that rows with negative
p-values can be recorded/reported separately in the `check_range_p_val` step.  
* `check_small_p_val` now reports the number of SNPs <= 5e-324 to console/logs. 
* Unit tests have been added for both `check_range_p_val` 
and `check_small_p_val`. 
* `parse_logs` can now extract information reported by `check_range_p_val` and 
`check_small_p_val`.  
* New internal function `logs_example` provides easy access to log file stored 
in *inst/extdata*, and includes documentation on how it was created.  
* Both `check_range_p_val` and `check_small_p_val` now use `#' @inheritParams format_sumstats` to improve consistency of documentation.  

### Bug fixes

* Reduced vignette sizes.
* Removed usage of `suppressWarnings` where possible.  
* Deleted old *.Rproj* file and hidden folder (contained large files). 
* Configured *.Rproj* so it doesn't store large data files. 
* Fix badger issues: https://github.com/GuangchuangYu/badger/issues/34 
* Prevent *test-index_tabix.R* from running due to errors (for now). 


## CHANGES IN VERSION 1.3.0

### New Features

* Version bump to align with Bioconductor release 3.14.

## CHANGES IN VERSION 1.1.27

### Bug fixes

* `validate_parameters` can now handle `ref_genome=NULL`  
* *.tsv.gz* no longer assigned suffix *.tsv*.   
* Made code width <80 characters.  
* Changed `to_GRanges`/`to_GRanges` functions to all-lowercase functions
(for consistency with other functions). 
* Set `nThread=1` in `data.table` test functions.

### New Features

* Added tests for `get_genome_builds`  
* Added early check for making sure the directory `save_path` is in was 
actually created (as opposed to finding out at the very end of the pipeline). 
* Tabix-indexing now available for tabular output data.
* `read_header` and `read_sumstats` now both work with .bgz files.  

## CHANGES IN VERSION 1.1.26

### New Features

* Extra mappings for FRQ column, see `data("sumstatsColHeaders")` for details  

## CHANGES IN VERSION 1.1.23

### New Features

* `format_sumstats(FRQ_filter)` added so SNPs can now be filtered by allele 
frequency 
* Mapping file now has mappings for allele frequency (AF) to FRQ
* VCF files with AF in INFO column e.g. 'AF=...' now converted to AF column
* `format_sumstats(frq_is_maf)` check added to infer if FRQ column values are
minor/effect allele frequencies or not. frq_is_maf allows users to rename the
FRQ column as MAJOR_ALLELE_FRQ if some values appear to be major allele 
frequencies

## CHANGES IN VERSION 1.1.19

### New Features

* `get_genome_builds()` can now be called to quickly get the genome build 
without running the whole reformatting.
* `format_sumstats(compute_n)` now has more methods to compute the effective 
sample size with "ldsc", "sum", "giant" or "metal". 
* `format_sumstats(convert_ref_genome)` now implemented which can perform 
liftover to GRCh38 from GRCh37 and vice-versa enabling better cohesion between
different study's summary statistics.

## CHANGES IN VERSION 1.1.11

### Bug fixes

* `check_no_rs_snp` can now handle extra information after an RS ID. So if you 
have `rs1234:A:G` that will be separated into two columns.
* `check_two_step_col` and `check_four_step_col`, the two checks for when 
multiple columns are in one, have been updated so if not all SNPs have multiple
columns or some have more than the expected number, this can now be handled.
* Extra mappings for the `FRQ` column have been added to the mapping file

## CHANGES IN VERSION 1.1.11

### New Features

* `check_multi_rs_snp` can now handle all punctuation with/without spaces. So if
a row contains `rs1234,rs5678` or `rs1234, rs5678` or any other punctuation 
character other than `,` these can be handled.
* `format_sumstats(path)` can now be passed a dataframe/datatable of the summary
statistics directly as well as a path to their saved location.
* Input summary statistics with `A0/A1` corresponding to ref/alt can now be 
handled by the mappign file as well as `A1/A2` corresponding to ref/alt.

## CHANGES IN VERSION 1.1.2

### New Features

*   `import_sumstats` reads GWAS sum stats directly from Open GWAS. Now 
parallelised and reports how long each dataset took to import/format in total. 
*   `find_sumstats` searches Open GWAS for datasets. 
*   `compute_z` computes Z-score from P. 
*   `compute_n` computes N for all SNPs from user defined smaple size.
*   `format_sumstats(ldsc_format=TRUE)` ensures sum stats can be fed directly 
into [LDSC](https://github.com/bulik/ldsc) without any additional munging. 
*   `read_sumstats`, `write_sumstas`, and `download_vcf` functions now exported.  
*   `format_sumstats(sort_coordinates=TRUE)` sorts results by their genomic 
coordinates. 
*   `format_sumstats(return_data=TRUE)` returns data directly to user. Can be 
returned in either `data.table` (default), `GRanges` or `VRanges` format using 
`format_sumstats(return_format="granges")`.  
*   `format_sumstats(N_dropNA=TRUE)` (default) drops rows where N is missing. 
*   `format_sumstats(snp_ids_are_rs_ids=TRUE)` (default) Should the SNP IDs 
inputted be inferred as RS IDs or some arbitrary ID.
*   `format_sumstats(write_vcf=TRUE)` writes a tabix-indexed VCF file instead of
tabular format. 
*   `format_sumstats(save_path=...)` lets users decide where their results are 
saved and what they're named. 
*   When the `save_path` indicates it's in `tempdir()`, message warns users that
these files will be deleted when R session ends.  
*   Summary of data is given at the beginning and the end of `format_sumstats` 
via `report_summary()`.  
*   Readability of `preview_sumstats()` messages improved.  
*   New checks standard error (SE) must >0 and BETA (and other effect columns) 
must not equal 0: `format_sumstats(pos_se=TRUE,effect_columns_nonzero=TRUE)`
*   Log directory containing all removed SNPs is now available and can be 
changed to a different directory by setting:
`format_sumstats(log_folder_ind=TRUE,log_folder=tempdir())`
*   All imputed data can now be identified with a column in the output using:
`format_sumstats(imputation_ind=TRUE)`
*   Users can now input their own mapping file to be used for the column header 
mapping in place of `data(sumstatsColHeaders)`. See 
`format_sumstats(mapping_file = mapping_file)`.


### Bug fixes 

*   CHR column now standardised (X and Y caps, no "chr" prefix).
*   Allele flipping done on a per-SNP basis (instead of whole-column). 
*   Allele flipping now includes FRQ column as well as effect columns.
*   The effect allele is now interpreted as the A2 allele consistent with [IEU GWAS VCF approach](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7805039/). A1 will always be the reference allele.
*   `read_vcf` upgraded to account for more VCF formats. 
*   `check_n_num` now accounts for situations where N is a character vector and converts to numeric.  


## CHANGES IN VERSION 1.1.1

### Bug fixes

*   Preprint publication citation added.


## CHANGES IN VERSION 1.0.0

### New Features

*   MungeSumstats released to Bioconductor.

