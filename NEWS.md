## CHANGES IN VERSION 1.17.1

### New features
* MungeSumstats can now handle local versions of dbSNP in tarball format 
(using the `dbSNP_tarball` parameter in `format_sumstats`). This was enabled to
help with dbSNP versions >=156, after the decision to no longer provide dbSNP 
releases as bioconductor packages. dbSNP 156 tarball is available 
here:http://149.165.171.124/SNPlocs/.

## CHANGES IN VERSION 1.15.3

### Bug fix
*Updated retrieval of IEU OpenGWAS to new approach requiring login. Also updated
to use the [IEU OpwnGWAS R package](https://github.com/MRCIEU/ieugwasr) as a 
dependency.

## CHANGES IN VERSION 1.15.2

### New features
*[FAQ Website](https://github.com/Al-Murphy/MungeSumstats/wiki/FAQ) updated.

## CHANGES IN VERSION 1.15.1

### New features
*[FAQ Website](https://github.com/Al-Murphy/MungeSumstats/wiki/FAQ) added.

## CHANGES IN VERSION 1.13.7

### Bug fix
* `infer_eff_direction` now includes A0 as an ambiguous case as well as A1/A2.

### New features
* `eff_on_minor_alleles` parameter added (off by default) - controls whether 
MungeSumstats should assume that the effects are majoritively measured on the 
minor alleles. Default is FALSE as this is an assumption that won't be 
appropriate in all cases. However, the benefit is that if we know the majority 
of SNPs have their effects based on the minor alleles, we can catch cases where 
the allele columns have been mislabelled. 

## CHANGES IN VERSION 1.13.6

### New features
* Mappings added to mapping file for risk and non risk allele.

## CHANGES IN VERSION 1.13.3

### Bug fix
* Bug fix for check 3 in infer effect column - previously A1 & A2 were swapped
when there were more matches for the ref genome in A1 rather than A2 which was
incorrect. Corrected now so it will only be flipped when A2 has more matches to
the reference genome.

## CHANGES IN VERSION 1.13.2

### New features
* Handling of -log10 p-values (outside of VCFs) added.

## CHANGES IN VERSION 1.13.1

### New features
* Mapping for OA (other Alllele) added to A1.

## CHANGES IN VERSION 1.11.10

### New features
* Can now pass local chain files for liftover (`local_chain` in
`format_sumstats()` and `liftover()`).

## CHANGES IN VERSION 1.11.9

### New features
* Can now control what columns are checked for missing data (`drop_na_cols` in
`format_sumstats()`). By default, SNP, effect columns and P/N columns are 
checked. Set to Null to check all columns or choose specific columns.

## CHANGES IN VERSION 1.11.7

### Bug fix
* Force no tab indexing when writing removed rows of SNPs. This avoids any 
issues where missing data causes sort errors.
* Issue fixed when sorting CHR column based on a format when CHR column is a 
factor.

## CHANGES IN VERSION 1.11.6

### Bug fix
* Catch for overflow when NA's in SNP col for `check_no_rs_snp()` check with 
`imputation_ind=TRUE`.

## CHANGES IN VERSION 1.11.4

### Bug fix
* Minor fix to `get_genome_builds()` to help with RAM & CPU usage during unit 
tests. No change in functionality for end user.

## CHANGES IN VERSION 1.11.3

### Bug fix
* For LDSC format, rename A1 and A2 as LDSC expects A1 to be the effect column 
rather than A2 (the opposite to MSS's default) - see more [here](https://groups.google.com/g/ldsc_users/c/S7FZK743w68).
Although, this didn't seem to make any difference to results in tests, see more
[here](https://github.com/neurogenomics/MungeSumstats/issues/160#issuecomment-1891899253).

## CHANGES IN VERSION 1.11.2

### Bug fix
* Remove unused argument `make_ordered` from `sort_coords()`
* Issue fixed with check ldsc format wehn compute_n type chosen 

## CHANGES IN VERSION 1.11.1

### Bug fix
* Speed up unit test timing for bioc checks (predominately for linux tests)

## CHANGES IN VERSION 1.9.19

### New features
* `infer_eff_direction` parameter added so user can decide whether to run the check

### Bug fix
* Typo in unit test for infer effect direction.
* IEU GWAS unit tests updated to account for server outages.

## CHANGES IN VERSION 1.9.18

### Bug fix
* Fixed column header mappings
  * Made all uncorrected header names uppercase and removed duplicates
  * "TOTALSAMPLESIZE" now maps to "N" instead of "NSTUDY"
  * "MAJORALLELE", "MAJOR_ALLELE", "MAJOR-ALLELE", and "MAJOR ALLELE" now map to 
    "A1" instead of "A2"
  * Removed the mappings for "OR-A1", "OR.A1", "OR_A1", and "BETA1" because MSS 
    assumes that A2 is the effect allele
  * Removed mappings for "A1FREQ", "A1FRQ", "AF1", "FREQ.A1.1000G.EUR", 
    "FREQ.A1.ESP.EUR", "FREQ.ALLELE1.HAPMAPCEU", "FREQ1", "FREQ1.HAPMAP", and 
    "FRQ_A1" because MSS defines "FRQ" to be the allele frequency of A2
  * Removed mappings for "CHR36", "BASE_GRCH36", "POSITION36", "POSGRCH36", 
    "BASEGRCH36", "POS36", "POS GRCH36", "POS.GRCH36", "POS-GRCH36", and 
    "POS_GRCH36" 
    because MSS does not support the GRCh36 genome build
  * Removed the ambiguous mapping "NMISS" -> "N" because "NMISS" can refer to 
    the number of samples with missing data
  * Removed the ambiguous mapping "WEIGHT" -> "N" because "WEIGHT" can refer to 
    coefficient weights
* Fixed inference of allele where ambiguous (A1, A2) naming used (see 
  infer_effect_column.R for code) but in short:
  * Three checks now made to infer which allele the effect/frequency information
    relates to. See infer_effect_column.R for further details.
  * See get_eff_frq_allele_combns.R for how effect/frequency columns that infer 
    the allele are captured in the mapping file

### New features
* New column header mappings:
  * "VARIANT_ID" and "RSIDS" --> "SNP"
  * "P_BOLT_LMM" --> "P"
  * "NCASES" --> "N_CAS"
  * "N_EFFECTIVE", "N_INFORMATIVE", and "TOTAL_N" --> "N"
  * "HET_P" --> "HETPVAL"
  * "HET_ISQ" --> "HETISQT"
  * "ALL_AF" --> "FRQ"
  * "DIRECT" --> "DIRECTION"
  * "ALT_EFFSIZE" --> "BETA"
  * "INFORMATIVE_ALT_AC" --> "AC"

## CHANGES IN VERSION 1.9.17

### Bug fix
* Cases checking ref genome where there are no indels would sometimes cause an 
error when joining. This resolved this issue.

## CHANGES IN VERSION 1.9.16

### New features
* flip_frq_as_biallelic parameter added enabling frequencies of non-bi-allelic
SNPs to be flipped as if they were bi-allelic (1 - frequency) i.e. ignoring the
frequencies of other alternative alleles (assuming these will be negligible). 
Note this will not be done as default as it is not fully correct but may be 
useful for some users.

## CHANGES IN VERSION 1.9.15

### Bug fix
* Fix for imputation column when imputing RS ID from CHR:BP. Avoids crash and 
ensures correct identification of imputed SNPs.
* Avoid running compute_nsize function when no imputation is wanted by user - 
also avoids message output in this situation.

## CHANGES IN VERSION 1.9.14

### Bug fix
* Fix reporting of genome-wide sign variants before formatting.

## CHANGES IN VERSION 1.9.13

### Bug fix
* In `check_bp_range` ensure that the BP column is numeric.

## CHANGES IN VERSION 1.9.12

### Bug fix
* In `check_no_rs_snp` the order of operations had to be reversed to ensure all 
values were present before sorting column headers when `imputation_ind=TRUE` and
imputing rsIDs.

## CHANGES IN VERSION 1.9.11

### New features
* The `rmv_chrPrefix` parameter in `format_sumstats()` has been replaced with
the new `chr_style` parameter, which allows users to specify their desired
chromosome name style. The supported chromosome styles are "NCBI", "UCSC", "dbSNP",
and "Ensembl" with "Ensembl" being the default.
* `check_chr()` now automatically removes all SNPs with nonstandard CHR entries
(anything other than 1-22, X, Y, and MT in the Ensembl naming style).

## CHANGES IN VERSION 1.9.10

### Bug fix
* Better method to detect vcf files - looks for vcf in extension not in name.

## CHANGES IN VERSION 1.9.9

### Bug fix
* Check ref genome change - if not match found for either genome build, an error
will now be thrown.
* Checks has been added so that if chrom col has chr as a prefix, this will be 
removed before testing genome build.

## CHANGES IN VERSION 1.9.8

### Bug fix
* Bug fix when using imputation_ind with NA in chr column.

## CHANGES IN VERSION 1.9.7

### New features
* `ignore_multi_trait` parameter added which will ignore any multi-trait 
p-values if set to TRUE. By default it is false to maintain the current default
running conditions for MSS.

## CHANGES IN VERSION 1.7.18

### New features
* Check added, ensure BP is between 1 - length of chromosome using reference
chromosome.

## CHANGES IN VERSION 1.7.17

### New features
* extra mapping for base-pair position (BP) column added

## CHANGES IN VERSION 1.7.14

### Bug fix

* Fix ensembl chain file retrieval so works on all environments

## CHANGES IN VERSION 1.7.13

### Bug fix

* `write_sumstats`:
    - Fix indexing issues due to incomplete genome coordinates sorting: 
        https://github.com/neurogenomics/MungeSumstats/issues/117
    - Add default `NULL` to `ref_genome`. 
    - Check `ref_genome` (only in conditions where its used). 
* `sort_coord`:
    - Renamed .R file from *sort_coordinates* to match current function name.
    - Add multiple `sort_methods`, 
        including improved/more robust `data.table`-native method.
    - Added dedicated unit tests within `test-index_tabular.R`.
* New helper function: `check_numeric`:
    - Ensures relevant sumstats cols are numeric.
    - Added internally to: `sort_coord`, `read_header`
* *rworkflows.yml*:
    - Omit Windows runner.
    - Turn on `run_biocheck`
* *to_GRanges.R* / *to_VRanges.R*:
    - Rename files to match current function names.
* Remove extra extdata files (I think these were created by accident):
   - ALSvcf.vcf.bgz
   - ALSvcf.vcf.bgz.bgz
   - ALSvcf.vcf.bgz.bgz.tbi
   - ALSvcf.vcf.bgz.tbi
   - ALSvcf.vcf.gz
* Remove *.DS_Store* files throughout.
* Don't check for duplicates based on RS ID with Indels, remove these first.

### New features

* Implement `rworkflows`. 
    - Removed old Dockerfile (not needed anymore) and workflow yaml.
* Add `drop_indels` parameter so a user can decide to remove indels from 
sumstats.    

## CHANGES IN VERSION 1.7.12

### Bug fix
* For downloading files use `sed -E` rather than `sed -r` as its compatible with 
mac which has issues with `sed -r`

### New features
* For instances where a single column contains CHR, BP, A1 and A2. The default 
order has been updated to CHR:BP:A1:A2 to align with  
[SPDI format](https://www.ncbi.nlm.nih.gov/variation/notation/). If your format 
differs and MSS doesn't pick up on it, update the column name to the true format
e.g. CHR:BP:A2:A1

## CHANGES IN VERSION 1.7.11

### New features
* Update to where SNP column is given by the four CHR, BP, A1, A2. Now, if A1 or
A2 is also a separate column, these will be used to infer the order.

## CHANGES IN VERSION 1.7.10

### Bug fix
* further fix for Latex issues when rendering PDF of examples.

## CHANGES IN VERSION 1.7.9

### Bug fix
* fix for Latex issues when rendering PDF of examples.

## CHANGES IN VERSION 1.7.3

### Bug fix
* fix for offline runs and accessing chain files from 1.7.2.

## CHANGES IN VERSION 1.7.2

### New features
* New chain files used for lifting over the genome build from Ensembl have now 
been added. These will now be set as the default chain file instead of UCSC due
to [licensing issues](https://github.com/neurogenomics/MungeSumstats/issues/128).
The choice to use UCSC files will still be there but the files will not be 
stored in the package themselves, they will instead be downloaded for use on the
fly.

## CHANGES IN VERSION 1.7.1

### New features
* The use of the `log_folder` parameter in `format_sumstats()` has been updated. 
It is still used to point to the directory for the log files and the log of 
MungeSumstats messages to be stored. And the default is still a temporary 
directory. However, now the name of the log files (log messages and log outputs)
are the same as the name of the file specified in the `save_path` parameter with
the extension '_log_msg.txt' and '_log_output.txt' respectively.

## CHANGES IN VERSION 1.5.18

### Bug fix
* GHA fix.

## CHANGES IN VERSION 1.5.17

### New features
* By default ES taken as BETA new parameter added so users can specify if this 
isn't the case (`es_is_beta`). If set to FALSE, mapping removed.
* Imputing BETA ordering has been changed so log(OR) will be sued before 
calculating from Z, SE.

## CHANGES IN VERSION 1.5.16

### New features
* A new method for computing the Z-score of a sumstats (`compute_z` input) has 
been added: BETA/SE. To use it set `compute_z = 'BETA'` to continue to use the
P-value calculation use `compute_z = 'P'`. Note the default is stil 
`compute_z = FALSE`.

### Bug fix
* Remove erroneous print statement. 

## CHANGES IN VERSION 1.5.15

### Bug fixes

* Fix NA representation for tabular outputs - By default, `data.table::fread()` 
leaves NAs blank instead of including a literal NA. That's fine for CSVs and if 
the output is read in by fread, but it breaks other tools for TSVs and is hard 
to read. Updated that and added a message when the table is switched to 
uncompressed for indexing.

## CHANGES IN VERSION 1.5.14

### New features

* `read_header`: 
    - Can now read entire files by setting `n=NULL`.
    - Improved reading in of VCF files (can read *.vcf.bgz* now). 
    - Now exported.
    - Added unit tests. 
* Remove `seqminer` from all code (too buggy). 
* Automatically remove residual .tsv files after tabix indexing. 
* `import_sumstats`:
    - Use `@inheritDotParams format_sumstats` for better documentation. 
* `parse_logs`: Added new fields. 
* `format_sumstats`: Added time report at the end (minutes taken total). 
    Since this is a message, will be included in the logs, 
    and is now parsed by `parse_logs` and put into the column "time".

### Bug fixes

* `index_tabular`: Fixed by replacing `seqminer` with `Rsamtools`. 
* When SNP ID's passed with format `1:123456789`, it will now be dealt with 
appropriately.
* `compute_n` can't handle SNP level N values for imputation only population 
level. An explanatory error message has now been added.

## CHANGES IN VERSION 1.5.13

### Bug fixes

* Special characters causing issues with find empty columns function. Now fixed.

## CHANGES IN VERSION 1.5.12

### Bug fixes

* Mitchondrial (MT) SNPs' chromosome value were being forced to NA by 
sort_coords function. This has been fixed.

## CHANGES IN VERSION 1.5.11

### Bug fixes

* Had to pass check_dups to other checks so they also wouldn't be run. Now 
independent of non-biallelic check.

## CHANGES IN VERSION 1.5.10

### New features

* check_dups parameter added so duplicates won't be removed if formatting QTL 
datasets

## CHANGES IN VERSION 1.5.9

### Bug fixes

* validate_parameters checks for incorrect version of dbSNP package, corrected.

## CHANGES IN VERSION 1.5.6

### Bug fixes

* MSS can now impute CHR, BP at a SNP level. For cases where CHR and/or BP are 
NA but the RS ID is present, these will now be imputed fromt he reference 
genome. Note previously, this imputation was done when the chr and/or bp column 
was missing.
* Print statement from liftover silenced when no liftover required
* check missing data function will no longer remove cases with NA's in SNP_INFO 
column. The SNP_INFO column is created by MSS for cases with RS ID and some 
other information in the same SNP column (like rs1234:.....). Rather than throw
out this info, it is stored in a new column - SNP_INFO. However, the remove 
missing data function was also looking in this column to remove SNPs. This has 
been corrected.
* `find_sumstats()`: 
    - Fix N column in metadata.

## CHANGES IN VERSION 1.5.5

### New features

* save_format parameter created for format_sumstats. This will replace
ldsc_format which is now deprecated. Use save_format="LDSC" instead. Other 
options for save_format are generic standardised (NULL) and IEU Open GWAS VCF
format ("openGWAS").
* dbSNP version 155 has now been added. Users can now control the version of 
dbSNP to be used for imputation (144 or 155). Note that with the 9x more SNPs
in dbSNP 155 vs 144, run times will increase.

### Bug fixes

* Change where sex chromosomes were made lower case removed to match UCSC

## CHANGES IN VERSION 1.5.4

### New features

* Further mappings added

### Bug fixes

* Duplication of non-bi-allelic and indels fixed
* Correct compute_nsize documentation

## CHANGES IN VERSION 1.5.1 

### New features

* Export `vcf2df`.
    - Move some post-processing function inside this function
    (e.g. drop duplicate cols/rows). 
* `read_vcf` can now be parallised: splits query into chunks, imports them, and (optionally) converts them to `data.table` before rbinding them back into one object. 
    - Added report of VCF size (variants x samples) before processing to give
    user an idea of long it will take to process. 
    - Added arg `mt_thresh` to avoid using parallelisation when VCFs are small, 
    due to the overhead outweighing the benefits in these cases.
* Added Linux installation instructions for *axel* downloader.
* Added 2nd `tryCatch` to `downloader` with different `download.file` 
parameters that may work better on certain machines. 
* Avoid using `file.path` to specify URL in:
    - `get_chain_file`
    - `import_sumstats` 
* Allow `download_vcf` to pass URLs directly (without downloading the files) 
when `vcf_download=FALSE`.  
* `download_vcf`:
    - Make timeout 10min instead of 30min.
    - Make axel verbose. 
* `load_ref_genome_data`:
    - Give more informative messages that 
        let user know which steps take a long time.
    - Speed up substring preprocessing. 
* `read_vcf_genome`: more robust way to get genome build from VCF. 
* `read_sumstats`: Speed up by using `remove_empty_cols(sampled_rows=)`, 
    and only run for tabular file (`read_vcf` already does this internally).   
    

### Bug fixes 

* `select_vcf_field`: Got rid of "REF col doesn't exists" warning by omitting `rowRanges`. 
* Ensured several unevaluated code chunks in `vignettes/MungeSumstats.Rmd` were
surrounding by ticks.  
* `vcf2df`: Accounted for scenarios where `writeVcf` accidentally converts `geno`
data into redundant 3D matrices. 
    - Use `data.table::rbindlist(fill=TRUE)` to bind chunks back together. 
* Remove unused functions after `read_vcf` upgrades:
    - `infer_vcf_sample_ids`
    - `is_vcf_parsed`
    - `check_tab_delimited`
    - `read_vcf_data`
    - `remove_nonstandard_vcf_cols`
* Remove redundant `dt_to_granges` by merging functionality into `to_granges`.
    - Adjusted `liftover` to accommodate the slight change. 
* Fix `is_tabix` (I had incorrectly made `path` all lowercase).  
* Let `index_vcf` recognize all compressed vcf suffixes. 
    - Add extra error handling when .gz is not actually bgz-compressed. 
* Set `BiocParallel` registered threads back to 1 after 
    `read_vcf_parallel` finishes, to avoid potential 
    conflicts with downstream steps. 

## CHANGES IN VERSION 1.5.0 

### New features

* Added "query" column to `find_sumstats` output to keep track of 
search parameters.
* `import_sumstats`: 
    - Check if formatted file (`save_path`) exists
    *before* downloading to save time.  
    - Pass up `force_new` in additional to `force_new_vcf`. 
* Updated *Description* tag in DESCRIPTION file to better reflect the scope 
of `MungeSumstats`. 
* Upgraded `read_vcf` to be more robust. 
* Edited Deps/Suggests
    - Elevate `IRanges` to Imports. 
    - Remove `stringr` (no longer used)
* Add new internal function `is_tabix` to check whether a file is already 
tabix-indexed. 
* `read_sumstats`: 
    - now takes `samples` as an arg. 
    - Parallises reading VCF using `GenomicFiles`. 
* `read_sumstats`: now takes `samples` as an arg.  
By default, only uses first sample (if multiple are present in file).  
* Remove `INFO_filter=` from ALS VCF examples in vignettes 
(no longer necessary now that INFO parsing has been corrected). 
* `download_vcf` can now handle situations with `vcf_url=`
is actually a local file (not remote).

### Bug fixes

* AF (allele frequency) was accidentally
being assigned as INFO column in VCFs where the INFO rows started with "AF". 
This caused a large number of SNPs to be incorrectly dropped 
during the `check_info_score` step.
* If INFO score is not available, INFO column is now dropped entirely 
(rather than assigning all 1s). 
    - Adjusted *test-vcf_formatting* to reflect this. 
This avoids ambiguity about whether the INFO score is real or not.
* `check_info_score`:
    - Added extra messages in various conditions where INFO is not 
    used for filtering,
    and don't add `log_files$info_filter` in these instances.   
    - Added unit tests.  
* `check_empty_cols` was accidentally dropping more columns than it should have.
* Fix GHA pkgdown building: 
    - The newest version of [git introduced bugs when building pkgdown sites](https://github.com/actions/checkout/issues/760) 
    from within Docker containers (e.g. via my Linux GHA workflow). 
    Adjusting GHA to fix this. 
* Fix `write_sumstats` when indexing VCF. 
* Ensure `read_sumstats` can read in any VCF files 
(local/remote, indexed/non-indexed). 
* Fix `test-vcf_formatting.R`
    - line 51: had wrong AF value in string
    - line 109: encountering error? due to duplicate SNPs?
* Fix `test-check_impute_se_beta`
    - lines 51/52: `setkey` on SNP 
    (now automatically renamed from ID by `read_vcf`). 
* Fix `test-read_sumstats`:
    - standardising of headers is now handled internally by `read_sumstats`. 
    - Ensure CHR is a character vector when being read in.
    - line 44: Ensure extra cols in `vcf_ss` are dropped. 
* `parse_logs`: Add lines to parsing subfunctions to allow handling of logs 
that don't contain certain info 
(thus avoid warnings when creating the final data.table).  
* '*Avoid the use of 'paste' in condition signals*' fixed: 
    - `check_pos_se`
    - `check_signed_col`
* Used to rely on *gunzip* to read bgz files, but apparently this functionality 
is no longer supported (possibly due to changes to how `Rsamtools::bgzip` does 
compression in Bioc 3.15. Switched to using `fread + readLines` in:
    - `read_header`
    - `read_sumstats` 
* `read_header`: wasn't reading in enough lines to get past the VCF header.
Increase to `readLines(n=1000)`.  
* `read_vcf`: Would sometimes induce duplicate rows. 
Now only unique rows are used (after sample and columns filtering). 
* Issue with mix of chr:bp:a1:a2 and chr:bp and rs id resolved
    

## CHANGES IN VERSION 1.3.19

### Bug fixes

- `format_sumstats` can now import remote files (other than OpenGWAS). 

### New features  

* New `sumstatsColHeaders` entries:
    - "PosGRCh37" --> "BP"
    - "testedAllele" --> "A1"

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

