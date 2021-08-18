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

