## CHANGES IN VERSION 1.1.1

### New Features

*   `import_sumstats` reads GWAS sum stats directly from Open GWAS. Now parallelised and reports how long each dataset took to import/format in total. 
*   `find_sumstats` searches Open GWAS for datasets. 
*   `compute_z` computes Z-score from P. 
*   `format_sumstats(ldsc_format=TRUE)` ensures sum stats can be fed directly into [LDSC](https://github.com/bulik/ldsc) without any additional munging. 
*   `read_sumstats`, `write_sumstas`, and `download_vcf` functions now exported.  
*   `format_sumstats(sort_coordinates=TRUE)` sorts results by genomic coordinates. 
*   `format_sumstats(return_data=TRUE)` returns data directly to user. Can be returned in either `data.table` (default), `GRanges` or `VRanges` format using `format_sumstats(return_format="granges")`.  
*   `format_sumstats(N_dropNA=TRUE)` (default) drops rows where N is missing.  
*   `format_sumstats(write_vcf=TRUE)` writes a tabix-indexed VCF file instead of tabular format. 
*   `format_sumstats(save_path=...)` lets users decide where their results are saved and what they're named. 
*   When the `save_path` indicates it's in `tempdir()`, message warns users that these files will be deleted when R session ends.  
*   Summary of data is given at the beginning and the end of `format_sumstats` via `report_summary()`.  
*   Readability of `preview_sumstats()` messages improved.  


### Bug fixes 

*   CHR column now standardised (X and Y caps, no "chr" prefix).
*   Allele flipping done on a per-SNP basis (instead of whole-column). 
*   Tests added for new functions. 
*   Existing tests made more robust by using `grep` searches for specific lines instead of 
*   New parameters added to `validate_parameters`.
*   `read_vcf` upgraded to account for more VCF formats. 
*   `check_n_num` now accounts for situations where N is a character vector and converts to numeric.  



## CHANGES IN VERSION 1.0.0

### New Features

*   MungeSumstats released to Bioconductor.

