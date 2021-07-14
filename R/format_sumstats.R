#' Check that summary statistics from GWAS are in a homogeneous format
#'
#' @return The address for the modified sumstats file
#'
#' @examples
#' #Pass path to Educational Attainment Okbay sumstat file to a temp directory
#' 
#' eduAttainOkbayPth <- system.file("extdata","eduAttainOkbay.txt", 
#'                                  package="MungeSumstats")
#' 
#' ## Call uses reference genome as default with more than 2GB of memory,
#' ## which is more than what 32-bit Windows can handle so remove certain checks
#' 
#' is_32bit_windows <- .Platform$OS.type == "windows" && .Platform$r_arch == "i386"
#' if (!is_32bit_windows) {
#' reformatted <- MungeSumstats::format_sumstats(path=eduAttainOkbayPth,
#'                                               ref_genome="GRCh37")
#' } else{
#' reformatted <- MungeSumstats::format_sumstats(path=eduAttainOkbayPth,
#'                                               ref_genome="GRCh37",
#'                                               on_ref_genome = FALSE,
#'                                               strand_ambig_filter=FALSE,
#'                                               bi_allelic_filter=FALSE,
#'                                               allele_flip_check=FALSE)
#' }
#' #returned location has the updated summary statistics file
#' @param path Filepath for the summary statistics file to be formatted
#' @param ref_genome name of the reference genome used for the GWAS ("GRCh37" or "GRCh38"). 
#' Argument is case-insensitive. Default is NULL which infers the reference genome from the data.
#' @param convert_small_p Binary, should p-values < 5e-324 be converted to 0? Small p-values pass the R limit and can cause errors with LDSC/MAGMA and should be converted. Default is TRUE.
#' @param convert_n_int Binary, if N (the number of samples) is not an integer, should this be rounded? Default is TRUE.
#' @param analysis_trait If multiple traits were studied, name of the trait for analysis from the GWAS. Default is NULL.
#' @param INFO_filter numeric The minimum value permissible of the imputation information score (if present in sumstatsfile). Default 0.9.
#' @param N_std numeric The number of standard deviations above the mean a SNP's N is needed to be removed. Default is 5.
#' @param rmv_chr vector or character The chromosomes on which the SNPs should be removed. Use NULL if no filtering necessary. Default is X, Y and mitochondrial. 
#' @param on_ref_genome Binary Should a check take place that all SNPs are on the reference genome by SNP ID. Default is TRUE.
#' @param strand_ambig_filter Binary Should SNPs with strand-ambiguous alleles be removed. Default is FALSE.
#' @param allele_flip_check Binary Should the allele columns be checked against reference genome to infer if flipping is necessary. Default is TRUE.
#' @param bi_allelic_filter Binary Should non-biallelic SNPs be removed. Default is TRUE.
#' @param sort_coordinates Whether to sort by coordinates.
#' @param nThread Number of threads to use for parallel processes. 
#' @param save_path File path to save formatted data. Defaults to \code{tempfile(fileext=".tsv.gz")}.
#' @param write_vcf Whether to write as VCF (TRUE) or tabular file (FALSE). 
#' @param tabix_index Index the formatted summary statistics with \href{http://www.htslib.org/doc/tabix.html}{tabix} for fast querying. 
#' @param return_data Return \code{data.table} directly to user. Otherwise, return the path to the save data. Default is FALSE.
#' @inheritParams convert_sumstats 
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom data.table setcolorder
#' @importFrom utils read.table
#' @importFrom utils data
#' @export
format_sumstats <- function(path,
                            ref_genome=NULL, 
                            convert_small_p=TRUE,
                            convert_n_int=TRUE, 
                            analysis_trait=NULL,
                            INFO_filter=0.9, 
                            N_std=5, 
                            N_dropNA=TRUE,
                            rmv_chr=c("X","Y","MT"),
                            on_ref_genome=TRUE,
                            strand_ambig_filter=FALSE, 
                            allele_flip_check=TRUE,
                            bi_allelic_filter=TRUE,
                            sort_coordinates=TRUE,
                            nThread=1,
                            save_path=tempfile(fileext=".tsv.gz"),
                            write_vcf=FALSE,
                            tabix_index=FALSE,
                            return_data=FALSE,
                            return_format="data.table",
                            force_new=FALSE
                            ){  
  #### Setup multi-threading ####
  data.table::setDTthreads(threads = nThread)
  #### Setup empty variables ####
  rsids = NULL
  orig_dims = NULL
  
  #### Check 1: Ensure save_path is correct.   #### 
  check_save_out <- check_save_path(save_path = save_path, 
                                    write_vcf = write_vcf)
  if(tabix_index && sort_coordinates==FALSE){
    message("Setting `sort_coordinates=TRUE` in order to tabix-index results.")
    sort_coordinates <- TRUE
  }
  
  #### Recognize previously formatted files ####
  if(file.exists(check_save_out$save_path) && force_new==FALSE){
    message("Importing previously formatted file. Set `force_new=TRUE` to override this.")
    message("    ",check_save_out$save_path)
  } else { 
    #Avoid reloading ref genome every time, save it to this parent environment
    #after being made once - speed up code
    
    #Check input parameters
    validate_parameters(path=path,
                        ref_genome=ref_genome, 
                        convert_small_p=convert_small_p, 
                        convert_n_int=convert_n_int, 
                        analysis_trait=analysis_trait, 
                        INFO_filter=INFO_filter, 
                        N_std=N_std,
                        N_dropNA=N_dropNA,
                        rmv_chr=rmv_chr,
                        on_ref_genome=on_ref_genome, 
                        strand_ambig_filter=strand_ambig_filter, 
                        allele_flip_check=allele_flip_check,
                        bi_allelic_filter=bi_allelic_filter,
                        write_vcf=write_vcf)
    # This almost surely modifies the file (since most sumstats from different
    # studies are differently formatted), so it makes more sense to just make a
    # temporary file <tmp>, and return the address of the temp 
    
    ####  Check 2: Check input format and import ####
    sumstats_return <- list()
    sumstats_return[["sumstats_dt"]] <- read_sumstats(path = path, 
                                                      nThread = nThread) 
    
    #### Check 3:Standardise headers for all OS ####
    sumstats_return <-
      standardise_sumstats_column_headers_crossplatform(sumstats_dt = sumstats_return$sumstats_dt,
                                                        path =  path) 
    
    ### Report the number of SNP/CHR/etc. before any filtering (but after header formatting)
    report_summary(sumstats_dt = sumstats_return$sumstats_dt)
    orig_dims <- dim(sumstats_return$sumstats_dt)
    
    #### Check 4: Check if multiple models used or multiple traits tested in GWAS ####
    sumstats_return <-  
      check_multi_gwas(sumstats_dt = sumstats_return$sumstats_dt,
                       path = path, 
                       analysis_trait = analysis_trait)
    
    #### Infer reference genome if necessary ####
    if(is.null(ref_genome))
      ref_genome <- get_genome_build(sumstats = sumstats_return$sumstats_dt, 
                                     standardise_headers = FALSE, ## Already done previously
                                     sampled_snps = 10000)
    
    #### Check 5: Check for uniformity in SNP col - no mix of rs/missing rs/chr:bp ####
    sumstats_return <- 
      check_no_rs_snp(sumstats_dt = sumstats_return$sumstats_dt,
                      path = path, 
                      ref_genome = ref_genome)
    
    #### Check 28: Check for combined allele column (A1 and A2) ####
    sumstats_return <- 
      check_allele_merge(sumstats_dt = sumstats_return$sumstats_dt, 
                         path = path)
    
    col_headers <- names(sumstats_return$sumstats_dt)
    
    # Series of checks if CHR or BP columns aren't present
    if(sum(c("CHR","BP") %in% col_headers)!=2){
      msg <- paste0("Summary statistics file does not have obvious CHR/BP colum",
                    "ns. Checking to see if they are joined in another column")
      message(msg)
      
      ####Check 6: check if CHR:BP:A2:A1 merged to 1 column
      sumstats_return <- check_four_step_col(sumstats_dt = sumstats_return$sumstats_dt, 
                                             path = path) 
      
      #### Check 7: check if there is a column of data with CHR:BP format ####
      sumstats_return <- check_two_step_col(sumstats_dt = sumstats_return$sumstats_dt,
                                            path = path)  
      #### Re-standardise in case the joined column headers were unusual ####
      sumstats_return <-
        standardise_sumstats_column_headers_crossplatform(sumstats_dt = sumstats_return$sumstats_dt, 
                                                          path = path)
    }
    
    #### Check 8: check if CHR and BP are missing but SNP is present ####
    sumstats_return <- 
      check_no_chr_bp(sumstats_dt = sumstats_return$sumstats_dt,
                      path = path, 
                      ref_genome = ref_genome,
                      rsids = rsids)
    rsids <- sumstats_return$rsids #update rsids
    sumstats_return$rsids <- NULL
    
    #### Check 9: check if CHR and BP are present but SNP is missing ####
    sumstats_return <- check_no_snp(sumstats_dt = sumstats_return$sumstats_dt, 
                                    path = path, 
                                    ref_genome = ref_genome)
    
    #### Check 25: check that all snps are present on reference genome ####
    sumstats_return <- check_on_ref_genome(sumstats_dt = sumstats_return$sumstats_dt, 
                                           path = path,
                                           ref_genome = ref_genome,
                                           on_ref_genome = on_ref_genome,
                                           rsids = rsids)
    rsids <- sumstats_return$rsids #update rsids
    sumstats_return$rsids <- NULL
    
    #### Check 10: check if SNP is present but A1 and/or A2 is missing ####
    sumstats_return <- 
      check_no_allele(sumstats_return$sumstats_dt, path, ref_genome, rsids)
    rsids <- sumstats_return$rsids #update rsids
    sumstats_return$rsids <- NULL
    
    #### Check 11: check that all the vital columns are present ###
    check_vital_col(sumstats_dt = sumstats_return$sumstats_dt)
    
    #### Check 12: check there is at least one signed sumstats column ###
    check_signed_col(sumstats_dt = sumstats_return$sumstats_dt)
    
    #### Check 13: check for allele flipping ####
    sumstats_return <- 
      check_allele_flip(sumstats_dt = sumstats_return$sumstats_dt,
                        path = path, 
                        ref_genome = ref_genome,
                        rsids = rsids,
                        allele_flip_check = allele_flip_check)
    rsids <- sumstats_return$rsids #update rsids
    sumstats_return$rsids <- NULL
    
    #### Check 14: check first three column headers are SNP, CHR, BP (in that order) ####
    sumstats_return <- check_col_order(sumstats_dt = sumstats_return$sumstats_dt,
                                       path = path)
    
    #### Check 15: Keep only rows which have the number of columns expected ####
    sumstats_return <- check_miss_data(sumstats_dt = sumstats_return$sumstats_dt, 
                                       path = path)
    
    #### Check 16: check for duplicated columns ####
    # The formatting process can (rarely) result in duplicated columns,
    # i.e. CHR, if CHR:BP is expanded and one already exists... delete duplicates
    sumstats_return <- check_dup_col(sumstats_dt = sumstats_return$sumstats_dt, 
                                     path = path)
    
    #### Check 17: check for small P-values (3e-400 or lower) ####
    sumstats_return <- 
      check_small_p_val(sumstats_dt = sumstats_return$sumstats_dt, 
                        path = path, 
                        convert_small_p = convert_small_p)
    
    #### Check 18: check is N column not all integers, if so round it up ####
    sumstats_return <- 
      check_n_int(sumstats_dt = sumstats_return$sumstats_dt, 
                  path = path, 
                  convert_n_int = convert_n_int)
    
    #### Check 19: check all rows have SNPs starting with SNP or rs, drop those don't ####
    sumstats_return <- check_row_snp(sumstats_dt = sumstats_return$sumstats_dt, 
                                     path = path)
    
    #### Check 20: check all rows for duplicated SNPs, remove any that are ####
    sumstats_return <- check_dup_snp(sumstats_dt = sumstats_return$sumstats_dt, 
                                     path = path)
    
    #### Check 21: check all rows for duplicated BPs, remove any that are ####
    sumstats_return <- check_dup_bp(sumstats_dt = sumstats_return$sumstats_dt, 
                                    path = path)
    
    #### Check 22: check for low INFO scores ####
    sumstats_return <- 
      check_info_score(sumstats_dt = sumstats_return$sumstats_dt, 
                       path = path, 
                       INFO_filter = INFO_filter)
    
    #### Check 23: check for N > X std dev above mean ####
    sumstats_return <- check_n_num(sumstats_dt = sumstats_return$sumstats_dt, 
                                   path = path, 
                                   N_std = N_std, 
                                   N_dropNA = N_dropNA)
    
    #### Check 24: check that no snps are on specific chromosomes ####
    sumstats_return <- check_chr(sumstats_dt = sumstats_return$sumstats_dt,
                                 path = path, 
                                 rmv_chr = rmv_chr)
    
    #### Check 26: check that all snps are not strand ambiguous ####
    sumstats_return <- check_strand_ambiguous(sumstats_dt = sumstats_return$sumstats_dt, 
                                              path = path,
                                              ref_genome = ref_genome, 
                                              strand_ambig_filter = strand_ambig_filter)
    
    #### Check 27: check for non-biallelic SNPS ####
    sumstats_return <- check_bi_allelic(sumstats_dt = sumstats_return$sumstats_dt, 
                                        path = path, 
                                        ref_genome = ref_genome, 
                                        bi_allelic_filter = bi_allelic_filter, 
                                        rsids = rsids)
    rsids <- sumstats_return$rsids #update rsids
    sumstats_return$rsids <- NULL
    
    #### Check 28: Sort rows by genomic coordinates ####
    sumstats_return$sumstats_dt <- sort_coords(sumstats_dt =  sumstats_return$sumstats_dt, 
                                               sort_coordinates = sort_coordinates)
    
    
    #### WRITE data.table TO PATH ####
    write_sumstats(sumstats_dt = sumstats_return$sumstats_dt, 
                   save_path=check_save_out$save_path,
                   sep=check_save_out$sep,
                   write_vcf = write_vcf,
                   tabix_index = tabix_index,
                   nThread = nThread)
    rm(rsids)#free up memory 
    
    #### Report summary #### 
    report_summary(sumstats_dt = sumstats_return$sumstats_dt,
                   orig_dims = orig_dims) 
  }
  
  

  #### Preview sumstats ####
  preview_sumstats(save_path = check_save_out$save_path,
                   nrows = 5L)
  
  if(return_data){
    message("Returning data directly.")
    #### Load data into memory when a pre-existing file is being used 
    if(!exists("sumstats_return")){ 
      sumstats_return <- list()
      sumstats_return[["sumstats_dt"]] <- read_sumstats(path = check_save_out$save_path, 
                                                        nThread = nThread)
    }
    out <- convert_sumstats(sumstats_dt = sumstats_return$sumstats_dt, 
                            return_format = return_format)
    return(out)
  } else {
    message("Returning path to saved data.")
    return(check_save_out$save_path) # Returns address of modified file
  } 
}




