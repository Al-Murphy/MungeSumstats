#' Infers the genome build of the summary statistics file (GRCh37 or GRCh38)
#' from the data. Uses SNP (RSID) & CHR & BP to get genome build.
#'
#' @param sumstats data table/data frame obj of the summary statistics file for
#' the GWAS ,or file path to summary statistics file.
#' @param nThread Number of threads to use for parallel processes.
#' @param sampled_snps Downsample the number of SNPs used when inferring genome
#' build to save time.
#' @param standardise_headers Run
#' @param standardise_headers Run
#' \code{standardise_sumstats_column_headers_crossplatform}.
#' @param mapping_file \pkg{MungeSumstats} has a pre-defined
#' column-name mapping file
#' which should cover the most common column headers and their interpretations.
#' However, if a column header that is in your file is missing of the mapping we
#' give is incorrect you can supply your own mapping file. Must be a 2 column
#' dataframe with column names "Uncorrected" and "Corrected". See
#' \code{data(sumstatsColHeaders)} for default mapping and necessary format.
#' @param dbSNP version of dbSNP to be used (144 or 155). Default is 155.
#' @param header_only Instead of reading in the entire \code{sumstats} file,
#' only read in the first N rows where N=\code{sampled_snps}.
#' This should help speed up cases where you have to read in \code{sumstats}
#' from disk each time.
#' @param allele_match_ref Instead of returning the genome_build this will 
#' return the proportion of matches to each genome build for each allele 
#' (A1,A2).
#' @inheritParams format_sumstats
#' @inheritParams get_genome_builds 
#' 
#' @return ref_genome the genome build of the data
#' @importFrom data.table setDT :=
#' @keywords internal
get_genome_build <- function(sumstats,
                             nThread = 1,
                             sampled_snps = 10000,
                             standardise_headers = TRUE,
                             mapping_file = sumstatsColHeaders,
                             dbSNP=155,
                             dbSNP_tarball = NULL,
                             header_only = FALSE,
                             allele_match_ref = FALSE,
                             ref_genome = NULL,
                             chr_filt = NULL) {
    ### Add this to avoid confusing BiocCheck
    seqnames <- CHR <- SNP <- BP <- alt_alleles <- NULL 
    if(isFALSE(allele_match_ref)){
      message("Inferring genome build.")
    }
    # if not a data.table, must be a path
    if (!is.data.frame(sumstats)) {
        # Checking if the file exists should happen first
        if (!file.exists(sumstats)) {
            stop("Path to GWAS sumstats is not valid")
        }
        if (header_only) {
            # Read in N lines of the data
            message("Reading in only the first ",
                    sampled_snps, " rows of sumstats.")
            sumstats <- read_sumstats(
                path = sumstats,
                nThread = nThread,
                nrows = sampled_snps
            )
        } else {
            # Read in full data
            sumstats <- read_sumstats(
                path = sumstats,
                nThread = nThread
            )
        }
    } else {
        # ensure data table obj
        sumstats <- data.table::setDT(sumstats)
    }
    # need SNP ID column (RS ID) CHR and BP (POS) to infer build 
    # - check these are present, considering all known names
    if (standardise_headers) {
        sumstats_return <-
            standardise_sumstats_column_headers_crossplatform(
                sumstats_dt = sumstats,
                mapping_file = mapping_file
            )
        sumstats <- sumstats_return$sumstats_dt
    }

    err_msg <-
        paste0(
            "SNP ID column (RS ID), CHR and BP (POSITION) columns are needed ",
            "to infer the genome build. These could not be\nfound in your ",
            "dataset. Please specify the genome build manually to run ",
            "format_sumstats()."
        )
    # Infer genome build using SNP & CHR & BP
    if (!all(c("SNP", "CHR", "BP") %in% colnames(sumstats))) {
      #want it returned rather than throwing an error
      if(isTRUE(allele_match_ref)){
        return(err_msg)
      }
      else{  
        stop(err_msg)
      }
    }

    #### Do some filtering first to avoid errors ####
    nrow_org <- nrow(sumstats)
    sumstats <- sumstats[complete.cases(SNP, BP, CHR)]
    err_msg2 <-
        paste0(
            "SNP ID column (RS ID), CHR and BP (POSITION)",
            "columns are needed to",
            " infer the genome build.",
            "These contain too many\nmissing values in",
            " your dataset to be used.",
            "Please specify the genome build manually",
            " to run format_sumstats()"
        )
    # also remove common incorrect formatting of SNP
    sumstats <- sumstats[grepl("^rs", SNP), ]
    sumstats <- sumstats[SNP != ".", ]
    #also deal with common misformatting of CHR
    # if chromosome col has chr prefix remove it
    sumstats[, CHR := gsub("chr", "", CHR)]
    
    #for internal testing - filter to specified chromosomes
    if(!is.null(chr_filt)){
      sumstats <- sumstats[CHR %in% chr_filt]
    }
    
    # if removing erroneous cases leads to <min(10k,50% org dataset) will fail -
    # NOT ENOUGH DATA TO INFER
    nrow_clean <- nrow(sumstats)
    size_okay <- FALSE
    if (nrow_clean > sampled_snps || (nrow_clean != 0 && 
                                      (nrow_clean / nrow_org) > .5)) {
        size_okay<- TRUE
    }
    if (!size_okay) {
      #want it returned rather than throwing an error
      if(isTRUE(allele_match_ref)){
        return(err_msg2)
      }
      else{  
        stop(err_msg2)
      }
    }
    #### Downsample SNPs to save time ####
    if ((nrow(sumstats) > sampled_snps) && !(is.null(sampled_snps))) {
        snps <- sample(sumstats$SNP, sampled_snps)
    } else { # nrow(sumstats)<10k
        snps <- sumstats$SNP
    }

    sumstats <- sumstats[SNP %in% snps, ]

    #now split into functions two roles
    if(isTRUE(allele_match_ref)){
      #1. checking for matches for A1/A2 to ref genomes
      if (is.null(ref_genome)){
        #have to check multiple
        snp_loc_data_37 <- load_ref_genome_data(
          snps = snps,
          ref_genome = "GRCH37",
          dbSNP = dbSNP,
          dbSNP_tarball = dbSNP_tarball
        )
        snp_loc_data_38 <- load_ref_genome_data(
          snps = snps,
          ref_genome = "GRCH38",
          dbSNP = dbSNP,
          dbSNP_tarball = dbSNP_tarball
        )
        # convert CHR filed in ref genomes to character not factor
        snp_loc_data_37[, seqnames := as.character(seqnames)]
        snp_loc_data_38[, seqnames := as.character(seqnames)]
        # convert CHR filed in data to character if not already
        sumstats[, CHR := as.character(CHR)]
        # Now check which genome build has more matches to data
        num_37 <-
          nrow(snp_loc_data_37[sumstats, ,
                               on = c("SNP"="SNP","pos"="BP","seqnames"="CHR"),
                               nomatch = FALSE
          ])
        num_38 <-
          nrow(snp_loc_data_38[sumstats, ,
                               on = c("SNP"="SNP","pos"="BP","seqnames"="CHR"),
                               nomatch = FALSE
          ])
        
        if (num_37 > num_38) {
          ref_gen_num <- num_37
          ref_genome <- "GRCH37"
          snp_loc_data <- snp_loc_data_37
          
        } else {
          ref_gen_num <- num_38
          ref_genome <- "GRCH38"
          snp_loc_data <- snp_loc_data_38
        }
      }
      else{
        #only check one chosen
        snp_loc_data <- load_ref_genome_data(
          snps = snps,
          ref_genome = ref_genome,
          dbSNP = dbSNP,
          dbSNP_tarball = dbSNP_tarball
        )
        # convert CHR filed in ref genomes to character not factor
        snp_loc_data[, seqnames := as.character(seqnames)]
        # convert CHR filed in data to character if not already
        sumstats[, CHR := as.character(CHR)]
      }
      # Now check which allele has more matches to data
      # want to match on ref and alt alleles too 
      # need to take first alt from list to do this
      snp_loc_data[,alt_alleles := unlist(lapply(alt_alleles, 
                                                 function(x) x[[1]]))]
      num_a1 <-
        nrow(snp_loc_data[sumstats, ,
                          on = c("SNP"="SNP","pos"="BP","seqnames"="CHR",
                                 "ref_allele"="A1","alt_alleles"="A2"),
                          nomatch = FALSE
        ])
      num_a2 <-
        nrow(snp_loc_data[sumstats, ,
                          on = c("SNP"="SNP","pos"="BP","seqnames"="CHR",
                                 "ref_allele"="A2","alt_alleles"="A1"),
                          nomatch = FALSE
        ])

      if(num_a1>=num_a2){
        message("Effect/frq column(s) relate to A2 in the inputted sumstats")
        #this is what MSS expects so no action required
        switch_req <- FALSE
      }else{#num_a1<num_a2
        message("Effect/frq column(s) relate to A1 in the inputted sumstats")
        switch_req <- TRUE
      }
      
      return(switch_req)
    }  
    else{
      #2. checking for ref genome
      # otherwise SNP, CHR, BP were all found and can infer
      snp_loc_data_37 <- load_ref_genome_data(
        snps = snps,
        ref_genome = "GRCH37",
        dbSNP = dbSNP,
        dbSNP_tarball = dbSNP_tarball,
        chr_filt = chr_filt
      )
      snp_loc_data_38 <- load_ref_genome_data(
        snps = snps,
        ref_genome = "GRCH38",
        dbSNP = dbSNP,
        dbSNP_tarball = dbSNP_tarball,
        chr_filt = chr_filt
      )
      # convert CHR filed in ref genomes to character not factor
      snp_loc_data_37[, seqnames := as.character(seqnames)]
      snp_loc_data_38[, seqnames := as.character(seqnames)]
      # convert CHR filed in data to character if not already
      sumstats[, CHR := as.character(CHR)]
      # Now check which genome build has more matches to data
      num_37 <-
        nrow(snp_loc_data_37[sumstats, ,
                             on = c("SNP"="SNP", "pos"="BP", "seqnames"="CHR"),
                             nomatch = FALSE
        ])
      num_38 <-
        nrow(snp_loc_data_38[sumstats, ,
                             on = c("SNP"="SNP", "pos"="BP", "seqnames"="CHR"),
                             nomatch = FALSE
        ])
      #if no matches throw error
      if(num_37==0 && num_38==0){
        msg_err <- 
          paste0("No matches found in either reference genome for your ",
                 "SNPs.\nPlease check their formatting (SNP, CHR and BP",
                 " columns) or supply the genome build.")
        stop(msg_err)
      }
      if (num_37 > num_38) {
        ref_gen_num <- num_37
        ref_genome <- "GRCH37"
      } else {
        ref_gen_num <- num_38
        ref_genome <- "GRCH38"
      }
      
      message("Inferred genome build: ", ref_genome)
      # add a warning if low proportion of matches found
      msg <- paste0(
        "WARNING: Less than 10% of your sampled SNPs matched that of ",
        "either reference genome, this may question the quality of ",
        "your summary statistics file."
      )
      if (ref_gen_num / length(snps) < 0.1) {
        message(msg)
      }
      
      return(ref_genome)
      
    }
}
