#' Read in VCF format
#'
#' @inheritParams format_sumstats
#' @param temp_save Save temprorary file before proceeding,
#' @param keep_extra_cols Keep non-standard columns.
#' @return The modified sumstats_file
#' @keywords internal
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom data.table tstrsplit
#' @importFrom data.table :=
#' @importFrom dplyr %>% rename 
#' @importFrom utils head
#' @importFrom stringr str_split
read_vcf <- function(path, 
                     nThread=1, 
                     temp_save=FALSE,
                     keep_extra_cols=FALSE, 
                     save_path=tempfile()){ 
    ########## OTHER VCF READERS/WRITERS ###########
    # 1. [vcfR](https://cran.r-project.org/web/packages/vcfR/vcfR.pdf)
    # 2. [VariantAnnotation](https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html)
    # 3. [seqminer](https://cran.r-project.org/web/packages/seqminer/index.html) 
    # 4. [Rsamtools](https://bioconductor.org/packages/release/bioc/html/Rsamtools.html) 
    ##################################
     
    message("Reading VCF file.")  
    # fileformat <- gsub("^##fileformat=","",header[1])
    # #First get the name of data column, held in the ##SAMPLE row 
    sample_id <- get_vcf_sample_ids(path = path) 
    
    #### Read in full data ####
    
    #### If the VCF is remotely stored, data.table automatically 
    # downloads it to a tmpdir without telling you where. This lets you know where.
    if(startsWith(path,"https://gwas.mrcieu.ac.uk")){
        vcf_suffixes <- supported_suffixes(tabular = FALSE, tabular_compressed = FALSE)
        tmpdir <- file.path(tempdir(),gsub(paste(vcf_suffixes,collapse = "|"),"",basename(path)),"")
        dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
        message("Downloading VCF ==> ",tmpdir) 
    }else {
        tmpdir <- tempdir()
    }
    
    sumstats_file <- data.table::fread(path, nThread = nThread, sep="\t", skip = "#CHR", 
                                       tmpdir = tmpdir) 
    sumstats_file <- sumstats_file %>% dplyr::rename(CHROM="#CHROM")
    
    #### Infer sample ID from data colnames if necessary #### 
    sample_id <- infer_vcf_sample_ids(sample_id = sample_id,
                                      sumstats_file = sumstats_file) 
    
    ## Get format of FORMAT col for parsing
    format <- sumstats_file$FORMAT[1]
    format_cols <- stringr::str_split(format, ":")[[1]]
    
    # where INFO="." then can be missing AF**
    # Remove unnecessary cols - need sample_id and FORMAT column
    sumstats_file <- remove_nonstandard_vcf_cols(sample_id = sample_id,
                                                 sumstats_file = sumstats_file,
                                                 keep_extra_cols = keep_extra_cols)
     
    #sample_id and FORMAT col will look like:
    #               FORMAT                                           IEU-a-2
    #1: ES:SE:LP:AF:SS:ID  -0.0067:0.0145:0.193006:0.9322:109823:rs12565286
    #Values stand for:
    #>       Number Description
    #>    ES A      Effect size estimate relative to the alternative allele
    #>    SE A      Standard error of effect size estimate
    #>    LP A      -log10 p-value for effect estimate
    #>    AF A      Alternate allele frequency in the association study
    #>    SS A      Sample size used to estimate genetic effect
    #>    EZ A      Z-score provided if it was used to derive the EFFECT and...
    #>    SI A      Accuracy score of summary data imputation
    #>    NC A      Number of cases used to estimate genetic effect
    #>    ID 1      Study variant identifier
    #
    
    #if sample_id col present, split it out into separate columns
    message("Parsing '",sample_id,"' data column.")
    if(length(sample_id)!=0){
        #split out format into separate values
        # format <- strsplit(format,":")[[1]]
        #First, there is an issue where INFO=".", AF from FORMAT is missing
        #Need to impute 0 for these values
        #See dataset for example of this: 
        #http://fileserve.mrcieu.ac.uk/vcf/IEU-a-2.vcf.gz 
        if("AF" %in% format && 
           "INFO" %in% names(sumstats_file) && 
           nrow(sumstats_file[INFO=="."])>0){
            #get positions of ":" separators for each SNP
            find_splits <- gregexpr(":", sumstats_file[,get(sample_id)])
            #get position of AF to be imputed
            AF_pos <- which("AF" == format)
            #Update sumstats_file where number of splits isn't correct
            #get row identifiers
            update_rows<-
                (lengths(regmatches(sumstats_file[,get(sample_id)],find_splits))!=
                     length(format)-1)
            #get char pos for imputation
            find_splits <- unlist(lapply(find_splits,function(x) x[AF_pos-1]))
            #add to dt
            sumstats_file[,find_splits:=find_splits]
            #Now update these by imputing 0 for AF
            sumstats_file[update_rows,(sample_id):=
                              paste0(substr(get(sample_id),0,find_splits),
                                     "0:",#AF impute
                                     substr(get(sample_id),find_splits+1,
                                            nchar(get(sample_id))))]
            sumstats_file[,find_splits:=NULL]
            #Then replace INFO="." with 0 at later stage
        }
        
        #check if any cols already present - ID likely will be, rename if so
        if(any(format %in% names(sumstats_file))){
            format_copy <- format[format %in% names(sumstats_file)]
            for(format_i in format_copy){
                #If it is ID just remove other ID column
                if(format_i=="ID"){
                    sumstats_file[,(format_i):=NULL]
                }
                else{#otherwise keep both columns just rename one
                    data.table::setnames(sumstats_file,format_i,paste0(format_i,"2"))
                }
            }  
        }  
    
        
        sumstats_file[, (format_cols) :=
                          data.table::tstrsplit(get(sample_id),
                                                split=":", fixed=TRUE)]
        #Now remove sample_id column
        sumstats_file[,(sample_id):=NULL]
    }
    #MungeSumstats::sumstatsColHeaders contains mappings for
    #ID to SNP
    #EZ to Z
    #NC to N_CAS
    #SS to N
    #ES ro BETA* -need confirmation 
    
    #### Check for empty columns #### 
    empty_cols <- check_empty_cols(sumstats_file = sumstats_file, 
                                   sampled_rows=1000)

    
    if("MarkerName" %in% colnames(sumstats_file)){ 
        if("ID" %in% empty_cols){
            message("Replacing empty ID col with MarkerName col.")
            sumstats_file$ID <- sumstats_file$MarkerName
            sumstats_file[,("MarkerName"):=NULL]
        }
    }
    
    #Need to convert P-value, currently -log10 
    if("Pval" %in% colnames(sumstats_file)){
        sumstats_file <- sumstats_file %>% dplyr::rename(P=Pval)
    } 
    if("LP" %in% names(sumstats_file)){ 
        if("LP" %in% empty_cols){
            message("LP column is empty. Cannot compute raw p-values.")
        } else {
            msg <- paste0("VCF file has -log10 P-values, these will be ",
                          "converted to unadjusted p-values in the 'P' column.")
            message(msg)
            sumstats_file[,P:=10^(-1*as.numeric(LP))]
        }
    }
    
    #Need to remove "AF=" at start of INFO column and replace any "." with 0
    if("INFO" %in% names(sumstats_file)){
        message("Formatting INFO column.")
        sumstats_file[,INFO:=gsub("^AF=","",INFO)]
        sumstats_file[INFO==".",INFO:=0]
        #update to numeric
        sumstats_file[,INFO:=as.numeric(INFO)]
    }
    if("INFO" %in% names(empty_cols)){
        message("WARNING: All INFO scores are empty. Replacing all with 1.")
        sumstats_file$INFO <- 1
    }
    
    #VCF format has dups of each row, get unique
    sumstats_file <- unique(sumstats_file)
    
    # #write new data
    if(temp_save){
        message("Storing intermediate file before proceeding ==> ",path)
        data.table::fwrite(x = sumstats_file,
                           file = save_path,
                           nThread = nThread,
                           sep="\t") 
    } 
    return(sumstats_file)
}

 
