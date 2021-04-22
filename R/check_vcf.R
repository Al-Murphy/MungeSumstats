#' Check if the inputted file is in VCF format, if so reformat as standard
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @return The modified sumstats_file
#' @keywords internal
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom data.table tstrsplit
#' @importFrom data.table :=
check_vcf <- function(sumstats_file, path){
  P = LP = INFO = NULL
  #if the file is a VCF, first line will look like: ##fileformat=VCFv4.2
  first_line <- sumstats_file[[1]]
  file_type <- gsub("^##fileformat=","",first_line)
  if(length(grep("^vcf",tolower(file_type)))==1){
    msg <- paste0("VCF format detected, this will be converted to a standard",
                  " summary statistics file format.")
    message(msg)
    #First get the name of data column, held in the ##SAMPLE row
    sample_id <- sumstats_file[grepl("^##SAMPLE",sumstats_file)]#gets ##SAMPLE
    sample_id <- gsub(",.*$", "", sample_id)#get rid of everything after ID
    sample_id <- substr(sample_id,10,nchar(sample_id))# get rid of ##SAMPLE=
    sample_id <- sub('.+=(.+)', '\\1', sample_id)# remove things before equals
    #Now remove all lines starting with ## to leave just the data rows
    sumstats_file <- sumstats_file[!grepl("^##",sumstats_file)]
    #Remove # at start of header line
    sumstats_file[1] <- gsub("^#","",sumstats_file[1])
    #Now write the data to the path and read in as DT
    writeLines(sumstats_file, con=path)
    sumstats_file <- data.table::fread(path)
    #Get format column values too
    format <- sumstats_file$FORMAT[1]#assume format same for all rows except 
    #where INFO="." then can be missing AF**
    #Remove unnecessary cols - need sample_id and FORMAT column
    colsToRemove <-
      names(sumstats_file)[!names(sumstats_file) %in% c("CHROM","POS","ID",
                                                          "REF","ALT","INFO",
                                                          sample_id)]
    sumstats_file[,(colsToRemove):=NULL]
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
    if(length(sample_id)!=0){
      #split out format into separate values
      format <- strsplit(format,":")[[1]]
      #First, there is an issue where INFO=".", AF from FORMAT is missing
      #Need to impute 0 for these values
      #See dataset for example of this: 
      #http://fileserve.mrcieu.ac.uk/vcf/IEU-a-2.vcf.gz
      if("AF" %in% format && "INFO" %in% names(sumstats_file) && 
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
      sumstats_file[, (format) :=
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

    #Need to convert P-value, currently -log10
    if("LP" %in% names(sumstats_file)){
      msg <- paste0("Inputted VCF format has -log10 P-values, these will be ",
                    "converted to unadjusted p-values in the 'P' column.")
      message(msg)
      sumstats_file[,P:=10^(-1*as.numeric(LP))]
    }
    
    #Need to remove "AF=" at start of INFO column and replace any "." with 0
    if("INFO" %in% names(sumstats_file)){
      sumstats_file[,INFO:=gsub("^AF=","",INFO)]
      sumstats_file[INFO==".",INFO:=0]
      #update to numeric
      sumstats_file[,INFO:=as.numeric(INFO)]
    }

    #VCF format has dups of each row, get unique
    sumstats_file <- unique(sumstats_file)

    #write new data
    data.table::fwrite(sumstats_file,file=path,sep="\t")
    sumstats_file <- readLines(path)

    return(sumstats_file)
  }
  else{
    return(sumstats_file)
  }
}
