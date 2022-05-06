#' VCF to DF
#' 
#' Function to convert a \pkg{VariantAnnotation}
#'  \code{CollapsedVCF}/\code{ExpandedVCF}
#'  object to a \code{data.frame}. 
#' @param vcf Variant Call Format (VCF) file imported into R 
#' as a \pkg{VariantAnnotation} 
#' \link[VariantAnnotation]{CollapsedVCF}/
#' \link[VariantAnnotation]{ExpandedVCF} object. 
#' @param add_sample_names Append sample names to column names 
#' (e.g. "EZ" --> "EZ_ubm-a-2929").
#' @param verbose Print messages.
#' 
#' @source https://gist.github.com/zhujack/849b75f5a8305edaeca1001dfb9c3fe9
#' @source 
#' \code{
#' #### VariantAnnotation ####
#' vcf_file <- system.file("extdata", "ALSvcf.vcf",
#'                         package = "MungeSumstats")
#' vcf <- VariantAnnotation::readVcf(file = vcf_file)
#' vcf_df <- MungeSumstats::vcf2df(vcf = vcf)
#' } 
#' @source {
#' #### vcfR ####
#' if(!require("pinfsc50")) install.packages("pinfsc50")
#' vcf_file <- system.file("extdata", "pinf_sc50.vcf.gz", package = "pinfsc50")
#' vcf <- read.vcfR( vcf_file, verbose = FALSE )
#' vcf_df_list <- vcfR::vcfR2tidy(vcf, single_frame=TRUE)
#' vcf_df <- data.table::data.table(vcf_df_list$dat)
#' }
#' @return data.frame version of VCF 
#' 
#' @keywords internal
#' @importFrom utils type.convert 
#' @importFrom data.table as.data.table
vcf2df <- function(vcf, 
                   add_sample_names=TRUE,
                   verbose=TRUE) {
    requireNamespace("VariantAnnotation")
    requireNamespace("MatrixGenerics")
  
    messager("Converting VCF to data.table.",v=verbose) 
    #### Automatically determine whether to add sample names ####
    if(is.null(add_sample_names)){
        header <- VariantAnnotation::header(x = vcf)
        samples <- VariantAnnotation::samples(header)
        add_sample_names <- length(samples)!=1
    }
    samples <- VariantAnnotation::samples(header)
    #### .anncols function ####
    .anncols = function(anncol,headerstring) {
        anncols <- strsplit(sub("Functional annotations: '",'',
                               headerstring),' \\| ')[[1]]
        dfannempty <- data.frame(matrix(vector(), 0, length(anncols),
                                       dimnames=list(c(), anncols)),
                                stringsAsFactors=FALSE)
        yy <- data.frame(
            suppressWarnings(
                do.call(
                    rbind,
                    c(dfannempty,lapply(lapply(anncol,`[`,1),
                                        function(x){strsplit(x,'\\|')[[1]]})
                      )
                    )
                ),
            stringsAsFactors=FALSE)
        yy = data.frame(lapply(yy,type.convert))
        colnames(yy) = paste("ANN",anncols,sep="_")
        return(yy)
    }
    #### v2df function ####
    v2df <- function(x) {
        t1 <- Sys.time()
        #### This step takes the longest ####
        ## as.data.table is better than as.data.frame bc it can handle duplicate
        ## row names.
        # path <- "https://gwas.mrcieu.ac.uk/files/ubm-a-2929/ubm-a-2929.vcf.gz"
        # vcf <- VariantAnnotation::readVcf(file = path)
        # x <- vcf[1:1000000,] 
        #### Get rowranges only if available ####
        if("rowRanges" %in% methods::slotNames(x)){
            gr <- MatrixGenerics::rowRanges(x)
        } else {
            gr <- NULL
        }
        df <- data.table::data.table(
            ID = names(gr),
            granges_to_dt(gr = gr),
            DF_to_dt(DF = VariantAnnotation::info(x)) 
        )
        remove(gr)
        if('ANN' %in% colnames(df)) {
            dfann <- .anncols(
                anncol = df$ANN,
                headerstring = VariantAnnotation::info(
                    VariantAnnotation::header(x)
                    )['ANN',]$Description
            )
            df <- df[,colnames(df)!="ANN"]
            df <- cbind(df,dfann)
        }
        ##### Convert geno data ##### 
        ## Works better for VCFHeader
        if(methods::is(x,"VCFHeader")){
            tmp <- data.table::as.data.table(VariantAnnotation::geno(x), 
                                             keep.rownames = "name")
            # tmp <- DF_to_dt(DF = VariantAnnotation::geno(x))    
        } else {
            ## Works better for VCF
            n <- names(VariantAnnotation::geno(x))
            #### Avoid parsing redundant columns ####
            if("ID" %in% colnames(df)) n <- n[n!="ID"]
            tmp <- lapply(n,function(col) {
                coldat <- VariantAnnotation::geno(x)[[col]]
                #### Drop empty columns within each matrix ####
                for(i in ncol(coldat)){
                    if(sum(!is.na(coldat[,i]))==0){
                        coldat <- coldat[,-i]
                    }
                }
                if(ncol(coldat)==0 | nrow(coldat)==0){
                    return(NULL)
                } else {
                    ## keeps colnames unchanged
                    data.table::as.data.table(coldat)
                } 
            })
            #### Remove NULL entries ####
            nulls <- unname(unlist(lapply(tmp,is.null)))
            tmp <- tmp[!nulls]
            n <- n[!nulls]
            ## Each element can potentially have >1 column 
            ncols <- unlist(lapply(tmp,ncol))
            tmp <- do.call(cbind, tmp)
            if(isTRUE(add_sample_names)){
                colnames(tmp) = paste(rep(n, times = ncols), 
                                      colnames(tmp),sep = "_") 
            } else {
                colnames(tmp) <- rep(n, times = ncols)
            } 
        } 
        df <- cbind(df, tmp) 
        if(verbose) methods::show(round(difftime(Sys.time(),t1),1))
        return(df)
    }
    #### Call functions ####
    if (methods::is(vcf,"CollapsedVCF")) {
        # message('Expanding VCF first, so number of rows may increase')
        return(v2df(VariantAnnotation::expand(x = vcf)))
    } else {
        # message('No VCF Expanding')
        return(v2df(x = vcf))
    } 
}
