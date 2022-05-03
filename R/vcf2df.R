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
#' 
#' @source https://gist.github.com/zhujack/849b75f5a8305edaeca1001dfb9c3fe9
#' @source 
#' \code{
#' vcf_file <- system.file("extdata", "ALSvcf.vcf",
#'                         package = "MungeSumstats")
#' vcf <- VariantAnnotation::readVcf(file = vcf_file)
#' vcf_df <- MungeSumstats::vcf2df(vcf = vcf)
#' } 
#' @return data.frame version of VCF 
#' 
#' @keywords internal
#' @importFrom utils type.convert 
#' @importFrom data.table as.data.table
vcf2df <- function(vcf, 
                   add_sample_names=TRUE) {
    requireNamespace("VariantAnnotation")
    requireNamespace("MatrixGenerics")
  
    messager("Converting VCF to data.table.") 
    #### .anncols function ####
    .anncols = function(anncol,headerstring) {
        anncols = strsplit(sub("Functional annotations: '",'',
                               headerstring),' \\| ')[[1]]
        dfannempty = data.frame(matrix(vector(), 0, length(anncols),
                                       dimnames=list(c(), anncols)),
                                stringsAsFactors=FALSE)
        yy = data.frame(
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
        df <- data.table::data.table(
            ID = names(MatrixGenerics::rowRanges(x)),
            granges_to_dt(gr = MatrixGenerics::rowRanges(x)),
            DF_to_dt(DF = VariantAnnotation::info(x)) 
        )
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
        # geno_dt <- DF_to_dt(DF = VariantAnnotation::geno(x))  
        n <- names(VariantAnnotation::geno(x))
        tmp <- lapply(n,function(col) { 
            ## keeps colnames unchanged
            data.table::as.data.table(
                VariantAnnotation::geno(x)[[col]]
            )
        })
        ## Each element can potentially have >1 column 
        ncols <- unlist(lapply(tmp,ncol))
        tmp <- do.call(cbind, tmp)
        if(isTRUE(add_sample_names)){
            colnames(tmp) = paste(rep(n, times = ncols), 
                                  colnames(tmp),sep = "_") 
        } else {
            colnames(tmp) <- rep(n, times = ncols)
        } 
        df <- cbind(df, tmp) 
        methods::show(round(difftime(Sys.time(),t1),1))
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
