#' VCF to DF
#' 
#' Function to convert a \pkg{VariantAnnotation}
#'  \code{CollapsedVCF}/\code{ExpandedVCF}
#'  object to a \code{data.frame}. 
#' @param vcf Variant Call Format (VCF) file imported into R 
#' as a \pkg{VariantAnnotation} 
#' \link[VariantAnnotation]{CollapsedVCF}/
#' \link[VariantAnnotation]{ExpandedVCF} object. 
#' @param expand Expand data into multiple columns using
#'  \code{VariantAnnotation::expand}.
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
#' @importFrom Biostrings strsplit
#' @importFrom utils type.convert 
vcf2df <- function(vcf, 
                   expand = TRUE) {
    requireNamespace("VariantAnnotation")
    requireNamespace("MatrixGenerics")
    v2df <- function(x, 
                     ...) {
        ## Function to parse ANN column in to a dataframe
        .anncols = function(anncol,
                            headerstring) {
            anncols = Biostrings::strsplit(
                sub("Functional annotations: '",'',
                    headerstring),' \\| ')[[1]]
            dfannempty = data.frame(matrix(vector(), 0, length(anncols),
                                           dimnames=list(c(), anncols)),
                                    stringsAsFactors=FALSE)
            dd <- lapply(lapply(anncol,`[`,1),
                         function(x){Biostrings::strsplit(x,'\\|')[[1]]})
            ncls <- max(unlist(lapply(dd, length)))
            
            yy = data.frame(suppressWarnings(
                do.call(rbind,
                        c(dfannempty[seq(1,ncls)], dd))),
                            stringsAsFactors=FALSE)
            yy = data.frame(lapply(yy,utils::type.convert))
            colnames(yy) = paste("ANN",anncols[seq(1,ncls)],sep="..")
            return(yy)
        }
        
        
        df = as.data.frame(MatrixGenerics::rowRanges(x)) 
        df = cbind(df, 
                   as.data.frame(VariantAnnotation::info(x))
                   ) 
        if ( any(c('ANN', 'EFF') %in% names(VariantAnnotation::info(x))) ) {
            ann = c('ANN', 'EFF')[ c('ANN', 'EFF') %in% names(
                VariantAnnotation::info(x)) ][1]
            dfann = .anncols(
                df$ANN, 
                VariantAnnotation::info(
                    VariantAnnotation::header(x)
                )[ann, ]$Description)
            df = df[, colnames(df) != ann]
            df = cbind(df, dfann)
        }
        geno_data <- VariantAnnotation::geno(x)
        SNPs <- rownames(geno_data[[1]])
        n  = names(geno_data)
        tmp = lapply(n, function(col) {
            return(as.data.frame(geno_data[[col]]))
        })
        ncols = unlist(lapply(tmp, FUN = ncol))
        tmp = do.call(cbind, tmp)
        rownames(tmp) <- NULL 
        colnames(tmp) = paste(rep(n, times = ncols), colnames(tmp),
                              sep = "_")
        df = cbind(df, tmp)
        df[vapply(df, is.list, FUN.VALUE = logical(1))] <- 
            apply(df[vapply(df, is.list, FUN.VALUE = logical(1))], 2, 
                  function(x) { 
                      unlist(lapply(x, paste, sep=",", collapse=";"))  } ) 
        #### Add SNPs back in ####
        df <- cbind(SNP=SNPs, df)
        #### Remove duplicate rows ####
        df <- unique(df)
        return(df)
    }
    
    if (expand) {
        # message('Expanding VCF first, so number of rows may increase')
        return(v2df(VariantAnnotation::expand(x = vcf)))
    } else {
        # message('No VCF Expanding')
        return(v2df(x = vcf))
    }
}
