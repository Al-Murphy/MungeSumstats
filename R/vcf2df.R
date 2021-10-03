#' Function to convert VariantAnnotation CollapsedVCF/ExpandedVCF
#'  objects to a data frame
#'  
#' @source https://gist.github.com/zhujack/849b75f5a8305edaeca1001dfb9c3fe9
#' 
#' @keywords internal
#' @importFrom Biostrings strsplit
vcf2df = function(v, 
                  expand = TRUE) {
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
                                    stringsAsFactors=F)
            dd <- lapply(lapply(anncol,`[`,1),
                         function(x){Biostrings:::strsplit(x,'\\|')[[1]]})
            ncls <- max(unlist(lapply(dd, length)))
            
            yy = data.frame(suppressWarnings(
                do.call(rbind,
                        c(dfannempty[1:ncls], dd))),
                            stringsAsFactors=FALSE)
            yy = data.frame(lapply(yy,type.convert))
            colnames(yy) = paste("ANN",anncols[1:ncls],sep="..")
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
                VariantAnnotation::info(header(x))[ann, ]$Description)
            df = df[, colnames(df) != ann]
            df = cbind(df, dfann)
        }
        n  = names(VariantAnnotation::geno(x))
        tmp = lapply(n, function(col) {
            return(as.data.frame(VariantAnnotation::geno(x)[[col]]))
        })
        ncols = sapply(tmp, ncol)
        tmp = do.call(cbind, tmp)
        colnames(tmp) = paste(rep(n, times = ncols), colnames(tmp),
                              sep = "_")
        df = cbind(df, tmp)
        df[sapply(df, is.list)] <- apply(df[sapply(df, is.list)], 2, 
                                         function(x) { 
                                             unlist(lapply(x, paste,  
                                                           sep=",", 
                                                           collapse=";"))
                                         } ) 
        return(df)
    }
    
    if (expand) {
        # message('Expanding VCF first, so number of rows may increase')
        return(v2df(VariantAnnotation::expand(v)))
    } else {
        # message('No VCF Expanding')
        return(v2df(v))
    }
}
