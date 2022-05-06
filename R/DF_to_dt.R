#' DataFrame to data.table
#' 
#' Efficiently convert \link[S4Vectors]{DataFrame} to 
#' \link[data.table]{data.table}. 
#' @source \href{https://support.bioconductor.org/p/66874/}{
#' Solution from Bioc forum}
#' 
#' @param DF \link[S4Vectors]{DataFrame} object.
#' 
#' @keywords internal
#' @importFrom data.table as.data.table data.table
#' @importFrom methods is
#' @importFrom Biostrings unstrsplit
#' @importFrom IRanges CharacterList
#' @returns VCF data in data.table format.
DF_to_dt <- function(DF){ 
    #### Check if DF is empty ####
    # DF <- cbind(DF,dummy=NA)
    if(nrow(DF)==0 | ncol(DF)==0) return(data.table::data.table())
    m <- mapply(DF, 
                FUN=function(s){ 
                    # s <- DF[["REF"]]
                    # s <- DF[["ALT"]]
                    # s <- DF[[1]]
                    if(methods::is(s,"DNAStringSet") ){
                        s <- as.character(s)
                    } else if(methods::is(s,"DNAStringSetList")){
                        s <- IRanges::CharacterList(s)
                        s <- Biostrings::unstrsplit(s, sep=",")
                    } else if(methods::is(s,"NumericList")){
                        s <- as.numeric(s)
                    } else if(methods::is(s,"list")){
                        s <- unlist(s)
                    } else {
                        s <- as.vector(s)
                    }
                    #### Check if empty ####
                    if(all(is.na(s)) || 
                       all(s==".") || 
                       all(s=="")){
                        return(NULL)
                    } else {
                        return(s)
                    } 
            }) 
    ## Turn named list or matrix into a data.table
    data.table::as.data.table(m)
}
