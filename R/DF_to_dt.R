#' DataFrame to data.table
#' 
#' Efficiently convert \link[S4Vectors]{DataFrame} to 
#' \link[data.table]{data.table}. 
#' @param DF \link[S4Vectors]{DataFrame} object.
#' @keywords internal
#' @source \href{https://support.bioconductor.org/p/66874/}{
#' Solution from Bioc forum}
#' @importFrom methods is
#' @importFrom Biostrings unstrsplit
#' @importFrom IRanges CharacterList
DF_to_dt <- function(DF){
    data.table::data.table(
        mapply(DF, 
               FUN=function(s){ 
                   # s <- DF[["REF"]]
                   # s <- DF[["ALT"]]
                   # s <- DF[[1]]
                   if(methods::is(s,"DNAStringSet") ){
                       s <- as.character(s)
                   } else if(methods::is(s,"DNAStringSetList")){
                       s <- IRanges::CharacterList(s)
                       s <- Biostrings::unstrsplit(s, sep=",")
                   } else {
                       s <- as.vector(s)
                   }
                   return(s)
               })
    )
}
