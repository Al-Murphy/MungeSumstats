#' Register cores
#' 
#' Register a multi-threaded instances using \pkg{BiocParallel}.
#' @inheritParams BiocParallel::SnowParam
#' @returns Null output.
register_cores <- function(workers=1,
                           progressbar=TRUE){
    requireNamespace("BiocParallel")
    BPPARAM <- if(.Platform$OS.type=="windows"){
        BiocParallel::SnowParam(workers = workers, 
                                progressbar = progressbar)
    } else {
        BiocParallel::MulticoreParam(workers = workers, 
                                     progressbar = progressbar)
    }
    BiocParallel::register(BPPARAM = BPPARAM)
}