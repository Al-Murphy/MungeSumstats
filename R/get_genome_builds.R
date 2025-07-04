#' Infer genome builds
#'
#' Infers the genome build of  summary statistics files (GRCh37 or GRCh38)
#' from the data. Uses SNP (RSID) & CHR & BP to get genome build.
#'
#' Iterative version of \code{get_genome_build}.
#'
#' @param sumstats_list A named list of paths to summary statistics,
#' or a named list of \code{data.table} objects.
#' @param header_only Instead of reading in the entire \code{sumstats} file,
#' only read in the first N rows where N=\code{sampled_snps}.
#' This should help speed up cases where you have to read in \code{sumstats}
#' from disk each time.
#' @param sampled_snps Downsample the number of SNPs used when inferring genome
#' build to save time.
#' @param names_from_paths Infer the name of each item in \code{sumstats_list}
#' from its respective file path.
#' Only works if \code{sumstats_list} is a list of paths.
#' @param dbSNP version of dbSNP to be used (144 or 155). Default is 155.
#' @param dbSNP_tarball Pass local versions of dbSNP in tarball format. Default 
#' of NULL uses the dbSNP version passed in `dbSNP` parmeter. `dbSNP_tarball` 
#' was enabled to help with dbSNP versions >=156, after the decision to no 
#' longer provide dbSNP releases as bioconductor packages. dbSNP 156 tarball is 
#' available here: http://149.165.171.124/SNPlocs/.
#' @param nThread Number of threads to use for parallel processes.
#' @param chr_filt Internal for testing - filter reference genomes and sumstats
#' to specific chromosomes for testing. Pass a list of chroms in format: 
#' c("1","2"). Default is NULL i.e. no filtering
#'
#' @return ref_genome the genome build of the data
#'
#' @examples
#' # Pass path to Educational Attainment Okbay sumstat file to a temp directory
#'
#' eduAttainOkbayPth <- system.file("extdata", "eduAttainOkbay.txt",
#'     package = "MungeSumstats"
#' )
#' sumstats_list <- list(ss1 = eduAttainOkbayPth, ss2 = eduAttainOkbayPth)
#'
#' ## Call uses reference genome as default with more than 2GB of memory,
#' ## which is more than what 32-bit Windows can handle so remove certain checks
#' is_32bit_windows <-
#'     .Platform$OS.type == "windows" && .Platform$r_arch == "i386"
#' if (!is_32bit_windows) {
#'     
#'     #multiple sumstats can be passed at once to get all their genome builds:
#'     #ref_genomes <- get_genome_builds(sumstats_list = sumstats_list)
#'     #just passing first here for speed
#'     sumstats_list_quick <- list(ss1 = eduAttainOkbayPth)
#'     ref_genomes <- get_genome_builds(sumstats_list = sumstats_list_quick,
#'                                      dbSNP=144)
#' }
#' @export
#' @importFrom parallel mclapply
#' @importFrom utils capture.output
#' @importFrom dplyr %>%
#' @importFrom methods is
get_genome_builds <- function(sumstats_list,
                              header_only = TRUE,
                              sampled_snps = 10000,
                              names_from_paths = FALSE,
                              dbSNP=155,
                              dbSNP_tarball=NULL,
                              nThread = 1,
                              chr_filt = NULL) {
    start <- Sys.time()
    #### Convert to list if it isn't already ####
    if (!methods::is(sumstats_list, "list")) {
        sumstats_list <- list(sumstats_list)
    }
    message(
        "Inferring genome build of ", length(sumstats_list),
        " sumstats file(s)."
    )
    #### Infer names from path ####
    if (names_from_paths &
        (methods::is(sumstats_list[[1]], "character"))) {
        message("Inferring names of sumstats_list items from paths.")
        names(sumstats_list) <- gsub(
            paste(supported_suffixes(), collapse = "|"),
            "", basename(unlist(sumstats_list))
        )
    }
    #### Make sure names are unique ###
    names(sumstats_list) <- make.unique(names(sumstats_list))
    #### Assign new name ####
    if (is.null(names(sumstats_list))) {
        message(
            "sumstats_list has no names. ",
            "Assigning names using ss# format."
        )
        names(sumstats_list) <- paste0("ss", seq(1, length(sumstats_list)))
    }
    #### Infer builds ####
    #Weirdly more efficient to not use parallel::mcapply() if only one
    if(length(names(sumstats_list))==1){
      build_ <- get_genome_build(
        sumstats = sumstats_list[[1]],
        sampled_snps = sampled_snps,
        dbSNP = dbSNP, 
        dbSNP_tarball = dbSNP_tarball,
        header_only = header_only,
        nThread = nThread,
        chr_filt = chr_filt
      )
      builds <- list(build_)
      names(builds) <- c(names(sumstats_list))
    }else{
      builds <- parallel::mclapply(names(sumstats_list),
                                   function(x,
                                            .sampled_snps = sampled_snps,
                                            .dbSNP=dbSNP,
                                            .dbSNP_tarball = dbSNP_tarball,
                                            .header_only = header_only) {
                                     message_parallel(x)
                                     get_genome_build(
                                       sumstats = sumstats_list[[x]],
                                       sampled_snps = .sampled_snps,
                                       dbSNP = .dbSNP, 
                                       dbSNP_tarball = .dbSNP_tarball,
                                       header_only = .header_only,
                                       nThread = 1,
                                       chr_filt = chr_filt
                                     )
                                   },
                                   mc.cores = nThread) %>%
        `names<-`(names(sumstats_list)) 
    }
    #### Report time ####
    message(utils::capture.output(difftime(Sys.time(), start)))
    #### Report build counts ####
    build_counts <- table(unlist(builds))
    for (i in seq(1, length(build_counts))) {
        message(names(build_counts)[i], ": ", build_counts[i], " file(s)")
    }
    return(builds)
}
