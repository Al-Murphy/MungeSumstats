#' Search Open GWAS for datasets matching criteria
#'
#' For each  argument, searches for any datasets matching
#' a case-insensitive substring search in the respective metadata column.
#' Users can supply a single character string or a
#' list/vector of character strings.
#' 
#' To authenticate, you need to generate a token from the OpenGWAS website. The 
#' token behaves like a password, and it will be used to authorise the requests 
#' you make to the OpenGWAS API. Here are the steps to generate the token and 
#' then have `ieugwasr` automatically use it for your queries:
#' 1. Login to https://api.opengwas.io/profile/
#' 2. Generate a new token
#' 3. Add `OPENGWAS_JWT=<token>` to your .Renviron file, thi can be edited in R
#' by running `usethis::edit_r_environ()`
#' 4. Restart your R session
#' 5. To check that your token is being recognised, run 
#' `ieugwasr::get_opengwas_jwt()`. If it returns a long random string then you 
#' are authenticated.
#' 6. To check that your token is working, run `ieugwasr::user()`. It will make a 
#' a request to the API for your user information using your token. It should 
#' return a list with your user information. If it returns an error, then your 
#' token is not working.
#' 7. Make sure you have submitted use
#'
#' By default, returns metadata for all studies currently in Open GWAS database.
#'
#' @return (Filtered) GWAS metadata table.
#'
#' @param ids List of Open GWAS study IDs
#' (e.g. \code{c("prot-a-664", "ieu-b-4760")}).
#' @param traits List of traits
#' (e.g. \code{c("parkinson", "Alzheimer")}).
#' @param years List of years
#' (e.g. \code{seq(2015,2021)} or \code{c(2010, 2012, 2021)}).
#' @param consortia List of consortia
#'  (e.g. \code{c("MRC-IEU","Neale Lab")}.
#' @param authors List of authors
#' (e.g. \code{c("Elsworth","Kunkle","Neale")}).
#' @param populations List of populations
#' (e.g. \code{c("European","Asian")}).
#' @param categories List of categories
#' (e.g. \code{c("Binary","Continuous","Disease","Risk factor"))}).
#' @param subcategories List of categories
#' (e.g. \code{c("neurological","Immune","cardio"))}).
#' @param builds List of genome builds
#' (e.g. \code{c("hg19","grch37")}).
#' @param pmids List of PubMed ID (exact matches only)
#' (e.g. \code{c(29875488, 30305740, 28240269)}).
#' @param min_sample_size Minimum total number of study participants
#' (e.g. \code{5000}).
#' @param min_ncase Minimum number of case participants
#' (e.g. \code{1000}).
#' @param min_ncontrol Minimum number of control participants
#' (e.g. \code{1000}).
#' @param min_nsnp Minimum number of SNPs
#' (e.g. \code{200000}).
#' @param include_NAs Include datasets with missing metadata for size criteria
#' (i.e. \code{min_sample_size}, \code{min_ncase}, or \code{min_ncontrol}).
#' 
#' @export
#' @importFrom dplyr %>% arrange desc mutate rowwise
#' @importFrom data.table setorderv
#' @importFrom ieugwasr gwasinfo
#' @examples
#' # Only run the examples if user has internet access
#' # and if access token has been added
#' if(try(is.character(getURL("www.google.com")))==TRUE && ieugwasr::get_opengwas_jwt()!=""){
#' ### By ID
#' metagwas <- find_sumstats(ids = c(
#'     "ieu-b-4760",
#'     "prot-a-1725",
#'     "prot-a-664"
#' ))
#' ### By ID and sample size
#' metagwas <- find_sumstats(
#'     ids = c("ieu-b-4760", "prot-a-1725", "prot-a-664"),
#'     min_sample_size = 5000
#' )
#' ### By criteria
#' metagwas <- find_sumstats(
#'     traits = c("alzheimer", "parkinson"),
#'     years = seq(2015, 2021)
#' )
#' }
find_sumstats <- function(ids = NULL,
                          traits = NULL,
                          years = NULL,
                          consortia = NULL,
                          authors = NULL,
                          populations = NULL,
                          categories = NULL,
                          subcategories = NULL,
                          builds = NULL,
                          pmids = NULL,
                          min_sample_size = NULL,
                          min_ncase = NULL,
                          min_ncontrol = NULL,
                          min_nsnp = NULL,
                          include_NAs = FALSE
                          ) {
    
    ## Set up fake empty variables to avoid confusing BiocCheck
    sample_size <- ncase <- ncontrol <- nsnp <- N <- NULL;
    
    message("Collecting metadata from Open GWAS.")
    if (!is.null(ids)) {
        metagwas <- ieugwasr::gwasinfo(
            id = ids
        )
        ## gwasinfo() doesn't always return all columns for some reason
        missing_cols <- c("ncase", "ncontrol")
        missing_cols <- missing_cols[!missing_cols %in% colnames(metagwas)]
        if (length(missing_cols) > 0) {
            for (x in missing_cols) {
                metagwas[[x]] <- NA
            }
        }
    } else {
        metagwas <- ieugwasr::gwasinfo()
    }
    message("Filtering metadata by substring criteria.")
    if (!is.null(traits)) {
        metagwas <- metagwas[grepl(paste(traits, collapse = "|"),
                                   metagwas$trait,
                                   ignore.case = TRUE
        ), ]
    }
    if (!is.null(years)) {
        metagwas <- metagwas[grepl(paste(years, collapse = "|"),
                                   metagwas$year,
                                   ignore.case = TRUE
        ), ]
    }
    if (!is.null(consortia)) {
        metagwas <- metagwas[grepl(paste(consortia, collapse = "|"),
                                   metagwas$consortium,
                                   ignore.case = TRUE
        ), ]
    }
    if (!is.null(authors)) {
        metagwas <- metagwas[grepl(paste(authors, collapse = "|"),
                                   metagwas$author,
                                   ignore.case = TRUE
        ), ]
    }
    if (!is.null(populations)) {
        metagwas <- metagwas[grepl(paste(populations, collapse = "|"),
                                   metagwas$population,
                                   ignore.case = TRUE
        ), ]
    }
    if (!is.null(categories)) {
        metagwas <- metagwas[grepl(paste(categories, collapse = "|"),
                                   metagwas$category,
                                   ignore.case = TRUE
        ), ]
    }
    if (!is.null(subcategories)) {
        metagwas <- metagwas[grepl(paste(subcategories, collapse = "|"),
                                   metagwas$subcategory,
                                   ignore.case = TRUE
        ), ]
    }
    if (!is.null(builds)) {
        metagwas <- metagwas[grepl(paste(builds, collapse = "|"),
                                   metagwas$build,
                                   ignore.case = TRUE
        ), ]
    }
    if (!is.null(pmids)) metagwas <- metagwas[metagwas$pmid %in% pmids, ]
    
    if (any(
        !is.null(min_sample_size),
        !is.null(min_ncase),
        !is.null(min_ncontrol),
        !is.null(min_nsnp)
    )) {
        message("Filtering metadata by sample/case/control/SNP size criteria.")
        if (include_NAs) {
            message("Including sample/case/control size with NAs.")
            if (!is.null(min_sample_size)) {
                metagwas <- subset(
                    metagwas,
                    sample_size >= min_sample_size |
                        is.na(sample_size) |
                        sample_size == "NA"
                )
            }
            if (!is.null(min_ncase)) {
                metagwas <- subset(
                    metagwas,
                    ncase >= min_ncase |
                        is.na(ncase) |
                        ncase == "NA"
                )
            }
            if (!is.null(min_ncontrol)) {
                metagwas <- subset(
                    metagwas,
                    ncontrol >= min_ncontrol |
                        is.na(ncontrol) |
                        ncontrol == "NA"
                )
            }
            if (!is.null(min_nsnp)) {
                metagwas <- subset(
                    metagwas,
                    nsnp >= min_nsnp |
                        is.na(nsnp) |
                        nsnp == "NA"
                )
            }
        } else {
            message("Excluding sample/case/control size with NAs.")
            if (!is.null(min_sample_size)) {
                metagwas <- subset(metagwas, sample_size >= min_sample_size)
            }
            if (!is.null(min_ncase)) {
                metagwas <- subset(metagwas, ncase >= min_ncase)
            }
            if (!is.null(min_ncontrol)) {
                metagwas <- subset(metagwas, ncontrol >= min_ncontrol)
            }
            if (!is.null(min_nsnp)) {
                metagwas <- subset(metagwas, nsnp >= min_nsnp)
            }
        }
    }
    #### Add N col #### 
    metagwas <- metagwas |>
        dplyr::rowwise() |>
        dplyr::mutate(N = ifelse(is.na(sample_size),
                                 sum(ncase, ncontrol, na.rm = TRUE),
                                 sample_size
        )) |>
        data.table::data.table()
    metagwas[N %in% c(NA,0,"","NA",NULL),]$N <- NA
    #### Ensure data.table format ####
    metagwas <- data.table::data.table(metagwas)
    #### Add query col to keep track of groups ####
    param_list <- as.list(match.call())  
    metagwas$query <- list(param_list[names(param_list)!=""] )
    #### Sort results  #### 
    data.table::setorderv(metagwas,
                          cols = c(
                              "trait", "N", "sample_size",
                              "ncase", "ncontrol", "year"
                          ),
                          order = c(1, -1, -1, -1, -1, -1)
    )
    message(
        "Found ", formatC(nrow(metagwas), big.mark = ","),
        " GWAS datasets matching search criteria across:",
        "\n   - ", formatC(length(unique(metagwas$trait)),
                           big.mark = ","
        ), " trait(s)",
        "\n   - ", formatC(length(unique(metagwas$population)),
                           big.mark = ","
        ), " population(s)",
        "\n   - ", formatC(length(unique(metagwas$category)),
                           big.mark = ","
        ), " category(ies)",
        "\n   - ", formatC(length(unique(metagwas$subcategory)),
                           big.mark = ","
        ), " subcategory(ies)",
        "\n   - ", formatC(length(unique(metagwas$pmid)),
                           big.mark = ","
        ), " publication(s)",
        "\n   - ", formatC(length(unique(metagwas$consortium)),
                           big.mark = ","
        ), " consortia(ium)",
        "\n   - ", formatC(length(unique(metagwas$build)),
                           big.mark = ","
        ), " genome build(s)"
    ) 
    return(metagwas)
}