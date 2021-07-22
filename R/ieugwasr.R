
################################################################
## Functions extracted from ieugwasr because Bioconductor doesn't allow any remotes. 
## All code in this file was written by the ieugwasr team and is available on GitHub: 
## https://github.com/MRCIEU/ieugwasr
################################################################


#' Get list of studies with available GWAS summary statistics through API
#'
#' @param id List of MR-Base IDs to retrieve. If NULL (default) retrieves all available datasets
#' @param access_token Google OAuth2 access token. Used to authenticate level of access to data
#'
#' @importFrom magrittr %>%
#' @keywords internal
#' @return Dataframe of details for all available studies
gwasinfo <- function(id=NULL, access_token = check_access_token())
{
    id <- legacy_ids(id)
    if(!is.null(id))
    {
        stopifnot(is.vector(id))
        out <- api_query('gwasinfo', query = list(id=id), access_token=access_token) %>% get_query_content()
    } else {
        out <- api_query('gwasinfo', access_token=access_token) %>% get_query_content()
    }
    if(length(out) == 0)
    {
        return(dplyr::tibble())
    }
    out <- dplyr::bind_rows(out) %>%
        dplyr::select("id", "trait", dplyr::everything())
    class(out) <- c("GwasInfo", class(out))
    return(out)
}



#' Get access token for OAuth2 access to MR Base
#'
#'
#' @keywords internal
#' @return access token string
#' @importFrom googleAuthR gar_auth
get_access_token <- function()
{
    message("Using access token. For info on how this is used see logging_info()")
    tf <- basename(tempfile())
    check <- suppressWarnings(file.create(tf))
    if(!check)
    {
        stop("You are currently in a directory which doesn't have write access.\n",
             "  In order to authenticate we need to store the credentials in a file called '.httr-oauth'.\n",
             "  Please setwd() to a different directory where you have write access.")
    } else {
        unlink(tf)
    }
    a <- googleAuthR::gar_auth(email=TRUE)
    if(! a$validate())
    {
        a$refresh()
    }
    return(a$credentials$access_token)
}


#' Check if authentication has been made
#'
#' If a call to get_access_token() has been made then it will have generated mrbase.oauth. Pass the token if it is present, if not, return NULL and do not authenticate.
#'
#' @keywords internal
#' @return NULL or access_token depending on current authentication state
check_access_token <- function()
{
    if(file.exists("ieugwasr_oauth"))
    {
        return(get_access_token())
    } else {
        return(NULL)
    }
}


#' Convert current IDs to legacy IDs
#'
#' @param x Vector of ids
#'
#' @keywords internal
#' @return vector of back compatible ids
#' @importFrom dplyr tibble
legacy_ids <- function(x)
{
    if(is.null(x)) return(NULL)
    changes <- dplyr::tibble(
        old = c("UKB-a:", "UKB-b:", "UKB-c:", "IEU-a:", "\\D"),
        new = c("ukb-a-", "ukb-b-", "ukb-c-", "ieu-a-", "ieu-a-")
    )
    
    y <- x
    for(i in seq(1,nrow(changes)))
    {
        index <- grepl(changes$old[i], x)
        if(changes$old[i] == "\\D")
        {
            index <- !grepl(changes$old[i], x)
        }
        if(any(index))
        {
            if(changes$old[i] == "\\D")
            {
                x[index] <- paste0(changes$new[i], x[index])
            } else {
                x[index] <- gsub(changes$old[i], changes$new[i], x[index])
            }
        }
    }
    
    # met datasets
    index <- x %in% paste0("ieu-a-", 303:754)
    x[index] <- gsub("ieu-a-", "met-a-", x[index])
    index <- x %in% paste0("ieu-a-", 119:269)
    x[index] <- gsub("ieu-a-", "met-b-", x[index])
    index <- x %in% paste0("ieu-a-", 838:960)
    x[index] <- gsub("ieu-a-", "met-c-", x[index])
    
    overallindex <- y != x
    if(any(overallindex))
    {
        message("Deprecated IDs being used? Detected numeric IDs. Trying to fix, but please note the changes below for future.")
        message(paste(y[overallindex], " -> ", x[overallindex], collapse="\n"))
    }
    return(x)
}





#' Wrapper for sending queries and payloads to API
#'
#' There are a number of different GET and POST endpoints in the GWAS database API. This is a generic way to access them
#'
#' @param path Either a full query path (e.g. for get) or an endpoint (e.g. for post) queries
#' @param query If post query, provide a list of arguments as the payload. NULL by default
#' @param access_token Google OAuth2 access token. Used to authenticate level of access to data. By default, checks if already authenticated through \code{get_access_token} and if not then does not perform authentication
#' @param method GET (default) or POST, DELETE etc
#' @param silent TRUE/FALSE to be passed to httr call. TRUE by default
#' @param encode Default = json, see httr::POST for options
#' @param timeout Default = 300, avoid increasing this, preferentially simplify the query first.
#'
#' @keywords internal
#' @return httr response object
#' @importFrom httr add_headers timeout DELETE GET POST 
api_query <- function(path, query=NULL, access_token=check_access_token(), 
                        method="GET", silent=TRUE, encode="json", timeout=300)
{
	ntry <- 0
	ntries <- 5
	headers <- httr::add_headers(
		# 'Content-Type'='application/json; charset=UTF-8',
		'X-Api-Token'=access_token,
		'X-Api-Source'=ifelse(is.null(options()$mrbase.environment), 'R/TwoSampleMR', 'mr-base-shiny')
	)

	retry_flag <- FALSE

	while(ntry <= ntries)
	{
		if(method == "DELETE")
		{
			r <- try(
				httr::DELETE(
					paste0(options()$ieugwasr_api, path),
					headers,
					httr::timeout(timeout)
				),
				silent=TRUE
			)
		} else if(!is.null(query)) {
			r <- try(
				httr::POST(
					paste0(options()$ieugwasr_api, path),
					body = query, 
					headers,
					encode=encode,
					httr::timeout(timeout)
				),
				silent=TRUE
			)
		} else {
			r <- try(
				httr::GET(
					paste0(options()$ieugwasr_api, path),
					headers,
					httr::timeout(timeout)
				),
				silent=TRUE
			)			
		}
		if('try-error' %in% class(r))
		{
			if(grepl("Timeout", as.character(attributes(r)$condition)))
			{
				stop("The query to MR-Base exceeded ", timeout, 
				        " seconds and timed out. Please simplify the query")
			}
		}
		if(! 'try-error' %in% class(r))
		{
			if(r$status_code >= 500 & r$status_code < 600)
			{
				message("Server code: ", r$status_code, "; Server is possibly experiencing traffic, trying again...")
				retry_flag <- TRUE
				Sys.sleep(1)
			} else {
				if(retry_flag)
				{
					message("Retry succeeded!")
				}
				break
			}
		}
		ntry <- ntry + 1
	}

	if(r$status_code >= 500 & r$status_code < 600)
	{
		message("Server issue: ", r$status_code)
		message("Unable to retrieve results from server. See status in",
		        " the returned object and contact the developers if the ",
		        "problem persists.")
		return(r)
	}
	if('try-error' %in% class(r))
	{
		if(grepl("Could not resolve host", 
		            as.character(attributes(r)$condition)))
		{
			stop("The MR-Base server appears to be down, the following issue", 
			        " was received:\n", as.character(attributes(r)$condition))
		} else {
			stop("The following issue was encountered in trying to query the ",
			     "MR-Base server:\n",as.character(attributes(r)$condition)
			)
		}
	}

	return(r)
}


#' Parse out json response from httr object
#'
#' @param response Output from httr
#'
#' @keywords internal
#' @return Parsed json output from query, often in form of data frame. 
#' If status code is not successful then return the actual response.
#' @importFrom jsonlite fromJSON
#' @importFrom httr content status_code
get_query_content <- function(response)
{
	if(httr::status_code(response) >= 200 & httr::status_code(response) < 300)
	{
		o <- jsonlite::fromJSON(httr::content(response, "text", encoding='UTF-8'))
		if('eaf' %in% names(o)) 
		{
			o[["eaf"]] <- as.numeric(o[["eaf"]])
		}
		return(o)
	} else {
		return(response)
		# stop("error code: ", httr::status_code(response), "\n  message: ", jsonlite::fromJSON(httr::content(response, "text", encoding='UTF-8')))
	}
}




#' Toggle API address between development and release 
#' 
#' From \code{ieugwasr}.
#'
#' @return No return
#'
#' @param where Which API to use. Choice between "local", "release", "test". Default = "local"
#'
#' @keywords internal
select_api <- function(where="public")
{
    url <- switch(where,
                  public = "http://gwas-api.mrcieu.ac.uk/",
                  private = "http://ieu-db-interface.epi.bris.ac.uk:8082/",
                  dev1 = "http://localhost:8019/",
                  dev2 = "http://127.0.0.1:5000/"
    )
    if(is.null(url))
    {
        url <- options()$ieugwasr_api
        warning("A valid API was not selected. No change")
    }
    
    options(ieugwasr_api=url)
    message("API: ", where, ": ", url)
}



