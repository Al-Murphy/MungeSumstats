#' SEnd messages to console even from within parallel processes
#' @return A message
#' @keywords internal 
message_parallel <- function(...){
    system(sprintf('echo "%s"', paste0(..., collapse="")))
}