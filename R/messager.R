#' Print messages
#'
#' Print messages with option to silence.
#'
#' @param ... Message input.
#' @param v Whether to print messages.
#' 
#' @return Null output.
#'
#' @keywords internal
messager <- function(..., v = TRUE) {
    if (v) {
        msg <- paste(...)
        message(msg)
    }
}