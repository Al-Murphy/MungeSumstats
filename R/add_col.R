add_col <- function(dt,
                    col,
                    force_new = FALSE) {
    col_exists <- col %in% names(dt)
    if (col_exists) {
        if (force_new) {
            message(col, " already exists (replacing).")
            return(TRUE)
        } else {
            message(col, " already exists (keeping original).")
            return(FALSE)
        }
    } else {
        return(TRUE)
    }
}
