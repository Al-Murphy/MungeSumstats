formatted_example <- function(path=system.file("extdata", "eduAttainOkbay.txt",
                                               package = "MungeSumstats")){
    sumstats_dt <- suppressMessages(
        read_sumstats(path = path)
    )
    sumstats_dt <-
       standardise_sumstats_column_headers_crossplatform(
            sumstats_dt = sumstats_dt)$sumstats_dt
    sumstats_dt <- sort_coords(sumstats_dt = sumstats_dt)
    return(sumstats_dt)
}