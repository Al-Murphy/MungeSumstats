#' Convert sumstats to desired object type
#'
#' @inheritParams format_sumstats
#' @param return_format Object type to convert to 
#' @keywords internal 
convert_sumstats <- function(sumstats_dt,
                             return_format=c("data.table","vranges","granges")){
    return_format <- return_format[1]
    if(tolower(return_format) %in% c("vr","vranges")){ 
        out <- to_VRanges(sumstats_dt = sumstats_dt)
        return(out)
    }else if(tolower(return_format) %in% c("gr","granges","genomicranges")){
        out <- to_GRanges(sumstats_dt = sumstats_dt)
        return(out)
    } else {
        return(sumstats_dt)
    }
}
