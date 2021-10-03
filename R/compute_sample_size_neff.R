#' Compute Neff/N
#'
#' @inheritParams compute_sample_size
#' @return No return
#' @keywords internal
compute_sample_size_neff <- function(sumstats_dt,
                                     method,
                                     force_new = FALSE,
                                     append_method_name = FALSE) {

    # Avoid confusing Biocheck
    Neff <- N <- N_ldsc <- P_ldsc <- N_CAS <- N_CON <- NULL

    if (is.character(method) &
        all(c("N_CAS", "N_CON") %in% colnames(sumstats_dt))) {

        #### Standardize method ####
        method <- tolower(method)[1]
        # get N col name
        N_col <- "Neff"
        if ((method %in% c("ldsc", TRUE))) {
            if (add_col(sumstats_dt, N_col, force_new)) {
                if (append_method_name) {
                      N_col <- "Neff_ldsc"
                  }
                #### Compute sample size: ldsc ####
                message(
                    "Computing effective sample size using the LDSC method:\n ",
                    "Neff = (N_CAS+N_CON) * (N_CAS/(N_CAS+N_CON)) /",
                    " mean((N_CAS/(N_CAS+N_CON))",
                    "[(N_CAS+N_CON)==max(N_CAS+N_CON)]))"
                )
                ### Create tmp cols (formula gets too confusing otherwise)
                #### Original python code in LDSC (converted below to R code)
                # N = dat.N_CAS + dat.N_CON
                # P = dat.N_CAS / N
                # dat['N'] = N * P / P[N == N.max()].mean()

                ## N_ldsc = Sum of cases vs. controls
                sumstats_dt[, N_ldsc := N_CAS + N_CON]
                ## P_ldsc = Proportion of cases vs. controls
                sumstats_dt[, P_ldsc := N_CAS / N_ldsc]

                sumstats_dt[, eval(N_col) :=
                    as.integer(N_ldsc * P_ldsc / 
                                   mean(P_ldsc[N_ldsc == max(
                                       N_ldsc, na.rm = TRUE)],
                                        na.rm=TRUE ))]
                ### Drop tmp cols
                sumstats_dt[, c("N_ldsc", "P_ldsc") := NULL]
                ### All in one go *but super confusing a liable to typos!
                # sumstats_dt[,Neff:=as.integer((N_CAS+N_CON) *
                # (N_CAS/(N_CAS+N_CON)) / 
                # mean((N_CAS/(N_CAS+N_CON))[(N_CAS+N_CON)
                # =max(N_CAS+N_CON)]))]
            }
        } else if (method == "metal") {
            if (add_col(sumstats_dt, N_col, force_new)) {
                if (append_method_name) {
                      N_col <- "Neff_metal"
                  }
                #### Compute sample size: metal ####
                message(
                    "Computing effective sample size",
                    " using the METAL method:\n ",
                    "Neff = 4 / (1/N_CAS + 1/N_CON)"
                )
                sumstats_dt[, eval(N_col) :=
                    as.integer(4L / (1L / N_CAS + 1L / N_CON))]
            }
        } else if (method == "giant") {
            if (add_col(sumstats_dt, N_col, force_new)) {
                if (append_method_name) {
                      N_col <- "Neff_giant"
                  }
                #### Compute sample size: giant ####
                message(
                    "Computing effective sample size",
                    " using the GIANT method:\n ",
                    "Neff = 2 / (1/N_CAS + 1/N_CON)"
                )
                sumstats_dt[, eval(N_col) :=
                    as.integer(2L / (1L / N_CAS + 1L / N_CON))]
            }
        } else if (method == "sum") {
            if (add_col(sumstats_dt, "N", force_new)) {
                #### Compute sample size: giant ####
                message(
                    "Computing sample size using the sum method:\n ",
                    "N = N_CAS + N_CON"
                )
                sumstats_dt[, N := as.integer(N_CAS + N_CON)]
            }
        } else {
            all_methods <- c("ldsc", "giant", "metal", "sum")
            msg <- paste0(
                "method must be one of:\n",
                paste0(" - ", all_methods, collapse = "\n")
            )
            message(msg)
        }
    } else {
        msg_err <- paste0(
            "WARNING: Neff column could not be calculated as the",
            " columns N_CAS & N_CON were not found in the datset"
        )
        message(msg_err)
    }
}
