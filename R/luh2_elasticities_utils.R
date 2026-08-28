# Internal helpers shared by luh2_load_faoland(), luh2_compute_gdppc(),
# and luh2_fit_elasticities(). Not exported.

.luh2_read_ctyset <- function(sets_xlsx) {
    ctyset <- openxlsx2::read_xlsx(file = sets_xlsx, sheet = "Regions",
                                   start_col = 5, start_row = 4)[, 1:2]
    ctyset <- ctyset[!is.na(ctyset$FCTY), ]
    names(ctyset) <- tolower(names(ctyset))
    ctyset$fcty <- as.numeric(ctyset$fcty)
    ctyset
}

# Linear interpolation of interior NAs, leaving leading/trailing NAs
# untouched (mirrors zoo::na.approx(x, na.rm = FALSE) without the dependency).
.luh2_approx_na <- function(x) {
    idx <- seq_along(x)
    ok  <- !is.na(x)
    if (sum(ok) < 2) return(x)
    x[!ok] <- stats::approx(idx[ok], x[ok], xout = idx[!ok])$y
    x
}

.luh2_safe_coef <- function(coef_tbl, term, stat, default = NA_real_) {
    if (term %in% rownames(coef_tbl)) coef_tbl[term, stat] else default
}
