#' Compute per-pixel linear trend for a land use pool
#'
#' Sums the specified LUH2 layers for each year in the range and fits a
#' per-pixel linear trend, returning the slope (change per year).
#'
#' @param luh A `SpatRaster` returned by [luh2_load()].
#' @param pool_vars Character vector of LUH2 layer name prefixes to sum
#'   (e.g. `c("c3ann", "c4ann", "c3per", "c4per", "c3nfx")` for crop,
#'   `c("primf", "secdf")` for natural forest,
#'   `c("primn", "secdn")` for other).
#' @param year_start Integer. First year of the trend period (e.g. 1990).
#' @param year_end Integer. Last year of the trend period (e.g. 2015).
#'
#' @return A single-layer `SpatRaster` named `"pool_slope"`.
#' @author Abhijeet Mishra, Claude Code
#' @export
luh2_pool_trend <- function(luh, pool_vars, year_start, year_end) {
    luh_start <- as.integer(terra::time(luh)[1])

    suffix_start <- year_start - luh_start + 1
    suffix_end   <- year_end   - luh_start + 1
    suffixes     <- seq(suffix_start, suffix_end)

    pool_stack <- terra::rast(lapply(suffixes, function(s) {
        layers    <- paste0(pool_vars, "_", s)
        available <- layers[layers %in% names(luh)]
        if (length(available) == 0) stop("No layers found for suffix ", s,
                                         " (pool_vars: ", paste(pool_vars, collapse = ", "), ")")
        terra::app(luh[[available]], fun = "sum", na.rm = TRUE)
    }))

    time_idx <- seq_along(suffixes)

    pool_trend <- terra::app(pool_stack, function(x) {
        if (all(is.na(x))) return(NA_real_)
        stats::coef(stats::lm(x ~ time_idx))[2]
    })

    names(pool_trend) <- "pool_slope"
    pool_trend
}
