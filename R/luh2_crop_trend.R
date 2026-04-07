#' Compute per-pixel crop expansion trend
#'
#' Fits a linear trend over time for the sum of crop land use shares across
#' five crop types (`c3ann`, `c4ann`, `c3per`, `c4per`, `c3nfx`), returning
#' the slope (change per year) for each pixel.
#'
#' @param luh A `SpatRaster` returned by [luh2_load()].
#' @param year_start Integer. First year of the trend period (e.g. 1990).
#' @param year_end Integer. Last year of the trend period (e.g. 2015).
#'
#' @return A single-layer `SpatRaster` named `"crop_slope"`.
#' @author Abhijeet Mishra, Claude Code
#' @export
luh2_crop_trend <- function(luh, year_start, year_end) {
    crop_vars <- c("c3ann", "c4ann", "c3per", "c4per", "c3nfx")

    luh_start <- as.integer(terra::time(luh)[1])

    suffix_start <- year_start - luh_start + 1
    suffix_end   <- year_end   - luh_start + 1
    suffixes     <- seq(suffix_start, suffix_end)

    crop_stack <- terra::rast(lapply(suffixes, function(s) {
        layers    <- paste0(crop_vars, "_", s)
        available <- layers[layers %in% names(luh)]
        if (length(available) == 0) stop("No crop layers found for suffix ", s)
        terra::app(luh[[available]], fun = "sum", na.rm = TRUE)
    }))

    time_idx <- seq_along(suffixes) # nolint: object_usage_linter

    crop_trend <- terra::app(crop_stack, function(x) {
        if (all(is.na(x))) return(NA_real_)
        stats::coef(stats::lm(x ~ time_idx))[2]
    })

    names(crop_trend) <- "crop_slope"
    crop_trend
}
