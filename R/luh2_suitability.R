#' Calculate land use suitability index from LUH2 historical data
#'
#' Computes a suitability index for each land pool by averaging historical
#' land use shares over a specified time window. Pixels with consistently
#' high shares of a given land type over centuries are considered more
#' suitable for that land type in future projections.
#'
#' @param luh SpatRaster from [luh2_load()].
#' @param year_start Integer. Start year for suitability window (e.g. 1700).
#' @param year_end Integer. End year for suitability window (e.g. 2015).
#' @param subsample Integer. Step size for subsampling years to speed up
#'   computation. Default is 10 (every 10th year).
#'
#' @return SpatRaster with 4 layers: \code{suit_crop}, \code{suit_natfor},
#'   \code{suit_past}, \code{suit_other}. Values range from 0 to 1,
#'   representing the mean share of each land pool over the specified window.
#'
#' @details
#' Land pool definitions follow LUH2 variable groupings:
#' \itemize{
#'   \item \code{crop}: c3ann, c4ann, c3per, c4per, c3nfx
#'   \item \code{natfor}: primf, secdf
#'   \item \code{past}: pastr, range
#'   \item \code{other}: primn, secdn
#' }
#' Within each year, variables are summed across pool members before
#' averaging across years.
#'
#' @examples
#' \dontrun{
#' luh  <- luh2_load("states.nc")
#' suit <- luh2_suitability(luh, year_start = 1700, year_end = 2015)
#' plot(suit)
#' }
#'
#' @export
luh2_suitability <- function(luh, year_start = 1715, year_end = 2015, subsample = 10) {

    pools <- list(
        crop   = c("c3ann", "c4ann", "c3per", "c4per", "c3nfx"),
        natfor = c("primf", "secdf"),
        past   = c("pastr", "range"),
        other  = c("primn", "secdn")
    )

    # LUH2 layer index = year - 849
    idx_subsample <- seq(year_start - 849, year_end - 849, by = subsample)

    suit <- terra::rast(lapply(names(pools), function(pool) {

        vars <- pools[[pool]]

        pool_layers <- terra::rast(lapply(idx_subsample, function(i) {
            pattern <- paste0("^(", paste(vars, collapse = "|"), ")_", i, "$")
            idx <- grep(pattern, names(luh))
            terra::app(luh[[idx]], sum, na.rm = TRUE)
        }))

        r <- terra::mean(pool_layers, na.rm = TRUE)
        names(r) <- paste0("suit_", pool)
        message("Done: ", pool)
        r
    }))

    suit
}
