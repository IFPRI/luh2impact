#' Extract land use areas for a target year
#'
#' Subsets the LUH2 raster to a single year, multiplies shares by cell area,
#' and returns both the raw area raster and the `cellarea` layer for downstream use.
#'
#' @param luh A `SpatRaster` returned by [luh2_load()].
#' @param static_nc Path to the LUH2 `staticData_quarterdeg.nc` file.
#' @param year Integer. Target year (e.g. 2015).
#'
#' @return A named list:
#' \describe{
#'   \item{`luarea`}{`SpatRaster` of land use areas in km² per pool.}
#'   \item{`cellarea`}{Single-layer `SpatRaster` of cell areas in km².}
#'   \item{`icwtr`}{Single-layer `SpatRaster` of ice/water fraction.}
#' }
#' @author Abhijeet Mishra, Claude Code
#' @export
luh2_extract_year <- function(luh, static_nc, year) {
    luh_start <- as.integer(terra::time(luh)[1])
    suffix    <- year - luh_start + 1

    year_layers <- grep(paste0("_", suffix, "$"), names(luh), value = TRUE)
    if (length(year_layers) == 0) {
        stop("Year ", year, " (suffix _", suffix, ") not found in LUH2 raster.")
    }

    luh_yr <- luh[[year_layers]]
    names(luh_yr) <- gsub("_\\d+$", "", names(luh_yr))

    static <- terra::rast(static_nc)
    cellarea <- static[["carea"]]
    icwtr    <- static[["icwtr"]]
    ccode    <- static[["ccode"]]

    luarea <- luh_yr * cellarea
    names(luarea) <- names(luh_yr)

    list(
        luarea   = luarea,
        cellarea = cellarea,
        icwtr    = icwtr,
        ccode    = ccode
    )
}
