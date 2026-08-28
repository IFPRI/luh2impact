#' Aggregate a pre-processed WDPA raster to IMPACT country totals
#'
#' Aggregates an already-processed protected area raster (e.g. the
#' pixel-level output of [wdpa_process()], km^2 per pixel) to IMPACT country
#' totals. Use this instead of [wdpa_process()] when a PA raster has already
#' been built, to skip the (slow) geodatabase intersection step.
#'
#' @param pa_tif Path to a protected area raster (km^2 per pixel) on the
#'   LUH2 grid, e.g. as produced by `wdpa_process(..., output = "pixel")`.
#' @param cty_shp Path to the IMPACT regions shapefile. Must contain a
#'   `NEW_REGION` field used as the country identifier.
#'
#' @return A data frame with columns `cty` and `value` (protected area in
#'   1000 ha).
#' @author Abhijeet Mishra, Claude Code
#' @export
luh2_pa_area_cty <- function(pa_tif, cty_shp) {
    pa_rast <- terra::rast(pa_tif)

    cty <- terra::vect(cty_shp)
    names(cty)[names(cty) %in% "NEW_REGION"] <- "cty"
    cty <- terra::project(cty, "OGC:CRS84")
    cty_rast <- terra::rasterize(cty, pa_rast, field = "cty", touches = TRUE)

    df <- as.data.frame(c(cty_rast, pa_rast), na.rm = TRUE)
    names(df) <- c("cty", "pa_area_km2")

    df |>
        dplyr::group_by(cty) |>
        dplyr::summarise(value = sum(pa_area_km2, na.rm = TRUE) * 100 / 1000, .groups = "drop") |>
        as.data.frame()
}
