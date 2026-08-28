#' Compute IMPACT-country urban area from LUH2 for a target year
#'
#' Aggregates the LUH2 `urban` land-use layer to IMPACT country totals for a
#' single year, for use as a fixed baseline land-use pool (see
#' [luh2_export_elasticities_gdx()]) alongside cropland, pasture, forest, and
#' other land.
#'
#' @param states_nc Path to the LUH2 `states.nc` file.
#' @param static_nc Path to the LUH2 `staticData_quarterdeg.nc` file.
#' @param cty_shp Path to the IMPACT regions shapefile. Must contain a
#'   `NEW_REGION` field used as the country identifier.
#' @param year Integer. Target year. Default `2021`.
#'
#' @return A data frame with columns `cty` and `value` (urban area in 1000 ha).
#' @author Abhijeet Mishra, Claude Code
#' @export
luh2_urban_area_cty <- function(states_nc, static_nc, cty_shp, year = 2021) {
    luh       <- luh2_load(states_nc)
    extracted <- luh2_extract_year(luh, static_nc, year)
    urban_km2 <- extracted$luarea[["urban"]]

    cty <- terra::vect(cty_shp)
    names(cty)[names(cty) %in% "NEW_REGION"] <- "cty"
    cty <- terra::project(cty, "OGC:CRS84")
    cty_rast <- terra::rasterize(cty, urban_km2, field = "cty", touches = TRUE)

    df <- as.data.frame(c(cty_rast, urban_km2), na.rm = TRUE)
    names(df) <- c("cty", "urban_area_km2")

    df |>
        dplyr::group_by(cty) |>
        dplyr::summarise(value = sum(urban_area_km2, na.rm = TRUE) * 100 / 1000, .groups = "drop") |>
        as.data.frame()
}
