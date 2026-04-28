#' Merge planted forest layer into LUH2 shares
#'
#' Loads a planted forest raster (e.g. from SDPT), carves out area first from
#' `other`, then from `natfor`, and caps any residual planted share at the
#' remaining space after all LUH2 categories. Returns adjusted share rasters.
#'
#' @param luarea A `SpatRaster` of land use areas (km²) returned inside the
#'   list from [luh2_extract_year()].
#' @param cellarea A single-layer `SpatRaster` of cell areas in km² (also from
#'   [luh2_extract_year()]).
#' @param planted_tif Path to the planted forest GeoTIFF (e.g. `sdpt_global.tif`).
#' @param fstnf Flag if a pixel is forested or not.
#' @param cty_shp Path to the IMPACT regions shapefile. Must contain a `NEW_REGION`
#'   field used as the country identifier.
#' @return A `SpatRaster` of adjusted land use shares (one layer per pool,
#'   including `planted_forest_share`).
#' @author Abhijeet Mishra, Claude Code
#' @export
luh2_merge_planted <- function(luarea, cellarea, planted_tif, fstnf, cty_shp) {
    cty <- terra::vect(cty_shp)
    names(cty)[names(cty) %in% "NEW_REGION"] <- "cty"

    # Project to match LUH2 CRS exactly
    cty <- terra::project(cty, "OGC:CRS84")
    # touches = TRUE ensures border pixels are assigned to a country rather than
    # dropped as NA, preventing holes along country edges
    cty_rast <- terra::rasterize(cty, luarea[[1]], field = "cty", touches = TRUE)

    # Fir all natural land is in one prim and secd (f/n)
    pool_map <- list(
        crop   = c("c3ann", "c4ann", "c3per", "c4per", "c3nfx"),
        natfor = c("primf", "secdf", "primn", "secdn"),
        range  = c("range", "pastr"),
        urban  = "urban",
        other = c("primn", "secdf", "secdn")
    )

    pool_areas <- terra::rast(lapply(names(pool_map), function(pool) {
        layers <- pool_map[[pool]][pool_map[[pool]] %in% names(luarea)]
        if (length(layers) == 0) stop("No layers found for pool: ", pool)
        r <- terra::app(luarea[[layers]], fun = "sum", na.rm = TRUE)
        names(r) <- pool
        r
    }))

    pool_areas$natfor <- terra::ifel(cty_rast == "ZAF",
                                     pool_areas$natfor,
                                     pool_areas$natfor * fstnf)
    pool_areas$other <- terra::ifel(cty_rast == "ZAF",
                                    pool_areas$other,
                                    pool_areas$other * (1 - fstnf))

    plant_total <- terra::rast(planted_tif)

    all_layers <- c(pool_areas, plant_total)
    names(all_layers)[terra::nlyr(all_layers)] <- "planted_forest"

    shares <- all_layers / cellarea
    names(shares) <- paste0(names(all_layers), "_share")

    planted <- shares[["planted_forest_share"]]

    # Carve from other first
    pool <- "other"
    pooln <- paste0(pool, "_share")
    take_from_pool <- terra::lapp(
        c(planted, shares[[pooln]]),
        fun = function(a, b) pmin(a, b)
    )
    shares[[pooln]] <- shares[[pooln]] - take_from_pool
    planted_remaining <- planted - take_from_pool

    # Then carve from past
    pool <- "range"
    pooln <- paste0(pool, "_share")
    take_from_pool <- terra::lapp(
        c(planted_remaining, shares[[pooln]]),
        fun = function(a, b) pmin(a, b)
    )
    shares[[pooln]] <- shares[[pooln]] - take_from_pool

    # Cap residual planted at remaining space
    luh2_pools <- c("crop_share", "natfor_share", "range_share", "other_share", "urban_share")
    luh2_total <- terra::app(shares[[luh2_pools]], fun = "sum", na.rm = TRUE)

    shares[["planted_forest_share"]] <- terra::lapp(
        c(luh2_total, shares[["planted_forest_share"]]),
        fun = function(a, b) pmin(b, pmax(1 - a, 0))
    )

    shares
}
