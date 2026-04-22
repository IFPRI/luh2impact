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
#'
#' @return A `SpatRaster` of adjusted land use shares (one layer per pool,
#'   including `planted_forest_share`).
#' @author Abhijeet Mishra, Claude Code
#' @export
luh2_merge_planted <- function(luarea, cellarea, planted_tif) {
    pool_map <- list(
        crop   = c("c3ann", "c4ann", "c3per", "c4per", "c3nfx"),
        natfor = c("primf", "primn", "secdf"),
        range  = c("range", "pastr"),
        other  = c("secdn"),
        urban  = "urban"
    )

    pool_areas <- terra::rast(lapply(names(pool_map), function(pool) {
        layers <- pool_map[[pool]][pool_map[[pool]] %in% names(luarea)]
        if (length(layers) == 0) stop("No layers found for pool: ", pool)
        r <- terra::app(luarea[[layers]], fun = "sum", na.rm = TRUE)
        names(r) <- pool
        r
    }))

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
