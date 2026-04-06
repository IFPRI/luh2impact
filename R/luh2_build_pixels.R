#' Build pixel-level data frame with weights
#'
#' Rasterizes the IMPACT country shapefile to match the shares grid, merges
#' all raster data into a pixel-level data frame, computes land pool areas,
#' available land, and spatial weights for each pool.
#'
#' @param shares A `SpatRaster` of adjusted land use shares returned by
#'   [luh2_merge_planted()].
#' @param crop_trend A single-layer `SpatRaster` of crop slope returned by
#'   [luh2_pool_trend()].
#' @param natfor_trend A single-layer `SpatRaster` of natural forest slope
#'   returned by [luh2_pool_trend()].
#' @param other_trend A single-layer `SpatRaster` of other land slope
#'   returned by [luh2_pool_trend()].
#' @param cellarea A single-layer `SpatRaster` of cell areas in km².
#' @param icwtr A single-layer `SpatRaster` of ice/water fraction.
#' @param cty_shp Path to the IMPACT regions shapefile. Must contain a `NEW_REGION`
#'   field used as the country identifier.
#' @param landx0 A data frame of IMPACT land use results with columns
#'   `cty`, `fland`, `yrs`, `value`.
#'
#' @return A data frame with one row per pixel containing areas, shares,
#'   weights, and country assignments.
#' @author Abhijeet Mishra, Claude Code
#' @export
luh2_build_pixels <- function(shares, crop_trend, natfor_trend, other_trend,
                              cellarea, icwtr, cty_shp, landx0) {
    cty <- terra::vect(cty_shp)
    cty$cty <- cty$NEW_REGION

    # Project to match LUH2 CRS exactly
    cty <- terra::project(cty, "OGC:CRS84")
    # touches = TRUE ensures border pixels are assigned to a country rather than
    # dropped as NA, preventing holes along country edges
    cty_rast <- terra::rasterize(cty, shares[[1]], field = "cty", touches = TRUE)

    # Rename range to past for consistency with IMPACT
    names(shares)[names(shares) == "range_share"] <- "past_share"
    names(shares)[names(shares) == "planted_forest_share"] <- "plant_share"

    df <- as.data.frame(c(cty_rast, shares), xy = TRUE, na.rm = TRUE)

    # Normalize shares to sum to 1
    share_cols <- c("crop_share", "natfor_share", "past_share", "other_share", "urban_share", "plant_share")
    df$sharetotal <- rowSums(df[, share_cols])
    df[, share_cols] <- pmin(pmax(df[, share_cols], 0), 1)
    df[, share_cols] <- df[, share_cols] / df$sharetotal
    df$sharetotal    <- rowSums(df[, share_cols])

    # Drop SSD (not in IMPACT)
    df <- df[df$cty != "SSD", ]

    pts <- terra::vect(df, geom = c("x", "y"), crs = terra::crs(cty_rast))

    # Cell area and available land
    df$carea_km2 <- terra::extract(cellarea, pts)[, 2]
    df$carea_kha <- df$carea_km2 * 100 / 1000

    df$icwtr <- terra::extract(icwtr, pts)[, 2]
    df$icwtr[is.na(df$icwtr)] <- 0
    df$land_frac <- 1 - df$icwtr

    # Pool areas
    landpools <- c("crop", "natfor", "past", "other", "urban", "plant")
    pool_cols  <- paste0(landpools, "_share")
    area_cols  <- paste0(landpools, "_area")
    df[, area_cols] <- df[, pool_cols] * df$carea_kha

    df$avail_kha <- pmax(df$carea_kha * df$land_frac - df$urban_area, 0)

    # Pool trends
    df$crop_slope   <- terra::extract(crop_trend,   pts)[, 2]
    df$natfor_slope <- terra::extract(natfor_trend, pts)[, 2]
    df$other_slope  <- terra::extract(other_trend,  pts)[, 2]
    df$crop_slope[is.na(df$crop_slope)]     <- 0
    df$natfor_slope[is.na(df$natfor_slope)] <- 0
    df$other_slope[is.na(df$other_slope)]   <- 0

    # Planted forest proximity weight
    pts <- terra::vect(df, geom = c("x", "y"), crs = terra::crs(cty_rast))
    plant_rast  <- terra::rasterize(pts, cty_rast, field = "plant_area", fun = "sum")
    plant_prox  <- terra::focal(plant_rast, w = 3, fun = "mean", na.rm = TRUE)
    names(plant_prox) <- "plant_prox"
    df$plant_prox <- terra::extract(plant_prox, pts)[, 2]
    df$plant_prox[is.na(df$plant_prox)] <- 0

    # Pixel weights per country
    df <- df |>
        dplyr::group_by(cty) |>
        dplyr::mutate(
            w_crop_exp    = (crop_area   * pmax(crop_slope,    0)) / sum(crop_area   * pmax(crop_slope,    0), na.rm = TRUE),
            w_crop_con    = (crop_area   * pmax(-crop_slope,   0)) / sum(crop_area   * pmax(-crop_slope,   0), na.rm = TRUE),
            w_natfor_loss = (natfor_area * pmax(-natfor_slope, 0)) / sum(natfor_area * pmax(-natfor_slope, 0), na.rm = TRUE),
            w_natfor_gain = (natfor_area * pmax( natfor_slope, 0)) / sum(natfor_area * pmax( natfor_slope, 0), na.rm = TRUE),
            w_other_loss  = (other_area  * pmax(-other_slope,  0)) / sum(other_area  * pmax(-other_slope,  0), na.rm = TRUE),
            w_other_gain  = (other_area  * pmax( other_slope,  0)) / sum(other_area  * pmax( other_slope,  0), na.rm = TRUE),
            w_past        = past_area  / sum(past_area,  na.rm = TRUE),
            w_plant       = (plant_area * plant_prox) / sum(plant_area * plant_prox, na.rm = TRUE)
        ) |>
        dplyr::ungroup()

    # Clean up NaN weights
    w_cols <- c("w_crop_exp", "w_crop_con", "w_natfor_loss", "w_natfor_gain",
                "w_other_loss", "w_other_gain", "w_past", "w_plant")
    df[, w_cols] <- lapply(df[, w_cols], function(x) ifelse(is.nan(x) | is.na(x), 0, x))
    df[is.na(df)] <- 0

    df$pid <- paste0("p", seq_len(nrow(df)))
    df
}
