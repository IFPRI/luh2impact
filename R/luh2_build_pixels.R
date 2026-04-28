# nolint start: object_usage_linter
#' Build pixel-level data frame with weights
#'
#' Rasterizes the IMPACT country shapefile to match the shares grid, merges
#' all raster data into a pixel-level data frame, computes land pool areas,
#' available land, and spatial weights for each pool.
#'
#' @param luh A `SpatRaster` of historical land use shares returned by
#'  [luh2_load()].
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
#' @importFrom dplyr summarise left_join mutate select group_by
#' @importFrom tidyr pivot_wider
#' @importFrom magrittr %>%
#' @export
luh2_build_pixels <- function(luh, shares, crop_trend, natfor_trend, other_trend,
                              cellarea, icwtr, cty_shp, landx0) {
    cty <- terra::vect(cty_shp)
    names(cty)[names(cty) %in% "NEW_REGION"] <- "cty"

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
    share_cols <- paste0(c("crop", "natfor", "past", "other", "urban", "plant"), "_share")
    df$sharetotal <- rowSums(df[, share_cols], na.rm = TRUE)

    # Drop SSD (not in IMPACT)
    df <- df[df$cty != "SSD", ]

    pts <- terra::vect(df, geom = c("x", "y"), crs = terra::crs(cty_rast))

    # Cell area and available land
    df$carea_km2 <- terra::extract(cellarea, pts)[, 2]
    df$carea_kha <- df$carea_km2 * 100 / 1000

    df$icwtr <- terra::extract(icwtr, pts)[, 2]
    df$land_frac <- 1 - df$icwtr

    # Pool areas
    landpools <- c("crop", "natfor", "past", "other", "urban", "plant")
    pool_cols  <- paste0(landpools, "_share")
    area_cols  <- paste0(landpools, "_area")
    df[, area_cols] <- df[, pool_cols] * df$carea_kha

    df$avail_kha <- pmax(df$carea_kha * df$land_frac - df$urban_area, 0)

    ##########################

    # IPF adjustment (itrerative proportional fitting)
    n_iter <- 20
    tol <- 1
    for (i in 1:n_iter) {

        cty_totals <- df %>%
            group_by(cty) %>%
            summarise(
                crop_tot   = sum(crop_area),
                natfor_tot = sum(natfor_area),
                past_tot   = sum(past_area),
                other_tot  = sum(other_area),
                plant_tot  = sum(plant_area)
            )

        deltas <- cty_totals %>%
            left_join(landx0 %>%
                          filter(yrs == 2021) %>%
                          mutate(yrs = NULL) %>%
                          pivot_wider(names_from = fland, values_from = value) %>%
                          mutate(cty = as.character(cty)), by = "cty") %>%
            mutate(
                d_crop   = crop   - crop_tot,
                d_natfor = natfor - natfor_tot,
                d_past   = past   - past_tot,
                d_other  = other  - other_tot,
                d_plant  = plant  - plant_tot
            ) %>%
            select(cty, starts_with("d_"))

        df <- df %>%
            left_join(deltas, by = "cty") %>%
            group_by(cty) %>%
            mutate(
                crop_area   = pmin(crop_area   + (crop_area   / sum(crop_area))   * d_crop,   avail_kha),
                natfor_area = pmin(natfor_area + (natfor_area / sum(natfor_area)) * d_natfor, avail_kha),
                past_area   = pmin(past_area   + (past_area   / sum(past_area))   * d_past,   avail_kha),
                other_area  = pmin(other_area  + (other_area  / sum(other_area))  * d_other,  avail_kha),
                plant_area  = pmin(plant_area  + (plant_area  / sum(plant_area))  * d_plant,  avail_kha)
            ) %>%
            ungroup() %>%
            select(-starts_with("d_"))

        df <- df %>%
            mutate(
                pool_total  = crop_area + natfor_area + past_area + other_area + plant_area,
                scale_down  = ifelse(pool_total > avail_kha, avail_kha / pool_total, 1),
                crop_area   = crop_area   * scale_down,
                natfor_area = natfor_area * scale_down,
                past_area   = past_area   * scale_down,
                other_area  = other_area  * scale_down,
                plant_area  = plant_area  * scale_down
            ) %>%
            select(-pool_total, -scale_down)

        # convergence check
        max_diff <- df %>%
            group_by(cty) %>%
            summarise(across(ends_with("_area"), sum, .names = "{.col}_tot")) %>%
            left_join(landx0 %>%
                          pivot_wider(names_from = fland, values_from = value) %>%
                          mutate(cty = as.character(cty)), by = "cty") %>%
            mutate(diff = abs(crop_area_tot - crop) + abs(natfor_area_tot - natfor) +
                       abs(past_area_tot - past)  + abs(other_area_tot - other) +
                       abs(plant_area_tot - plant)) %>%
            pull(diff) %>% max(na.rm = TRUE)

        cat("Iteration", i, "- max diff:", max_diff, "\n")
        if(max_diff < tol) break
    }

    # recalc shares ----

    df <- df %>%
        mutate(
            crop_share   = crop_area   / carea_kha,
            natfor_share = natfor_area / carea_kha,
            past_share   = past_area   / carea_kha,
            other_share  = other_area  / carea_kha,
            plant_share  = plant_area  / carea_kha,
            urban_share  = urban_area  / carea_kha,
            sharetotal   = crop_share + natfor_share + past_share +
                other_share + plant_share + urban_share
        )

    ##########################

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

    # Generate suitability numbers

    suit <- luh2_suitability(luh = luh)

    pts <- terra::vect(df, geom = c("x", "y"), crs = terra::crs(cty_rast))

    df$suit_crop   <- terra::extract(suit[["suit_crop"]],   pts)[, 2]
    df$suit_natfor <- terra::extract(suit[["suit_natfor"]], pts)[, 2]
    df$suit_past   <- terra::extract(suit[["suit_past"]],   pts)[, 2]
    df$suit_other  <- terra::extract(suit[["suit_other"]],  pts)[, 2]

    df[, c("suit_crop","suit_natfor","suit_past","suit_other")] <-
        lapply(df[, c("suit_crop","suit_natfor","suit_past","suit_other")],
               function(x) ifelse(is.na(x), 0, x))

    # Pixel weights per country
    df <- df |>
        dplyr::group_by(cty) |>
        dplyr::mutate(
            w_crop_exp    = (crop_area   * pmax(crop_slope,    0) * suit_crop)   / sum(crop_area   * pmax(crop_slope,    0) * suit_crop,   na.rm = TRUE),
            w_crop_con    = (crop_area   * pmax(-crop_slope,   0) * suit_crop)   / sum(crop_area   * pmax(-crop_slope,   0) * suit_crop,   na.rm = TRUE),
            w_natfor_loss = (natfor_area * pmax(-natfor_slope, 0) * suit_natfor) / sum(natfor_area * pmax(-natfor_slope, 0) * suit_natfor, na.rm = TRUE),
            w_natfor_gain = (natfor_area * pmax(natfor_slope,  0) * suit_natfor) / sum(natfor_area * pmax(natfor_slope,  0) * suit_natfor, na.rm = TRUE),
            w_other_loss  = (other_area  * pmax(-other_slope,  0) * suit_other)  / sum(other_area  * pmax(-other_slope,  0) * suit_other,  na.rm = TRUE),
            w_other_gain  = (other_area  * pmax(other_slope,   0) * suit_other)  / sum(other_area  * pmax(other_slope,   0) * suit_other,  na.rm = TRUE),
            w_past        = (past_area   * suit_past)                             / sum(past_area   * suit_past,                           na.rm = TRUE),
            w_plant       = (plant_area  * plant_prox)                            / sum(plant_area  * plant_prox,                          na.rm = TRUE)
        ) |>
        dplyr::ungroup()

    # Clean up NaN weights
    w_cols <- c("w_crop_exp", "w_crop_con", "w_natfor_loss", "w_natfor_gain",
                "w_other_loss", "w_other_gain", "w_past", "w_plant")
    df[, w_cols] <- lapply(df[, w_cols], function(x) ifelse(is.nan(x) | is.na(x), 0, x))
    df[is.na(df)] <- 0

    df$pid <- paste0("p", seq_len(nrow(df)))
    df
    # nolint end
}
