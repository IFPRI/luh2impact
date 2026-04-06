#' Export pixel data and IMPACT totals to GDX
#'
#' Writes all sets and parameters required by the LUMEN GAMS model into a
#' single GDX file using `gamstransfer`.
#'
#' @param df A pixel-level data frame returned by [luh2_build_pixels()].
#' @param landx0 A data frame of IMPACT land use results with columns
#'   `cty`, `fland`, `yrs`, `value`.
#' @param output_gdx Path to write the output GDX file.
#'
#' @return Invisibly returns `output_gdx`. Called for its side effect.
#' @author Abhijeet Mishra, Claude Code
#' @export
luh2_export_gdx <- function(df, landx0, output_gdx) {
    m <- gamstransfer::Container$new()

    # --- Sets ---
    pid_set   <- m$addSet("pid",   records = df$pid,                                         description = "pixel IDs")
    cty_set   <- m$addSet("cty",   records = levels(factor(df$cty)),                         description = "country IDs")
    fland_set <- m$addSet("fland", records = c("crop", "natfor", "past", "other", "plant"),  description = "land pools")
    yrs_set   <- m$addSet("yrs",   records = as.numeric(as.character(unique(landx0$yrs))),   description = "years")
    coord_set <- m$addSet("coord", records = c("x", "y"),                                    description = "coordinates")

    m$addSet("pid_cty", c(pid_set, cty_set),
             records = df[, c("pid", "cty")],
             description = "pixel to country mapping")

    # --- Parameters ---

    # luh2_base: baseline areas by pool
    luh2_base_df <- df |>
        dplyr::select(pid,
                      crop   = crop_area,
                      natfor = natfor_area,
                      past   = past_area,
                      other  = other_area,
                      plant  = plant_area) |>
        tidyr::pivot_longer(-pid, names_to = "fland", values_to = "value")
    luh2_base_df$value[is.na(luh2_base_df$value)] <- 0

    m$addParameter("luh2_base", c(pid_set, fland_set),
                   records = luh2_base_df,
                   description = "LUH2 baseline pixel areas by pool (1000 ha)")

    m$addParameter("urban_area", pid_set,
                   records = df[, c("pid", "urban_area")],
                   description = "fixed urban area per pixel (1000 ha)")

    m$addParameter("carea", pid_set,
                   records = df[, c("pid", "carea_kha")],
                   description = "cell area per pixel (1000 ha)")

    m$addParameter("avail_kha", pid_set,
                   records = df[, c("pid", "avail_kha")],
                   description = "available land per pixel excl. urban and ice/water (1000 ha)")

    # IMPACT control totals
    impact_total_df <- landx0 |>
        dplyr::select(cty, fland, yrs, value) |>
        dplyr::mutate(yrs = as.character(yrs))

    m$addParameter("impact_total", c(cty_set, fland_set, yrs_set),
                   records = impact_total_df,
                   description = "IMPACT control totals by country pool year (1000 ha)")

    # Pixel coordinates
    coords_df <- df |>
        dplyr::select(pid, x, y) |>
        tidyr::pivot_longer(-pid, names_to = "coord", values_to = "value")

    m$addParameter("pid_coords", c(pid_set, coord_set),
                   records = coords_df,
                   description = "pixel centroid coordinates")

    # Weights
    m$addParameter("w_crop_exp", pid_set,
                   records = df[, c("pid", "w_crop_exp")],
                   description = "crop expansion weight per pixel")

    m$addParameter("w_crop_con", pid_set,
                   records = df[, c("pid", "w_crop_con")],
                   description = "crop contraction weight per pixel")

    m$addParameter("w_natfor_loss", pid_set,
                   records = df[, c("pid", "w_natfor_loss")],
                   description = "natural forest loss weight per pixel")

    m$addParameter("w_natfor_gain", pid_set,
                   records = df[, c("pid", "w_natfor_gain")],
                   description = "natural forest gain weight per pixel")

    m$addParameter("w_other_loss", pid_set,
                   records = df[, c("pid", "w_other_loss")],
                   description = "other land loss weight per pixel")

    m$addParameter("w_other_gain", pid_set,
                   records = df[, c("pid", "w_other_gain")],
                   description = "other land gain weight per pixel")

    m$addParameter("w_past", pid_set,
                   records = df[, c("pid", "w_past")],
                   description = "pasture weight per pixel (reserved)")

    m$addParameter("w_plant", pid_set,
                   records = df[, c("pid", "w_plant")],
                   description = "planted forest proximity weight per pixel")

    m$write(output_gdx)
    invisible(output_gdx)
}
