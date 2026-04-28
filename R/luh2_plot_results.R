#' Plot LUMEN land use results
#'
#' Reads GeoTIFFs written by [luh2_write_tifs()] and produces three plots:
#' \describe{
#'   \item{`changes`}{Faceted map of land use change vs 2021 for years 2035 and 2050,
#'     one panel per pool. Diverging red-blue palette, symmetric scale.}
#'   \item{`shares`}{Faceted map of land use shares for 2021 and 2050,
#'     faceted by pool and year.}
#'   \item{`dominant`}{Map of dominant land use pool per pixel for 2021 and 2050,
#'     coloured by pool with alpha proportional to share.}
#' }
#'
#' @param output_dir Path to directory containing `lu_<pool>.tif` files and the
#'   solution GDX. Country name is derived via `basename(output_dir)`.
#' @param cty_shp Path to the IMPACT regions shapefile.
#' @param save_png Logical. If `TRUE` (default), saves each plot as a PNG file
#'   in `output_dir`.
#'
#' @return A named list of `ggplot` objects: `changes`, `shares`, `dominant`.
#' @author Abhijeet Mishra, Claude Code
#' @export
luh2_plot_results <- function(output_dir, cty_shp, save_png = TRUE) {
    cty_name <- basename(output_dir)
    pools <- c("crop", "natfor", "past", "other", "plant")

    pool_colors <- c(
        crop   = "#f5c842",
        natfor = "#1a7a2e",
        past   = "#a8d98f",
        other  = "#c2a06e",
        plant  = "#2d6e9e"
    )

    # Load country shapefile and filter to current country
    cty_vect <- terra::vect(cty_shp)
    ctyx     <- cty_vect[cty_vect$NEW_REGION %in% cty_name, ]

    # Clip CONUS for USA to avoid Alaska distortion
    if (cty_name == "USA") {
        ctyx    <- terra::crop(ctyx, terra::ext(-125, -66, 24, 50))
    }

    # Load all pool rasters
    lu_all <- stats::setNames(
        lapply(pools, function(p) terra::rast(file.path(output_dir, paste0("lu_", p, ".tif")))),
        pools
    )

    yrs_all <- names(lu_all[[1]])
    ext_rast  <- terra::ext(lu_all[[1]][[1]])
    ctyx      <- terra::crop(ctyx, ext_rast)

    # Load available area for share calculation
    m_out    <- gamstransfer::Container$new(file.path(output_dir, "solution_lu.gdx"))
    avail_df <- m_out["avail_kha"]$records
    names(avail_df) <- c("pid", "avail_kha")

    pid_coords <- m_out["pid_coords"]$records |>
        tidyr::pivot_wider(names_from = "coord", values_from = "value")

    avail_df <- avail_df |>
        dplyr::left_join(pid_coords[, c("pid", "x", "y")], by = "pid")

    # Precalc share stuff

    share_yrs <- intersect(c("2021", "2050"), yrs_all)

    share_df <- dplyr::bind_rows(lapply(share_yrs, function(yr) {
        pool_stack <- terra::rast(lapply(pools, function(p) {
            r <- lu_all[[p]][[yr]]
            names(r) <- p
            r
        }))

        as.data.frame(pool_stack, xy = TRUE) |>
            tidyr::pivot_longer(-c(x, y), names_to = "pool", values_to = "area") |>
            dplyr::left_join(avail_df[, c("x", "y", "avail_kha")], by = c("x", "y")) |>
            dplyr::mutate(share = area / avail_kha, year = yr)
    }))

    # ----------------------------------------------------------------
    # Plot 2: Shares faceted by pool and year (2021 and 2050)
    # ----------------------------------------------------------------

    ext_run <- c(range(share_df$x, na.rm = TRUE), range(share_df$y, na.rm = TRUE))

    p_shares <- ggplot2::ggplot() +
        tidyterra::geom_spatvector(data = ctyx, fill = NA, color = "black") +
        ggplot2::geom_tile(
            data = share_df |> dplyr::filter(share >= 0.005, !is.na(share), !is.infinite(share)),
            ggplot2::aes(x = x, y = y, fill = pool, alpha = share)
        ) +
        ggplot2::scale_fill_manual(values = pool_colors, name = "pool") +
        ggplot2::scale_alpha_continuous(
            name   = "share",
            range  = c(0, 1),
            limits = c(0.005, 1),
            breaks = seq(0.1, 1, 0.15)
        ) +
        ggplot2::coord_sf(
            xlim = c(ext_run[1] - 1, ext_run[2] + 1),
            ylim = c(ext_run[3] - 1, ext_run[4] + 1)
        ) +
        ggplot2::facet_grid(year ~ pool) +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "Land use share")

    # ----------------------------------------------------------------
    # Plot 3: Area (kha) faceted by pool, years 2021 and 2050
    # ----------------------------------------------------------------
    area_filtered <- share_df |>
        dplyr::filter(!is.na(area), area > 0, share >= 0.005) |>
        dplyr::pull(area)

    area_breaks <- seq(0, ceiling(max(area_filtered)), by = ceiling(max(area_filtered)) / 6)


    p_area <- ggplot2::ggplot() +
        tidyterra::geom_spatvector(data = ctyx, fill = NA, color = "black") +
        ggplot2::geom_tile(
            data = share_df |> dplyr::filter(!is.na(area), area > 0, share >= 0.005),
            ggplot2::aes(x = x, y = y, fill = share, alpha = area)
        ) +
        ggplot2::scale_fill_distiller(
            palette   = "YlOrRd",
            direction = 1,
            name      = "share",
            na.value  = "transparent"
        ) +
        ggplot2::scale_alpha_continuous(
            name   = "area (kha)",
            range  = c(0.2, 1),
            trans  = "log",
            breaks = area_breaks
        ) +
        ggplot2::coord_sf(
            xlim = c(ext_run[1] - 1, ext_run[2] + 1),
            ylim = c(ext_run[3] - 1, ext_run[4] + 1)
        ) +
        ggplot2::facet_grid(year ~ pool) +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "Land use area (kha)") +
        ggplot2::scale_alpha_binned(
            name   = "area (kha)",
            range  = c(0, 1),
            breaks = area_breaks,
            limits = c(0, ceiling(max(area_filtered)))
        ) +
    ggplot2::guides(
        alpha = ggplot2::guide_bins(override.aes = list(fill = "black"))
    )

    # ----------------------------------------------------------------
    # Plot 3a: Area (kha) faceted by pool, years 2050 vs 2021
    # ----------------------------------------------------------------
    area_change <- share_df |>
        filter(year %in% c(2021, 2050), share >= 0.005) |>
        group_by(x, y, pool) |>
        mutate(
            area_change = if (all(c(2021, 2050) %in% year)) area - area[year == 2021] else NA_real_
        ) |>
        filter(!is.na(area_change), area_change != 0)

    area_filtered <- area_change |> pull(area_change)

    area_breaks <- seq(-ceiling(max(area_filtered)),
                       ceiling(max(area_filtered)),
                       by = ceiling(max(area_filtered)) / 3)


    p_area_change <- ggplot2::ggplot() +
        ggplot2::geom_tile(
            data = area_change |> dplyr::filter(!is.na(area), area > 0, share >= 0.005),
            ggplot2::aes(x = x, y = y, fill = area_change)
        ) +
        tidyterra::geom_spatvector(data = ctyx, fill = NA, color = "black") +
        ggplot2::scale_fill_gradient2(
            low      = "red",
            mid      = "white",
            high     = "green",
            midpoint = 0,
            name     = "area change (kha)",
            na.value = "transparent"
        ) +
        ggplot2::coord_sf(
            xlim = c(ext_run[1] - 1, ext_run[2] + 1),
            ylim = c(ext_run[3] - 1, ext_run[4] + 1)
        ) +
        ggplot2::facet_wrap(. ~ pool) +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "Land use area change 2021-2050 (Kha)")

    # ----------------------------------------------------------------
    # Plot 4: Dominant land use per pixel
    # ----------------------------------------------------------------
    share_breaks <- seq(0, 1, by = 0.2)

    plot_df_dominant <- share_df |>
        dplyr::filter(share >= 0.05, !is.na(share), !is.infinite(share)) |>
        dplyr::group_by(x, y, year) |>
        dplyr::slice_max(share, n = 1, with_ties = FALSE) |>
        dplyr::ungroup()

    p_dominant <- ggplot2::ggplot() +
        tidyterra::geom_spatvector(data = ctyx, fill = NA, color = "black") +
        ggplot2::geom_tile(
            data = plot_df_dominant,
            ggplot2::aes(x = x, y = y, fill = pool, alpha = share)
        ) +
        ggplot2::scale_fill_manual(values = pool_colors, name = "pool") +
        ggplot2::scale_alpha_binned(
            name   = "share",
            range  = c(0, 1),
            breaks = share_breaks,
            limits = c(0, 1)
        ) +
        ggplot2::coord_sf(
            xlim = c(ext_run[1] - 1, ext_run[2] + 1),
            ylim = c(ext_run[3] - 1, ext_run[4] + 1)
        ) +
        ggplot2::facet_wrap(year ~ .) +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "Dominant land use share")

    # ----------------------------------------------------------------
    # Save PNGs
    # ----------------------------------------------------------------
    if (save_png) {
        ggplot2::ggsave(
            plot     = p_shares,
            filename = file.path(output_dir, paste0(cty_name, "-shares.png")),
            width    = 16, height = 7, bg = "white"
        )
        ggplot2::ggsave(
            plot     = p_dominant,
            filename = file.path(output_dir, paste0(cty_name, "-dominant.png")),
            width    = 16, height = 8, bg = "white"
        )
        ggplot2::ggsave(
            plot     = p_area,
            filename = file.path(output_dir, paste0(cty_name, "-area.png")),
            width    = 16, height = 7, bg = "white"
        )
        ggplot2::ggsave(
            plot     = p_area_change,
            filename = file.path(output_dir, paste0(cty_name, "-areachange.png")),
            width    = 12, height = 8, bg = "white"
        )
        message("Plots saved to ", output_dir)
    }

    invisible(list(
        shares   = p_shares,
        area     = p_area,
        dominant = p_dominant
    ))
}
