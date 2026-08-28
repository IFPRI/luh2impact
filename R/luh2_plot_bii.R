#' Plot BII results from LUMEN output
#'
#' Reads `bii_pixel` from `solution_lu.gdx` and produces a three-panel figure:
#' \describe{
#'   \item{Left}{BII score in 2021 (baseline), sequential green palette.}
#'   \item{Centre}{BII score in 2050, same scale as 2021.}
#'   \item{Right}{Change in BII 2021-2050, diverging red-blue palette centred at 0.}
#' }
#'
#' @param output_dir Path to directory containing `solution_lu.gdx`.
#'   Country name is derived via `basename(output_dir)`.
#' @param cty_shp Path to the IMPACT regions shapefile.
#' @param save_png Logical. If `TRUE` (default), saves the combined figure as a
#'   PNG in `output_dir`.
#'
#' @return A `ggplot` object (the combined `ggarrange` figure), invisibly.
#' @author Abhijeet Mishra, Claude Code
#' @export
luh2_plot_bii <- function(output_dir, cty_shp, save_png = TRUE) {

    cty_name <- basename(output_dir)

    # ----------------------------------------------------------------
    # Load GDX
    # ----------------------------------------------------------------
    m_out <- gamstransfer::Container$new(
        file.path(output_dir, "solution_lu.gdx")
    )

    bii_df <- m_out["bii_pixel"]$records
    names(bii_df) <- c("pid", "year", "bii")

    # pid_coords for lat/lon
    pid_coords <- m_out["pid_coords"]$records |>
        tidyr::pivot_wider(names_from = "coord", values_from = "value")

    # carea for net hectare calculation (kept for future use / display)
    carea_df <- m_out["carea"]$records
    names(carea_df) <- c("pid", "carea_kha")

    # Join spatial coordinates
    bii_df <- bii_df |>
        dplyr::left_join(pid_coords[, c("pid", "x", "y")], by = "pid") |>
        dplyr::left_join(carea_df, by = "pid")

    # ----------------------------------------------------------------
    # Country boundary
    # ----------------------------------------------------------------
    cty_vect <- terra::vect(cty_shp)
    ctyx     <- cty_vect[cty_vect$NEW_REGION %in% cty_name, ]

    if (cty_name == "USA") {
        ctyx <- terra::crop(ctyx, terra::ext(-125, -66, 24, 50))
    }

    ext_run <- c(
        range(bii_df$x, na.rm = TRUE),
        range(bii_df$y, na.rm = TRUE)
    )

    # ----------------------------------------------------------------
    # Subset years
    # ----------------------------------------------------------------
    bii_2021 <- bii_df |> dplyr::filter(year == "2021", !is.na(bii))
    bii_2050 <- bii_df |> dplyr::filter(year == "2050", !is.na(bii))

    # Shared absolute scale across 2021 and 2050 panels
    bii_lim <- range(c(bii_2021$bii, bii_2050$bii), na.rm = TRUE)

    # Change: 2050 minus 2021
    bii_change <- bii_2021 |>
        dplyr::select(pid, x, y, bii_base = bii) |>
        dplyr::left_join(
            bii_2050 |> dplyr::select(pid, bii_2050 = bii),
            by = "pid"
        ) |>
        dplyr::mutate(delta = bii_2050 - bii_base) |>
        dplyr::filter(!is.na(delta))

    delta_lim <- max(abs(bii_change$delta), na.rm = TRUE)

    # ----------------------------------------------------------------
    # Shared map theme
    # ----------------------------------------------------------------
    map_theme <- ggplot2::theme_minimal(base_size = 15) +
        ggplot2::theme(
            axis.title        = ggplot2::element_blank(),
            axis.text         = ggplot2::element_text(size = 5),
            legend.position   = "bottom",
            legend.key.width  = ggplot2::unit(1.2, "cm"),
            legend.key.height = ggplot2::unit(0.4, "cm")
        )

    coord_run <- ggplot2::coord_sf(
        xlim = c(ext_run[1] - 1, ext_run[2] + 1),
        ylim = c(ext_run[3] - 1, ext_run[4] + 1)
    )

    # ----------------------------------------------------------------
    # Panel 1: BII 2021
    # ----------------------------------------------------------------
    p_2021 <- ggplot2::ggplot() +
        ggplot2::geom_tile(
            data = bii_2021,
            ggplot2::aes(x = x, y = y, fill = bii)
        ) +
        tidyterra::geom_spatvector(
            data = ctyx, fill = NA, color = "black", linewidth = 0.3
        ) +
        ggplot2::scale_fill_distiller(
            palette   = "YlGn",
            direction = 1,
            limits    = bii_lim,
            name      = "BII score",
            na.value  = "transparent"
        ) +
        coord_run +
        map_theme +
        ggplot2::labs(title = "BII 2021 (baseline)")

    # ----------------------------------------------------------------
    # Panel 2: BII 2050
    # ----------------------------------------------------------------
    p_2050 <- ggplot2::ggplot() +
        ggplot2::geom_tile(
            data = bii_2050,
            ggplot2::aes(x = x, y = y, fill = bii)
        ) +
        tidyterra::geom_spatvector(
            data = ctyx, fill = NA, color = "black", linewidth = 0.3
        ) +
        ggplot2::scale_fill_distiller(
            palette   = "YlGn",
            direction = 1,
            limits    = bii_lim,
            name      = "BII score",
            na.value  = "transparent"
        ) +
        coord_run +
        map_theme +
        ggplot2::labs(title = "BII 2050")

    # ----------------------------------------------------------------
    # Panel 3: ΔBII (2050 - 2021)
    # ----------------------------------------------------------------
    p_change <- ggplot2::ggplot() +
        ggplot2::geom_tile(
            data = bii_change,
            ggplot2::aes(x = x, y = y, fill = delta)
        ) +
        tidyterra::geom_spatvector(
            data = ctyx, fill = NA, color = "black", linewidth = 0.3
        ) +
        ggplot2::scale_fill_gradient2(
            low      = "#d73027",
            mid      = "white",
            high     = "#1a9850",
            midpoint = 0,
            limits   = c(-delta_lim, delta_lim),
            name     = "\u0394BII",
            na.value = "transparent"
        ) +
        coord_run +
        map_theme +
        ggplot2::labs(title = "BII change 2021\u20132050")

    # ----------------------------------------------------------------
    # Combine
    # ----------------------------------------------------------------
    p_combined <- ggpubr::ggarrange(
        p_2021, p_2050, p_change,
        ncol   = 3,
        nrow   = 1,
        labels = c("A", "B", "C"),
        align  = "hv"
    )

    # ----------------------------------------------------------------
    # Save
    # ----------------------------------------------------------------
    if (save_png) {
        out_path <- file.path(output_dir, paste0(cty_name, "-bii.png"))
        ggplot2::ggsave(
            plot     = p_combined,
            filename = out_path,
            width    = 18,
            height   = 7,
            bg       = "white"
        )
        message("BII plot saved to ", out_path)
    }

    invisible(p_combined)
}
