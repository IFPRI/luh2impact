#' Plot land-use elasticity distributions
#'
#' Visualizes the country-level elasticities returned by
#' [luh2_fit_elasticities()]: a boxplot of elasticities by driver and land
#' use with outlier countries labelled, and a histogram of the same.
#'
#' @param results The `results` data frame from [luh2_fit_elasticities()].
#' @param region_map_xlsx Optional path to a region-mapping workbook with
#'   columns `cty` and `region`, used to restrict the plot to a single
#'   region (e.g. the world aggregate). If `NULL` (default), all countries
#'   in `results` are plotted.
#' @param region Character. Region to filter to when `region_map_xlsx` is
#'   supplied. Default `"WLD"`.
#' @param output_dir Optional directory to save `elas_box.png` and
#'   `elas_hist.png`. If `NULL` (default), plots are not written to disk.
#'
#' @return A list with `box` and `hist` `ggplot` objects, invisibly.
#' @author Abhijeet Mishra, Claude Code
#' @export
luh2_plot_elasticities <- function(results, region_map_xlsx = NULL, region = "WLD", output_dir = NULL) {
    vis <- results

    if (!is.null(region_map_xlsx)) {
        region_map <- openxlsx2::read_xlsx(file = region_map_xlsx)
        vis <- dplyr::left_join(vis, region_map, by = "cty")
        vis <- vis[vis$region %in% region, ]
    }

    vis_long <- vis |>
        dplyr::filter(flag == "ok") |>
        dplyr::select(cty, land_use, elast_crop, elast_gdppc) |>
        tidyr::pivot_longer(cols = c(elast_crop, elast_gdppc), names_to = "driver", values_to = "elasticity") |>
        dplyr::filter(!is.na(elasticity), elasticity != 0)

    outlier_data <- vis_long |>
        dplyr::group_by(land_use, driver) |>
        dplyr::mutate(
            q1  = stats::quantile(elasticity, 0.25, na.rm = TRUE),
            q3  = stats::quantile(elasticity, 0.75, na.rm = TRUE),
            iqr = q3 - q1,
            is_outlier = elasticity < q1 - 1.5 * iqr | elasticity > q3 + 1.5 * iqr
        ) |>
        dplyr::filter(is_outlier)

    driver_colors <- c(elast_crop = "#2196F3", elast_gdppc = "#FF9800")

    p_box <- ggplot2::ggplot(vis_long, ggplot2::aes(x = driver, y = elasticity, fill = driver)) +
        ggplot2::geom_boxplot(outlier.shape = 21, alpha = 0.4, staplewidth = 0.2) +
        ggrepel::geom_text_repel(data = outlier_data, ggplot2::aes(label = cty), size = 3) +
        ggplot2::scale_fill_manual(values = driver_colors) +
        ggplot2::labs(x = NULL, y = "Elasticity", fill = "Driver") +
        ggplot2::facet_wrap(~land_use, scales = "free_y") +
        ggplot2::theme_bw(base_size = 15) +
        ggplot2::theme(legend.position = "none")

    p_hist <- ggplot2::ggplot(vis_long, ggplot2::aes(x = elasticity, fill = driver)) +
        ggplot2::geom_histogram(alpha = 0.6, bins = 30, position = "identity") +
        ggplot2::scale_fill_manual(values = driver_colors) +
        ggplot2::labs(x = "Elasticity", y = "Count", fill = "Driver") +
        ggplot2::facet_wrap(~land_use, scales = "free") +
        ggplot2::theme_bw(base_size = 10)

    if (!is.null(output_dir)) {
        ggplot2::ggsave(file.path(output_dir, "elas_box.png"),  plot = p_box,  width = 8, height = 8)
        ggplot2::ggsave(file.path(output_dir, "elas_hist.png"), plot = p_hist, width = 7, height = 6)
        message("Plots saved to ", output_dir)
    }

    invisible(list(box = p_box, hist = p_hist))
}
