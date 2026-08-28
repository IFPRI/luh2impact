#' Fit land-use elasticities to cropland and GDP per capita
#'
#' Builds a per-country annual panel of land-use areas and GDP per capita,
#' imputes missing years by linear interpolation, first-differences the logs,
#' and fits a per-country OLS regression of each non-cropland land use's
#' log-change on cropland and GDP-per-capita log-changes.
#'
#' @param faoland A data frame returned by [luh2_load_faoland()].
#' @param gdppc A data frame returned by [luh2_compute_gdppc()].
#' @param year_max Integer. Latest year to include in the panel. Default `2023`.
#' @param min_obs Integer. Minimum number of complete observations required
#'   to fit a country's regression; below this the elasticity is flagged
#'   `"insufficient_obs"` and set to `0`. Default `10`.
#'
#' @return A list with:
#' \describe{
#'   \item{results}{One row per country x land use, with `elast_crop`,
#'     `elast_gdppc` (each clipped to `[-4, 4]` / `[-2, 2]`), standard errors,
#'     p-values, `n_obs`, and a `flag` (`"ok"`, `"insufficient_obs"`,
#'     `"model_failed"`, `"dropped_term"`).}
#'   \item{wide_imputed}{The imputed wide panel (`cty`, `year`, `land`,
#'     `crop`, `past`, `natfor`, `plant`, `other`, `gdppc`) used to fit the
#'     regressions.}
#' }
#' @author Abhijeet Mishra, Claude Code
#' @export
luh2_fit_elasticities <- function(faoland, gdppc, year_max = 2023, min_obs = 10) {
    non_crop_uses <- c("past", "natfor", "other", "plant")

    df <- rbind(faoland, gdppc)
    df <- df[order(df$cty, df$year), ]
    df <- df[df$year %in% intersect(unique(faoland$year), unique(gdppc$year)), ]

    wide <- df |>
        dplyr::filter(year <= year_max) |>
        dplyr::select(cty, item, year, value) |>
        tidyr::pivot_wider(names_from = item, values_from = value) |>
        dplyr::arrange(cty, year)

    wide_imputed <- wide |>
        dplyr::group_by(cty) |>
        tidyr::complete(year = seq(min(year), max(year))) |>
        dplyr::mutate(dplyr::across(c(crop, dplyr::all_of(non_crop_uses), gdppc, land), .luh2_approx_na)) |>
        dplyr::ungroup()

    fd <- wide_imputed |>
        dplyr::group_by(cty) |>
        dplyr::mutate(dplyr::across(c(crop, dplyr::all_of(non_crop_uses), gdppc),
                                    \(x) c(NA, diff(log(ifelse(x <= 0, NA, x)))),
                                    .names = "dln_{.col}")) |>
        dplyr::ungroup() |>
        dplyr::filter(!is.na(dln_crop))

    fit_elasticity <- function(data, use) {
        dv  <- paste0("dln_", use)
        frm <- stats::as.formula(paste(dv, "~ dln_crop + dln_gdppc"))

        data |>
            dplyr::group_by(cty) |>
            dplyr::group_map(\(d, g) {
                d_clean <- d |> dplyr::select(dplyr::all_of(c(dv, "dln_crop", "dln_gdppc"))) |> tidyr::drop_na()
                n_obs <- nrow(d_clean)

                if (n_obs < min_obs) {
                    data.frame(cty = g$cty, land_use = use,
                              elast_crop = 0, se_crop = NA_real_, pval_crop = NA_real_,
                              elast_gdppc = 0, se_gdppc = NA_real_, pval_gdppc = NA_real_,
                              n_obs = n_obs, flag = "insufficient_obs")
                } else {
                    fit <- tryCatch(stats::lm(frm, data = d_clean), error = \(e) NULL)
                    if (is.null(fit)) {
                        data.frame(cty = g$cty, land_use = use,
                                  elast_crop = 0, se_crop = NA_real_, pval_crop = NA_real_,
                                  elast_gdppc = 0, se_gdppc = NA_real_, pval_gdppc = NA_real_,
                                  n_obs = n_obs, flag = "model_failed")
                    } else {
                        coef_tbl <- summary(fit)$coefficients
                        data.frame(
                            cty = g$cty, land_use = use,
                            elast_crop  = .luh2_safe_coef(coef_tbl, "dln_crop",  "Estimate", 0),
                            se_crop     = .luh2_safe_coef(coef_tbl, "dln_crop",  "Std. Error"),
                            pval_crop   = .luh2_safe_coef(coef_tbl, "dln_crop",  "Pr(>|t|)"),
                            elast_gdppc = .luh2_safe_coef(coef_tbl, "dln_gdppc", "Estimate", 0),
                            se_gdppc    = .luh2_safe_coef(coef_tbl, "dln_gdppc", "Std. Error"),
                            pval_gdppc  = .luh2_safe_coef(coef_tbl, "dln_gdppc", "Pr(>|t|)"),
                            n_obs = n_obs,
                            flag = if (any(!c("dln_crop", "dln_gdppc") %in% rownames(coef_tbl))) "dropped_term" else "ok"
                        )
                    }
                }
            }, .keep = TRUE) |>
            dplyr::bind_rows()
    }

    results <- lapply(non_crop_uses, function(u) fit_elasticity(fd, u)) |>
        dplyr::bind_rows() |>
        dplyr::arrange(land_use, cty) |>
        dplyr::mutate(
            elast_crop  = pmax(-4, pmin(4, elast_crop)),
            elast_gdppc = pmax(-2, pmin(2, elast_gdppc))
        )

    list(results = results, wide_imputed = wide_imputed)
}
