#' Compute per-country BII coefficients per land pool from pixel imagery
#'
#' Reads the pixel-level `BII_COEFF(cty,fland)` table derived from Trong
#' Nguyen's 2020 High-Resolution Harmonized Land Use (HHLU, 0.004deg) and
#' Biodiversity Intactness Index (BII, 0.02deg) rasters, as an alternative
#' to the literature-based [luh2_bii_coeff_cty()]. Unlike that function,
#' every pool (`crop`, `urban`, `past`, `plant`, `natfor`, `other`) gets a
#' country-specific value derived directly from imagery, not just `natfor`.
#'
#' @details
#' `crop`/`urban`/`past`/`plant` are direct area-weighted mean BII scores
#' over HHLU pixels of the matching land-use class within each country.
#'
#' `natfor`/`other` cannot be split from HHLU pixel values alone - HHLU's
#' three "natural vegetation" classes (primary minimal-use / primary /
#' secondary) are organized by successional/management-intensity stage, not
#' by cover type, and no grouping of them correlated well with this
#' package's existing FAOSTAT-derived natfor/other area split (best
#' correlation tested across 6 groupings was ~0.35). Instead, the combined
#' natural-land pixel score was split back into `natfor`/`other` assuming
#' `other_coeff = 0.90 * natfor_coeff`, solved per country so the
#' `LAND00`-area-weighted average of the two still matches the
#' pixel-observed combined value. In countries where `other` area vastly
#' exceeds `natfor` area (arid/desert countries: DZA, GRL, LBY, MLI, MRT,
#' NER, TCD), this algebra would push `natfor` above the physical BII
#' ceiling of 1.0; those 7 countries have `natfor` capped at 1.0 with
#' `other` re-solved from the mass-balance equation alone.
#'
#' See `pixel_bii_coeff_r090.csv`'s sibling `README.txt`, shipped alongside
#' it under `inst/extdata`, for the full two-script pipeline that produced
#' it from the raw rasters plus this package's `LAND00`.
#'
#' @param pixel_csv Path to the pixel-derived BII coefficient CSV. Defaults
#'   to the copy shipped with the package
#'   (`inst/extdata/pixel_bii_coeff_r090.csv`).
#'
#' A handful of countries have no HHLU pixels at all for `past` (1 country)
#' or `plant` (3 countries) - those rows are dropped (`value` would
#' otherwise be `NA`), matching how `LAND00` already omits a `fland` pool
#' entirely for a country with none of it.
#'
#' @return A data frame with columns `cty`, `fland`, `value` (BII
#'   coefficient, 0-1), one row per country per pool with data - same shape
#'   as [luh2_bii_coeff_cty()]'s output, ready to pass as the `bii_coeff`
#'   argument to [luh2_export_elasticities_gdx()].
#' @author Abhijeet Mishra, Claude Code
#' @export
luh2_bii_coeff_cty_pixel <- function(pixel_csv = system.file("extdata", "pixel_bii_coeff_r090.csv",
                                                              package = "luh2impact")) {
    raw <- utils::read.csv(pixel_csv, stringsAsFactors = FALSE)

    pool_cols <- c(crop = "bii_crop", urban = "bii_urban", past = "bii_past",
                   plant = "bii_plant", natfor = "bii_natfor", other = "bii_other")

    out <- do.call(rbind, lapply(names(pool_cols), function(fl) {
        data.frame(cty = raw$cty, fland = fl, value = raw[[pool_cols[[fl]]]])
    }))
    out[!is.na(out$value), ]
}
