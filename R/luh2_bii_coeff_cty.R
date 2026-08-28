#' Compute per-country BII coefficients per land pool
#'
#' Blends [luh2_natfor_primary_share_cty()]'s FAOSTAT-derived primary-forest
#' share with primary/secondary forest BII coefficients from Leclere et al.
#' (2018), Figure 3, to give `natfor` a country-specific BII coefficient
#' rather than one flat global value. All other pools (`crop`, `past`,
#' `other`, `plant`, `urban`) use the same flat coefficients as the existing
#' LUMEN pixel-level implementation (`lumen.gms`'s `bii_coeff`) - no
#' validated country-level split exists yet for those pools (see
#' [luh2_natfor_primary_share_cty()] for what was tried and abandoned for
#' `crop`/`past`/`other`).
#'
#' @details
#' `bii_secondary_forested` (default `0.85`) is a visual read off Figure 3's
#' forested-panel secondary-vegetation bars (excluding Timber, which maps to
#' `plant`, not `natfor`) - an approximation, not a precisely digitized
#' value. Adjust and re-run if a more precise source becomes available.
#' `bii_urban` (default `0.65`) is likewise a visual read of Figure 3's
#' `Urban` bars - `lumen.gms` has no precedent for it since urban is held
#' fixed there rather than entering the pixel BII calc.
#'
#' Countries with no reported "Naturally regenerating forest" (so
#' `natfor_primary_share` is `NA`) fall back to the midpoint of
#' `bii_primary_forested`/`bii_secondary_forested`, rather than dropping the
#' country or assuming all-primary/all-secondary.
#'
#' @param fao_land_csv Path to the FAOSTAT `Inputs_LandUse_E_All_Data_(Normalized).csv` file.
#' @param sets_xlsx Path to the IMPACT `Sets.xlsx` workbook.
#' @param year Integer. Year to compute the `natfor` primary-forest share for. Default `2021`.
#' @param bii_primary_forested Numeric. BII coefficient for primary forest. Default `1.0`.
#' @param bii_secondary_forested Numeric. BII coefficient for secondary forest (approximate, see Details). Default `0.85`.
#' @param bii_other,bii_past,bii_plant,bii_crop Numeric. Flat BII
#'   coefficients matching `lumen.gms`'s existing `bii_coeff` table.
#'   Defaults `0.888, 0.791, 0.672, 0.576`.
#' @param bii_urban Numeric. Flat BII coefficient for urban land
#'   (approximate, see Details). Default `0.65`.
#'
#' @return A data frame with columns `cty`, `fland`, `value` (BII
#'   coefficient, 0-1), one row per country per pool - ready to export as a
#'   GDX parameter alongside `LAND00`.
#' @author Abhijeet Mishra, Claude Code
#' @export
luh2_bii_coeff_cty <- function(fao_land_csv, sets_xlsx, year = 2021,
                               bii_primary_forested = 1.0, bii_secondary_forested = 0.85,
                               bii_other = 0.888, bii_past = 0.791, bii_plant = 0.672,
                               bii_crop = 0.576, bii_urban = 0.65) {
    shares <- luh2_natfor_primary_share_cty(fao_land_csv, sets_xlsx, year = year)

    natfor_coeff <- shares |>
        dplyr::mutate(
            fland = "natfor",
            value = ifelse(is.na(natfor_primary_share),
                           (bii_primary_forested + bii_secondary_forested) / 2,
                           natfor_primary_share * bii_primary_forested +
                               (1 - natfor_primary_share) * bii_secondary_forested)
        ) |>
        dplyr::select(cty, fland, value)

    flat_pools <- data.frame(
        fland = c("other", "past", "plant", "crop", "urban"),
        value = c(bii_other, bii_past, bii_plant, bii_crop, bii_urban)
    )
    flat_coeff <- do.call(rbind, lapply(unique(shares$cty), function(c) {
        data.frame(cty = c, fland = flat_pools$fland, value = flat_pools$value)
    }))

    rbind(as.data.frame(natfor_coeff), flat_coeff)
}
