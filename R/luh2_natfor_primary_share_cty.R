#' Compute FAOSTAT-native primary-forest share of the `natfor` pool
#'
#' For each IMPACT country, computes the share of `natfor` (naturally
#' regenerating forest) that is primary forest, from FAOSTAT items already
#' present in the bulk land-use file - no LUH2 data needed. Follows the
#' FAO/FRA definitional hierarchy: naturally regenerating forest is primary
#' forest plus other naturally regenerated ("secondary") forest, so
#' secondary forest = `Naturally regenerating forest - Primary Forest`.
#'
#' Spot-checked against countries with known large tracts of intact forest
#' (Brazil, DRC, Indonesia, Russia, Canada), and the shares come out
#' sensibly high relative to the global median (~1%). Equivalent splits for
#' `past`/`other`/`crop` were explored and abandoned: FAOSTAT's pasture
#' sub-items are not populated for most countries, and LUH2's
#' primary/secondary non-forest classification did not produce a
#' trustworthy signal for `other` (e.g. Saudi Arabia's near-zero denominator).
#' This function stays intentionally narrow to what's been validated.
#'
#' Returned as a share, not an absolute area - meant to be used as a
#' blending weight against `natfor`'s existing FAOSTAT-sourced total (e.g.
#' in `LAND00`), not as a substitute for it.
#'
#' @param fao_land_csv Path to the FAOSTAT `Inputs_LandUse_E_All_Data_(Normalized).csv` file.
#' @param sets_xlsx Path to the IMPACT `Sets.xlsx` workbook (see [luh2_load_faoland()]).
#' @param year Integer. Year to compute the share for. Default `2021`.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{cty}{IMPACT country code.}
#'   \item{natfor_primary_share}{Share of `natfor` that is primary forest (0-1, `NA` if no naturally regenerating forest reported).}
#'   \item{flag_natfor}{`"ok"` or `"primary_exceeds_natregen"` if `Primary Forest > Naturally
#'     regenerating forest` (data inconsistency; share is clipped to 1).}
#' }
#' @author Abhijeet Mishra, Claude Code
#' @export
luh2_natfor_primary_share_cty <- function(fao_land_csv, sets_xlsx, year = 2021) {
    items_needed <- c("Primary Forest", "Naturally regenerating forest")

    raw <- as.data.frame(data.table::fread(file = fao_land_csv, header = TRUE))
    names(raw) <- tolower(names(raw))

    raw <- raw[raw$unit %in% "1000 ha" & raw$item %in% items_needed & raw$year == year,
              c("area code", "item", "value")]
    names(raw)[1] <- "fcty"
    raw$fcty <- as.numeric(raw$fcty)

    ctyset <- .luh2_read_ctyset(sets_xlsx)
    raw <- dplyr::left_join(raw, ctyset, by = "fcty")
    raw <- raw[!is.na(raw$cty), ]

    agg <- raw |>
        dplyr::group_by(cty, item) |>
        dplyr::summarise(value = sum(value, na.rm = TRUE), .groups = "drop") |>
        tidyr::pivot_wider(names_from = item, values_from = value, values_fill = 0)

    for (col in items_needed) if (!col %in% names(agg)) agg[[col]] <- 0

    agg$flag_natfor <- ifelse(agg[["Primary Forest"]] > agg[["Naturally regenerating forest"]] + 1e-6,
                              "primary_exceeds_natregen", "ok")
    agg$natfor_primary_share <- ifelse(agg[["Naturally regenerating forest"]] > 0,
                                       agg[["Primary Forest"]] / agg[["Naturally regenerating forest"]], NA_real_)
    agg$natfor_primary_share <- pmin(pmax(agg$natfor_primary_share, 0), 1)

    as.data.frame(agg[, c("cty", "natfor_primary_share", "flag_natfor")])
}
