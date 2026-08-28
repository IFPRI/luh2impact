#' Compute IMPACT-country GDP per capita from SSP data with WB fallback
#'
#' Builds annual GDP per capita (thousand USD) by IMPACT country from SSP2
#' GDP and population series in an SSP scenario GDX, aggregated from the
#' SSPdb country level to IMPACT regions. Countries missing from the SSP GDX
#' mapping are backfilled from World Bank GDP per capita (PPP), with manual
#' extrapolations for Afghanistan, Syria, and Venezuela - see Details.
#'
#' @details
#' Manual CAGR-based backfills/forward-fills applied after the World Bank
#' merge, replicating known analyst adjustments used in the original
#' elasticity workbook because World Bank coverage does not reach back to
#' `year_min` (Afghanistan, Syria) or capture recent years cleanly
#' (Venezuela, due to hyperinflation-era swings):
#' \itemize{
#'   \item AFG: years `year_min`-1999 back-filled from the 1999-2024 CAGR
#'     `(2201/813)^(1/25)-1`.
#'   \item SYR: years `year_min`-2016 back-filled from the 2016-2022 CAGR
#'     `(4772/3265)^(1/6)-1`; 2023 manually set to 4800 (current USD, pre-scaling).
#'   \item VEN: years 2012-2023 forward-filled from a 55-year CAGR
#'     `(21241/13889)^(1/55)-1`.
#' }
#'
#' @param ssp_gdx Path to the SSP scenario GDX file (must contain `OECD_GDP`
#'   and `POP` symbols).
#' @param fao_country_map_csv Path to a FAOSTAT country-code-to-ISO3 mapping
#'   CSV (columns 4 and 8 must be the FAO country code and ISO3 code).
#' @param sets_xlsx Path to the IMPACT `Sets.xlsx` workbook (see [luh2_load_faoland()]).
#' @param wb_gdppc_csv Path to the World Bank `NY.GDP.PCAP.PP.CD` GDP per
#'   capita (PPP) CSV, used to backfill IMPACT countries missing from the SSP GDX.
#' @param ssp Character. SSP scenario to use. Default `"SSP2"`.
#' @param year_min,year_max Integer. Year range to retain from the World Bank series.
#'
#' @return A data frame with columns `cty`, `year`, `value` (thousand USD
#'   per capita), `item = "gdppc"`, `unit = "000USD"`.
#' @author Abhijeet Mishra, Claude Code
#' @export
luh2_compute_gdppc <- function(ssp_gdx, fao_country_map_csv, sets_xlsx, wb_gdppc_csv,
                               ssp = "SSP2", year_min = 1995, year_max = 2023) {

    ctyset <- .luh2_read_ctyset(sets_xlsx)

    fao_master_list <- utils::read.csv(fao_country_map_csv)[, c(4, 8)]
    names(fao_master_list) <- c("fcty", "iso")
    fao_master_list$fcty[fao_master_list$iso %in% c("CHN", "TWN", "HKG", "MAC")] <- 351

    gdp <- DOORMAT::readGDX(gdx = ssp_gdx, name = "OECD_GDP", quick_df = TRUE, use_model_name = FALSE)$data
    names(gdp)[1] <- "ssp"
    names(gdp)[3] <- "year"
    gdp <- gdp[gdp$ssp %in% ssp & gdp$var %in% "GDP|PPP", -4]

    gdp_impact <- dplyr::left_join(gdp, fao_master_list, by = "iso")
    gdp_impact <- dplyr::left_join(gdp_impact, ctyset, by = "fcty")
    gdp_impact <- gdp_impact[!is.na(gdp_impact$cty), ]
    gdp_impact <- gdp_impact |>
        dplyr::group_by(cty, year) |>
        dplyr::summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

    pop <- DOORMAT::readGDX(gdx = ssp_gdx, name = "POP", quick_df = TRUE, use_model_name = FALSE)$data
    names(pop)[1] <- "ssp"
    names(pop)[5] <- "year"
    pop <- pop[pop$ssp %in% ssp & pop$gender %in% "BOTH" & pop$tranche %in% "PTOTL", c(-1, -3, -4)]

    pop_impact <- dplyr::left_join(pop, fao_master_list, by = "iso")
    pop_impact <- dplyr::left_join(pop_impact, ctyset, by = "fcty")
    pop_impact <- pop_impact[!is.na(pop_impact$cty), ]
    pop_impact <- pop_impact |>
        dplyr::group_by(cty, year) |>
        dplyr::summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

    gdppc <- dplyr::left_join(gdp_impact, pop_impact, by = c("cty", "year"), suffix = c(".gdp", ".pop"))
    gdppc$value <- (gdppc$value.gdp / gdppc$value.pop) / 1e3
    gdppc <- gdppc[, c("cty", "year", "value")]
    gdppc$item <- "gdppc"
    gdppc$unit <- "000USD"
    gdppc$year <- as.integer(as.character(gdppc$year))

    # Backfill countries missing from the SSP GDX mapping from World Bank data
    mcty <- setdiff(unique(ctyset$cty), unique(gdppc$cty))

    wbgdppc <- utils::read.csv(wb_gdppc_csv, header = TRUE, skip = 4)
    wbgdppc <- wbgdppc[wbgdppc$Country.Code %in% mcty,
                       !(names(wbgdppc) %in% c("Country.Name", "Indicator.Name", "Indicator.Code"))]
    wbgdppc <- tidyr::pivot_longer(wbgdppc, cols = dplyr::starts_with("X"),
                                   names_to = "year", values_to = "value")
    wbgdppc$year <- as.numeric(gsub("X", "", wbgdppc$year))
    wbgdppc <- wbgdppc[!is.na(wbgdppc$year) & wbgdppc$year %in% year_min:year_max, ]

    # AFG: back-fill from the 1999-2024 CAGR
    cagr <- ((2201 / 813)^(1 / 25)) - 1
    for (i in seq(1999, year_min, -1)) {
        oldval <- wbgdppc$value[wbgdppc$Country.Code %in% "AFG" & wbgdppc$year == i + 1]
        wbgdppc$value[wbgdppc$Country.Code %in% "AFG" & wbgdppc$year == i] <- oldval * (1 - cagr)
    }

    # SYR: back-fill from the 2016-2022 CAGR; manual 2023 estimate
    cagr <- ((4772 / 3265)^(1 / 6)) - 1
    for (i in seq(2016, year_min, -1)) {
        oldval <- wbgdppc$value[wbgdppc$Country.Code %in% "SYR" & wbgdppc$year == i + 1]
        wbgdppc$value[wbgdppc$Country.Code %in% "SYR" & wbgdppc$year == i] <- oldval * (1 - cagr)
    }
    wbgdppc$value[wbgdppc$Country.Code %in% "SYR" & wbgdppc$year == 2023] <- 4800

    # VEN: forward-fill from a 55-year CAGR to smooth hyperinflation-era swings
    cagr <- ((21241 / 13889)^(1 / 55)) - 1
    for (i in seq(2012, 2023, 1)) {
        oldval <- wbgdppc$value[wbgdppc$Country.Code %in% "VEN" & wbgdppc$year == i - 1]
        wbgdppc$value[wbgdppc$Country.Code %in% "VEN" & wbgdppc$year == i] <- oldval * (1 + cagr)
    }

    names(wbgdppc)[names(wbgdppc) == "Country.Code"] <- "cty"
    wbgdppc <- wbgdppc[, c("cty", "year", "value")]
    wbgdppc$value <- wbgdppc$value / 1e3
    wbgdppc$item <- "gdppc"
    wbgdppc$unit <- "000USD"

    rbind(gdppc, wbgdppc)
}
