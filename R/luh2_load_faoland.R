#' Load and map FAOSTAT land use data to IMPACT countries
#'
#' Reads the FAOSTAT "Inputs: Land Use" bulk download, keeps the land-use
#' items used for elasticity estimation, maps FAO countries to IMPACT
#' regions via the `Regions` sheet of the IMPACT sets workbook, and
#' aggregates to IMPACT country level.
#'
#' @param fao_land_csv Path to the FAOSTAT `Inputs_LandUse_E_All_Data_(Normalized).csv` file.
#' @param sets_xlsx Path to the IMPACT `Sets.xlsx` workbook. Must contain a
#'   `Regions` sheet mapping FAO country codes (`FCTY`) to IMPACT regions (`CTY`).
#' @param year_min Integer. Earliest year to retain. Default `1995`.
#'
#' @return A data frame with columns `cty`, `item`, `year`, `unit`, `value`,
#'   where `item` is one of `land`, `crop`, `past`, `natfor`, `plant`, `other`.
#' @author Abhijeet Mishra, Claude Code
#' @export
luh2_load_faoland <- function(fao_land_csv, sets_xlsx, year_min = 1995) {
    faoland <- as.data.frame(data.table::fread(file = fao_land_csv, header = TRUE))
    names(faoland) <- tolower(names(faoland))

    itemlist <- c("Land area", "Cropland", "Naturally regenerating forest",
                  "Other land", "Permanent meadows and pastures", "Planted Forest")

    faoland <- faoland[faoland$unit %in% "1000 ha" & faoland$item %in% itemlist,
                       c("area code", "item", "year", "unit", "value")]
    names(faoland)[1] <- "fcty"
    faoland$fcty <- as.numeric(faoland$fcty)

    faoland$item <- factor(faoland$item,
                           levels = c("Land area",
                                      "Cropland", "Permanent meadows and pastures",
                                      "Naturally regenerating forest", "Planted Forest",
                                      "Other land"),
                           labels = c("land",
                                      "crop", "past",
                                      "natfor", "plant",
                                      "other"))

    ctyset  <- .luh2_read_ctyset(sets_xlsx)
    faoland <- dplyr::left_join(faoland, ctyset, by = "fcty")
    faoland <- faoland[!is.na(faoland$cty), ]

    faoland <- faoland |>
        dplyr::group_by(cty, item, year, unit) |>
        dplyr::summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

    faoland <- as.data.frame(faoland[faoland$year >= year_min, ])
    faoland
}
