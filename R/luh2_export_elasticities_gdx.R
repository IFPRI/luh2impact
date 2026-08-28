#' Export land-use elasticities and initial areas to GDX
#'
#' Writes elasticity estimates (raw and "ok"-flagged), observation counts,
#' and baseline FAOSTAT land areas to a GDX file for use as LUMEN model inputs.
#'
#' @param results The `results` data frame from [luh2_fit_elasticities()].
#' @param faoland A data frame returned by [luh2_load_faoland()], used to
#'   pull baseline land areas for `base_year`.
#' @param output_gdx Path to write the output GDX file.
#' @param base_year Integer. Year to use for the `LAND00` baseline land
#'   parameter. Default `2021`.
#' @param extra_land00 Optional data frame with columns `cty`, `item`,
#'   `value` to append to `LAND00` as additional baseline land-use pools
#'   that are genuinely distinct land covers (e.g. `item = "urban"` from
#'   [luh2_urban_area_cty()]). These pools have no corresponding
#'   `LANDELAS` entry - they are fixed baseline stocks, not fitted against
#'   cropland/GDP dynamics. Protected area does NOT belong here - it is not
#'   a distinct land cover, it overlays crop/natfor/past/other, so it is
#'   written separately as `LANDPA00` via the `pa00` argument. Default
#'   `NULL` (no extra pools).
#' @param pa00 Optional data frame with columns `cty`, `value` (aggregate
#'   protected non-crop land area, 1000 ha), e.g. from [luh2_pa_area_cty()]
#'   or [wdpa_process()]. Written as a separate `LANDPA00(cty)` GDX
#'   parameter - a floor on protected land, not a `LAND00` pool. Default
#'   `NULL` (skip).
#' @param bii_coeff Optional data frame with columns `cty`, `fland`, `value`
#'   (BII coefficient, 0-1), e.g. from [luh2_bii_coeff_cty()]. Written as a
#'   `BII_COEFF(cty,fland)` GDX parameter for computing country-level BII in
#'   IMPACT. Default `NULL` (skip).
#' @param author Character. Recorded in the GDX as provenance metadata
#'   (see `gdx_metadata` in the Return section). Defaults to the OS user
#'   running the current R session (`Sys.info()[["user"]]`) if `NULL` - not
#'   the package maintainer - since the point is tracking who actually
#'   generated this specific file, which matters once it's version
#'   controlled and multiple people may run this function.
#'
#' @return Invisibly returns `output_gdx`. Called for its side effect. The
#'   GDX also contains a `gdx_metadata` set recording `author=...`,
#'   `created_at=...` (local timestamp), and `package=...` (installed
#'   `luh2impact` version) - GDX has no native file-level metadata, so this
#'   is stored as a small set whose elements carry the information as text.
#' @author Abhijeet Mishra, Claude Code
#' @export
luh2_export_elasticities_gdx <- function(results, faoland, output_gdx, base_year = 2021,
                                         extra_land00 = NULL, pa00 = NULL, bii_coeff = NULL,
                                         author = NULL) {
    if (is.null(author)) {
        author <- tryCatch(Sys.info()[["user"]], error = function(e) NA_character_)
        if (is.na(author) || !nzchar(author)) author <- "unknown"
    }
    pkg_version <- tryCatch(as.character(utils::packageVersion("luh2impact")),
                            error = function(e) "unknown")
    land00 <- faoland[faoland$year == base_year & faoland$item %in% c("crop", "natfor", "other", "past", "plant"),
                      c("cty", "item", "value")]
    land00$item <- as.character(land00$item)

    if (!is.null(extra_land00)) {
        land00 <- rbind(land00, extra_land00[, c("cty", "item", "value")])
    }

    raw_long <- results |>
        tidyr::pivot_longer(
            cols = c(elast_crop, se_crop, pval_crop, elast_gdppc, se_gdppc, pval_gdppc),
            names_to = c(".value", "driver"),
            names_pattern = "^(elast|se|pval)_(crop|gdppc)$"
        ) |>
        dplyr::rename(elasticity = elast)

    landelas_long <- raw_long |>
        dplyr::filter(flag == "ok", !is.na(elasticity)) |>
        dplyr::select(cty, land_use, driver, elasticity)

    m <- gamstransfer::Container$new()

    add_stat <- function(stat_col, param_name, desc) {
        m$addParameter(
            name = param_name,
            domain = list("*", "*", "*"),
            records = raw_long |>
                dplyr::select(cty, land_use, driver, value = dplyr::all_of(stat_col)) |>
                dplyr::filter(!is.na(value)) |>
                as.data.frame(),
            description = desc
        )
    }

    add_stat("elasticity", "LANDELAS_RAW_ELAST", "Elasticity estimates, all flags")
    add_stat("se",         "LANDELAS_RAW_SE",    "Standard errors, all flags")
    add_stat("pval",       "LANDELAS_RAW_PVAL",  "P-values, all flags")

    m$addParameter(
        name = "LANDELAS_RAW_NOBS",
        domain = list("*", "*"),
        records = results |> dplyr::select(cty, land_use, n_obs) |> dplyr::distinct() |> as.data.frame(),
        description = "Observation count per country x land use"
    )

    m$addParameter(
        name = "LANDELAS",
        domain = list("*", "*", "*"),
        records = as.data.frame(landelas_long),
        description = "Land-use elasticities wrt cropland and gdppc, ok-flagged only"
    )

    m$addParameter(
        name = "LAND00",
        domain = list("*", "*"),
        records = as.data.frame(land00),
        description = "Initial FAOSTAT land"
    )

    if (!is.null(pa00)) {
        m$addParameter(
            name = "LANDPA00",
            domain = list("*"),
            records = as.data.frame(pa00[, c("cty", "value")]),
            description = "Protected non-crop land area, aggregate floor (1000 ha)"
        )
    }

    if (!is.null(bii_coeff)) {
        m$addParameter(
            name = "BII_COEFF",
            domain = list("*", "*"),
            records = as.data.frame(bii_coeff[, c("cty", "fland", "value")]),
            description = "BII coefficient per country per land pool (0-1)"
        )
    }

    m$addSet(
        "gdx_metadata",
        records = c(
            paste0("author=", author),
            paste0("created_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
            paste0("package=luh2impact ", pkg_version)
        ),
        description = "Provenance metadata for this GDX file"
    )

    m$write(output_gdx)
    invisible(output_gdx)
}
