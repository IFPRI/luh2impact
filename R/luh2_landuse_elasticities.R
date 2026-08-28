#' Run the full FAOSTAT land-use elasticity pipeline
#'
#' Convenience wrapper that loads FAOSTAT land-use data and GDP per capita,
#' fits per-country land-use elasticities, optionally plots them, and writes
#' the result to a GDX file. Each step can be run and tested individually
#' using the underlying functions.
#'
#' Urban land (`item = "urban"` in `LAND00`, from [luh2_urban_area_cty()])
#' is a fixed baseline land-cover pool. Protected area is NOT a land cover -
#' it overlays crop/natfor/past/other - so it is written separately as an
#' aggregate `LANDPA00(cty)` floor (from [luh2_pa_area_cty()] or
#' [wdpa_process()]), not as a `LAND00` item. Neither gets a fitted
#' `LANDELAS` entry yet. BII coefficients (`BII_COEFF(cty,fland)`, from
#' [luh2_bii_coeff_cty()]) are always computed and written, since they only
#' need `fao_land_csv`/`sets_xlsx` (already required). Pass `bii_pixel_csv`
#' to use [luh2_bii_coeff_cty_pixel()]'s imagery-derived coefficients
#' instead.
#'
#' @inheritParams luh2_load_faoland
#' @inheritParams luh2_compute_gdppc
#' @param output_gdx Path to write the output GDX file.
#' @param plot_dir Optional directory to save elasticity distribution plots.
#'   If `NULL` (default), plots are not generated.
#' @param region_map_xlsx,region Passed to [luh2_plot_elasticities()].
#' @param year_max Integer. Latest year to include when fitting elasticities. Default `2023`.
#' @param base_year Integer. Year to use for baseline land areas in the GDX. Default `2021`.
#' @param min_obs Integer. Minimum observations required per country regression. Default `10`.
#' @param pa_tif Optional path to an already-processed WDPA raster (km^2 per
#'   pixel, e.g. `wdpa_process(..., output = "pixel")` output). If supplied
#'   (along with `cty_shp`), protected area is aggregated to IMPACT
#'   countries via [luh2_pa_area_cty()] and written as `LANDPA00`. Takes
#'   precedence over `gdb_path`. Default `NULL` (skip).
#' @param gdb_path Optional path to the WDPA geodatabase (`.gdb` folder),
#'   used only when `pa_tif` is not supplied. If supplied (along with
#'   `static_nc` and `cty_shp`), protected area is processed from scratch via
#'   [wdpa_process()] and written as `LANDPA00`. Default `NULL` (skip).
#' @param states_nc Optional path to the LUH2 `states.nc` file. If supplied
#'   (along with `static_nc` and `cty_shp`), urban area is aggregated to
#'   IMPACT countries via [luh2_urban_area_cty()] and added to `LAND00` as
#'   `item = "urban"`. Default `NULL` (skip).
#' @param static_nc Path to the LUH2 `staticData_quarterdeg.nc` file.
#'   Required when `gdb_path` (without `pa_tif`) or `states_nc` is supplied.
#' @param cty_shp Path to the IMPACT regions shapefile. Required when
#'   `pa_tif`, `gdb_path`, or `states_nc` is supplied.
#' @param wdpa_year Integer. Only PAs designated by this year are included
#'   when processing from `gdb_path` (see [wdpa_process()]). Default `2021`.
#' @param urban_year Integer. Target year for the urban area aggregation.
#'   Defaults to `base_year` when `NULL`.
#' @param bii_year Integer. Year to compute the `natfor` primary-forest share
#'   for (see [luh2_bii_coeff_cty()]). Defaults to `base_year` when `NULL`.
#'   Ignored when `bii_pixel_csv` is supplied.
#' @param bii_pixel_csv Optional path to a pixel-derived BII coefficient CSV
#'   (see [luh2_bii_coeff_cty_pixel()]). If supplied, `BII_COEFF` is built
#'   from this imagery-derived table instead of the literature-based.
#'   use bii_pixel_csv <- system.file("extdata", "pixel_bii_coeff_fraction_r090.csv", package = "luh2impact")
#'   [luh2_bii_coeff_cty()]. Default `NULL` (use the literature-based
#'   coefficients).
#'
#' @return Invisibly returns a list with `results` (elasticity estimates)
#'   and `output_gdx` (path written).
#' @author Abhijeet Mishra, Claude Code
#' @export
luh2_landuse_elasticities <- function(fao_land_csv, sets_xlsx, ssp_gdx, fao_country_map_csv, wb_gdppc_csv,
                                      output_gdx, ssp = "SSP2", year_min = 1995, year_max = 2023,
                                      base_year = 2021, min_obs = 10,
                                      plot_dir = NULL, region_map_xlsx = NULL, region = "WLD",
                                      pa_tif = NULL, gdb_path = NULL,
                                      states_nc = NULL, static_nc = NULL, cty_shp = NULL,
                                      wdpa_year = 2021, urban_year = NULL, bii_year = NULL,
                                      bii_pixel_csv = NULL) {

    if ((!is.null(pa_tif) || !is.null(gdb_path) || !is.null(states_nc)) && is.null(cty_shp)) {
        stop("cty_shp is required when pa_tif, gdb_path, or states_nc is supplied")
    }
    if (is.null(pa_tif) && !is.null(gdb_path) && is.null(static_nc)) {
        stop("static_nc is required when gdb_path is supplied without pa_tif")
    }
    if (!is.null(states_nc) && is.null(static_nc)) {
        stop("static_nc is required when states_nc is supplied")
    }

    message("Loading FAOSTAT land use data...")
    faoland <- luh2_load_faoland(fao_land_csv, sets_xlsx, year_min = year_min)

    message("Computing GDP per capita by IMPACT country...")
    gdppc <- luh2_compute_gdppc(ssp_gdx, fao_country_map_csv, sets_xlsx, wb_gdppc_csv,
                                ssp = ssp, year_min = year_min, year_max = year_max)

    message("Fitting land-use elasticities...")
    fit <- luh2_fit_elasticities(faoland, gdppc, year_max = year_max, min_obs = min_obs)

    if (!is.null(plot_dir)) {
        message("Plotting elasticities...")
        luh2_plot_elasticities(fit$results, region_map_xlsx = region_map_xlsx, region = region, output_dir = plot_dir)
    }

    pa00 <- NULL

    if (!is.null(pa_tif)) {
        message("Aggregating pre-processed WDPA raster to IMPACT country level...")
        pa00 <- luh2_pa_area_cty(pa_tif, cty_shp)
    } else if (!is.null(gdb_path)) {
        message("Processing WDPA protected areas...")
        pa_cty <- wdpa_process(gdb_path, static_nc, output = "cty", cty_shp = cty_shp, year = wdpa_year)
        pa00 <- data.frame(cty = pa_cty$cty, value = pa_cty$pa_area_km2 * 100 / 1000)
    }

    extra_land00 <- NULL

    if (!is.null(states_nc)) {
        message("Computing urban area by IMPACT country...")
        urban_cty <- luh2_urban_area_cty(states_nc, static_nc, cty_shp,
                                         year = if (is.null(urban_year)) base_year else urban_year)
        extra_land00 <- rbind(extra_land00,
                              data.frame(cty = urban_cty$cty, item = "urban", value = urban_cty$value))
    }

    if (!is.null(bii_pixel_csv)) {
        message("Computing BII coefficients from pixel imagery (", bii_pixel_csv, ")...")
        bii_coeff <- luh2_bii_coeff_cty_pixel(bii_pixel_csv)
    } else {
        message("Computing BII coefficients...")
        bii_coeff <- luh2_bii_coeff_cty(fao_land_csv, sets_xlsx,
                                        year = if (is.null(bii_year)) base_year else bii_year)
    }

    message("Exporting GDX to ", output_gdx, "...")
    luh2_export_elasticities_gdx(fit$results, faoland, output_gdx, base_year = base_year,
                                 extra_land00 = extra_land00, pa00 = pa00, bii_coeff = bii_coeff)

    message("Done.")
    invisible(list(results = fit$results, output_gdx = output_gdx))
}
