#' Run the full LUH2-to-GDX pipeline
#'
#' Convenience wrapper that calls all processing steps in sequence and writes
#' the final GDX file. Each intermediate step can be run and tested individually
#' using the underlying functions.
#'
#' @param states_nc Path to the LUH2 `states.nc` file.
#' @param static_nc Path to the LUH2 `staticData_quarterdeg.nc` file.
#' @param planted_tif Path to the planted forest GeoTIFF (e.g. `sdpt_global.tif`).
#' @param cty_shp Path to the IMPACT regions shapefile.
#' @param impact_gdx Path to the IMPACT scenario GDX file.
#' @param output_gdx Path to write the output GDX file.
#' @param year Integer. Target baseline year to extract from LUH2 (e.g. 2015).
#' @param year_start Integer. Start year for trend estimation (e.g. 1990).
#' @param year_end Integer. End year for trend estimation (e.g. 2015).
#'
#' @return Invisibly returns the path to the written GDX file.
#' @author Abhijeet Mishra, Claude Code
#' @importFrom dplyr filter select mutate group_by ungroup summarise left_join across pull starts_with ends_with
#' @importFrom tidyr pivot_wider
#' @export
luh2gdx <- function(states_nc,
                    static_nc,
                    planted_tif,
                    cty_shp,
                    impact_gdx,
                    output_gdx,
                    year = 2015,
                    year_start = 1990,
                    year_end = 2015) {

    message("Loading LUH2 states...")
    luh <- luh2_load(states_nc)

    message("Computing crop trend (", year_start, "-", year_end, ")...")
    crop_trend <- luh2_pool_trend(luh,
                                  pool_vars  = c("c3ann", "c4ann", "c3per", "c4per", "c3nfx"),
                                  year_start = year_start,
                                  year_end   = year_end)

    message("Computing natural forest trend (", year_start, "-", year_end, ")...")
    natfor_trend <- luh2_pool_trend(luh,
                                    pool_vars  = c("primf", "primn", "secdf"),
                                    year_start = year_start,
                                    year_end   = year_end)

    message("Computing other land trend (", year_start, "-", year_end, ")...")
    other_trend <- luh2_pool_trend(luh,
                                   pool_vars  = c("secdn"),
                                   year_start = year_start,
                                   year_end   = year_end)

    message("Extracting year ", year, "...")
    extracted <- luh2_extract_year(luh, static_nc, year)

    message("Merging planted forest...")
    shares <- luh2_merge_planted(luarea = extracted$luarea,
                                 cellarea = extracted$cellarea,
                                 planted_tif = planted_tif,
                                 fstnf = extracted$fstnf,
                                 cty_shp = cty_shp)

    message("Reading IMPACT land use data...")
    landx0 <- DOORMAT::readGDX(gdx = impact_gdx, name = "LANDX0", quick_df = TRUE)$data

    message("Building pixel data frame...")
    df <- luh2_build_pixels(
        luh          = luh,
        shares       = shares,
        crop_trend   = crop_trend,
        natfor_trend = natfor_trend,
        other_trend  = other_trend,
        cellarea     = extracted$cellarea,
        icwtr        = extracted$icwtr,
        cty_shp      = cty_shp,
        landx0       = landx0
    )

    message("Exporting GDX to ", output_gdx, "...")
    luh2_export_gdx(df, landx0, output_gdx)

    message("Done.")
    invisible(output_gdx)
}
