#' Write GeoTIFFs from LUMEN solution GDX
#'
#' Reads the LUMEN solution GDX and writes one multi-layer GeoTIFF per land
#' pool, with one layer per year. Output files are named `lu_<pool>.tif` and
#' written to `output_dir`.
#'
#' @param gdx_path Path to the LUMEN solution GDX file
#'   (e.g. `"outputs/CHM/solution_lu.gdx"`).
#' @param output_dir Path to the directory where GeoTIFFs will be written.
#'   Also used to derive the country name via `basename(output_dir)`.
#'
#' @return Invisibly returns a named list of `SpatRaster` objects, one per pool.
#' @author Abhijeet Mishra, Claude Code
#' @export
luh2_write_tifs <- function(gdx_path, output_dir) {
    m_out <- gamstransfer::Container$new(gdx_path)

    lu_out <- m_out$getSymbols("lu_out")[[1]]$records

    pid_coords <- m_out["pid_coords"]$records |>
        tidyr::pivot_wider(names_from = "coord", values_from = "value")

    lu_out <- lu_out |>
        dplyr::left_join(pid_coords[, c("pid", "x", "y")], by = "pid")

    template <- terra::rast(
        resolution = 0.25,
        extent     = terra::ext(-180, 180, -90, 90),
        crs        = "OGC:CRS84"
    )

    pools   <- c("crop", "natfor", "past", "other", "plant")
    yrs_seq <- as.character(sort(unique(lu_out$yrs)))

    result <- list()

    for (pool in pools) {
        pool_stack <- terra::rast(lapply(yrs_seq, function(yr) {
            sub <- lu_out |>
                dplyr::filter(fland == pool, yrs == yr) |>
                dplyr::select(x, y, value)

            # Build raster directly from xyz coordinates — avoids point snapping
            # issues that cause NA holes when rasterizing from SpatVector points
            r <- terra::rast(sub, type = "xyz", crs = "OGC:CRS84")

            # Resample to standard LUH2 template to ensure consistent grid alignment
            terra::resample(r, template, method = "near")
        }))

        names(pool_stack) <- yrs_seq
        out_path <- file.path(output_dir, paste0("lu_", pool, ".tif"))
        terra::writeRaster(pool_stack, out_path, overwrite = TRUE)
        message("Written: ", basename(out_path))
        result[[pool]] <- pool_stack
    }

    invisible(result)
}
