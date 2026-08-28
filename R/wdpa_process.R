#' Process WDPA protected areas to pixel or country level
#'
#' Loads the WDPA geodatabase, filters to ecologically meaningful IUCN
#' categories and designated status, and either computes protected area
#' coverage in km^2 per LUH2 pixel, or aggregates to IMPACT country totals.
#'
#' IUCN categories retained (all except III):
#' \itemize{
#'   \item Ia - Strict Nature Reserve
#'   \item Ib - Wilderness Area
#'   \item II - National Park
#'   \item IV - Habitat/Species Management Area
#'   \item V  - Protected Landscape/Seascape
#'   \item VI - Protected Area with Sustainable Use
#' }
#' Category III (Natural Monument) is excluded as features are too small
#' for 0.25 degree resolution.
#'
#' @param gdb_path Path to the WDPA geodatabase (`.gdb` folder).
#' @param static_nc A `SpatRaster` used to match resolution, extent, and
#'   CRS - use static_quarterdeg.nc from LUh2.
#' @param output_tif Path to write the output GeoTIFF (pixel output only).
#' @param year Integer. Only PAs designated by this year (`STATUS_YR <= year`)
#'   are included. Default is 2021.
#' @param output Character. Either `"pixel"` (default) for a per-pixel PA area
#'   raster, or `"cty"` for a country-level summary data frame.
#' @param cty_shp Path to the IMPACT regions shapefile. Required when
#'   `output = "cty"`.
#'
#' @return If `output = "pixel"`: a single-layer `SpatRaster` of PA area in
#'   km^2 per pixel, also written to `output_tif`.
#'   If `output = "cty"`: a data frame with columns `cty` and `pa_area_km2`.
#' @export
wdpa_process <- function(gdb_path, static_nc, output_tif = NULL,
                         year = 2021, output = "pixel", cty_shp = NULL) {
    if (output == "cty" && is.null(cty_shp)) {
        stop("cty_shp must be provided when output = 'cty'")
    }

    iucn_keep   <- c("Ia", "Ib", "II", "IV", "V", "VI")
    status_keep <- c("Designated", "Inscribed", "Established", "Adopted")

    message("Loading WDPA polygons (write rds save here eventually so that i can load it easily later)...")
    wdpa <- terra::vect(gdb_path, layer = "WDPA_poly_Jul2026")

    message("Filtering...")
    wdpa_df  <- as.data.frame(wdpa)
    keep_idx <- which(
        wdpa_df$IUCN_CAT %in% iucn_keep &
            wdpa_df$REALM %in% c("Terrestrial", "Coastal") &
            wdpa_df$STATUS   %in% status_keep &
            !is.na(wdpa_df$STATUS_YR) &
            wdpa_df$STATUS_YR > 0 &
            wdpa_df$STATUS_YR <= year
    )
    message("PAs retained: ", length(keep_idx))
    wdpa <- wdpa[keep_idx, ]

    message("Projecting to OGC:CRS84...")
    wdpa <- terra::project(wdpa, "OGC:CRS84")

    # Compute actual PA area per pixel in km^2 - avoids binary presence/absence
    # problem where small PAs incorrectly fill entire 0.25 degree pixels (~700 km^2) so cannot use "touch" in rasterizing
    message("Computing PA area per pixel (km^2) - this one take quite some time [oh so fun..fml]")

    message("only process land data water we dont care for now ")
    static    <- terra::rast(static_nc)
    template_land <- static[["ccode"]] # only pickup countries with ccod we need to knwo land thats it
    template_land <- terra::mask(template_land, template_land == 0, maskvalue = TRUE)

    # Convert template to pixel polygons once
    template_unique <- template_land
    terra::values(template_unique) <- 1:(terra::ncell(template_unique))
    pixel_grid <- terra::as.polygons(template_unique)
    names(pixel_grid) <- "pid_wdpa"

    # Break rasterization ----

    out_dir    <- "C:/Local/Work/IFPRI/Landuse/WDPA/Processed Polygons"
    lat_breaks <- seq(-90, 90, by = 18)

    for (i in seq_len(length(lat_breaks) - 1)) {
        lat_min  <- lat_breaks[i]
        lat_max  <- lat_breaks[i + 1]
        out_file <- file.path(out_dir, paste0("pa_band_", lat_min, "_", lat_max, ".tif"))

        if (file.exists(out_file)) {
            message("Band ", lat_min, " to ", lat_max, " already exists, skipping")
            next
        }

        message("Processing lat band: ", lat_min, " to ", lat_max)

        ext_band   <- terra::ext(-180, 180, lat_min, lat_max)
        wdpa_chunk <- terra::crop(wdpa, ext_band)

        if (is.null(wdpa_chunk) || nrow(wdpa_chunk) == 0) {
            message("No PAs in this band, skipping")
            next
        }

        grid_chunk <- terra::crop(pixel_grid, ext_band)

        message("Intersecting ", nrow(wdpa_chunk), " PAs with ", nrow(grid_chunk), " pixels...")
        wdpa_chunk <- terra::intersect(wdpa_chunk, terra::as.polygons(terra::ext(grid_chunk)))

        isect <- tryCatch(
            terra::intersect(wdpa_chunk, grid_chunk),
            error = function(e) {
                message("Intersection failed: ", e$message)
                NULL
            }
        )

        if (is.null(isect) || nrow(isect) == 0) {
            message("No intersection in this band, skipping")
            next
        }

        # Compute area of each intersection piece in km^2
        isect$area_km2 <- terra::expanse(isect, unit = "km")

        # Aggregate to pixel level
        isect_df <- as.data.frame(isect)[, c("pid_wdpa", "area_km2")]
        pix_area <- isect_df |>
            dplyr::group_by(pid_wdpa) |>
            dplyr::summarise(pa_area_km2 = sum(area_km2, na.rm = TRUE)) |>
            dplyr::ungroup()

        # Join back to grid and rasterize
        grid_chunk$pa_area_km2 <- 0
        idx <- match(grid_chunk$pid_wdpa, pix_area$pid_wdpa)
        grid_chunk$pa_area_km2[!is.na(idx)] <- pix_area$pa_area_km2[idx[!is.na(idx)]]

        template_band <- terra::crop(template_land, ext_band)
        r <- terra::rasterize(grid_chunk, template_band, field = "pa_area_km2", background = 0)
        cell_size <- terra::cellSize(r, unit = "km")
        r <- min(r, cell_size)
        names(r) <- "pa_area_km2_2021"

        terra::writeRaster(r, out_file, overwrite = TRUE)
        message("Saved: ", basename(out_file))
        gc()
    }

    # Merge all bands
    # Merge all bands
    tifs    <- list.files(out_dir, pattern = "\\.tif$", full.names = TRUE)
    chunks  <- lapply(tifs, terra::rast)
    pa_rast <- do.call(terra::merge, chunks)
    pa_rast <- terra::classify(pa_rast, cbind(NA, 0))
    names(pa_rast) <- "pa_area_km2_2021"

    if (output == "pixel") {
        if (is.null(output_tif)) stop("output_tif must be provided when output = 'pixel'")
        message("Writing to ", output_tif)
        terra::writeRaster(pa_rast, output_tif, overwrite = TRUE)
        message("Done. Total PA area: ",
                round(sum(terra::values(pa_rast), na.rm = TRUE), 0), " km^2")
        return(pa_rast)
    }

    if (output == "cty") {
        message("Aggregating to IMPACT country level...")
        cty       <- terra::vect(cty_shp)
        cty$cty   <- cty$NEW_REGION
        cty       <- terra::project(cty, "OGC:CRS84")
        cty_rast  <- terra::rasterize(cty, template_land, field = "cty", touches = TRUE)

        df <- as.data.frame(c(cty_rast, pa_rast), na.rm = TRUE)
        names(df) <- c("cty", "pa_area_km2")

        cty_summary <- df |>
            dplyr::group_by(cty) |>
            dplyr::summarise(pa_area_km2 = sum(pa_area_km2, na.rm = TRUE)) |>
            dplyr::ungroup()

        message("Done. Countries with PA coverage: ",
                sum(cty_summary$pa_area_km2 > 0))
        return(cty_summary)
    }
}
