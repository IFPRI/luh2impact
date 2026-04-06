#' Load and clean LUH2 states NetCDF
#'
#' Reads the LUH2 `states.nc` file and removes non-share layers (`secma`, `secmb`).
#'
#' @param states_nc Path to the LUH2 `states.nc` file.
#'
#' @return A `SpatRaster` with all land use state share layers.
#' @author Abhijeet Mishra, Claude Code
#' @export
luh2_load <- function(states_nc) {
    luh <- terra::rast(states_nc)
    drop <- grep("secmb|secma", names(luh), value = TRUE)
    luh <- luh[[!names(luh) %in% drop]]
    luh
}
