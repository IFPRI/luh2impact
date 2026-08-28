utils::globalVariables(c(
    # luh2_build_pixels
    "cty", "crop_area", "crop_slope", "natfor_area", "natfor_slope",
    "other_area", "other_slope", "past_area", "plant_area", "plant_prox",
    # luh2_export_gdx
    "pid", "fland", "yrs", "value", "coord",
    "crop_area", "natfor_area", "past_area", "other_area", "plant_area",
    # luh2_write_tifs
    "lu_out", "pid_coords", "x", "y",
    # luh2_plot_results
    "pool", "area", "avail_kha", "share", "lyr", "year",
    "NEW_REGION",
    # luh2_crop_trend and pool_trend
    "time_idx",
    # variables
    "d_crop",
    "d_natfor",
    "d_past",
    "d_other",
    "d_plant",
    "crop", "natfor", "past", "other", "plant",
    "crop_tot", "natfor_tot", "past_tot", "other_tot", "plant_tot",
    "pool_total", "scale_down", "crop_area_tot", "natfor_area_tot", "past_area_tot",
    "other_area_tot", "plant_area_tot", "carea_kha", "urban_area",
    "crop_share", "natfor_share", "past_share", "other_share", "plant_share", "urban_share",
    "suit_crop", "suit_natfor", "suit_past", "suit_other",
    # luh2_load_faoland, luh2_compute_gdppc, luh2_fit_elasticities
    "item", "unit", "gdppc", "land", "dln_crop", "land_use",
    "elast_crop", "elast_gdppc",
    # luh2_plot_elasticities, luh2_export_elasticities_gdx
    "elasticity", "driver", "q1", "q3", "iqr", "is_outlier", "elast", "n_obs",
    "flag", "se_crop", "pval_crop", "se_gdppc", "pval_gdppc",
    # luh2_urban_area_cty
    "urban_area_km2",
    # luh2_pa_area_cty
    "pa_area_km2",
    # luh2_bii_coeff_cty
    "natfor_primary_share",
    # luh2_export_gdx (pa_area aggregation)
    "pa_area",
    # luh2_plot_bii
    "bii", "bii_base", "delta",
    # luh2_plot_results (dominant-pool RGB blend)
    "total", "r", "g", "b", "r_mix", "g_mix", "b_mix", "max_ch", "hex",
    # wdpa_process
    "pid_wdpa", "area_km2"
))
