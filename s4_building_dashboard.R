# s4_building_dashboard.R
#
# Builds a static, self-contained local dashboard (dashboard/) for browsing,
# per species x model (base/anthro) x climate scenario (RCP4.5/RCP8.5):
#   - the tracking-slope plot: route trends by SDM class 1-7 with the species'
#     tracking slope (v2_species_tracking_slope.R)
#   - the range-change-category plot: route trends pooled into Decreasing
#     suitability / Neutral / Increasing suitability, with the category means and
#     the steps between them (v2_species_tracking_slope.R)
#   - the range-shift map
#   - the route-trend-distribution raincloud plot (v3_species_trend_distributions.R)
#   - a Decreasing suitability/Neutral/Increasing suitability table (n, mean, sd) for the selected
#     model + scenario
#
# Scope: only species with model output (those in
# output/species_routes_covariates/per_species_sdm/), not the full
# data/spp_names_codes_group_aou.csv list.
#
# No fetch()/server needed: manifest + table data are written as plain .js
# files (`const X = {...}`) loaded via <script src>, so dashboard/index.html
# works by double-clicking it, same as opening any other local file.
#
# Reads:
#   output/species_routes_covariates/per_species_sdm/*_base_route_trends_sdm.csv
#     (species list: code, name, group, slug -- one row read per file)
#   output/visualizations/species_tracking_slope/<model>/<rcp>/*.png
#     (tracking-slope plots, v2_species_tracking_slope.R)
#   output/visualizations/species_category_trend/<model>/<rcp>/*.png
#     (range-change-category plots, v2_species_tracking_slope.R; these replace the
#     old *_trend_by_group_*_cld.png violins of the former
#     4c_statistical_analysis_and_visualization_covariates.R)
#   output/visualizations/species_maps/<model>/<rcp>/*.png
#     (range-shift maps)
#   output/visualizations/species_trend_distributions/<model>/<rcp>/*.png
#     (raincloud plots; a species/model/scenario with <5 categorized routes
#     has no file here -- not an error, see v3_species_trend_distributions.R)
#   output/files/all_species_habitat_contraction_stable_expansion_2010_2025_{base,anthro}.xlsx
#     (Decreasing suitability/Neutral/Increasing suitability n/mean/sd table, both RCP sheets)
#
# Writes:
#   dashboard/data/manifest.js, dashboard/data/tables.js
#   dashboard/img/tracking_slope/<code>_<model>_<rcp>.webp (skipped where no source exists)
#   dashboard/img/violin_group/<code>_<model>_<rcp>.webp  (skipped where no source exists)
#   dashboard/img/maps/<code>_<model>_<rcp>.webp
#   dashboard/img/trend_dist/<code>_<model>_<rcp>.webp    (skipped where no source exists)
#
# Missing-source handling: maps are required (every species/model/scenario
# has one; a gap is a real error). Violins and trend-distributions can
# legitimately be absent for a given species/model/scenario (no routes in
# any Decreasing suitability/Neutral/Increasing suitability category, or fewer than v3's
# min_routes_per_plot) -- those are skipped during conversion, and
# dashboard/app.js shows a placeholder when the webp is missing rather than
# needing a precomputed "is this available" flag.
#
# Idempotent: re-running rebuilds all generated files.

library(here)
library(tidyverse)
library(openxlsx)
library(magick)
library(jsonlite)

here::i_am("s4_building_dashboard.R")

# Settings -------------------------------------------------------------------
models    <- c("base", "anthro")
scenarios <- c("rcp45", "rcp85")
scenario_sheet <- c(rcp45 = "RCP4.5", rcp85 = "RCP8.5")
grp_levels <- c("Decreasing suitability", "Neutral", "Increasing suitability")

violin_width  <- 700   # px, width of the tracking-slope and category plots
map_width     <- 1400  # px, maps are all fig_width=10in wide -> uniform scale factor
webp_quality  <- 82

sdm_dir        <- here::here("output", "species_routes_covariates", "per_species_sdm")
category_dir   <- here::here("output", "visualizations", "species_category_trend")
tracking_dir   <- here::here("output", "visualizations", "species_tracking_slope")
maps_src_dir   <- here::here("output", "visualizations", "species_maps")
trend_dist_dir <- here::here("output", "visualizations", "species_trend_distributions")
files_dir      <- here::here("output", "files")

dash_dir        <- here::here("dashboard")
data_dir        <- file.path(dash_dir, "data")
tracking_out    <- file.path(dash_dir, "img", "tracking_slope")
violin_grp_out  <- file.path(dash_dir, "img", "violin_group")
map_out         <- file.path(dash_dir, "img", "maps")
trend_dist_out  <- file.path(dash_dir, "img", "trend_dist")

# Wipe and recreate each image output dir so files from an old naming
# scheme don't linger.
old_violin_rcp <- file.path(dash_dir, "img", "violin_rcp")   # replaced by tracking_slope
for (d in c(old_violin_rcp, tracking_out, violin_grp_out, map_out, trend_dist_out)) {
  if (dir.exists(d)) unlink(d, recursive = TRUE)
}
for (d in c(data_dir, tracking_out, violin_grp_out, map_out, trend_dist_out)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ==========================================================================
# 1. Species master list -- one row per species, read straight off the
#    per_species_sdm CSVs (the "base" file; species/code/group are identical
#    in the "anthro" file for the same species).
# ==========================================================================
cat("=== Building species list from per_species_sdm CSVs ===\n")

sdm_files <- list.files(sdm_dir, pattern = "_base_route_trends_sdm\\.csv$", full.names = TRUE)
species_meta <- map_dfr(sdm_files, function(f) {
  d <- read.csv(f, nrows = 1, stringsAsFactors = FALSE)
  tibble(code = d$species_code, name = d$species, group = d$group,
         slug = sub("_base_route_trends_sdm\\.csv$", "", basename(f)))
}) %>% arrange(name)

cat("Species found:", nrow(species_meta), "\n")
stopifnot(!any(duplicated(species_meta$code)))
stopifnot(!any(duplicated(species_meta$slug)))

# Map filenames use the species' common name with apostrophes stripped and
# spaces -> underscores.
species_meta <- species_meta %>%
  mutate(map_slug = str_replace_all(str_remove_all(name, "'"), " ", "_"))

# ==========================================================================
# 2. Coverage check -- verify every expected source file exists before
#    copying anything, so a silent gap doesn't turn into a silently-broken
#    dashboard panel.
# ==========================================================================
cat("\n=== Checking source file coverage ===\n")

check_missing <- function(paths, label, required = TRUE) {
  missing <- paths[!file.exists(paths)]
  cat(label, "-- missing:", length(missing), "of", length(paths),
      if (required) "(required)" else "(legitimately optional -- see header note)", "\n")
  if (length(missing) > 0) print(head(missing, 10))
  missing
}

# Maps: every species/model/scenario has one -- any gap here is a real error.
map_paths <- expand_grid(i = seq_len(nrow(species_meta)), model = models, scenario = scenarios) %>%
  mutate(path = file.path(maps_src_dir, model, scenario,
                          paste0(species_meta$map_slug[i], "_", model, "_", scenario, "_2010_2025.png")))
missing_maps <- check_missing(map_paths$path, "Range-shift maps", required = TRUE)
stopifnot("Missing map files -- run v1_species_map.R first" = length(missing_maps) == 0)

# Violins and trend-distributions: can legitimately be absent for a given
# species/model/scenario (see header note) -- just reported, not asserted.
tracking_paths <- expand_grid(i = seq_len(nrow(species_meta)), model = models, scenario = scenarios) %>%
  mutate(path = file.path(tracking_dir, model, scenario,
                          paste0(species_meta$slug[i], "_", model, "_", scenario, "_2010_2025.png")))
invisible(check_missing(tracking_paths$path, "Tracking-slope plots", required = FALSE))

category_paths <- expand_grid(i = seq_len(nrow(species_meta)), model = models, scenario = scenarios) %>%
  mutate(path = file.path(category_dir, model, scenario,
                          paste0(species_meta$slug[i], "_", model, "_", scenario, "_2010_2025.png")))
invisible(check_missing(category_paths$path, "Range-change-category plots", required = FALSE))

trend_dist_paths <- expand_grid(i = seq_len(nrow(species_meta)), model = models, scenario = scenarios) %>%
  mutate(path = file.path(trend_dist_dir, model, scenario,
                          paste0(species_meta$slug[i], "_", model, "_", scenario, "_2010_2025.png")))
invisible(check_missing(trend_dist_paths$path, "Trend-distribution plots", required = FALSE))

# ==========================================================================
# 3. Table data from the two xlsx workbooks (RCP4.5 + RCP8.5 sheets each)
# ==========================================================================
cat("\n=== Reading xlsx tables ===\n")

read_grp_sheet <- function(model, scenario) {
  path <- file.path(files_dir, paste0("all_species_habitat_contraction_stable_expansion_2010_2025_", model, ".xlsx"))
  d <- openxlsx::read.xlsx(path, sheet = scenario_sheet[[scenario]])
  # The workbooks still carry the old "Contraction"/"Stable"/"Expansion"
  # headers -- relabel to grp_levels so "Contraction.n" becomes "Decreasing suitability.n".
  names(d) <- sub("^Stable\\.", "Neutral.",
                  sub("^Expansion\\.", "Increasing suitability.",
                      sub("^Contraction\\.", "Decreasing suitability.", names(d))))
  cat(" ", basename(path), scenario_sheet[[scenario]], "--", nrow(d), "rows\n")
  d
}

xlsx_raw <- expand_grid(model = models, scenario = scenarios) %>%
  mutate(data = map2(model, scenario, read_grp_sheet))

# tables[[code]][[model]][[scenario]] = list(Decreasing suitability=list(n,mean,sd)|NULL, Neutral=..., Increasing suitability=...)
# openxlsx::read.xlsx() turns header "Contraction n" into "Contraction.n";
# read_grp_sheet() relabels that to "Decreasing suitability.n".
build_species_row <- function(row) {
  out <- list()
  for (g in grp_levels) {
    n <- row[[paste0(g, ".n")]]
    if (is.null(n) || is.na(n)) {
      out[[g]] <- NULL
    } else {
      out[[g]] <- list(n = unname(n),
                        mean = unname(row[[paste0(g, ".mean")]]),
                        sd = unname(row[[paste0(g, ".sd")]]))
    }
  }
  out
}

tables <- setNames(vector("list", nrow(species_meta)), species_meta$code)
for (code in species_meta$code) tables[[code]] <- setNames(vector("list", length(models)), models)

missing_from_xlsx <- character(0)
for (r in seq_len(nrow(xlsx_raw))) {
  model <- xlsx_raw$model[r]; scenario <- xlsx_raw$scenario[r]
  d <- xlsx_raw$data[[r]]
  d_by_name <- split(d, d$Species)
  for (i in seq_len(nrow(species_meta))) {
    code <- species_meta$code[i]; name <- species_meta$name[i]
    if (is.null(tables[[code]][[model]])) tables[[code]][[model]] <- list()
    if (!name %in% names(d_by_name)) {
      missing_from_xlsx <- c(missing_from_xlsx, paste(name, model, scenario))
      tables[[code]][[model]][[scenario]] <- list()  # empty -> "no data" in the UI
      next
    }
    tables[[code]][[model]][[scenario]] <- build_species_row(d_by_name[[name]][1, ])
  }
}
cat("\nSpecies/model/scenario combos with no xlsx row (expected: RBME x2 models x2 scenarios):",
    length(missing_from_xlsx), "\n")
if (length(missing_from_xlsx) > 0) print(missing_from_xlsx)

# ==========================================================================
# 4. Copy + downsize images to WebP
# ==========================================================================
cat("\n=== Converting images (this takes a few minutes) ===\n")

convert_one <- function(src, dst, width) {
  img <- magick::image_read(src)
  img <- magick::image_resize(img, paste0(width, "x"))
  magick::image_write(img, dst, format = "webp", quality = webp_quality)
}

t0 <- Sys.time()
n_done <- 0
report_every <- 200

convert_if_exists <- function(src, dst, width) {
  if (!file.exists(src)) return(FALSE)
  convert_one(src, dst, width)
  TRUE
}

for (i in seq_len(nrow(species_meta))) {
  code <- species_meta$code[i]
  for (model in models) {
    for (scenario in scenarios) {
      src <- file.path(tracking_dir, model, scenario,
                       paste0(species_meta$slug[i], "_", model, "_", scenario, "_2010_2025.png"))
      convert_if_exists(src, file.path(tracking_out, paste0(code, "_", model, "_", scenario, ".webp")), violin_width)
      n_done <- n_done + 1

      # range-change-category plot; the dashboard folder keeps its old name (img/violin_group)
      src <- file.path(category_dir, model, scenario,
                       paste0(species_meta$slug[i], "_", model, "_", scenario, "_2010_2025.png"))
      convert_if_exists(src, file.path(violin_grp_out, paste0(code, "_", model, "_", scenario, ".webp")), violin_width)
      n_done <- n_done + 1

      src <- file.path(maps_src_dir, model, scenario,
                       paste0(species_meta$map_slug[i], "_", model, "_", scenario, "_2010_2025.png"))
      convert_one(src, file.path(map_out, paste0(code, "_", model, "_", scenario, ".webp")), map_width)
      n_done <- n_done + 1

      src <- file.path(trend_dist_dir, model, scenario,
                       paste0(species_meta$slug[i], "_", model, "_", scenario, "_2010_2025.png"))
      convert_if_exists(src, file.path(trend_dist_out, paste0(code, "_", model, "_", scenario, ".webp")), violin_width * 1.3)
      n_done <- n_done + 1
    }
  }
  if (i %% report_every == 0) {
    cat(sprintf("  %d/%d species (%s elapsed)\n", i, nrow(species_meta),
               format(Sys.time() - t0, digits = 3)))
  }
}
cat(sprintf("Converted images for %d species in %s\n", nrow(species_meta), format(Sys.time() - t0, digits = 3)))

# ==========================================================================
# 5. Write manifest.js + tables.js
# ==========================================================================
cat("\n=== Writing manifest.js + tables.js ===\n")

manifest <- species_meta %>%
  select(code, name, group) %>%
  arrange(name)

groups <- sort(unique(manifest$group))

manifest_json <- jsonlite::toJSON(
  list(species = manifest, groups = groups),
  auto_unbox = TRUE, pretty = FALSE, na = "null"
)
writeLines(paste0("const MANIFEST = ", manifest_json, ";"), file.path(data_dir, "manifest.js"))

tables_json <- jsonlite::toJSON(tables, auto_unbox = TRUE, pretty = FALSE, na = "null")
writeLines(paste0("const TABLES = ", tables_json, ";"), file.path(data_dir, "tables.js"))

cat("\n=== Dashboard data built ===\n")
cat("Species:", nrow(manifest), "\n")
cat("Open dashboard/index.html in a browser to view.\n")
