# v2_species_tracking_slope.R
#
# One "tracking slope" figure per species x model (base / anthro) x climate
# scenario (RCP 4.5 / RCP 8.5), in the style of 4j's concept figure. Replaces
# 4c's per-species violin plots.
#
#   - violins: route trends (posterior means, %/yr) in each SDM class 1-7
#     (classes with fewer than min_violin_n routes are drawn as points)
#   - black ticks: class means; numbers under the axis: routes per class
#   - line: the species' tracking slope (change in route trend per SDM class
#     step), drawn through (mean score, mean trend) of its routes, as in 4j
#
# Where the slope comes from:
#   - If stage 1 (4a) and stage 2 (4b) outputs exist for that model x scenario,
#     the line is the stage-1 posterior-mean slope, the shaded fan covers
#     slopes within its 95% credible interval, and the colour is the species'
#     class from the stage-2 Student-t model (Tracks / Inconclusive /
#     Contrary). Currently only anthro x rcp45 has these.
#   - Otherwise the line is the least-squares slope through the posterior-mean
#     route trends, drawn in grey with no fan. For anthro x rcp45 this matches
#     the stage-1 slope closely (r = 0.9996, median difference 0.006), but it
#     carries no uncertainty. Run 4a/4b for that model x scenario to get the
#     full version.
#
# Reads:  output/species_routes_covariates/per_species_sdm/<Species>_<model>_route_trends_sdm.csv  (3)
#         output/files/sdm_tracking/sdm_tracking_stage1_species_<model>_<scenario>_<years>.csv     (4a, optional)
#         output/files/sdm_tracking/sdm_tracking_stage2_species_<model>_<scenario>_<years>.csv     (4b, optional)
# Writes: output/species_tracking_slope/<model>/<scenario>/<Species>_<model>_<scenario>_<years>.png
#         output/species_tracking_slope/manifest_<years>.csv  (one row per figure attempted)
#
# Resumable like v1: a figure that already exists is skipped unless
# overwrite = TRUE.

library(here)
library(tidyverse)

here::i_am("v2_species_tracking_slope.R")

# Settings -------------------------------------------------------------------
firstYear <- 2010
lastYear  <- 2025
models    <- c("base", "anthro")
scenarios <- c(rcp45 = "RCP 4.5", rcp85 = "RCP 8.5")

only_species <- NULL     # e.g. c("Eastern_Towhee") (file-name slugs) for a test run
overwrite    <- FALSE
require_route_converged <- TRUE   # same rule as 4c / v1

min_violin_n <- 5        # classes with fewer routes are drawn as points, not a violin
fig_width    <- 7
fig_height   <- 5.6
fig_dpi      <- 200

class_colours <- c(Tracks = "#2a78d6", Inconclusive = "#4c5852", Contrary = "#d9423f",
                   `No credible interval` = "#8a8a8a")
class_labels  <- c("1\nExtirpation", "2\nWorsening", "3\nSlightly\nworsening", "4\nNeutral",
                   "5\nSlightly\nimproving", "6\nImproving", "7\nColonization")

sdm_dir   <- here::here("output", "species_routes_covariates", "per_species_sdm")
track_dir <- here::here("output", "files", "sdm_tracking")
out_dir   <- here::here("output", "species_tracking_slope")

# Stage-1 / stage-2 results, where available -------------------------------------
read_stage <- function(stage, model, scenario) {
  f <- file.path(track_dir, paste0("sdm_tracking_", stage, "_species_", model, "_", scenario, "_",
                                   firstYear, "_", lastYear, ".csv"))
  if (file.exists(f)) read_csv(f, show_col_types = FALSE) else NULL
}

bayes <- list()
for (m in models) for (sc in names(scenarios)) {
  s1 <- read_stage("stage1", m, sc)
  s2 <- read_stage("stage2", m, sc)
  if (!is.null(s1) && !is.null(s2)) {
    bayes[[paste(m, sc)]] <- s1 %>%
      filter(status == "ok") %>%
      select(species, b = slope_mean, lo = slope_q025, hi = slope_q975) %>%
      left_join(s2 %>% select(species, label = class_studentt), by = "species")
  }
}
cat("Stage-1/2 results available for:", if (length(bayes)) paste(names(bayes), collapse = ", ") else "none", "\n")

# Plot one species x model x scenario ------------------------------------------------
line_at <- function(ybar, xbar, b, cls) ybar + b * ((cls - 4) - xbar)
fmt <- function(v) sprintf("%+.2f", v)

plot_one <- function(routes, sp_name, model, scenario_label, sl) {
  class_stats <- routes %>% group_by(cls) %>% summarise(n = n(), mean = mean(trend), .groups = "drop")
  big <- class_stats$cls[class_stats$n >= min_violin_n]
  xbar <- mean(routes$score); ybar <- mean(routes$trend)

  has_line <- !is.null(sl)
  has_fan  <- has_line && !is.na(sl$lo)
  if (has_line) {
    line_df <- tibble(cls = c(1, 7), y = line_at(ybar, xbar, sl$b, cls))
    fan_df  <- if (has_fan) tibble(cls = c(1, 7, 7, 1),
                                   y = line_at(ybar, xbar, c(sl$lo, sl$lo, sl$hi, sl$hi), cls)) else NULL
  }

  # y window: central 98% of routes plus class means and the line, so a few
  # extreme routes don't flatten everything else
  ys <- c(quantile(routes$trend, c(0.01, 0.99), names = FALSE), class_stats$mean, 0,
          if (has_line) line_df$y, if (has_fan) fan_df$y)
  y_lim <- range(ys); y_lim <- y_lim + c(-1, 1) * 0.06 * diff(y_lim)
  n_cut <- sum(routes$trend < y_lim[1] | routes$trend > y_lim[2])
  strip_y <- y_lim[1] - 0.10 * diff(y_lim)

  title2 <- if (!has_line) "fewer than 2 SDM classes · no slope" else if (has_fan)
    paste0("slope ", fmt(sl$b), " (", fmt(sl$lo), " to ", fmt(sl$hi), ") · ", sl$label) else
    paste0("slope ", fmt(sl$b), " · posterior-mean fit, no credible interval")
  colour <- if (!has_line) NA else class_colours[[sl$label]]

  p <- ggplot() +
    geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.5) +
    geom_violin(data = routes %>% filter(cls %in% big), aes(x = cls, y = trend, group = cls),
                fill = "grey70", colour = NA, alpha = 0.55, scale = "width", width = 0.8, trim = TRUE) +
    geom_point(data = routes %>% filter(!cls %in% big), aes(x = cls, y = trend), colour = "grey55", size = 1.2)
  if (has_fan)  p <- p + geom_polygon(data = fan_df, aes(x = cls, y = y), fill = colour, alpha = 0.18)
  if (has_line) p <- p + geom_line(data = line_df, aes(x = cls, y = y), colour = colour,
                                   linewidth = 1.2, lineend = "round")
  p +
    geom_segment(data = class_stats, aes(x = cls - 0.22, xend = cls + 0.22, y = mean, yend = mean),
                 linewidth = 0.9, colour = "black", lineend = "round") +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = y_lim[1] - 0.5 * diff(y_lim),
             ymax = y_lim[1] - 0.03 * diff(y_lim), fill = "white", colour = "white") +
    geom_text(data = class_stats, aes(x = cls, y = strip_y, label = n), size = 3.4, colour = "grey30") +
    scale_x_continuous(breaks = 1:7, labels = class_labels, limits = c(0.5, 7.5), expand = c(0, 0)) +
    scale_y_continuous(breaks = function(l) { b <- scales::breaks_extended(5)(y_lim); b[b >= y_lim[1] & b <= y_lim[2]] }) +
    coord_cartesian(ylim = c(y_lim[1] - 0.15 * diff(y_lim), y_lim[2]), expand = FALSE) +
    labs(title = sp_name, subtitle = title2,
         x = paste0("SDM class (", scenario_label, ")"), y = "Route trend (% per year)",
         caption = paste0("Model: ", model, ". Violins: route trends (posterior means); black ticks: class means; ",
                          "numbers: routes per class.\n",
                          if (has_fan) "Line: tracking slope (stage 1); shaded: slopes within its 95% credible interval; colour: stage-2 class.\n"
                          else if (has_line) "Line: least-squares slope through the posterior-mean route trends.\n" else "",
                          if (n_cut > 0) paste0(n_cut, " route(s) outside the plotted range are not shown.") else "")) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          plot.title = element_text(face = "bold", size = 15),
          plot.subtitle = element_text(face = "bold", size = 12, colour = if (has_line) colour else "grey35"),
          plot.caption = element_text(hjust = 0, size = 7.5, colour = "grey30"),
          axis.title = element_text(face = "bold", size = 12),
          axis.text.x = element_text(size = 8.5, colour = "black", lineheight = 0.9),
          axis.text.y = element_text(size = 11))
}

# Main loop ---------------------------------------------------------------------
files <- list.files(sdm_dir, pattern = "_route_trends_sdm\\.csv$", full.names = TRUE)
files <- files[str_detect(basename(files), paste0("_(", paste(models, collapse = "|"), ")_route_trends_sdm\\.csv$"))]
if (length(files) == 0) stop("No *_route_trends_sdm.csv files in ", sdm_dir, " -- run 3_add_SDM_covariates.R first.")

manifest <- list()
n_done <- 0
for (f in files) {
  model <- str_match(basename(f), "_(base|anthro)_route_trends_sdm\\.csv$")[, 2]
  slug  <- sub(paste0("_", model, "_route_trends_sdm\\.csv$"), "", basename(f))
  if (!is.null(only_species) && !slug %in% only_species) next

  d <- NULL
  for (sc in names(scenarios)) {
    out_file <- file.path(out_dir, model, sc, paste0(slug, "_", model, "_", sc, "_", firstYear, "_", lastYear, ".png"))
    if (file.exists(out_file) && !overwrite) next
    if (is.null(d)) {
      d <- read_csv(f, show_col_types = FALSE)
      if (require_route_converged) d <- d %>% filter(route_converged %in% TRUE)
    }
    sp_name <- d$species[1]
    routes <- d %>%
      filter(.data[[sc]] %in% 1:7, !is.na(trend)) %>%
      transmute(cls = as.integer(.data[[sc]]), score = cls - 4, trend)

    status <- "ok"
    if (nrow(routes) == 0) {
      status <- "no_classified_routes"
    } else {
      sl <- NULL
      key <- paste(model, sc)
      if (key %in% names(bayes) && sp_name %in% bayes[[key]]$species) {
        sl <- bayes[[key]] %>% filter(species == sp_name) %>% slice(1)
      } else if (n_distinct(routes$score) >= 2) {
        x <- routes$score
        sl <- tibble(b = sum((x - mean(x)) * routes$trend) / sum((x - mean(x))^2),
                     lo = NA_real_, hi = NA_real_, label = "No credible interval")
        status <- "ok_posterior_mean_fit"
      } else {
        status <- "ok_single_class"
      }
      p <- plot_one(routes, sp_name, model, scenarios[[sc]], sl)
      dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
      ggsave(out_file, p, width = fig_width, height = fig_height, dpi = fig_dpi,
             device = ragg::agg_png, bg = "white")
      n_done <- n_done + 1
      if (n_done %% 100 == 0) cat("  ", n_done, "figures written\n")
    }
    manifest[[length(manifest) + 1]] <- tibble(species = sp_name, slug = slug, model = model,
                                               scenario = sc, n_routes = nrow(routes), status = status)
  }
}

# Summary ------------------------------------------------------------------------
if (length(manifest) > 0) {
  manifest_df <- bind_rows(manifest)
  write_csv(manifest_df, file.path(out_dir, paste0("manifest_", firstYear, "_", lastYear, ".csv")))
  cat("\n=== Tracking-slope figures done ===\n")
  print(manifest_df %>% count(model, scenario, status))
  cat("Written to:", out_dir, "\n")
} else {
  cat("\nNothing to do -- every figure already exists (set overwrite <- TRUE to redraw).\n")
}
