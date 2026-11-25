# v2_species_tracking_slope.R
#
# Two figure sets per species x model (base / anthro) x climate scenario (RCP 4.5 /
# RCP 8.5), styled like Figure 3 (vf3_sdm_tracking_slope_concept_figure.R):
#
#   1. Tracking slope (as Figure 3, top row):
#      - violins: route trends (posterior means, %/yr) in each SDM class 1-7
#        (classes with fewer than min_violin_n routes are drawn as points)
#      - black ticks: class means; numbers above the violins: routes per class
#      - line: the species' tracking slope (change in route trend per SDM class
#        step), drawn through (mean score, mean trend) of its routes
#
#   2. Range-change categories (as Figure 3, bottom row):
#      - violins: route trends pooled into Decreasing suitability (classes 1-3),
#        Neutral (4) and Increasing suitability (5-7); numbers above: routes
#      - black ticks (and bars): category means (with 95% credible intervals)
#      - lines: the two steps between categories; the subtitle reports them relative to
#        Neutral, as in 4c: Decreasing - Neutral and Increasing - Neutral. Line colour
#        follows the line's direction (a rise from Decreasing to Neutral is blue). No caption.
#
# Where the slope, category means and colours come from:
#   - If stage 1 (4b) and stage 2 (4c) outputs exist for that model x scenario:
#     slope = stage-1 posterior-mean slope, fan = slopes within its 95% credible
#     interval, colour = the species' class from the stage-2 Student-t model
#     (Tracks / Inconclusive / Contrary). Category means and their 95% intervals
#     are the stage-1 ones; each step is coloured by its own P(> 0) with the same
#     class_threshold as 4c, using the normal approximation of vf3 / 4c
#     (difference of the two means, sd = sqrt(sd1^2 + sd2^2)), and coloured by the
#     line's direction: blue if it rises with P >= class_threshold. Currently only
#     anthro x rcp45 has these.
#   - Otherwise: the slope is the least-squares slope through the posterior-mean
#     route trends and the category means are the means of those route trends,
#     all drawn in grey with no interval. For anthro x rcp45 the slope matches the
#     stage-1 slope closely (r = 0.9996, median difference 0.006), but it carries
#     no uncertainty. Run 4b/4c for that model x scenario to get the full version.
#
# Reads:  output/species_routes_covariates/per_species_sdm/<Species>_<model>_route_trends_sdm.csv  (3)
#         output/sdm_tracking/sdm_tracking_stage1_species_<model>_<scenario>_<years>.csv     (4b, optional)
#         output/sdm_tracking/sdm_tracking_stage2_species_<model>_<scenario>_<years>.csv     (4c, optional;
#           the all-species file, so non-native species also get a class)
# Writes: output/visualizations/species_tracking_slope/<model>/<scenario>/<Species>_<model>_<scenario>_<years>.png
#         output/visualizations/species_tracking_slope/manifest_<years>.csv        (one row per slope figure attempted)
#         output/visualizations/species_category_trend/<model>/<scenario>/<Species>_<model>_<scenario>_<years>.png
#         output/visualizations/species_category_trend/manifest_<years>.csv        (one row per category figure attempted)
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
require_route_converged <- TRUE   # same rule as v1
class_threshold <- 0.90           # same as 4c / vf3

min_violin_n <- 5        # classes / categories with fewer routes are drawn as points, not a violin
fig_width    <- 7        # slope figures
cat_width    <- 6        # category figures
fig_height   <- 5.8
fig_dpi      <- 200

# same colours as Figure 3 (vf3); grey for fits without a credible interval
class_colours <- c(Tracks = "#2a78d6", Inconclusive = "#2b2b2b", Contrary = "#d9423f",
                   `No credible interval` = "#8a8a8a")
class_labels  <- c("1\nExtirpation", "2\nWorsening", "3\nSlightly\nworsening", "4\nNeutral",
                   "5\nSlightly\nimproving", "6\nImproving", "7\nColonization")
categories    <- tibble(cat   = c("decreasing", "neutral", "increasing"),
                        x     = 1:3,
                        label = c("Decreasing\nsuitability\n(1–3)", "Neutral\n(4)", "Increasing\nsuitability\n(5–7)"))

sdm_dir     <- here::here("output", "species_routes_covariates", "per_species_sdm")
track_dir   <- here::here("output", "sdm_tracking")
out_dir     <- here::here("output", "visualizations", "species_tracking_slope")
out_dir_cat <- here::here("output", "visualizations", "species_category_trend")

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
      select(species, b = slope_mean, lo = slope_q025, hi = slope_q975,
             starts_with("mean_decreasing_"), starts_with("mean_neutral_"), starts_with("mean_increasing_")) %>%
      left_join(s2 %>% select(species, label = class_studentt), by = "species")
  }
}
cat("Stage-1/2 results available for:", if (length(bayes)) paste(names(bayes), collapse = ", ") else "none", "\n")

# Shared helpers ------------------------------------------------------------------
line_at  <- function(ybar, xbar, b, cls) ybar + b * ((cls - 4) - xbar)
fmt      <- function(v) sprintf("%+.2f", v)
ci_txt   <- function(m, lo, hi) paste0(fmt(m), " (", fmt(lo), " to ", fmt(hi), ")")
classify <- function(p_pos) case_when(p_pos >= class_threshold     ~ "Tracks",
                                      1 - p_pos >= class_threshold ~ "Contrary",
                                      TRUE                         ~ "Inconclusive")

# y window: central 98% of routes plus the summary values, so a few extreme routes
# don't flatten everything else
window_of <- function(trend, extra) {
  y_lim <- range(c(quantile(trend, c(0.01, 0.99), names = FALSE), extra, 0), na.rm = TRUE)
  y_lim + c(-1, 1) * 0.06 * diff(y_lim)
}

# Pieces both figure types share: violins / points, the white strip with route counts
# above the violins (as in Figure 3), y scale, coord and theme.
add_base <- function(p, routes, stats, y_lim) {
  big <- stats$x[stats$n >= min_violin_n]
  p +
    geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.5) +
    geom_violin(data = routes %>% filter(x %in% big), aes(x = x, y = trend, group = x),
                fill = "grey70", colour = NA, alpha = 0.55, scale = "width", width = 0.8, trim = TRUE) +
    geom_point(data = routes %>% filter(!x %in% big), aes(x = x, y = trend), colour = "grey55", size = 1.2)
}
add_finish <- function(p, stats, y_lim, subtitle_colour) {
  p +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = y_lim[2] + 0.03 * diff(y_lim), ymax = Inf,
             fill = "white", colour = "white") +
    geom_text(data = stats, aes(x = x, y = y_lim[2] + 0.085 * diff(y_lim),
                                label = format(n, big.mark = ",", trim = TRUE)),
              size = 3.4, colour = "grey30") +
    scale_y_continuous(breaks = function(l) { b <- scales::breaks_extended(5)(y_lim); b[b >= y_lim[1] & b <= y_lim[2]] }) +
    coord_cartesian(ylim = c(y_lim[1] - 0.02 * diff(y_lim), y_lim[2] + 0.13 * diff(y_lim)), expand = FALSE) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          plot.title = element_text(face = "bold", size = 15),
          plot.subtitle = element_text(face = "bold", size = 11, colour = subtitle_colour, lineheight = 1.1),
          plot.caption = element_text(hjust = 0, size = 7.5, colour = "grey30"),
          axis.title = element_text(face = "bold", size = 12),
          axis.text.x = element_text(size = 8.5, colour = "black", lineheight = 0.9),
          axis.text.y = element_text(size = 11))
}

# Figure 1: tracking slope ---------------------------------------------------------
plot_slope <- function(routes, sp_name, model, scenario_label, sl) {
  routes <- routes %>% mutate(x = cls)
  stats  <- routes %>% group_by(x) %>% summarise(n = n(), mean = mean(trend), .groups = "drop")
  xbar <- mean(routes$score); ybar <- mean(routes$trend)

  has_line <- !is.null(sl)
  has_fan  <- has_line && !is.na(sl$lo)
  if (has_line) {
    line_df <- tibble(x = c(1, 7), y = line_at(ybar, xbar, sl$b, x))
    fan_df  <- if (has_fan) tibble(x = c(1, 7, 7, 1),
                                   y = line_at(ybar, xbar, c(sl$lo, sl$lo, sl$hi, sl$hi), x)) else NULL
  }
  y_lim <- window_of(routes$trend, c(stats$mean, if (has_line) line_df$y, if (has_fan) fan_df$y))
  n_cut <- sum(routes$trend < y_lim[1] | routes$trend > y_lim[2])

  subtitle <- if (!has_line) "fewer than 2 SDM classes · no slope" else if (has_fan)
    paste0("slope ", ci_txt(sl$b, sl$lo, sl$hi), " · ", sl$label) else
    paste0("slope ", fmt(sl$b), " · posterior-mean fit, no credible interval")
  colour <- if (!has_line) "grey35" else class_colours[[sl$label]]

  p <- add_base(ggplot(), routes, stats, y_lim)
  if (has_fan)  p <- p + geom_polygon(data = fan_df, aes(x = x, y = y), fill = colour, alpha = 0.18)
  if (has_line) p <- p + geom_line(data = line_df, aes(x = x, y = y), colour = colour,
                                   linewidth = 1.2, lineend = "round")
  p <- p +
    geom_segment(data = stats, aes(x = x - 0.22, xend = x + 0.22, y = mean, yend = mean),
                 linewidth = 0.9, colour = "black", lineend = "round") +
    scale_x_continuous(breaks = 1:7, labels = class_labels, limits = c(0.5, 7.5), expand = c(0, 0)) +
    labs(title = sp_name, subtitle = subtitle,
         x = paste0("SDM class (", scenario_label, ")"), y = "Route trend (% per year)",
         caption = paste0("Model: ", model, ". Violins: route trends (posterior means); black ticks: class means; ",
                          "numbers above the violins: routes per class.\n",
                          if (has_fan) "Line: tracking slope (stage 1); shaded: slopes within its 95% credible interval; colour: stage-2 class.\n"
                          else if (has_line) "Line: least-squares slope through the posterior-mean route trends.\n" else "",
                          if (n_cut > 0) paste0(n_cut, " route(s) outside the plotted range are not shown.") else ""))
  add_finish(p, stats, y_lim, colour)
}

# Figure 2: range-change categories ------------------------------------------------
# cm: one row per category present, with mean and (if Bayesian) q025 / q975 / sd
plot_category <- function(routes, sp_name, model, scenario_label, cm, bayesian) {
  routes <- routes %>%
    mutate(cat = case_when(cls <= 3 ~ "decreasing", cls == 4 ~ "neutral", TRUE ~ "increasing")) %>%
    left_join(categories %>% select(cat, x), by = "cat")
  stats <- routes %>% group_by(x) %>% summarise(n = n(), .groups = "drop")
  cm <- cm %>% left_join(categories %>% select(cat, x), by = "cat")

  # the two steps (only where both categories exist)
  # both differences are reported relative to Neutral (as in 4c): Decreasing - Neutral and
  # Increasing - Neutral. The line colour follows the line's direction left to right
  # (d = right - left), so a rise from Decreasing to Neutral (Decreasing - Neutral < 0) is blue.
  step_def <- tibble(from = c("decreasing", "neutral"), to = c("neutral", "increasing"),
                     name = c("Decreasing − Neutral", "Increasing − Neutral"),
                     sign = c(-1, 1))                  # reported value = sign * d
  steps <- step_def %>%
    inner_join(cm %>% select(from = cat, x0 = x, y0 = mean, sd0 = sd), by = "from") %>%
    inner_join(cm %>% select(to = cat,   x1 = x, y1 = mean, sd1 = sd), by = "to") %>%
    mutate(d = y1 - y0, v = sign * d)
  if (bayesian) {
    steps <- steps %>% mutate(s = sqrt(sd0^2 + sd1^2), p = pnorm(d / s), line_class = classify(p),
                              txt = paste0(name, "  ", ci_txt(v, v - 1.96 * s, v + 1.96 * s)))
  } else {
    steps <- steps %>% mutate(line_class = "No credible interval", txt = paste0(name, "  ", fmt(v)))
  }
  steps <- steps %>% mutate(dashed = line_class == "Inconclusive")

  y_lim <- window_of(routes$trend, c(cm$mean, if (bayesian) c(cm$q025, cm$q975)))

  subtitle <- if (nrow(steps) == 0) "fewer than 2 adjacent categories · no step" else
    paste0(paste(steps$txt, collapse = "\n"), if (!bayesian) "\nposterior-mean fit, no credible interval" else "")
  # one subtitle colour for two steps: the steps' class if they agree, else charcoal
  sub_colour <- if (!bayesian) "grey35" else if (n_distinct(steps$line_class) == 1)
    class_colours[[steps$line_class[1]]] else class_colours[["Inconclusive"]]

  p <- add_base(ggplot(), routes, stats, y_lim)
  if (nrow(steps) > 0) {
    p <- p + geom_segment(data = steps, aes(x = x0, xend = x1, y = y0, yend = y1, colour = line_class, linetype = dashed),
                          linewidth = 1.2, lineend = "round")
  }
  if (bayesian) {
    p <- p + geom_errorbar(data = cm, aes(x = x, ymin = q025, ymax = q975), width = 0, linewidth = 0.7, colour = "black")
  }
  p <- p +
    geom_segment(data = cm, aes(x = x - 0.22, xend = x + 0.22, y = mean, yend = mean),
                 linewidth = 0.9, colour = "black", lineend = "round") +
    scale_colour_manual(values = class_colours, guide = "none") +
    scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22"), guide = "none") +
    scale_x_continuous(breaks = categories$x, labels = categories$label, limits = c(0.5, 3.5), expand = c(0, 0)) +
    labs(title = sp_name, subtitle = subtitle,
         x = paste0("Range-change category (", scenario_label, ")"), y = "Route trend (% per year)")   # no caption
  add_finish(p, stats, y_lim, sub_colour)
}

save_fig <- function(p, out_file, width) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_file, p, width = width, height = fig_height, dpi = fig_dpi, device = ragg::agg_png, bg = "white")
}

# Main loop ---------------------------------------------------------------------
files <- list.files(sdm_dir, pattern = "_route_trends_sdm\\.csv$", full.names = TRUE)
files <- files[str_detect(basename(files), paste0("_(", paste(models, collapse = "|"), ")_route_trends_sdm\\.csv$"))]
if (length(files) == 0) stop("No *_route_trends_sdm.csv files in ", sdm_dir, " -- run 3_add_SDM_covariates.R first.")

manifest <- list(); manifest_cat <- list()
n_done <- 0
for (f in files) {
  model <- str_match(basename(f), "_(base|anthro)_route_trends_sdm\\.csv$")[, 2]
  slug  <- sub(paste0("_", model, "_route_trends_sdm\\.csv$"), "", basename(f))
  if (!is.null(only_species) && !slug %in% only_species) next

  d <- NULL
  for (sc in names(scenarios)) {
    stem     <- paste0(slug, "_", model, "_", sc, "_", firstYear, "_", lastYear, ".png")
    out_file <- file.path(out_dir,     model, sc, stem)
    out_cat  <- file.path(out_dir_cat, model, sc, stem)
    do_slope <- overwrite || !file.exists(out_file)
    do_cat   <- overwrite || !file.exists(out_cat)
    if (!do_slope && !do_cat) next
    if (is.null(d)) {
      d <- read_csv(f, show_col_types = FALSE)
      if (require_route_converged) d <- d %>% filter(route_converged %in% TRUE)
    }
    sp_name <- d$species[1]
    routes <- d %>%
      filter(.data[[sc]] %in% 1:7, !is.na(trend)) %>%
      transmute(cls = as.integer(.data[[sc]]), score = cls - 4, trend)
    key <- paste(model, sc)
    bs  <- if (key %in% names(bayes)) bayes[[key]] %>% filter(species == sp_name) %>% slice(1) else NULL
    if (!is.null(bs) && nrow(bs) == 0) bs <- NULL

    # -- slope figure
    if (do_slope) {
      status <- "ok"
      if (nrow(routes) == 0) {
        status <- "no_classified_routes"
      } else {
        sl <- NULL
        if (!is.null(bs)) {
          sl <- bs %>% select(b, lo, hi, label)
        } else if (n_distinct(routes$score) >= 2) {
          x <- routes$score
          sl <- tibble(b = sum((x - mean(x)) * routes$trend) / sum((x - mean(x))^2),
                       lo = NA_real_, hi = NA_real_, label = "No credible interval")
          status <- "ok_posterior_mean_fit"
        } else {
          status <- "ok_single_class"
        }
        save_fig(plot_slope(routes, sp_name, model, scenarios[[sc]], sl), out_file, fig_width)
        n_done <- n_done + 1
      }
      manifest[[length(manifest) + 1]] <- tibble(species = sp_name, slug = slug, model = model,
                                                 scenario = sc, n_routes = nrow(routes), status = status)
    }

    # -- category figure
    if (do_cat) {
      status <- "ok"
      if (nrow(routes) == 0) {
        status <- "no_classified_routes"
      } else {
        bayesian <- !is.null(bs) && any(is.finite(c(bs$mean_decreasing_mean, bs$mean_neutral_mean, bs$mean_increasing_mean)))
        if (bayesian) {
          cm <- bs %>%
            select(starts_with("mean_decreasing_"), starts_with("mean_neutral_"), starts_with("mean_increasing_")) %>%
            pivot_longer(everything(), names_to = c("cat", ".value"), names_pattern = "mean_(decreasing|neutral|increasing)_(.*)") %>%
            filter(is.finite(mean))
        } else {
          cm <- routes %>%
            mutate(cat = case_when(cls <= 3 ~ "decreasing", cls == 4 ~ "neutral", TRUE ~ "increasing")) %>%
            group_by(cat) %>% summarise(mean = mean(trend), .groups = "drop") %>%
            mutate(sd = NA_real_, q025 = NA_real_, q975 = NA_real_)
          status <- "ok_posterior_mean_fit"
        }
        if (n_distinct(cm$cat) < 2) status <- paste0(status, "_single_category")
        save_fig(plot_category(routes, sp_name, model, scenarios[[sc]], cm, bayesian), out_cat, cat_width)
        n_done <- n_done + 1
      }
      manifest_cat[[length(manifest_cat) + 1]] <- tibble(species = sp_name, slug = slug, model = model,
                                                         scenario = sc, n_routes = nrow(routes), status = status)
    }
    if (n_done %% 200 == 0 && n_done > 0) cat("  ", n_done, "figures written\n")
  }
}

# Summary ------------------------------------------------------------------------
write_manifest <- function(rows, dir, what) {
  if (length(rows) == 0) return(invisible(NULL))
  df <- bind_rows(rows)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(df, file.path(dir, paste0("manifest_", firstYear, "_", lastYear, ".csv")))
  cat("\n===", what, "figures ===\n"); print(df %>% count(model, scenario, status)); cat("Written to:", dir, "\n")
}
write_manifest(manifest,     out_dir,     "Tracking-slope")
write_manifest(manifest_cat, out_dir_cat, "Category")
if (length(manifest) + length(manifest_cat) == 0) {
  cat("\nNothing to do -- every figure already exists (set overwrite <- TRUE to redraw).\n")
}
