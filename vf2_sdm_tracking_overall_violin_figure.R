# vf2_sdm_tracking_overall_violin_figure.R
#
# Figure 2. All-species version of vf3's concept figure: one panel pooling the
# route-level trends of every species in the analysis (native species with a
# tracking slope from 4b, status == "ok"; as in 4c).
#
#   - violins: route trends (posterior means, %/yr) of all species, pooled, in each SDM class 1-7
#   - black ticks and bars: the average species' trend in each class with its 95%
#     credible interval, from the 4c meta-analysis of the species' class means
#     ("Class k" rows of the category-levels file); until 4c has produced them, the
#     plain average of the species' class means, without bars. Either way each species
#     counts once, so common species don't dominate -- the same "average species" scale
#     as the line
#   - line: the stage-2 overall average tracking slope mu ("Overall" model in 4c),
#     shaded fan = slopes within mu's 95% credible interval
#
# Where the line sits: each species' stage-1 line passes through that species'
# (mean score, mean trend). Averaging those lines over species gives a line through
# (mean of species' xbar, mean of species' ybar), drawn here with slope mu.
#
# Inputs:
#   output/species_routes_covariates/per_species_sdm/<Species>_<model_tag>_route_trends_sdm.csv
#   output/sdm_tracking/sdm_tracking_stage1_species_<suffix>.csv
#   output/sdm_tracking/sdm_tracking_stage2_models_<suffix>_native.csv   (4c)
#   output/sdm_tracking/sdm_tracking_category_levels_<suffix>_native.csv (4c; class estimates)
#   data/spp_names_codes_group_aou_taxa_redlist_diet_migration_invasive.csv  (non-native list)
# Output:
#   output/figures/figure2.{png,pdf}

library(here)
library(tidyverse)

here::i_am("vf2_sdm_tracking_overall_violin_figure.R")

# ---- settings ----
model_tag <- "anthro"
scenario  <- "rcp45"
firstYear <- 2010
lastYear  <- 2025
suffix    <- paste(model_tag, scenario, firstYear, lastYear, sep = "_")

# Row of the 4c models table whose mu is drawn: the Student-t overall model. In the
# current 4c it is "Overall"; output from the older 4c (Normal "Overall") has it as
# "Overall (Student-t)", which is used when present.
stage2_model <- "Overall"
y_lim        <- c(-5, 5)    # y-axis window (%/yr); routes beyond it are cut off visually
line_colour  <- "#2a78d6"     # the "Tracks" blue 

route_dir  <- here::here("output", "species_routes_covariates", "per_species_sdm")
track_dir  <- here::here("output", "sdm_tracking")
traits_csv <- here::here("data", "spp_names_codes_group_aou_taxa_redlist_diet_migration_invasive.csv")
fig_dir    <- here::here("output", "figures")
fig_stem   <- "figure2"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# ---- data ----
stage1 <- read_csv(file.path(track_dir, paste0("sdm_tracking_stage1_species_", suffix, ".csv")), show_col_types = FALSE)
models <- read_csv(file.path(track_dir, paste0("sdm_tracking_stage2_models_", suffix, "_native.csv")),
                   show_col_types = FALSE)

if ("Overall (Student-t)" %in% models$model) stage2_model <- "Overall (Student-t)"
mu <- models %>% filter(model == stage2_model, response == "Tracking slope")
cat("Slope line from model:", stage2_model, sprintf("(mu %.3f [%.3f, %.3f])\n", mu$mu_mean, mu$mu_lo, mu$mu_hi))
stopifnot(nrow(mu) == 1)

# native species only (US-RIIS non-native species dropped, as in 4c)
non_native <- read_csv(traits_csv, show_col_types = FALSE) %>%
  filter(in_bbs %in% TRUE, invasive != "Not listed") %>% pull(Code)
species_ok <- stage1 %>% filter(status == "ok", !species_code %in% non_native) %>% pull(species)

route_file <- function(sp) {
  file.path(route_dir, paste0(gsub(" ", "_", gsub("'", "", sp)), "_", model_tag, "_route_trends_sdm.csv"))
}
missing <- species_ok[!file.exists(route_file(species_ok))]
if (length(missing)) warning(length(missing), " species have no route file and are skipped: ",
                             paste(head(missing, 5), collapse = ", "))

routes <- map_dfr(setdiff(species_ok, missing), function(sp) {
  read_csv(route_file(sp), show_col_types = FALSE,
           col_types = cols(.default = col_guess(), route = col_character())) %>%
    filter(route_converged %in% TRUE, .data[[scenario]] %in% 1:7) %>%
    transmute(species, cls = as.integer(.data[[scenario]]), score = cls - 4, trend)
})

n_species <- n_distinct(routes$species)
cat(nrow(routes), "route x species rows from", n_species, "species\n")

# per-species centroids -> average-species line
centroids <- routes %>% group_by(species) %>%
  summarise(xbar = mean(score), ybar = mean(trend), .groups = "drop")
xbar <- mean(centroids$xbar)
ybar <- mean(centroids$ybar)
line_at <- function(b, cls) ybar + b * ((cls - 4) - xbar)

fitted_line <- tibble(cls = c(1, 7), y = line_at(mu$mu_mean, cls))
fan <- tibble(cls = c(1, 7, 7, 1), b = c(mu$mu_lo, mu$mu_lo, mu$mu_hi, mu$mu_hi)) %>%
  mutate(y = line_at(b, cls))

class_stats <- routes %>%
  group_by(cls, species) %>% summarise(m = mean(trend), n = n(), .groups = "drop") %>%
  group_by(cls) %>%
  summarise(mean = mean(m), n_routes = sum(n), n_sp = n(), .groups = "drop") %>%
  mutate(lo = NA_real_, hi = NA_real_)

# Per-class average species from the 4c meta-analysis (Student-t; "Class k" rows of the
# category-levels file): replaces the simple species average and adds its 95% CrI.
# Falls back to the simple average, without intervals, if 4c has not produced them yet.
levels_csv <- file.path(track_dir, paste0("sdm_tracking_category_levels_", suffix, "_native.csv"))
class_levels <- if (file.exists(levels_csv)) {
  read_csv(levels_csv, show_col_types = FALSE) %>%
    filter(species_set == "main", str_detect(quantity, "^Class [1-7]")) %>%
    transmute(cls = as.integer(str_extract(quantity, "[1-7]")), mu, mu_lo, mu_hi)
} else tibble()
has_ci <- nrow(class_levels) > 0
if (has_ci) {
  class_stats <- class_stats %>% left_join(class_levels, by = "cls") %>%
    mutate(mean = mu, lo = mu_lo, hi = mu_hi) %>% select(-mu, -mu_lo, -mu_hi)
  cat("Class ticks: meta-analysis estimates with 95% CrI (4c)\n")
} else {
  cat("Class ticks: simple species averages (no 'Class k' rows in", basename(levels_csv), "yet; rerun 4b and 4c)\n")
}

n_cut <- sum(routes$trend < y_lim[1] | routes$trend > y_lim[2])
fmt <- function(v) sprintf("%+.2f", v)

# ---- plot ----
p <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.5) +
  geom_violin(data = routes, aes(x = cls, y = trend, group = cls),
              fill = "grey70", colour = NA, alpha = 0.55, scale = "width", width = 0.8, trim = TRUE) +
  geom_polygon(data = fan, aes(x = cls, y = y), fill = line_colour, alpha = 0.18) +
  geom_line(data = fitted_line, aes(x = cls, y = y), colour = line_colour, linewidth = 1.2, lineend = "round") +
  { if (has_ci) geom_errorbar(data = class_stats, aes(x = cls, ymin = lo, ymax = hi),
                              width = 0, linewidth = 0.7, colour = "black") } +
  geom_segment(data = class_stats, aes(x = cls - 0.22, xend = cls + 0.22, y = mean, yend = mean),
               linewidth = 0.9, colour = "black", lineend = "round") +
  scale_x_continuous(breaks = 1:7, limits = c(0.5, 7.5), expand = c(0, 0),
                     # routes / species per class on a top axis, so the bottom axis is just the class number
                     sec.axis = dup_axis(name = NULL, breaks = class_stats$cls,
                                         labels = paste0(format(class_stats$n_routes, big.mark = ",", trim = TRUE),
                                                         " routes
", class_stats$n_sp, " spp"))) +
  coord_cartesian(ylim = y_lim) +
  labs(
    # title = paste0("All species (n = ", n_species, ")"),
    subtitle = paste0("Average tracking slope ", fmt(mu$mu_mean), " (", fmt(mu$mu_lo), " to ", fmt(mu$mu_hi),
                      ") % per year per SDM class step"),
    x = "SDM class",
    y = "Route trend (% per year)",
    # caption = paste0("Violins: route trends (posterior means) of all species pooled, by SDM class; black ticks: class means averaged over species.\n",
    #                  "Numbers above each class: routes and species. Line: average species' tracking slope (stage-2 '", stage2_model,
    #                  "' model, mu);\nshaded: 95% credible interval. Routes outside ", y_lim[1], " to +", y_lim[2],
    #                  " %/yr are cut off (", format(n_cut, big.mark = ","), " of ", format(nrow(routes), big.mark = ","), ").")
    ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold", size = 18),
        plot.subtitle = element_text(size = 14),
        axis.text.y = element_text(size = 13),
        plot.caption = element_text(hjust = 0, size = 8, colour = "grey30"),
        axis.title = element_text(face = "bold", size = 15),
        axis.text.x = element_text(size = 14, face = "bold", colour = "black"),
        axis.text.x.top = element_text(size = 12, face = "plain", colour = "grey30", lineheight = 0.9))

for (ext in c("png", "pdf")) {
  out_file <- file.path(fig_dir, paste0(fig_stem, ".", ext))
  ggsave(out_file, p, width = 10, height = 6, dpi = 300, bg = "white")
  cat("Saved:", out_file, "\n")
}
