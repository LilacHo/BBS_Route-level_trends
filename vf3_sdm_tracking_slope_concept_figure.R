# vf3_sdm_tracking_slope_concept_figure.R
#
# Figure 3. Concept figure for the SDM-tracking analysis (4b stage 1, 4c stage 2):
# "what is a tracking slope, and what hides inside Inconclusive?"
#
#   Row 1 - overall tracking slope (classes 1-7), one example species per stage-2 class:
#     A  Tracks          B  Inconclusive          C  Contrary
#   Row 2 - the Inconclusive species split by their range-change CATEGORY means
#           (Decreasing suitability = classes 1-3, Neutral = 4, Increasing suitability = 5-7):
#     D  step up from Decreasing to Neutral      (trends lower where suitability declines)
#     E  step up from Neutral to Increasing      (trends higher where suitability improves)
#     F  no clear step up                        (optional, show_panel_f)
#
# Row 1: violins of route-level trends (posterior means, %/yr) per SDM class, class means
# (black ticks), and the stage-1 tracking slope with a shaded fan of lines whose slope lies
# within the 95% credible interval. Stage 1 fits, in each posterior draw, a least-squares
# line of route trend on SDM score (-3 ... +3; class 4 = 0), weighting every route equally;
# such a line passes through (mean score, mean trend), so it is drawn as
# trend = ybar + slope * (score - xbar)  with the stage-1 posterior mean slope.
#
# Row 2 (layout as the old per-species "_trend_by_group_..._cld" plots): the routes pooled
# into the three categories, one violin each; black ticks and error bars = stage-1 category
# means with 95% credible intervals; lines joining the means = the two steps.
#
# Colours, in every panel: the classification of THAT line with the same rule as the
# stage-2 classes (class_threshold): blue if P(> 0) >= threshold, red if P(< 0) >= threshold,
# otherwise charcoal (dashed for the row-2 steps). Row 1 uses the stage-2 Student-t class of the
# species. The row-2 steps use  Neutral - Decreasing  and  Increasing - Neutral,  with mean =
# difference of the stage-1 category means and sd = sqrt(sd1^2 + sd2^2) (normal
# approximation treating the two means as independent, as in 4c's "... - Neutral" models;
# the joint draws are not saved by 4b). The step that defines the subgroup is shaded.
#
# Subgroups in row 2 (counted over the native Inconclusive species; titles show n):
#   D: P(Neutral - Decreasing > 0) >= class_threshold
#   E: P(Increasing - Neutral > 0) >= class_threshold
#   F: neither (includes flat species, species with a step down, and species missing a category)
#
# Both rows share one facet_wrap with free x scales: row 1 sits at x = 1..7 (classes) and
# row 2 at x = 11..13 (categories), so each panel keeps only its own axis labels.
#
# Inputs:
#   output/species_routes_covariates/per_species_sdm/<Species>_<model_tag>_route_trends_sdm.csv
#   output/sdm_tracking/sdm_tracking_stage1_species_<suffix>.csv
#   output/sdm_tracking/sdm_tracking_stage2_species_<suffix>_native.csv
# Output:
#   output/figures/figure3.{png,pdf}

library(here)
library(tidyverse)

here::i_am("vf3_sdm_tracking_slope_concept_figure.R")

# ---- settings ----
model_tag <- "anthro"
scenario  <- "rcp45"
firstYear <- 2010
lastYear  <- 2025
suffix    <- paste(model_tag, scenario, firstYear, lastYear, sep = "_")

class_threshold <- 0.90   # same as 4c

# Example species; names must match the 'species' column. Change freely.
#   Tracks:       Eastern Towhee (class means rise almost steadily, Spearman 0.93)
#   Inconclusive: Baltimore Oriole (flat overall and between categories: -2.20 / -2.19 / -2.36;
#                 118 / 855 / 298 routes). Alternatives: Scarlet Tanager, Cooper's Hawk.
#   Contrary:     White-breasted Nuthatch
#   D (Neutral above Decreasing, Increasing flat): Black-and-white Warbler (166 / 296 / 130
#                 routes; step +1.0, P = 0.95). Alternatives: Ring-billed Gull (step +6.6 but
#                 very wide spread), Acadian Flycatcher (only 9 routes in classes 1-2).
#   E (Increasing above Neutral, Decreasing flat): Cedar Waxwing (456 / 1,028 / 25 routes).
#                 Alternatives: Hooded Warbler, Double-crested Cormorant; Wood Thrush and
#                 Common Yellowthroat also step up but are valleys (Decreasing above Neutral).
#   F (no clear step up): same species as panel B, so F reads as "B split by category"
panels <- tribble(
  ~id, ~species,                  ~kind,
  "A", "Eastern Towhee",          "overall",
  "B", "Baltimore Oriole",        "overall",
  "C", "White-breasted Nuthatch", "overall",
  "D", "Black-and-white Warbler", "step_dn",
  "E", "Cedar Waxwing",           "step_in",
  "F", "Baltimore Oriole",        "step_none"
)
show_panel_f <- TRUE

y_lim         <- c(-12, 14)   # y-axis window (%/yr); routes beyond it are cut off visually
min_violin_n  <- 5            # classes / categories with fewer routes are drawn as points, not a violin
class_colours <- c(Tracks = "#2a78d6", Inconclusive = "#2b2b2b", Contrary = "#d9423f")   # charcoal: stays distinct from the red under protanopia
shade_colour  <- "#2a78d6"    # background of the step that defines the row-2 subgroup

# row-2 categories: x position (offset so the two rows never share axis breaks) and axis label
categories <- tibble(cat   = c("decreasing", "neutral", "increasing"),
                     x     = c(11, 12, 13),
                     label = c("Decreasing\n(1–3)", "Neutral\n(4)", "Increasing\n(5–7)"))

route_dir <- here::here("output", "species_routes_covariates", "per_species_sdm")
track_dir <- here::here("output", "sdm_tracking")
fig_dir   <- here::here("output", "figures")
fig_stem  <- "figure3"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

if (!show_panel_f) panels <- filter(panels, id != "F")

# ---- data ----
stage1 <- read_csv(file.path(track_dir, paste0("sdm_tracking_stage1_species_", suffix, ".csv")), show_col_types = FALSE)
stage2 <- read_csv(file.path(track_dir, paste0("sdm_tracking_stage2_species_", suffix, "_native.csv")), show_col_types = FALSE)

classify <- function(p_pos) case_when(p_pos >= class_threshold     ~ "Tracks",
                                      1 - p_pos >= class_threshold ~ "Contrary",
                                      TRUE                         ~ "Inconclusive")

# category steps per species (normal approximation, see header)
stage1 <- stage1 %>%
  mutate(dn_mean = mean_neutral_mean - mean_decreasing_mean,
         dn_sd   = sqrt(mean_neutral_sd^2 + mean_decreasing_sd^2),
         in_mean = mean_increasing_mean - mean_neutral_mean,
         in_sd   = sqrt(mean_increasing_sd^2 + mean_neutral_sd^2),
         dn_p_pos = pnorm(dn_mean / dn_sd),
         in_p_pos = pnorm(in_mean / in_sd))

# species counts per class and per Inconclusive subgroup (native species, as in 4c)
spp <- stage2 %>% select(species, label = class_studentt) %>%
  left_join(stage1 %>% select(species, dn_p_pos, in_p_pos), by = "species") %>%
  mutate(dn_up = coalesce(dn_p_pos >= class_threshold, FALSE),
         in_up = coalesce(in_p_pos >= class_threshold, FALSE))
inc <- filter(spp, label == "Inconclusive")
n_group <- c(Tracks       = sum(spp$label == "Tracks"),
             Inconclusive = nrow(inc),
             Contrary     = sum(spp$label == "Contrary"),
             step_dn      = sum(inc$dn_up),
             step_in      = sum(inc$in_up),
             step_none    = sum(!inc$dn_up & !inc$in_up))
n_both <- sum(inc$dn_up & inc$in_up)
cat(sprintf("Native species: %d Tracks, %d Inconclusive, %d Contrary\n",
            n_group[["Tracks"]], n_group[["Inconclusive"]], n_group[["Contrary"]]))
cat(sprintf("Inconclusive: %d Neutral > Decreasing, %d Increasing > Neutral (%d both), %d neither\n",
            n_group[["step_dn"]], n_group[["step_in"]], n_both, n_group[["step_none"]]))

# file names use the species name with spaces -> "_" and apostrophes dropped
route_file <- function(sp) {
  file.path(route_dir, paste0(gsub(" ", "_", gsub("'", "", sp)), "_", model_tag, "_route_trends_sdm.csv"))
}
route_data <- map_dfr(unique(panels$species), function(sp) {
  read_csv(route_file(sp), show_col_types = FALSE) %>%
    filter(route_converged %in% TRUE, .data[[scenario]] %in% 1:7) %>%
    transmute(species, cls = as.integer(.data[[scenario]]), score = cls - 4, trend)
})

# per-panel info and titles
fmt  <- function(v) sprintf("%+.2f", v)
ci   <- function(m, lo, hi) paste0(fmt(m), " (", fmt(lo), " to ", fmt(hi), ")")
ci_n <- function(m, s) ci(m, m - 1.96 * s, m + 1.96 * s)
info <- panels %>%
  left_join(stage1, by = "species") %>%
  left_join(stage2 %>% select(species, label = class_studentt), by = "species")
if (any(is.na(info$label))) stop("Not in the stage-2 (native) output: ", paste(info$species[is.na(info$label)], collapse = ", "))

info <- info %>%
  mutate(head = case_when(
           kind == "overall"   ~ paste0(id, "  ", label, " (", n_group[label], " species)"),
           kind == "step_dn"   ~ paste0(id, "  Inconclusive: Neutral above Decreasing (", n_group[["step_dn"]], ")"),
           kind == "step_in"   ~ paste0(id, "  Inconclusive: Increasing above Neutral (", n_group[["step_in"]], ")"),
           kind == "step_none" ~ paste0(id, "  Inconclusive: no clear step up (", n_group[["step_none"]], ")")),
         body = if_else(kind == "overall",
                        paste0(species, "\nslope 1–7  ", ci(slope_mean, slope_q025, slope_q975)),
                        paste0(species, "\nNeutral − Decreasing  ", ci_n(dn_mean, dn_sd),
                               "\nIncreasing − Neutral  ", ci_n(in_mean, in_sd))),
         panel = paste0(head, "\n", body))
info$panel <- factor(info$panel, levels = info$panel)

# every route once per panel; row-1 panels use the class (x = 1..7), row-2 panels the category (x = 11..13)
routes <- panels %>% select(id, species, kind) %>%
  left_join(info %>% select(id, panel), by = "id") %>%
  inner_join(route_data, by = "species", relationship = "many-to-many") %>%
  mutate(cat = case_when(cls <= 3 ~ "decreasing", cls == 4 ~ "neutral", TRUE ~ "increasing")) %>%
  left_join(categories %>% select(cat, cat_x = x), by = "cat") %>%
  mutate(x = if_else(kind == "overall", as.numeric(cls), cat_x))

# ---- row 1: slope line and fan ----
line_at <- function(ybar, xbar, b, cls) ybar + b * ((cls - 4) - xbar)
slope_def <- info %>% filter(kind == "overall") %>%
  transmute(id, panel, b = slope_mean, lo = slope_q025, hi = slope_q975, line_class = label) %>%
  left_join(routes %>% group_by(id) %>% summarise(xbar = mean(score), ybar = mean(trend), .groups = "drop"), by = "id")
slope_segs <- slope_def %>% mutate(x0 = 1, x1 = 7, y0 = line_at(ybar, xbar, b, 1), y1 = line_at(ybar, xbar, b, 7))
slope_fans <- slope_def %>%
  crossing(tibble(x = c(1, 7, 7, 1), which = c("lo", "lo", "hi", "hi"))) %>%
  mutate(y = line_at(ybar, xbar, if_else(which == "lo", lo, hi), x))

# ---- row 2: category means and steps ----
step_info <- info %>% filter(kind != "overall")
cat_means <- step_info %>%
  select(id, panel, starts_with("mean_decreasing_"), starts_with("mean_neutral_"), starts_with("mean_increasing_")) %>%
  pivot_longer(-c(id, panel), names_to = c("cat", ".value"), names_pattern = "mean_(decreasing|neutral|increasing)_(.*)") %>%
  filter(is.finite(mean)) %>%
  left_join(categories %>% select(cat, x), by = "cat")
steps <- bind_rows(
  step_info %>% transmute(id, panel, from = "decreasing", to = "neutral",    p = dn_p_pos),
  step_info %>% transmute(id, panel, from = "neutral",    to = "increasing", p = in_p_pos)
) %>%
  filter(is.finite(p)) %>%
  left_join(cat_means %>% select(id, from = cat, x0 = x, y0 = mean), by = c("id", "from")) %>%
  left_join(cat_means %>% select(id, to = cat,   x1 = x, y1 = mean), by = c("id", "to")) %>%
  mutate(line_class = classify(p), dashed = line_class == "Inconclusive")

shade <- info %>% filter(kind %in% c("step_dn", "step_in")) %>%
  transmute(panel, xmin = if_else(kind == "step_dn", 10.5, 11.5), xmax = if_else(kind == "step_dn", 12.5, 13.5))

# ---- violins, class means, counts ----
group_stats <- routes %>%
  group_by(id, kind, panel, x) %>%
  summarise(n = n(), mean = mean(trend), .groups = "drop")
big_groups    <- group_stats %>% filter(n >= min_violin_n) %>% select(panel, x)
violin_data   <- routes %>% semi_join(big_groups, by = c("panel", "x"))
sparse_points <- routes %>% anti_join(big_groups, by = c("panel", "x"))
class_ticks   <- group_stats %>% filter(kind == "overall")   # row 1 only; row 2 uses the stage-1 category means
n_cut <- routes %>% group_by(id, species) %>%
  summarise(n = sum(trend < y_lim[1] | trend > y_lim[2]), .groups = "drop")

# x range of each panel (free x scales), and one axis carrying both rows' labels
x_range <- info %>% select(panel, kind) %>%
  crossing(end = 1:2) %>%
  mutate(x = case_when(kind == "overall" & end == 1 ~ 0.5, kind == "overall" ~ 7.5,
                       end == 1 ~ 10.5, TRUE ~ 13.5))
x_breaks <- c(1:7, categories$x)
x_labels <- c(as.character(1:7), categories$label)


# ---- plot ----
make_plot <- function(ids, nrow, title_size = 12.5) {
  keep <- function(d) d %>% filter(panel %in% info$panel[info$id %in% ids])
  cut  <- n_cut %>% filter(id %in% ids) %>% distinct(species, n)
  cut_txt <- paste0("Routes outside ", y_lim[1], " to +", y_lim[2], " %/yr are cut off (",
                    paste0(cut$species, " ", cut$n, collapse = ", "), ").")
  caption <- paste0(
    "Violins: route trends (posterior means), by SDM class (A–C) or pooled by range-change category (D–F: Decreasing suitability = classes 1–3, ",
    "Neutral = 4, Increasing suitability = 5–7); numbers above the violins: routes.\n",
    "A–C: black ticks = class means; line = stage-1 tracking slope over classes 1–7, shaded fan = slopes within its 95% credible interval; ",
    "colour = species class from the stage-2 Student-t model.\n",
    "D–F: black ticks and bars = category means with 95% credible intervals (stage 1); lines = steps between categories.\n",
    "Colour: blue = increase, red = decrease, charcoal (dashed in D–F) = inconclusive, at posterior probability >= ", class_threshold,
    ". Shaded background: the step that defines the subgroup. Species counts: native species.\n", cut_txt)
  ggplot() +
    geom_blank(data = keep(x_range), aes(x = x)) +
    geom_rect(data = keep(shade), aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = shade_colour, alpha = 0.06) +
    geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.5) +
    geom_violin(data = keep(violin_data), aes(x = x, y = trend, group = x),
                fill = "grey70", colour = NA, alpha = 0.55, scale = "width", width = 0.8, trim = TRUE) +   # trim: no tails beyond the data
    geom_point(data = keep(sparse_points), aes(x = x, y = trend), colour = "grey55", size = 1) +
    # row 1: slope fan, slope line, class means
    geom_polygon(data = keep(slope_fans), aes(x = x, y = y, fill = line_class, group = id), alpha = 0.18) +
    geom_segment(data = keep(slope_segs), aes(x = x0, xend = x1, y = y0, yend = y1, colour = line_class),
                 linewidth = 1.2, lineend = "round") +
    geom_segment(data = keep(class_ticks), aes(x = x - 0.22, xend = x + 0.22, y = mean, yend = mean),
                 linewidth = 0.9, colour = "black", lineend = "round") +
    # row 2: steps between category means, then the means with their 95% CI
    geom_segment(data = keep(steps), aes(x = x0, xend = x1, y = y0, yend = y1, colour = line_class, linetype = dashed),
                 linewidth = 1.2, lineend = "round") +
    geom_errorbar(data = keep(cat_means), aes(x = x, ymin = q025, ymax = q975), width = 0, linewidth = 0.7, colour = "black") +
    geom_segment(data = keep(cat_means), aes(x = x - 0.22, xend = x + 0.22, y = mean, yend = mean),
                 linewidth = 0.9, colour = "black", lineend = "round") +
    # routes per class / category on top of the violins (as in Figure 2), in a white strip above
    # y_lim[2] so they never overlap data; a top sec.axis like Figure 2's can't differ between facets
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = y_lim[2] + 0.03 * diff(y_lim), ymax = Inf, fill = "white", colour = "white") +
    geom_text(data = keep(group_stats), aes(x = x, y = y_lim[2] + 0.085 * diff(y_lim),
                                            label = format(n, big.mark = ",", trim = TRUE)),
              size = 3.8, colour = "grey30") +
    facet_wrap(~ panel, nrow = nrow, scales = "free_x") +
    scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0, 0)) +
    scale_colour_manual(values = class_colours, guide = "none") +
    scale_fill_manual(values = class_colours, guide = "none") +
    scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22"), guide = "none") +
    scale_y_continuous(breaks = function(l) { b <- scales::breaks_extended(5)(y_lim); b[b >= y_lim[1] & b <= y_lim[2]] }) +   # ticks only inside y_lim, not in the count strip
    coord_cartesian(ylim = c(y_lim[1] - 0.02 * diff(y_lim), y_lim[2] + 0.13 * diff(y_lim)), expand = FALSE) +
    labs(x = "SDM class (A–C)  ·  range-change category (D–F)",
         y = "Route trend (% per year)",
         caption = caption) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.text = element_text(hjust = 0, face = "bold", size = title_size, lineheight = 1.1),
          plot.caption = element_text(hjust = 0, size = 8.5, colour = "grey30"),
          axis.title = element_text(face = "bold", size = 15),
          axis.text.x = element_text(size = 13, face = "bold", colour = "black", lineheight = 0.9),
          axis.text.y = element_text(size = 12),
          panel.spacing.x = unit(1.4, "lines"),
          panel.spacing.y = unit(1.6, "lines"))
}

save_plot <- function(p, stem, width, height) {
  for (ext in c("png", "pdf")) {
    out_file <- file.path(fig_dir, paste0(stem, ".", ext))
    if (ext == "pdf") {
      ggsave(out_file, p, width = width, height = height, bg = "white", device = cairo_pdf)   # cairo: keeps the en dashes
    } else {
      ggsave(out_file, p, width = width, height = height, dpi = 300, bg = "white")
    }
    cat("Saved:", out_file, "\n")
  }
}

save_plot(make_plot(panels$id, nrow = 2, title_size = 11.5), fig_stem, width = 15, height = 12)
