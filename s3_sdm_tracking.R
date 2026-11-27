# s3_sdm_tracking.R
#
# Supplementary Note 3 analyses of the SDM-tracking pipeline. Two parts, each switched on
# or off below:
#
# Part 1 (run_normal_vs_t; needs Stan). Why the meta-analysis uses a Student-t rather than
#   a normal distribution of species' true values. Refits the main models of 4c with a
#   Normal species distribution and compares the two by approximate leave-one-species-out
#   cross-validation (PSIS-LOO; positive elpd_diff = the Student-t predicts held-out
#   species better).
#     1. Tracking slope: mu, tau, % of species with a negative true slope and elpd under
#        each distribution; species classes (Tracks / Inconclusive / Contrary) under the
#        Normal, cross-tabulated against the Student-t classes from 4c
#     2. Category means (main species set: Decreasing suitability / Neutral / Increasing
#        suitability): mu, tau and elpd under each distribution
#   Uses the same data, priors, sampler settings and seeds as 4c and 4d
#   (helper/sdm_tracking_meta_common.R), so the Student-t fits here reproduce 4c's.
#
# Part 2 (run_step_plane; no Stan). Supplementary Fig. 3, the "step plane": for each species
#   with routes in all three range-change categories, where does it respond to suitability
#   change, at the trailing edge (losing suitability), at the leading edge (gaining
#   suitability), both, or neither?
#     x = trailing-edge step = mean(Neutral) - mean(Decreasing suitability)   (%/yr)
#     y = leading-edge step  = mean(Increasing suitability) - mean(Neutral)   (%/yr)
#   The category means and SDs come from stage 1 (4b). Steps use the same normal
#   approximation as vf3 and the 4c "... - Neutral" models: mean = difference of means,
#   sd = sqrt(sd1^2 + sd2^2). Both steps share the Neutral mean, so their errors are slightly
#   negatively correlated (cov = -var(Neutral)); Neutral usually has the most routes, so
#   this is small. Also prints the numbers of species with a clear step in each direction
#   on each axis, against the number expected by chance, as quoted in the Results.
#   Encoding
#     colour   stage-2 class of the overall tracking slope (Tracks / Inconclusive / Contrary, Student-t)
#     fill     filled = at least one step is credible (P(step > 0) >= class_threshold or
#              P(step < 0) >= class_threshold); open = neither step is credible
#     axes     pseudo-log (asinh) scale, so the crowd near 0 spreads out but large steps stay visible
#     y = x    dashed line: the same step at both edges (symmetric response)
#     corners  per quadrant (sign of the posterior mean steps): % of plotted species, and in brackets
#              the number whose two steps are both credible in that quadrant's direction
#     labels   the five example species of main-text Figure 3, labelled in the empty left / right
#              edge band with a leader line, so the two figures link
#     insets   one schematic per quadrant: the three category means (D, N, I) giving that
#              quadrant's shape (rising, valley, hump, falling)
#   Species missing the Decreasing or Increasing category have no step on that axis and are
#   left out; the console reports how many.
#
# Reads:  output/sdm_tracking/sdm_tracking_stage1_species_<suffix>.csv   (4b)
#         output/sdm_tracking/sdm_tracking_stage2_species_<suffix>.csv   (4c; Student-t species classes)
#         data/spp_names_codes_group_aou_taxa_redlist_diet_migration_invasive.csv   (Part 1)
#         models/sdm_tracking_meta.stan   (Part 1)
# Writes ("_native" added when native_only):
#   Part 1 (output/sdm_tracking/):
#     sdm_tracking_s3_normal_vs_studentt_SUMMARY_<suffix>.txt
#     sdm_tracking_s3_normal_vs_studentt_<suffix>.csv          one row per quantity x distribution
#     sdm_tracking_s3_species_normal_<suffix>.csv              species slopes and classes under the Normal
#   Part 2: output/figures/supp_figure3.{png,pdf}
#
# Needs: tidyverse; Part 1 also cmdstanr (with CmdStan installed), posterior, loo.

library(here)
library(tidyverse)

here::i_am("s3_sdm_tracking.R")

# Settings (must match 4a-4d) --------------------------------------------------
model_tag <- "anthro"; scenario <- "rcp45"; firstYear <- 2010; lastYear <- 2025
exclude_classes <- integer(0)
native_only <- TRUE

run_normal_vs_t <- TRUE   # Part 1: Normal vs Student-t (Stan; run on the fit machine)
run_step_plane  <- TRUE   # Part 2: Supplementary Fig. 3 (no Stan; runs anywhere)

# =============================================================================
# Part 1. Normal vs Student-t species distribution (Stan)
# =============================================================================
if (run_normal_vs_t) {
  source(here::here("helper", "sdm_tracking_meta_common.R"))   # settings, data d, fit functions

  species_csv <- file.path(out_dir, paste0("sdm_tracking_stage2_species_", suffix, ".csv"))
  out_txt     <- file.path(out_dir, paste0("sdm_tracking_s3_normal_vs_studentt_SUMMARY_", suffix, ".txt"))
  out_csv     <- file.path(out_dir, paste0("sdm_tracking_s3_normal_vs_studentt_", suffix, ".csv"))
  out_species <- file.path(out_dir, paste0("sdm_tracking_s3_species_normal_", suffix, ".csv"))

  # Fit one quantity under both distributions and compare them --------------------
  compare_dist <- function(quantity, dd, ycol, secol, prior_scale, seed_use) {
    cat("Fitting:", quantity, "\n")
    ms <- list(`Student-t` = fit_meta(dd, ycol, secol, student = 1, prior_scale = prior_scale, seed_use = seed_use),
               Normal      = fit_meta(dd, ycol, secol, student = 0, prior_scale = prior_scale, seed_use = seed_use))
    rs <- imap(ms, ~ summarise_model(.x, .y, quantity))
    d_i <- rs$`Student-t`$loo$pointwise[, "elpd_loo"] - rs$Normal$loo$pointwise[, "elpd_loo"]
    out <- bind_rows(map(rs, "row")) %>%
      rename(distribution = model, quantity = response) %>%
      mutate(elpd_t_minus_normal = sum(d_i), elpd_t_minus_normal_se = sqrt(length(d_i)) * sd(d_i))
    list(table = out, normal = ms$Normal)
  }

  # 1. Tracking slope
  slope <- compare_dist("Tracking slope", d, "slope_mean", "slope_sd", prior_scale_slope, seed)

  # Species classes under the Normal (in the Normal model theta is exact; no species-level
  # refit needed), against the Student-t classes from 4c
  sp_normal <- species_theta(slope$normal) %>%
    rename_with(~ paste0(.x, "_normal"), -species_code) %>%
    mutate(class_normal = classify(theta_p_pos_normal))
  if (file.exists(species_csv)) {
    sp_normal <- read.csv(species_csv, stringsAsFactors = FALSE) %>%
      select(species, species_code, group, theta_mean_studentt, class_studentt) %>%
      right_join(sp_normal, by = "species_code")
  } else warning(basename(species_csv), " not found (run 4c): no comparison with the Student-t classes")

  # 2. Category means (main species set)
  cats <- list(c("Decreasing suitability (classes 1-3)", "mean_decreasing"),
               c("Neutral (class 4)",                    "mean_neutral"),
               c("Increasing suitability (classes 5-7)", "mean_increasing"))
  cat_tabs <- map(cats, ~ compare_dist(.x[1], d, paste0(.x[2], "_mean"), paste0(.x[2], "_sd"),
                                       prior_scale_mean, seed_mean)$table)

  results <- bind_rows(slope$table, cat_tabs) %>%
    select(quantity, distribution, n_species, mu_mean, mu_lo, mu_hi, tau_mean, tau_lo, tau_hi, nu_mean,
           pct_species_negative, elpd_loo, elpd_loo_se, elpd_t_minus_normal, elpd_t_minus_normal_se,
           pareto_k_over_0.7, max_rhat, min_ess_bulk, divergences, ppc_p_sd)
  write.csv(results, out_csv, row.names = FALSE)
  write.csv(sp_normal, out_species, row.names = FALSE)

  # Plain-text summary ---------------------------------------------------------
  con <- file(out_txt, open = "wt", encoding = "UTF-8")
  w <- function(...) { cat(..., "\n", sep = "", file = con); cat(..., "\n", sep = "") }
  w("NORMAL vs STUDENT-T SPECIES DISTRIBUTION (", suffix, ")   ", format(Sys.time(), "%Y-%m-%d %H:%M"))
  w("Species: ", nrow(d), ". elpd(t - Normal) > 0: the Student-t predicts held-out species better ",
    "(|diff| < ~2 se = no clear difference). tau: Student-t scale / Normal SD. Intervals: 95% CrI.")
  w("")
  for (qn in unique(results$quantity)) {
    x <- results %>% filter(quantity == qn)
    w("== ", qn, " (n = ", x$n_species[1], ") ==")
    for (r in seq_len(nrow(x))) {
      w(sprintf("  %-9s mu %s  tau %.3f  nu %s  %% species negative %.1f%%", x$distribution[r],
                ci(x$mu_mean[r], x$mu_lo[r], x$mu_hi[r]), x$tau_mean[r],
                if (is.na(x$nu_mean[r])) "-" else f2(x$nu_mean[r], 1), x$pct_species_negative[r]))
    }
    w(sprintf("  elpd(t - Normal) = %.1f (se %.1f)", x$elpd_t_minus_normal[1], x$elpd_t_minus_normal_se[1]))
    w("")
  }
  w("== Species classes, tracking slope ==")
  tb_n <- table(factor(sp_normal$class_normal, levels = c("Tracks", "Inconclusive", "Contrary")))
  w(sprintf("Normal:     Tracks %3d | Inconclusive %3d | Contrary %3d", tb_n[1], tb_n[2], tb_n[3]))
  if ("class_studentt" %in% names(sp_normal)) {
    tb_t <- table(factor(sp_normal$class_studentt, levels = c("Tracks", "Inconclusive", "Contrary")))
    w(sprintf("Student-t:  Tracks %3d | Inconclusive %3d | Contrary %3d  (4c)", tb_t[1], tb_t[2], tb_t[3]))
    w("Rows: Student-t class; columns: Normal class")
    xt <- table(Student_t = factor(sp_normal$class_studentt, levels = c("Tracks", "Inconclusive", "Contrary")),
                Normal    = factor(sp_normal$class_normal,   levels = c("Tracks", "Inconclusive", "Contrary")))
    for (l in capture.output(print(xt))) w("  ", l)
    chg <- sp_normal %>% filter(class_studentt != class_normal) %>% arrange(class_studentt, species)
    w("Species classified differently (Student-t -> Normal): ", nrow(chg))
    for (i in seq_len(nrow(chg))) w(sprintf("  %-30s %-12s -> %s", chg$species[i], chg$class_studentt[i], chg$class_normal[i]))
  }
  w("")
  w("== Diagnostics ==")
  for (r in seq_len(nrow(results))) {
    x <- results[r, ]
    w(sprintf("  %-38s %-9s max Rhat %.3f  min ESS %5.0f  divergences %d  Pareto k>0.7: %d  PPC p(sd) %.2f",
              x$quantity, x$distribution, x$max_rhat, x$min_ess_bulk, x$divergences, x$pareto_k_over_0.7, x$ppc_p_sd))
  }
  close(con)
  cat("\nWrote:\n ", out_txt, "\n ", out_csv, "\n ", out_species, "\n")
}

# =============================================================================
# Part 2. Step plane: Supplementary Fig. 3 (no Stan)
# =============================================================================
if (run_step_plane) {
  sfx1 <- paste0(paste(model_tag, scenario, firstYear, lastYear, sep = "_"),
                 if (length(exclude_classes) > 0) paste0("_no", paste(sort(exclude_classes), collapse = "")) else "")
  sfx2 <- paste0(sfx1, if (native_only) "_native" else "")
  class_threshold <- 0.90   # same as 4c and vf3 (also set by the helper in Part 1)
  class_colours   <- c(Tracks = "#2a78d6", Inconclusive = "#2b2b2b", Contrary = "#d9423f")   # same as vf3
  asinh_sigma     <- 1      # pseudo-log scale: linear within about +-sigma %/yr, log-like beyond
  axis_breaks     <- c(-10, -5, -2, -1, 0, 1, 2, 5, 10)

  # Main-text Figure 3 example species, labelled here (panel letter in Figure 3 in brackets)
  examples <- tribble(
    ~species,                  ~fig3,
    "Eastern Towhee",          "3A",
    "Baltimore Oriole",        "3B, 3F",
    "White-breasted Nuthatch", "3C",
    "Black-and-white Warbler", "3D",
    "Wood Thrush",             "3E"
  )

  track_dir <- here::here("output", "sdm_tracking")
  fig_dir   <- here::here("output", "figures")
  fig_stem  <- "supp_figure3"
  if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

  # ---- data ----
  stage1 <- read_csv(file.path(track_dir, paste0("sdm_tracking_stage1_species_", sfx1, ".csv")), show_col_types = FALSE)
  stage2 <- read_csv(file.path(track_dir, paste0("sdm_tracking_stage2_species_", sfx2, ".csv")), show_col_types = FALSE)

  steps <- stage2 %>%
    select(species, class = class_studentt) %>%
    left_join(stage1, by = "species") %>%
    transmute(species, class,
              n_routes,
              dn_mean = mean_neutral_mean - mean_decreasing_mean,
              dn_sd   = sqrt(mean_neutral_sd^2 + mean_decreasing_sd^2),
              in_mean = mean_increasing_mean - mean_neutral_mean,
              in_sd   = sqrt(mean_increasing_sd^2 + mean_neutral_sd^2),
              dn_p_pos = pnorm(dn_mean / dn_sd),
              in_p_pos = pnorm(in_mean / in_sd))

  n_all     <- nrow(steps)
  n_no_dn   <- sum(is.na(steps$dn_mean))
  n_no_in   <- sum(is.na(steps$in_mean))
  plot_data <- steps %>% filter(!is.na(dn_mean), !is.na(in_mean))
  n_plot    <- nrow(plot_data)
  cat(sprintf("Native species: %d; plotted (both steps available): %d; no Decreasing step: %d; no Increasing step: %d\n",
              n_all, n_plot, n_no_dn, n_no_in))

  credible <- function(p) case_when(p >= class_threshold ~ 1L, 1 - p >= class_threshold ~ -1L, TRUE ~ 0L)
  plot_data <- plot_data %>%
    mutate(dn_cred = credible(dn_p_pos),
           in_cred = credible(in_p_pos),
           any_cred = if_else(dn_cred != 0 | in_cred != 0, "At least one step credible", "Neither step credible"),
           class = factor(class, levels = c("Inconclusive", "Contrary", "Tracks"))) %>%   # draw order: Tracks on top
    arrange(class, desc(any_cred))

  # ---- quadrant summaries ----
  quad <- plot_data %>%
    mutate(sx = sign(dn_mean), sy = sign(in_mean)) %>%
    group_by(sx, sy) %>%
    summarise(n = n(),
              n_both = sum(dn_cred == sx & in_cred == sy),
              n_tracks = sum(class == "Tracks"),
              .groups = "drop") %>%
    mutate(pct = 100 * n / n_plot)
  print(quad)

  # Species with a clear step on each axis (P >= class_threshold in one direction), among the
  # plotted species (both steps available; the numbers quoted in the Results); with no
  # association, about (1 - class_threshold) x n would be clear in each direction by chance
  for (ax in c("dn", "in")) {
    p  <- plot_data[[paste0(ax, "_p_pos")]]
    cat(sprintf("%s step: %d species; clearly positive %d, clearly negative %d; expected by chance %.0f each way\n",
                if (ax == "dn") "Trailing-edge" else "Leading-edge", length(p),
                sum(p >= class_threshold), sum(1 - p >= class_threshold), (1 - class_threshold) * length(p)))
  }

  # Axis window: the data reach about +-lim_data; the band between lim_data and lim_plot (on both
  # axes) is empty, and the four corners of that band hold the quadrant labels and schematic insets.
  lim_data <- max(abs(c(plot_data$dn_mean, plot_data$in_mean))) * 1.08
  lim_plot <- 34
  asinh_trans <- scales::pseudo_log_trans(sigma = asinh_sigma, base = exp(1))
  tr  <- asinh_trans$transform  # layout is worked out on the transformed scale, then mapped back
  inv <- asinh_trans$inverse    # (annotation_custom and geoms take positions in data units)

  # Quadrant names and the shape of the three category means (Decreasing, Neutral, Increasing)
  quad_names <- tribble(
    ~sx, ~sy, ~name,                        ~shape,
     1,   1,  "Tracks at both edges",        list(c(-1, 0,  1)),
    -1,   1,  "Valley: Neutral lowest",      list(c( 1, 0,  1)),
     1,  -1,  "Hump: Neutral highest",       list(c(-1, 0, -1)),
    -1,  -1,  "Reversed at both edges",      list(c( 1, 0, -1))
  ) %>% mutate(shape = map(shape, 1))

  # inset band, in transformed units: from just outside the data to the panel edge
  band_in  <- tr(lim_data) + 0.04
  band_out <- tr(lim_plot)
  band_w   <- band_out - band_in
  quad_lab <- quad %>%
    left_join(quad_names, by = c("sx", "sy")) %>%
    mutate(x = asinh_trans$inverse(sx * (band_in - 0.06)),      # text sits beside the inset, toward the centre
           y = asinh_trans$inverse(sy * (band_in + band_w / 2)),
           hjust = if_else(sx > 0, 1, 0), vjust = 0.5,
           label = sprintf("%s\n%.0f%% of species [%d credible]", name, pct, n_both))

  # one schematic inset per quadrant: three category means joined by steps, blue = up, red = down
  make_inset <- function(yv) {
    pts  <- tibble(x = 1:3, y = yv)
    segs <- tibble(x = 1:2, xend = 2:3, y = yv[1:2], yend = yv[2:3],
                   dir = if_else(yv[2:3] > yv[1:2], "up", "down"))
    g <- ggplot() +
      geom_segment(data = segs, aes(x, y, xend = xend, yend = yend, colour = dir), linewidth = 1.1, lineend = "round") +
      geom_point(data = pts, aes(x, y), size = 1.8) +
      annotate("text", x = 1:3, y = -1.75, label = c("D", "N", "I"), size = 2.6, colour = "grey30") +
      scale_colour_manual(values = c(up = class_colours[["Tracks"]], down = class_colours[["Contrary"]]), guide = "none") +
      scale_x_continuous(limits = c(0.6, 3.4)) +
      scale_y_continuous(limits = c(-2.1, 1.4)) +
      theme_void() +
      theme(plot.background = element_rect(fill = "white", colour = "grey60", linewidth = 0.4),
            plot.margin = margin(2, 2, 2, 2))
    ggplotGrob(g)
  }
  insets <- quad_names %>%
    mutate(grob = map(shape, make_inset),
           xmin = inv(if_else(sx > 0,  band_in, -band_out)), xmax = inv(if_else(sx > 0, band_out, -band_in)),
           ymin = inv(if_else(sy > 0,  band_in, -band_out)), ymax = inv(if_else(sy > 0, band_out, -band_in)))
  inset_layers <- pmap(insets, function(grob, xmin, xmax, ymin, ymax, ...)
    annotation_custom(grob, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax))

  # Example-species labels sit in the empty left / right edge band (left if the trailing step < 0),
  # at the species' own height, joined to the point by a thin leader line.
  lab_data <- plot_data %>% inner_join(examples, by = "species") %>%
    mutate(label = str_wrap(paste0(species, " (", fig3, ")"), width = 16),
           side  = if_else(dn_mean < 0, -1, 1),
           lab_x = inv(side * (band_in + 0.03)),
           hjust = if_else(side > 0, 0, 1))
  if (nrow(lab_data) < nrow(examples)) {
    warning("Example species not plotted: ", paste(setdiff(examples$species, lab_data$species), collapse = ", "))
  }

  # Spread labels on the same side so they are at least label_gap apart (transformed y units):
  # overlapping neighbours are pushed apart symmetrically until no pair is too close.
  label_gap <- 0.5
  spread <- function(y, gap) {
    o <- order(y); v <- y[o]
    for (iter in 1:100) {
      moved <- FALSE
      for (i in seq_len(length(v) - 1)) {
        d <- v[i + 1] - v[i]
        if (d < gap) { v[i] <- v[i] - (gap - d) / 2; v[i + 1] <- v[i + 1] + (gap - d) / 2; moved <- TRUE }
      }
      if (!moved) break
    }
    y[o] <- v; y
  }
  # Labels must also clear the horizontal y = 0 line, which runs through the label band: a label
  # closer than zero_gap to 0 is moved to +zero_gap or -zero_gap, whichever side has more room
  # from the other labels on its side, and the labels are then spread again.
  zero_gap <- 0.3
  off_zero <- function(y, gap, zgap) {
    y <- spread(y, gap)
    for (i in which(abs(y) < zgap)) {
      room <- function(z) min(c(Inf, abs(y[-i] - z)))
      y[i] <- if (room(zgap) >= room(-zgap)) zgap else -zgap
    }
    spread(y, gap)
  }
  lab_data <- lab_data %>%
    group_by(side) %>%
    mutate(lab_y = inv(off_zero(tr(in_mean), label_gap, zero_gap))) %>%
    ungroup()
  # leader lines: flat from the point to an elbow just inside the band, then angled to the label
  lab_lines <- lab_data %>%
    mutate(elbow_x = inv(side * (band_in - 0.15))) %>%
    transmute(species,
              x = map2(dn_mean, elbow_x, ~ c(.x, .y)) %>% map2(lab_x, ~ c(.x, .y)),
              y = map2(in_mean, in_mean, ~ c(.x, .y)) %>% map2(lab_y, ~ c(.x, .y))) %>%
    unnest(c(x, y))

  # ---- plot ----
  fig_step <- ggplot(plot_data, aes(dn_mean, in_mean)) +
    geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.5) +
    geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.5) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.5, linetype = "dashed") +
    geom_point(aes(colour = class, fill = interaction(class, any_cred)), shape = 21, size = 2.3, stroke = 0.7) +
    geom_text(data = quad_lab, aes(x, y, label = label, hjust = hjust, vjust = vjust),
              size = 3.4, lineheight = 0.95, colour = "grey20", fontface = "bold") +
    inset_layers +
    geom_path(data = lab_lines, aes(x, y, group = species), colour = "grey45", linewidth = 0.35) +
    geom_label(data = lab_data, aes(x = lab_x, y = lab_y, label = label, hjust = hjust), size = 3, lineheight = 0.9,
               linewidth = 0.25, label.padding = unit(0.18, "lines"), fill = "white") +
    scale_colour_manual(values = class_colours, breaks = names(class_colours),
                        labels = sprintf("%s (%d)", names(class_colours),
                                         as.integer(table(plot_data$class)[names(class_colours)])),
                        name = "Overall tracking slope") +
    scale_fill_manual(values = c(setNames(class_colours, paste0(names(class_colours), ".At least one step credible")),
                                 setNames(rep("white", 3), paste0(names(class_colours), ".Neither step credible"))),
                      guide = "none") +
    scale_x_continuous(trans = asinh_trans, breaks = axis_breaks, limits = c(-lim_plot, lim_plot), expand = c(0, 0)) +
    scale_y_continuous(trans = asinh_trans, breaks = axis_breaks, limits = c(-lim_plot, lim_plot), expand = c(0, 0)) +
    coord_equal() +
    labs(x = "Trailing-edge step: Neutral − Decreasing suitability (% per year)",
         y = "Leading-edge step: Increasing suitability − Neutral (% per year)",
         caption = sprintf(paste0("Insets: mean route trend in Decreasing (D), Neutral (N) and Increasing (I) suitability; blue = step up, red = step down.\n",
                                  "Filled: at least one step credible (P ≥ %.2f); open: neither. ",
                                  "%d of %d native species plotted (%d lack a Decreasing or Increasing category)."),
                           class_threshold, n_plot, n_all, n_all - n_plot)) +
    guides(colour = guide_legend(override.aes = list(fill = class_colours, size = 3))) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom",
          legend.title = element_text(face = "bold"),
          axis.title = element_text(face = "bold", size = 12),
          axis.text = element_text(size = 10.5),
          plot.caption = element_text(hjust = 0, size = 8.5, colour = "grey30"))

  # ---- save ----
  save_plot <- function(p, stem, width, height) {
    for (ext in c("png", "pdf")) {
      out_file <- file.path(fig_dir, paste0(stem, ".", ext))
      if (ext == "pdf") {
        ggsave(out_file, p, width = width, height = height, bg = "white", device = cairo_pdf)   # cairo: keeps the minus signs
      } else {
        ggsave(out_file, p, width = width, height = height, dpi = 300, bg = "white")
      }
      cat("Saved:", out_file, "\n")
    }
  }
  save_plot(fig_step, fig_stem, width = 8.5, height = 9)
}
