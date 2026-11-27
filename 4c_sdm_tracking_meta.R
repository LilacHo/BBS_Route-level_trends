# 4c_sdm_tracking_meta.R
#
# Step 3 of the SDM-tracking analysis (Methods, Combining species: hierarchical
# meta-analysis). Combines the species-level quantities from 4b across species.
# Fast (one row per species); needs Stan but not the trend-model fits.
# Species traits are in 4d_sdm_tracking_traits.R.
#
# The model, data preparation and analysis settings (priors, sampler, seeds, class
# threshold) are shared with 4d in helper/sdm_tracking_meta_common.R:
#   y[s] ~ Normal(theta[s], se[s]),   theta[s] ~ Student-t(nu, mu, tau)
# Student-t throughout. The comparison with a Normal species distribution
# (Supplementary Note 3) is in s3_sdm_tracking.R.
#
# Sections (as in the Methods):
#   1. Tracking slope, average species: mu, tau, nu, P(mu > 0), % of species with a negative
#      slope
#   2. Species classification: Tracks if P(theta > 0) >= class_threshold, Contrary if
#      P(theta < 0) >= class_threshold, Inconclusive otherwise
#   3. Category means (Decreasing suitability / Neutral / Increasing suitability), fitted
#      separately per category; species sets "main" (each category uses the species with
#      routes in it) and "common" (species with routes in all three); differences from
#      Neutral on the common set, with se = sqrt(se_a^2 + se_b^2) (the categories use
#      disjoint routes, so their draws are close to independent; an approximation);
#      and the same for each single SDM class 1-7 (Figure 2's per-class estimates)
#
# Reads:  output/sdm_tracking/sdm_tracking_stage1_species_<suffix>.csv   (4b)
#         data/spp_names_codes_group_aou_taxa_redlist_diet_migration_invasive.csv  (traits, non-native list)
#         models/sdm_tracking_meta.stan
# Writes (output/sdm_tracking/; "_native" added when native_only):
#   sdm_tracking_stage2_SUMMARY_<suffix>.txt   <- the file to read / share; all key numbers
#   sdm_tracking_stage2_models_<suffix>.csv    one row per slope model: mu, tau, nu, LOO, diagnostics
#   sdm_tracking_stage2_species_<suffix>.csv   one row per species: estimated and shrunken slopes, class
#   sdm_tracking_category_levels_<suffix>.csv  one row per category / class quantity x species set
#
# Next: 4d_sdm_tracking_traits.R (uses the species classes written here).
#
# Needs: tidyverse, cmdstanr (with CmdStan installed), posterior, loo.

library(here)
library(tidyverse)

here::i_am("4c_sdm_tracking_meta.R")

# Settings (must match 4a, 4b and 4d) ------------------------------------------
model_tag <- "anthro"; scenario <- "rcp45"; firstYear <- 2010; lastYear <- 2025
# integer(0) = main analysis, c(7) = colonization excluded ("_no7").
exclude_classes <- integer(0)
# TRUE = drop US-RIIS listed non-native species; outputs get a "_native" tag.
# FALSE = all species (sensitivity analysis, s2_non_native.R).
native_only <- TRUE

source(here::here("helper", "sdm_tracking_meta_common.R"))   # settings, data d, fit functions

out_txt     <- file.path(out_dir, paste0("sdm_tracking_stage2_SUMMARY_", suffix, ".txt"))
out_models  <- file.path(out_dir, paste0("sdm_tracking_stage2_models_", suffix, ".csv"))
out_species <- file.path(out_dir, paste0("sdm_tracking_stage2_species_", suffix, ".csv"))
out_levels  <- file.path(out_dir, paste0("sdm_tracking_category_levels_", suffix, ".csv"))

# =============================================================================
# 1. Tracking slope: the average species
# =============================================================================
t_name <- "Overall (species-level)"
fit_slope_model("Overall")                     # Student-t, marginal: mu, tau, nu

# =============================================================================
# 2. Species classification (species-level Student-t)
# =============================================================================
fit_slope_model(t_name, explicit = 1)

# =============================================================================
# 3. Category means
# =============================================================================
fit_level <- function(dd, ycol, secol) {
  m <- fit_meta(dd, ycol, secol, prior_scale = prior_scale_mean, seed_use = seed_mean)
  dr <- m$fit$draws(c("mu", "tau", "nu"), format = "draws_df")
  tibble(n_species = nrow(m$d), nu = median(dr$nu),
         mu = mean(dr$mu), mu_lo = q(dr$mu)[["lo"]], mu_hi = q(dr$mu)[["hi"]],
         p_mu_neg = mean(dr$mu < 0), tau = mean(dr$tau),
         max_rhat = max(m$fit$summary(c("mu", "tau"))$rhat),
         divergences = sum(m$fit$diagnostic_summary(quiet = TRUE)$num_divergent))
}
common <- d %>%
  filter(is.finite(mean_decreasing_mean), is.finite(mean_neutral_mean), is.finite(mean_increasing_mean)) %>%
  mutate(dn_mean = mean_decreasing_mean - mean_neutral_mean, dn_sd = sqrt(mean_decreasing_sd^2 + mean_neutral_sd^2),
         in_mean = mean_increasing_mean - mean_neutral_mean, in_sd = sqrt(mean_increasing_sd^2 + mean_neutral_sd^2))
cat_specs <- list(
  list("main",   "Decreasing suitability (classes 1-3)", d,      "mean_decreasing"),
  list("main",   "Neutral (class 4)",                    d,      "mean_neutral"),
  list("main",   "Increasing suitability (classes 5-7)", d,      "mean_increasing"),
  list("common", "Decreasing suitability (classes 1-3)", common, "mean_decreasing"),
  list("common", "Neutral (class 4)",                    common, "mean_neutral"),
  list("common", "Increasing suitability (classes 5-7)", common, "mean_increasing"),
  list("common", "Decreasing suitability - Neutral",     common, "dn"),
  list("common", "Increasing suitability - Neutral",     common, "in")
)
# Single SDM classes 1-7 (average species in each class; Figure 2's ticks and intervals),
# each using the species with routes in that class
class_labels <- c("extirpation", "worsening", "slightly worsening", "neutral",
                  "slightly improving", "improving", "colonization")
class_specs <- map(setdiff(1:7, exclude_classes), function(k)
  list("main", paste0("Class ", k, " (", class_labels[k], ")"), d, paste0("mean_class", k)))
levels_out <- map_dfr(c(cat_specs, class_specs), function(s) {
  cat("Fitting category mean:", s[[1]], s[[2]], "\n")
  fit_level(s[[3]], paste0(s[[4]], "_mean"), paste0(s[[4]], "_sd")) %>%
    mutate(species_set = s[[1]], quantity = s[[2]], .before = 1)
}) %>% mutate(across(where(is.double), unname))

# Collect results ------------------------------------------------------------
models_out <- bind_rows(unname(rows))

# Species' shrunken slopes from the species-level Student-t model (classification). (In the
# marginal Student-t model theta is only a Normal-model approximation, so it is not used.)
spt <- species_theta(fits[[t_name]]) %>% rename_with(~ paste0(.x, "_studentt"), -species_code)
species_out <- d %>%
  select(species, species_code, group, family, order, diet, migration, redlist, invasive, n_routes, n_classes,
         starts_with("n_class"), slope_mean, slope_sd, slope_p_pos) %>%
  left_join(spt, by = "species_code") %>%
  mutate(class_stage1 = classify(slope_p_pos),
         class_studentt = classify(theta_p_pos_studentt)) %>%
  arrange(theta_mean_studentt)

write.csv(models_out, out_models, row.names = FALSE)
write.csv(species_out, out_species, row.names = FALSE)
write.csv(as.data.frame(levels_out), out_levels, row.names = FALSE)

# Plain-text summary ---------------------------------------------------------
con <- file(out_txt, open = "wt", encoding = "UTF-8")
w <- function(...) { cat(..., "\n", sep = "", file = con); cat(..., "\n", sep = "") }
mr <- function(name) models_out %>% filter(model == name, response == "Tracking slope")

w("SDM TRACKING - META-ANALYSIS SUMMARY (", suffix, ")   ", format(Sys.time(), "%Y-%m-%d %H:%M"))
w("Tracking slope = change in route trend (%/yr) per SDM class step; positive = tracks the projection.")
w("Species in the analysis: ", nrow(d), " (4b rows: ", nrow(d_all), "; status counts: ",
  paste(names(table(d_all$status)), table(d_all$status), sep = " = ", collapse = "; "), ")")
w("Student-t species distribution throughout. Class threshold: P >= ", class_threshold, ". Intervals are 95% credible intervals.")
w("Species traits: see sdm_tracking_traits_SUMMARY_", suffix, ".txt (4d). Normal vs Student-t: s3_sdm_tracking.R.")
w("")
w("== 1. TRACKING SLOPE: THE AVERAGE SPECIES ==")
o <- mr("Overall")
w("mu (average species' slope): ", ci(o$mu_mean, o$mu_lo, o$mu_hi), "   P(mu > 0) = ", f2(o$p_mu_pos, 4))
w("tau (between-species scale of true slopes): ", ci(o$tau_mean, o$tau_lo, o$tau_hi), ";  nu = ", f2(o$nu_mean, 1))
w("Expected % of species with a negative true slope: ", f2(o$pct_species_negative, 1), "%")
w("Raw (unpooled) species slopes: median ", f2(median(d$slope_mean), 3))
w("")
w("== 2. SPECIES CLASSIFICATION (species-level Student-t) ==")
for (cl in c("class_studentt", "class_stage1")) {
  tb <- table(factor(species_out[[cl]], levels = c("Tracks", "Inconclusive", "Contrary")))
  w(sprintf("%-15s Tracks %3d | Inconclusive %3d | Contrary %3d", sub("class_", "", cl), tb[1], tb[2], tb[3]))
}
w("Contrary species:")
ctr <- species_out %>% filter(class_studentt == "Contrary") %>% arrange(desc(theta_mean_studentt))
if (nrow(ctr) == 0) w("  none") else for (i in seq_len(nrow(ctr))) {
  x <- ctr[i, ]
  w(sprintf("  %-30s %-20s slope %s  P(>0)=%.3f  routes %4d", x$species, x$group,
            ci(x$theta_mean_studentt, x$theta_lo_studentt, x$theta_hi_studentt), x$theta_p_pos_studentt,
            x$n_routes))
}
w("Strongest-tracking 10 species:")
top <- species_out %>% arrange(desc(theta_mean_studentt)) %>% head(10)
for (i in seq_len(nrow(top))) w(sprintf("  %-30s %-20s slope %.3f  P(>0)=%.3f", top$species[i], top$group[i], top$theta_mean_studentt[i], top$theta_p_pos_studentt[i]))
w("")
w("== 3. CATEGORY MEANS (average species' trend, %/yr) ==")
for (r in seq_len(nrow(levels_out))) {
  x <- levels_out[r, ]
  w(sprintf("%-7s %-38s n = %3d  mu %s  P(mu<0) = %.3f  tau %.2f  nu %.1f", x$species_set, x$quantity, x$n_species,
            ci(x$mu, x$mu_lo, x$mu_hi, 2), x$p_mu_neg, x$tau, x$nu))
}
w("")
w("== 4. DIAGNOSTICS ==")
diag_lines(models_out, w)
w("Category means: max Rhat ", f2(max(levels_out$max_rhat), 3), ", divergences ", sum(levels_out$divergences))
close(con)
cat("\nWrote:\n ", out_txt, "\n ", out_models, "\n ", out_species, "\n ", out_levels, "\n")
