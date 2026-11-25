# 4c_sdm_tracking_meta.R
#
# Step 3 of the SDM-tracking analysis (Methods, Combining species: hierarchical
# meta-analysis). Combines the species-level quantities from 4b across species.
# Fast (one row per species); needs Stan but not the trend-model fits.
#
# Model (models/sdm_tracking_meta.stan), one fit per quantity:
#   y[s] ~ Normal(theta[s], se[s]),   theta[s] ~ Student-t(nu, mu + a[g[s]], tau)
#   y[s], se[s]  species s's estimate and its standard error (posterior mean and SD from 4b)
#   theta[s]     species s's true value
#   mu           the average species
#   tau          between-species scale of true values, net of estimation error
#   nu           Student-t degrees of freedom (tail weight; estimated)
#   a[g]         optional group effects (traits), partially pooled with SD sigma
# The Student-t is used throughout because species' values are heavy-tailed: a few
# species track very strongly, and under a Normal they pull mu up and inflate tau.
# A Normal version of the overall tracking-slope model is fitted only for comparison.
# The marginal form (theta integrated out) gives mu, tau, nu, group effects and the
# leave-one-species-out comparisons (PSIS-LOO); the species-level form keeps theta[s]
# as parameters to classify species (not LOO-comparable).
#
# Sections (as in the Methods):
#   1. Tracking slope, average species: mu, tau, nu, P(mu > 0), % of species with a negative
#      slope; the Normal model for comparison
#   2. Species classification: Tracks if P(theta > 0) >= class_threshold, Contrary if
#      P(theta < 0) >= class_threshold, Inconclusive otherwise
#   3. Worsening half (classes 1-4) and improving half (classes 4-7) of the gradient
#   4. Category means (Decreasing suitability / Neutral / Increasing suitability), fitted
#      separately per category; species sets "main" (each category uses the species with
#      routes in it) and "common" (species with routes in all three); differences from
#      Neutral on the common set, with se = sqrt(se_a^2 + se_b^2) (the categories use
#      disjoint routes, so their draws are close to independent; an approximation);
#      and the same for each single SDM class 1-7 (Figure 2's per-class estimates)
#   5. Species traits: one model per trait (habitat group, family, order + family, diet;
#      migration, Red List and invasive status as fixed effects, too few levels to
#      estimate a between-group SD), each compared with the no-trait model by PSIS-LOO
#
# Reads:  output/sdm_tracking/sdm_tracking_stage1_species_<suffix>.csv   (4b)
#         data/spp_names_codes_group_aou_taxa_redlist_diet_migration_invasive.csv  (traits, non-native list)
#         models/sdm_tracking_meta.stan
# Writes (output/sdm_tracking/; "_native" added when native_only):
#   sdm_tracking_stage2_SUMMARY_<suffix>.txt   <- the file to read / share; all key numbers
#   sdm_tracking_stage2_models_<suffix>.csv    one row per slope model: mu, tau, sigmas, LOO, diagnostics
#   sdm_tracking_stage2_groups_<suffix>.csv    one row per trait model x group level
#   sdm_tracking_stage2_species_<suffix>.csv   one row per species: estimated and shrunken slopes, class
#   sdm_tracking_category_levels_<suffix>.csv  one row per category / class quantity x species set
#
# Needs: tidyverse, cmdstanr (with CmdStan installed), posterior, loo.

library(here)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(loo)

here::i_am("4c_sdm_tracking_meta.R")

# Settings -------------------------------------------------------------------
model_tag <- "anthro"; scenario <- "rcp45"; firstYear <- 2010; lastYear <- 2025
# Must match 4a and 4b: integer(0) = main analysis, c(7) = colonization excluded ("_no7").
exclude_classes <- integer(0)
# TRUE = drop US-RIIS listed non-native species (and the invasive-status trait model);
# outputs get a "_native" tag. FALSE = all species (sensitivity analysis, s2_non_native.R).
native_only <- TRUE

class_threshold   <- 0.90  # posterior probability needed to call a species/group Tracks or Contrary
prior_scale_slope <- 1     # %/yr per class step (observed scale ~0.2)
prior_scale_mean  <- 3     # %/yr, for the category means (observed scale ~1)
chains <- 4; iter_warmup <- 1000; iter_sampling <- 2000; adapt_delta <- 0.99
seed      <- 20261014      # slope models
seed_mean <- 20261112      # category-mean models

run_tag <- if (length(exclude_classes) > 0) paste0("_no", paste(sort(exclude_classes), collapse = "")) else ""
stage1_suffix <- paste0(paste(model_tag, scenario, firstYear, lastYear, sep = "_"), run_tag)
suffix <- paste0(stage1_suffix, if (native_only) "_native" else "")
stage1_csv <- here::here("output", "sdm_tracking", paste0("sdm_tracking_stage1_species_", stage1_suffix, ".csv"))
traits_csv <- here::here("data", "spp_names_codes_group_aou_taxa_redlist_diet_migration_invasive.csv")
stan_file  <- here::here("models", "sdm_tracking_meta.stan")
out_dir    <- here::here("output", "sdm_tracking")

out_txt     <- file.path(out_dir, paste0("sdm_tracking_stage2_SUMMARY_", suffix, ".txt"))
out_models  <- file.path(out_dir, paste0("sdm_tracking_stage2_models_", suffix, ".csv"))
out_groups  <- file.path(out_dir, paste0("sdm_tracking_stage2_groups_", suffix, ".csv"))
out_species <- file.path(out_dir, paste0("sdm_tracking_stage2_species_", suffix, ".csv"))
out_levels  <- file.path(out_dir, paste0("sdm_tracking_category_levels_", suffix, ".csv"))

# Data -----------------------------------------------------------------------
stage1 <- read.csv(stage1_csv, stringsAsFactors = FALSE)
if (!all(c("mean_neutral_mean", "mean_class1_mean") %in% names(stage1))) {
  stop(basename(stage1_csv), " has no category / class means: rerun 4b_sdm_tracking_species.R")
}
traits <- read.csv(traits_csv, stringsAsFactors = FALSE) %>%
  filter(in_bbs %in% TRUE) %>% distinct(Code, .keep_all = TRUE) %>%
  transmute(species_code = Code, group = Group, order, family,
            diet = case_when(diet %in% c("Seed", "Plant", "Herbivore") ~ "Plant-eater",
                             diet %in% c("Vertebrate", "Carnivore", "Scavenger") ~ "Meat-eater", TRUE ~ diet),
            migration,
            redlist = case_when(redlistCategory %in% c("Vulnerable", "Endangered", "Critically Endangered") ~ "Threatened (VU + EN)",
                                TRUE ~ redlistCategory),
            invasive = if_else(invasive == "Not listed", "Native (not listed)", invasive))
d_all <- stage1 %>% left_join(traits, by = "species_code")
if (native_only) {
  dropped <- d_all %>% filter(invasive != "Native (not listed)")
  cat("Dropping", nrow(dropped), "non-native species:", paste(sort(dropped$species), collapse = ", "), "\n")
  d_all <- d_all %>% filter(invasive == "Native (not listed)")
}
# Species of the analysis: those with a tracking slope (routes in >= min_classes classes, from 4b)
d <- d_all %>% filter(status == "ok", is.finite(slope_mean), slope_sd > 0)
cat("Species in the analysis:", nrow(d), "of", nrow(d_all), "\n")

mod <- cmdstan_model(stan_file)

# Helpers --------------------------------------------------------------------
idx <- function(v) { f <- factor(v); list(i = as.integer(f), K = nlevels(f), levels = levels(f)) }

# Fit the meta-analysis to column ycol (estimate) with standard error secol.
# student = 1: Student-t species distribution (0 = Normal, comparison only)
# explicit = 1: species-level Student-t (theta[s] as parameters; for classification, not LOO)
# fixed1 = 1: grouping 1 as separate effects without a between-group SD (2-3 level traits)
fit_meta <- function(d, ycol, secol, grp1 = NULL, grp2 = NULL, student = 1, prior_scale,
                     explicit = 0, fixed1 = 0, seed_use = seed) {
  d <- d %>% filter(is.finite(.data[[ycol]]), is.finite(.data[[secol]]), .data[[secol]] > 0)
  g1 <- if (is.null(grp1)) list(i = rep(0L, nrow(d)), K = 0L, levels = character()) else idx(d[[grp1]])
  g2 <- if (is.null(grp2)) list(i = rep(0L, nrow(d)), K = 0L, levels = character()) else idx(d[[grp2]])
  run <- function() mod$sample(data = list(N = nrow(d), y = d[[ycol]], se = d[[secol]],
                                           K1 = g1$K, g1 = g1$i, K2 = g2$K, g2 = g2$i,
                                           student = student, explicit_t = explicit, fixed1 = fixed1,
                                           prior_scale = prior_scale),
                               chains = chains, parallel_chains = chains, iter_warmup = iter_warmup,
                               iter_sampling = iter_sampling, adapt_delta = adapt_delta, seed = seed_use,
                               refresh = 0, show_messages = FALSE)
  fit <- run()
  # On Windows a freshly compiled model executable can be blocked on its first launch
  # (e.g. by antivirus scanning), so that no chain starts. Retry once after a pause.
  if (all(fit$return_codes() != 0)) {
    message("No chains finished; retrying once")
    Sys.sleep(10)
    fit <- run()
    if (all(fit$return_codes() != 0)) stop("No chains finished on retry; rerun with show_messages = TRUE to see why")
  }
  list(fit = fit, g1 = g1, g2 = g2, d = d, ycol = ycol, explicit = explicit, fixed1 = fixed1)
}

q <- function(v) c(mean = mean(v), lo = unname(quantile(v, 0.025)), hi = unname(quantile(v, 0.975)))
classify <- function(p_pos) case_when(p_pos >= class_threshold ~ "Tracks",
                                      1 - p_pos >= class_threshold ~ "Contrary",
                                      TRUE ~ "Inconclusive")

# One row of results for a fitted slope model (mu, tau, sigmas, nu, LOO, diagnostics)
summarise_model <- function(m, name, response) {
  fit <- m$fit
  dr <- fit$draws(c("mu", "tau", "sigma1", "sigma2", "nu"), format = "draws_df")
  pooled1 <- m$g1$K > 0 && m$fixed1 == 0
  student <- all(is.finite(dr$nu))
  sum_vars <- c("mu", "tau", if (m$g1$K > 0) "z1", if (pooled1) "sigma1", if (m$g2$K > 0) c("sigma2", "z2"),
                if (student) "nu")
  sm <- fit$summary(sum_vars) %>% filter(!is.na(rhat))
  diag <- fit$diagnostic_summary(quiet = TRUE)
  lo <- suppressWarnings(fit$loo(variables = "log_lik", cores = 1))
  # share of species whose true value is negative, under the fitted species distribution
  neg <- if (student) pt(-dr$mu / dr$tau, df = dr$nu) else pnorm(-dr$mu / dr$tau)
  row <- tibble(model = name, response = response, n_species = nrow(m$d),
                mu_mean = mean(dr$mu), mu_lo = q(dr$mu)[["lo"]], mu_hi = q(dr$mu)[["hi"]],
                p_mu_pos = mean(dr$mu > 0),
                tau_mean = mean(dr$tau), tau_lo = q(dr$tau)[["lo"]], tau_hi = q(dr$tau)[["hi"]],
                pct_species_negative = 100 * mean(neg),
                group_effects = if (m$g1$K == 0) "none" else if (m$fixed1 == 1) "fixed" else "partially pooled",
                sigma1_mean = if (pooled1) mean(dr$sigma1) else NA, sigma1_lo = if (pooled1) q(dr$sigma1)[["lo"]] else NA,
                sigma1_hi = if (pooled1) q(dr$sigma1)[["hi"]] else NA,
                sigma2_mean = if (m$g2$K > 0) mean(dr$sigma2) else NA, sigma2_lo = if (m$g2$K > 0) q(dr$sigma2)[["lo"]] else NA,
                sigma2_hi = if (m$g2$K > 0) q(dr$sigma2)[["hi"]] else NA,
                nu_mean = if (student) mean(dr$nu) else NA,
                elpd_loo = lo$estimates["elpd_loo", "Estimate"], elpd_loo_se = lo$estimates["elpd_loo", "SE"],
                pareto_k_over_0.7 = sum(lo$diagnostics$pareto_k > 0.7),
                max_rhat = max(sm$rhat, na.rm = TRUE), min_ess_bulk = min(sm$ess_bulk, na.rm = TRUE),
                divergences = sum(diag$num_divergent))
  # posterior predictive check: does the spread of simulated species values match the observed?
  yrep <- fit$draws("y_rep", format = "draws_matrix")
  row$ppc_p_sd <- mean(apply(yrep, 1, sd) >= sd(m$d[[m$ycol]]))
  if (m$explicit == 1) {   # conditional log_lik: not comparable in LOO
    row$elpd_loo <- NA; row$elpd_loo_se <- NA; row$pareto_k_over_0.7 <- NA
  }
  list(row = row, loo = if (m$explicit == 1) NULL else lo)
}

# Group-level values (mu + group effect) for each level of a trait
group_table <- function(m, name, which = 1) {
  g <- if (which == 1) m$g1 else m$g2
  if (g$K == 0) return(NULL)
  dr <- m$fit$draws(c("mu", paste0("a", which)), format = "draws_matrix")
  mu <- dr[, "mu"]
  map_dfr(seq_len(g$K), function(k) {
    v <- as.vector(mu + dr[, paste0("a", which, "[", k, "]")])
    tibble(model = name, grouping = if (which == 1) "g1" else "g2", level = g$levels[k],
           n_species = sum(g$i == k),
           mean = mean(v), lo = q(v)[["lo"]], hi = q(v)[["hi"]], p_pos = mean(v > 0),
           class = classify(mean(v > 0)))
  })
}

# Each species' shrunken value (theta) from a fitted model
species_theta <- function(m) {
  th <- m$fit$draws("theta", format = "draws_matrix")
  tibble(species_code = m$d$species_code,
         theta_mean = colMeans(th),
         theta_lo = apply(th, 2, quantile, 0.025), theta_hi = apply(th, 2, quantile, 0.975),
         theta_p_pos = colMeans(th > 0))
}

fits <- list(); rows <- list(); loos <- list(); groups <- list()
fit_slope_model <- function(name, grp1 = NULL, grp2 = NULL, student = 1, explicit = 0, fixed1 = 0) {
  cat("Fitting:", name, "\n")
  m <- fit_meta(d, "slope_mean", "slope_sd", grp1, grp2, student, prior_scale_slope,
                explicit = explicit, fixed1 = fixed1)
  r <- summarise_model(m, name, "Tracking slope")
  fits[[name]]   <<- m
  rows[[name]]   <<- r$row
  if (!is.null(r$loo)) loos[[name]] <<- r$loo
  groups[[name]] <<- bind_rows(group_table(m, name, 1), group_table(m, name, 2))
  invisible(m)
}

# =============================================================================
# 1. Tracking slope: the average species
# =============================================================================
t_name <- "Overall (species-level)"
n_name <- "Overall (Normal)"
fit_slope_model("Overall")                     # Student-t, marginal: mu, tau, nu; LOO reference
fit_slope_model(n_name, student = 0)           # Normal: comparison only

# =============================================================================
# 2. Species classification (species-level Student-t)
# =============================================================================
fit_slope_model(t_name, explicit = 1)

# =============================================================================
# 3. Worsening and improving halves of the gradient
# =============================================================================
side <- list()
top_class <- max(setdiff(1:7, exclude_classes))
for (resp in c("slope_low", "slope_up")) {
  lab <- c(slope_low = "Worsening-side slope (classes 1-4)",
           slope_up = paste0("Improving-side slope (classes 4-", top_class, ")"))[[resp]]
  cat("Fitting:", lab, "\n")
  m <- fit_meta(d, paste0(resp, "_mean"), paste0(resp, "_sd"), prior_scale = prior_scale_slope)
  side[[resp]] <- summarise_model(m, "Overall", lab)$row
}

# =============================================================================
# 4. Category means
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

# =============================================================================
# 5. Species traits
# =============================================================================
trait_specs <- list(   # name, grouping 1, grouping 2, fixed group effects
  list("Habitat group",      "group",     NULL,     0),
  list("Family",             "family",    NULL,     0),
  list("Order + family",     "order",     "family", 0),
  list("Diet (8)",           "diet",      NULL,     0),
  list("Migration",          "migration", NULL,     1),
  list("Red List (3)",       "redlist",   NULL,     1),
  list("Invasive (US-RIIS)", "invasive",  NULL,     1)
)
if (native_only) trait_specs <- Filter(function(s) s[[1]] != "Invasive (US-RIIS)", trait_specs)
for (s in trait_specs) fit_slope_model(s[[1]], s[[2]], s[[3]], fixed1 = s[[4]])

# Collect results ------------------------------------------------------------
side_tab   <- bind_rows(unname(side))
models_out <- bind_rows(c(unname(rows), unname(side)))
# Each model vs the no-trait Student-t "Overall" model: positive elpd_diff = the model
# predicts held-out species better (for "Overall (Normal)", negative = the Normal is worse).
cmp_tab <- map_dfr(names(loos), function(nm) {
  d_i <- loos[[nm]]$pointwise[, "elpd_loo"] - loos[["Overall"]]$pointwise[, "elpd_loo"]
  tibble(model = nm, response = "Tracking slope",
         elpd_diff = if (nm == "Overall") 0 else sum(d_i),
         se_diff   = if (nm == "Overall") 0 else sqrt(length(d_i)) * sd(d_i))
})
models_out <- models_out %>% left_join(cmp_tab, by = c("model", "response"))
groups_out <- bind_rows(unname(groups))

# Species' shrunken slopes: the species-level Student-t model (main classification) and
# the Normal model for comparison. (In the marginal Student-t models theta is only a
# Normal-model approximation, so it is not used for classification.)
spt <- species_theta(fits[[t_name]]) %>% rename_with(~ paste0(.x, "_studentt"), -species_code)
spn <- species_theta(fits[[n_name]]) %>% rename_with(~ paste0(.x, "_normal"), -species_code)
species_out <- d %>%
  select(species, species_code, group, family, order, diet, migration, redlist, invasive, n_routes, n_classes,
         starts_with("n_class"), slope_mean, slope_sd, slope_p_pos, slope_low_mean, slope_up_mean) %>%
  left_join(spt, by = "species_code") %>% left_join(spn, by = "species_code") %>%
  mutate(class_stage1 = classify(slope_p_pos),
         class_studentt = classify(theta_p_pos_studentt),
         class_normal = classify(theta_p_pos_normal)) %>%
  arrange(theta_mean_studentt)

write.csv(models_out, out_models, row.names = FALSE)
write.csv(groups_out, out_groups, row.names = FALSE)
write.csv(species_out, out_species, row.names = FALSE)
write.csv(as.data.frame(levels_out), out_levels, row.names = FALSE)

# Plain-text summary ---------------------------------------------------------
con <- file(out_txt, open = "wt", encoding = "UTF-8")
w <- function(...) { cat(..., "\n", sep = "", file = con); cat(..., "\n", sep = "") }
f2 <- function(x, d = 2) formatC(x, format = "f", digits = d)
ci <- function(m, lo, hi, d = 3) paste0(f2(m, d), " [", f2(lo, d), ", ", f2(hi, d), "]")
mr <- function(name) models_out %>% filter(model == name, response == "Tracking slope")

w("SDM TRACKING - META-ANALYSIS SUMMARY (", suffix, ")   ", format(Sys.time(), "%Y-%m-%d %H:%M"))
w("Tracking slope = change in route trend (%/yr) per SDM class step; positive = tracks the projection.")
w("Species in the analysis: ", nrow(d), " (4b rows: ", nrow(d_all), "; status counts: ",
  paste(names(table(d_all$status)), table(d_all$status), sep = " = ", collapse = "; "), ")")
w("Student-t species distribution throughout. Class threshold: P >= ", class_threshold, ". Intervals are 95% credible intervals.")
w("")
w("== 1. TRACKING SLOPE: THE AVERAGE SPECIES ==")
o <- mr("Overall")
w("mu (average species' slope): ", ci(o$mu_mean, o$mu_lo, o$mu_hi), "   P(mu > 0) = ", f2(o$p_mu_pos, 4))
w("tau (between-species scale of true slopes): ", ci(o$tau_mean, o$tau_lo, o$tau_hi), ";  nu = ", f2(o$nu_mean, 1))
w("Expected % of species with a negative true slope: ", f2(o$pct_species_negative, 1), "%")
w("Raw (unpooled) species slopes: median ", f2(median(d$slope_mean), 3))
on <- mr(n_name)
w("Normal model (comparison): mu ", ci(on$mu_mean, on$mu_lo, on$mu_hi), ", tau (SD) ", f2(on$tau_mean, 3),
  ", % species negative ", f2(on$pct_species_negative, 1), "%;  elpd Student-t minus Normal: ",
  f2(-on$elpd_diff, 1), " (se ", f2(on$se_diff, 1), ")")
w("")
w("== 2. SPECIES CLASSIFICATION (species-level Student-t) ==")
for (cl in c("class_studentt", "class_normal", "class_stage1")) {
  tb <- table(factor(species_out[[cl]], levels = c("Tracks", "Inconclusive", "Contrary")))
  w(sprintf("%-15s Tracks %3d | Inconclusive %3d | Contrary %3d", sub("class_", "", cl), tb[1], tb[2], tb[3]))
}
w("Contrary species (Normal-model class in brackets):")
ctr <- species_out %>% filter(class_studentt == "Contrary") %>% arrange(desc(theta_mean_studentt))
if (nrow(ctr) == 0) w("  none") else for (i in seq_len(nrow(ctr))) {
  x <- ctr[i, ]
  w(sprintf("  %-30s %-20s slope %s  P(>0)=%.3f  routes %4d  [%s]", x$species, x$group,
            ci(x$theta_mean_studentt, x$theta_lo_studentt, x$theta_hi_studentt), x$theta_p_pos_studentt,
            x$n_routes, x$class_normal))
}
w("Strongest-tracking 10 species:")
top <- species_out %>% arrange(desc(theta_mean_studentt)) %>% head(10)
for (i in seq_len(nrow(top))) w(sprintf("  %-30s %-20s slope %.3f  P(>0)=%.3f", top$species[i], top$group[i], top$theta_mean_studentt[i], top$theta_p_pos_studentt[i]))
w("")
w("== 3. WORSENING vs IMPROVING HALF OF THE GRADIENT ==")
for (r in seq_len(nrow(side_tab))) {
  x <- side_tab[r, ]
  w(sprintf("%-38s n = %3d  mu %s  P(mu>0) = %.3f  tau %.3f  nu %.1f", x$response, x$n_species,
            ci(x$mu_mean, x$mu_lo, x$mu_hi), x$p_mu_pos, x$tau_mean, x$nu_mean))
}
w("")
w("== 4. CATEGORY MEANS (average species' trend, %/yr) ==")
for (r in seq_len(nrow(levels_out))) {
  x <- levels_out[r, ]
  w(sprintf("%-7s %-38s n = %3d  mu %s  P(mu<0) = %.3f  tau %.2f  nu %.1f", x$species_set, x$quantity, x$n_species,
            ci(x$mu, x$mu_lo, x$mu_hi, 2), x$p_mu_neg, x$tau, x$nu))
}
w("")
w("== 5. SPECIES TRAITS (model comparison) ==")
w("sigma = SD of true group effects (near 0 = groups do not differ). elpd_diff: LOO improvement over the no-trait Overall model ",
  "(positive = trait helps; |diff| < ~2*se_diff = no clear difference).")
mt <- models_out %>% filter(response == "Tracking slope", !model %in% c(t_name, n_name))
for (r in seq_len(nrow(mt))) {
  x <- mt[r, ]
  w(sprintf("%-22s sigma1 %-22s sigma2 %-22s elpd_diff %s  tau %.3f",
            x$model,
            if (x$group_effects == "fixed") "(fixed effects)" else if (is.na(x$sigma1_mean)) "-" else sprintf("%.3f [%.3f,%.3f]", x$sigma1_mean, x$sigma1_lo, x$sigma1_hi),
            if (is.na(x$sigma2_mean)) "-" else sprintf("%.3f [%.3f,%.3f]", x$sigma2_mean, x$sigma2_lo, x$sigma2_hi),
            sprintf("%6.1f (se %4.1f)", x$elpd_diff, x$se_diff), x$tau_mean))
}
w("Group-level slopes (mu + group effect):")
for (nm in intersect(c("Habitat group", "Diet (8)", "Migration", "Red List (3)", "Invasive (US-RIIS)", "Family"), names(fits))) {
  gt <- groups_out %>% filter(model == nm, grouping == "g1") %>% arrange(desc(mean))
  if (nm == "Family") gt <- gt %>% filter(n_species >= 5)
  w("-- ", nm, if (nm == "Family") " (families with >= 5 species)" else "", " --")
  col <- c(`Habitat group` = "group", `Diet (8)` = "diet", Migration = "migration",
           `Red List (3)` = "redlist", `Invasive (US-RIIS)` = "invasive", Family = "family")[[nm]]
  for (i in seq_len(nrow(gt))) {
    sc <- species_out %>% filter(.data[[col]] == gt$level[i])
    tb <- table(factor(sc$class_studentt, levels = c("Tracks", "Inconclusive", "Contrary")))
    w(sprintf("  %-28s n=%3d  %s  P(>0)=%.3f  %-12s species T/I/C %d/%d/%d", gt$level[i], gt$n_species[i],
              ci(gt$mean[i], gt$lo[i], gt$hi[i]), gt$p_pos[i], gt$class[i], tb[1], tb[2], tb[3]))
  }
}
w("")
w("== 6. DIAGNOSTICS ==")
for (r in seq_len(nrow(models_out))) {
  x <- models_out[r, ]
  w(sprintf("%-26s %-38s max Rhat %.3f  min ESS %5.0f  divergences %d  Pareto k>0.7: %s  PPC p(sd) %.2f",
            x$model, x$response, x$max_rhat, x$min_ess_bulk, x$divergences,
            if (is.na(x$pareto_k_over_0.7)) "n/a" else x$pareto_k_over_0.7, x$ppc_p_sd))
}
w("Category means: max Rhat ", f2(max(levels_out$max_rhat), 3), ", divergences ", sum(levels_out$divergences))
close(con)
cat("\nWrote:\n ", out_txt, "\n ", out_models, "\n ", out_groups, "\n ", out_species, "\n ", out_levels, "\n")
