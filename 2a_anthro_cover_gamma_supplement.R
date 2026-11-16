# 2a_anthro_cover_gamma_supplement.R
#
# Species-level effect of anthropogenic land cover (gamma) from the anthro
# model, for the Supplementary Information.
#
#   log(lambda[r,t]) = alpha[r] + beta[r]*t + gamma * A[r,t] + ...
#
# A[r,t] is the proportion of route r's buffer in anthropogenic cover in year t
# (0-1, not centred), so exp(gamma) is the multiplicative change in expected
# count between a route with no anthropogenic cover and one fully covered, and
# 100 * (exp(0.1 * gamma) - 1) is the % change per 10 percentage points of
# cover. gamma has a Normal(0, 1) prior.
#
# Reports:
#   0. gamma posterior summaries: gamma's row (mean, median, sd, mad, q5, q95,
#      rhat, ess_bulk, ess_tail) from every species' anthro *_summ_fit.rds,
#      plus whole-model convergence (max Rhat / min bulk ESS across every
#      parameter; converged = Rhat < 1.01 and ESS > 400). The fit files are
#      large and only exist where the models were fitted, so when rds_dir is
#      absent the saved table (posterior_csv) is read instead.
#   1. Covariate variation: how much A differs among routes vs. changes within
#      routes over 2010-2025 (gamma is identified mostly by the former).
#   2. Per species: gamma posterior mean, sd, 90% CI (q5/q95), and whether the
#      CI excludes zero (credible positive / credible negative / not credible).
#   3. Overall and per habitat group: counts by sign, % credible, mean and
#      median gamma, mean-median gap (flags skewed posteriors), and
#      DerSimonian-Laird pooled means -- fixed-effect (one true gamma per
#      group) and random-effects (between-species heterogeneity; 90% CI, tau,
#      I2) -- using each species' posterior mean and sd as estimate and SE.
#   4. Effect on route trends: paired difference in route trend (anthro - base,
#      %/yr) for every route present in both models, and its relation to gamma.
#
# Reads:  output/rds/<species>_iCAR_anthro_2010_2025_summ_fit.rds   (1a/1b), if present;
#           otherwise output/gamma/gamma_posterior_summary_2010_2025.csv
#         data/spp_names_codes_group_aou.csv   (species list and group)
#         output/species_routes_covariates/all_route_trends_all_species_anthro_2010_2025.csv
#           (2_generate_route_trend_csvs_covariates.R)
#         data/Anthro.csv
# Writes: output/gamma/gamma_posterior_summary_2010_2025.csv  gamma posterior + convergence, one row per species
#           (only rewritten when the *_summ_fit.rds files are available)
#         output/gamma/gamma_species_2010_2025.csv        one row per species
#         output/gamma/gamma_group_summary_2010_2025.csv  one row per group + "All species" (Table S1 and more)
#         output/gamma/gamma_key_numbers_2010_2025.csv    name/value pairs quoted in the text
#         output/gamma/figures/gamma_by_group_2010_2025.png      (Fig. S1)
#         output/gamma/figures/gamma_trend_shift_2010_2025.png   (Fig. S2)

library(here)
library(tidyverse)

here::i_am("2a_anthro_cover_gamma_supplement.R")

# Settings -------------------------------------------------------------------
firstYear <- 2010
lastYear  <- 2025
require_route_converged <- TRUE   # same rule as 4c / v1

rds_dir    <- here::here("output", "rds")
spp_csv    <- here::here("data", "spp_names_codes_group_aou.csv")
routes_csv <- here::here("output", "species_routes_covariates",
                         paste0("all_route_trends_all_species_anthro_", firstYear, "_", lastYear, ".csv"))
anthro_csv <- here::here("data", "Anthro.csv")

out_dir <- here::here("output", "gamma")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
posterior_csv <- file.path(out_dir, paste0("gamma_posterior_summary_", firstYear, "_", lastYear, ".csv"))

col_pos <- "#2a78d6"   # credible positive
col_neg <- "#e34948"   # credible negative
col_ns  <- "#a8a6a0"   # not credible
sign_cols <- c("Credible positive" = col_pos, "Not credible" = col_ns, "Credible negative" = col_neg)

group_labels <- c(
  aridlands = "Aridlands", arctic = "Arctic", boreal_forests = "Boreal forests",
  coastal = "Coastal", eastern_forests = "Eastern forests", generalists = "Generalists",
  grasslands = "Grasslands", marshlands = "Marshlands", subtropical_forests = "Subtropical forests",
  urban_suburban = "Urban/suburban", waterbirds = "Waterbirds", western_forests = "Western forests"
)

key <- list()   # name -> value, written to gamma_key_numbers
add_key <- function(name, value) key[[name]] <<- value

pct_per_10pp <- function(g) 100 * (exp(0.1 * g) - 1)

#' DerSimonian & Laird (1986) method-of-moments pooling of estimates yi with
#' sampling variances vi. Returns the fixed-effect mean (one true value shared
#' by all) and the random-effects mean (adds between-estimate variance tau2,
#' estimated from Cochran's Q and floored at 0), with SEs, tau2 and I2 (% of
#' total variation that is between-estimate heterogeneity, 0-100).
dl_random_effects <- function(yi, vi) {
  keep <- is.finite(yi) & is.finite(vi) & vi > 0
  yi <- yi[keep]; vi <- vi[keep]
  k  <- length(yi)
  if (k == 0) return(tibble(k = 0L, fe_mean = NA_real_, fe_se = NA_real_, re_mean = NA_real_,
                            re_se = NA_real_, tau2 = NA_real_, i2 = NA_real_))
  if (k == 1) return(tibble(k = 1L, fe_mean = yi, fe_se = sqrt(vi), re_mean = yi,
                            re_se = sqrt(vi), tau2 = 0, i2 = 0))
  w_fe    <- 1 / vi
  fe_mean <- sum(w_fe * yi) / sum(w_fe)
  Q    <- sum(w_fe * (yi - fe_mean)^2)
  C    <- sum(w_fe) - sum(w_fe^2) / sum(w_fe)
  tau2 <- max(0, (Q - (k - 1)) / C)
  i2   <- if (Q > 0) max(0, 100 * (Q - (k - 1)) / Q) else 0
  w_re <- 1 / (vi + tau2)
  tibble(k = k, fe_mean = fe_mean, fe_se = sqrt(1 / sum(w_fe)),
         re_mean = sum(w_re * yi) / sum(w_re), re_se = sqrt(1 / sum(w_re)),
         tau2 = tau2, i2 = i2)
}

# 1. Covariate variation -----------------------------------------------------
anthro <- read.csv(anthro_csv) %>%
  filter(!is.na(Anthro), year >= firstYear, year <= lastYear) %>%
  mutate(route_key = paste(StateNum, Route, sep = "-"))

route_cov <- anthro %>%
  group_by(route_key) %>%
  summarise(n_years = n(), mean_A = mean(Anthro),
            change_A = Anthro[which.max(year)] - Anthro[which.min(year)],
            .groups = "drop") %>%
  filter(n_years >= 2)

add_key("anthro_n_routes", nrow(route_cov))
add_key("anthro_route_mean_median", median(route_cov$mean_A))
add_key("anthro_route_mean_q10", quantile(route_cov$mean_A, 0.10, names = FALSE))
add_key("anthro_route_mean_q90", quantile(route_cov$mean_A, 0.90, names = FALSE))
add_key("anthro_among_route_sd", sd(route_cov$mean_A))
add_key("anthro_within_route_change_median", median(route_cov$change_A))
add_key("anthro_within_route_change_abs_q90", quantile(abs(route_cov$change_A), 0.90, names = FALSE))
add_key("anthro_pct_routes_increased", 100 * mean(route_cov$change_A > 0))

# 2. Per-species gamma ---------------------------------------------------------
summ_files <- if (dir.exists(rds_dir)) {
  list.files(rds_dir, pattern = paste0("_iCAR_anthro_", firstYear, "_", lastYear, "_summ_fit\\.rds$"))
} else character(0)

if (length(summ_files) > 0) {
  # Species list and group; species listed under two groups keep the first.
  target_spp <- read.csv(spp_csv, stringsAsFactors = FALSE) %>%
    filter(in_bbs == TRUE) %>%
    distinct(Common.Name, Code, .keep_all = TRUE) %>%
    arrange(Common.Name)

  read_gamma <- function(sp, code, grp) {
    sp_f <- gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)   # same as species_to_f()
    f <- file.path(rds_dir, paste0(sp_f, "_iCAR_anthro_", firstYear, "_", lastYear, "_summ_fit.rds"))
    if (!file.exists(f)) return(NULL)
    summ <- readRDS(f)
    g <- summ[summ$variable == "gamma1", ]
    if (nrow(g) != 1) return(NULL)
    max_rhat <- suppressWarnings(max(summ$rhat, na.rm = TRUE))
    min_ess  <- suppressWarnings(min(summ$ess_bulk, na.rm = TRUE))
    g %>%
      select(-variable) %>%
      mutate(species = sp, species_code = code, model = "anthro", group = grp,
             model_max_rhat = max_rhat, model_min_ess_bulk = min_ess,
             model_converged = is.finite(max_rhat) && is.finite(min_ess) &&
               max_rhat < 1.01 && min_ess > 400,
             .before = 1)
  }

  cat("Reading gamma from", length(summ_files), "anthro summ_fit.rds files in", rds_dir, "\n")
  posterior <- bind_rows(pmap(list(target_spp$Common.Name, target_spp$Code, target_spp$Group), read_gamma)) %>%
    select(any_of(c("species", "species_code", "model", "group", "mean", "median", "sd", "mad",
                    "q5", "q95", "rhat", "ess_bulk", "ess_tail",
                    "model_max_rhat", "model_min_ess_bulk", "model_converged")), everything())
  write.csv(posterior, posterior_csv, row.names = FALSE)
  cat("Wrote", nrow(posterior), "species to", posterior_csv, "\n")
} else {
  if (!file.exists(posterior_csv)) {
    stop("No *_summ_fit.rds files in ", rds_dir, " and no saved ", posterior_csv,
         " -- run this script once where the model fits are stored.")
  }
  cat("No fit files in", rds_dir, "-- reading saved", posterior_csv, "\n")
  posterior <- read.csv(posterior_csv, stringsAsFactors = FALSE)
}

gamma_sp <- posterior %>%
  filter(model == "anthro") %>%
  distinct(species_code, .keep_all = TRUE) %>%
  transmute(species, species_code, group,
            gamma = mean, gamma_median = median, gamma_sd = sd,
            gamma_q5 = q5, gamma_q95 = q95, rhat, ess_bulk, model_converged,
            sign = case_when(q5 > 0  ~ "Credible positive",
                             q95 < 0 ~ "Credible negative",
                             TRUE    ~ "Not credible"),
            sign = factor(sign, levels = names(sign_cols)),
            pct_change_per_10pp = pct_per_10pp(gamma))

n_sp <- nrow(gamma_sp)
add_key("n_species", n_sp)
add_key("n_not_converged", sum(!gamma_sp$model_converged %in% TRUE))
add_key("n_credible_positive", sum(gamma_sp$sign == "Credible positive"))
add_key("n_credible_negative", sum(gamma_sp$sign == "Credible negative"))
add_key("n_not_credible", sum(gamma_sp$sign == "Not credible"))
add_key("pct_credible_any", 100 * mean(gamma_sp$sign != "Not credible"))
add_key("gamma_median", median(gamma_sp$gamma))
add_key("gamma_q25", quantile(gamma_sp$gamma, 0.25, names = FALSE))
add_key("gamma_q75", quantile(gamma_sp$gamma, 0.75, names = FALSE))
add_key("gamma_min", min(gamma_sp$gamma))
add_key("gamma_max", max(gamma_sp$gamma))
add_key("gamma_sd_median", median(gamma_sp$gamma_sd))

# 3. Pooled mean, overall and by group ------------------------------------------
pool <- function(d) {
  p <- dl_random_effects(d$gamma, d$gamma_sd^2)
  z <- qnorm(0.95)   # 90% CI, matching the species-level q5/q95
  tibble(
    n_species   = nrow(d),
    n_positive  = sum(d$sign == "Credible positive"),
    n_negative  = sum(d$sign == "Credible negative"),
    n_not_cred  = sum(d$sign == "Not credible"),
    pct_credible_90 = 100 * mean(d$sign != "Not credible"),
    mean_gamma   = mean(d$gamma),
    median_gamma = median(d$gamma),
    mean_median_gap = mean(abs(d$gamma - d$gamma_median)),
    re_mean     = p$re_mean,
    re_lo90     = p$re_mean - z * p$re_se,
    re_hi90     = p$re_mean + z * p$re_se,
    tau         = sqrt(p$tau2),
    i2          = p$i2,
    fe_mean     = p$fe_mean,
    fe_lo90     = p$fe_mean - z * p$fe_se,
    fe_hi90     = p$fe_mean + z * p$fe_se
  )
}

overall <- pool(gamma_sp) %>% mutate(group = "All species", .before = 1)
by_group <- gamma_sp %>%
  group_by(group) %>%
  group_modify(~ pool(.x)) %>%
  ungroup()

group_summary <- bind_rows(overall, by_group %>% arrange(desc(n_species))) %>%
  mutate(group_label = recode(group, !!!group_labels),
         re_pct_per_10pp = pct_per_10pp(re_mean), .after = group)

add_key("re_mean", overall$re_mean)
add_key("re_lo90", overall$re_lo90)
add_key("re_hi90", overall$re_hi90)
add_key("re_pct_per_10pp", pct_per_10pp(overall$re_mean))
add_key("tau", overall$tau)
add_key("i2", overall$i2)

# 4. Effect on route trends (anthro - base) --------------------------------------
routes <- read.csv(routes_csv, stringsAsFactors = FALSE)
if (require_route_converged) routes <- routes %>% filter(route_converged %in% TRUE)
# A few rows have no route ID (upstream join gap, also dropped in v1) and
# can't be paired across models.
add_key("routes_missing_id_dropped", sum(is.na(routes$route)))
routes <- routes %>% filter(!is.na(route))

paired <- routes %>%
  select(species_code, route, model, trend) %>%
  pivot_wider(names_from = model, values_from = trend) %>%
  filter(!is.na(base), !is.na(anthro)) %>%
  mutate(diff = anthro - base)

add_key("paired_n_routes", nrow(paired))
add_key("paired_n_species", n_distinct(paired$species_code))
add_key("diff_mean", mean(paired$diff))
add_key("diff_median", median(paired$diff))
add_key("diff_sd", sd(paired$diff))
add_key("diff_abs_q95", quantile(abs(paired$diff), 0.95, names = FALSE))
add_key("pct_abs_diff_lt_0.1", 100 * mean(abs(paired$diff) < 0.1))
add_key("pct_abs_diff_lt_0.5", 100 * mean(abs(paired$diff) < 0.5))
add_key("cor_trend_base_anthro", cor(paired$base, paired$anthro))
add_key("fig_s2_clip_limit", quantile(abs(paired$diff), 0.995, names = FALSE))
add_key("fig_s2_n_beyond_clip", sum(abs(paired$diff) > quantile(abs(paired$diff), 0.995, names = FALSE)))

sp_shift <- paired %>%
  group_by(species_code) %>%
  summarise(n_routes_paired = n(), mean_trend_base = mean(base),
            mean_trend_anthro = mean(anthro), mean_trend_shift = mean(diff),
            .groups = "drop")

gamma_sp <- gamma_sp %>% left_join(sp_shift, by = "species_code")
sh <- gamma_sp %>% filter(!is.na(mean_trend_shift))
add_key("species_shift_median", median(sh$mean_trend_shift))
add_key("species_shift_abs_max", max(abs(sh$mean_trend_shift)))
add_key("cor_gamma_species_shift", cor(sh$gamma, sh$mean_trend_shift, method = "spearman"))

# Write tables -----------------------------------------------------------------
write.csv(gamma_sp, file.path(out_dir, paste0("gamma_species_", firstYear, "_", lastYear, ".csv")), row.names = FALSE)
write.csv(group_summary, file.path(out_dir, paste0("gamma_group_summary_", firstYear, "_", lastYear, ".csv")), row.names = FALSE)
key_df <- tibble(name = names(key), value = unlist(key))
write.csv(key_df, file.path(out_dir, paste0("gamma_key_numbers_", firstYear, "_", lastYear, ".csv")), row.names = FALSE)

# Fig. S1: gamma by habitat group -------------------------------------------------
fig1_drop_groups <- c("arctic")   # single species; still in the tables
grp_order <- by_group %>% filter(!group %in% fig1_drop_groups) %>% arrange(re_mean) %>% pull(group)
set.seed(1)   # reproducible jitter
plot_sp <- gamma_sp %>%
  filter(group %in% grp_order) %>%
  mutate(group = factor(group, levels = grp_order),
         y = as.integer(group) + runif(n(), -0.28, 0.28))
plot_grp <- by_group %>%
  filter(group %in% grp_order) %>%
  mutate(group = factor(group, levels = grp_order), y = as.integer(group))

fig1 <- ggplot() +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_point(data = plot_sp, aes(x = gamma, y = y, fill = sign),
             shape = 21, size = 2, stroke = 0.25, colour = "white", alpha = 0.85) +
  geom_errorbar(data = plot_grp, aes(xmin = re_lo90, xmax = re_hi90, y = y),
                orientation = "y", width = 0, linewidth = 1.1, colour = "grey10") +
  geom_point(data = plot_grp, aes(x = re_mean, y = y), shape = 23, size = 3.2,
             fill = "white", colour = "grey10", stroke = 0.9) +
  scale_fill_manual(values = sign_cols, name = "Species γ (90% CI)") +
  scale_y_continuous(breaks = seq_along(grp_order),
                     labels = paste0(recode(grp_order, !!!group_labels), " (",
                                     plot_grp$n_species[match(grp_order, plot_grp$group)], ")"),
                     expand = expansion(add = 0.5)) +
  labs(x = "γ (effect of anthropogenic-cover proportion on log expected count)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        legend.position = "top")

ggsave(file.path(fig_dir, paste0("gamma_by_group_", firstYear, "_", lastYear, ".png")),
       fig1, width = 8, height = 6.2, dpi = 300, device = ragg::agg_png, bg = "white")

# Fig. S2: route-level trend shift (anthro - base) ------------------------------------
lim <- quantile(abs(paired$diff), 0.995, names = FALSE)
fig2 <- ggplot(paired %>% filter(abs(diff) <= lim), aes(x = diff)) +
  geom_histogram(binwidth = lim / 60, fill = "#2a78d6", colour = "white", linewidth = 0.15) +
  geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.4) +
  labs(x = "Route trend, anthropogenic-cover model minus base model (%/yr)",
       y = "Routes") +
  scale_y_continuous(labels = scales::label_comma()) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, paste0("gamma_trend_shift_", firstYear, "_", lastYear, ".png")),
       fig2, width = 7, height = 4, dpi = 300, device = ragg::agg_png, bg = "white")

# Console summary ---------------------------------------------------------------------
cat("\n=== gamma (anthro model) ===\n")
print(as.data.frame(key_df), digits = 4, row.names = FALSE)
cat("\n")
print(as.data.frame(group_summary %>% select(-group)), digits = 3, row.names = FALSE)
cat("\nWritten to:", out_dir, "\n")
