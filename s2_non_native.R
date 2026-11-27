# s2_non_native.R
#
# Supplementary Note 2: exclusion of non-native species. Collects the checks
# and summaries behind Note 2 and behind the native-only numbers in main.tex.
#
# Non-native = listed in the lower-48 US-RIIS v2.0 (Widespread invasive or
# Established; column `invasive` from 0d_prepare_invasive.R).
#
# The native-only analyses themselves are run by 4a_sdm_tracking_routes.R (route
# counts), 4c_sdm_tracking_meta.R (category means and tracking models) and
# 4d_sdm_tracking_traits.R (trait models), all with native_only <- TRUE (outputs tagged "_native"). Stage 1
# (4b_sdm_tracking_species.R) is computed species by species with a fixed seed,
# so its output still holds the non-native species and is filtered here.
#
# Sections:
#   1. Non-native species in the species list
#   2. Route-level trends of the non-native species (anthro fit, converged routes)
#   3. Why five non-native species have no trend estimates (fitting log)
#   4. Non-native species in the tracking analysis (share of data, slopes)
#   5. Category means, native only (and vs all species when that file exists)
#   6. Stage-2 tracking models, native only (and vs all species when present)
#   7. SDM-matched route counts and pooled means, all vs native (main.tex Results)
#   8. Trait counts for the native species in the trait analysis (main.tex Methods)

library(here)
library(tidyverse)

here::i_am("s2_non_native.R")

model_tag <- "anthro"; scenario <- "rcp45"; firstYear <- 2010; lastYear <- 2025
suffix <- paste(model_tag, scenario, firstYear, lastYear, sep = "_")
min_classes <- 2   # species need routes in at least this many SDM classes (as in 4a/4b)

traits_csv <- here::here("data", "spp_names_codes_group_aou_taxa_redlist_diet_migration_invasive.csv")
trends_csv <- here::here("output", "species_routes_covariates",
                         paste0("all_route_trends_all_species_", model_tag, "_", firstYear, "_", lastYear, ".csv"))
fit_log    <- here::here("output", "files", "1c_species_iCAR_covariates log 20260827.txt")
sdm_dir    <- here::here("output", "species_routes_covariates", "per_species_sdm")
track_dir  <- here::here("output", "sdm_tracking")

read_if <- function(f) if (file.exists(f)) read_csv(f, show_col_types = FALSE) else NULL

traits <- read_csv(traits_csv, show_col_types = FALSE)

# 1. Non-native species -------------------------------------------------------
non_native <- traits %>%
  filter(in_bbs %in% TRUE, invasive != "Not listed") %>%
  select(Common.Name, Code, Group, invasive) %>%
  arrange(invasive, Common.Name)
cat("\n== 1. Non-native species (US-RIIS lower 48) ==\n")
print(table(traits$invasive, useNA = "ifany"))
print(as.data.frame(non_native))
nn_codes <- non_native$Code

# 2. Route-level trends of the non-native species -----------------------------
trends <- read_csv(trends_csv, show_col_types = FALSE) %>%
  filter(model == model_tag, route_converged %in% TRUE) %>%
  left_join(traits %>% select(Code, invasive, redlistPopulationTrend), by = c(species_code = "Code"))

nn_trend_summary <- trends %>%
  filter(species_code %in% nn_codes) %>%
  group_by(species, species_code, invasive) %>%
  summarise(n_routes     = n(),
            median_trend = median(trend),
            mean_trend   = mean(trend),
            abund_wmean  = weighted.mean(trend, rel_abundance),
            pct_neg      = 100 * mean(trend < 0),
            pct_sig_dec  = 100 * mean(trend_uci < 0),   # 90% CrI entirely below 0
            pct_sig_inc  = 100 * mean(trend_lci > 0),   # 90% CrI entirely above 0
            iucn_trend   = first(redlistPopulationTrend),
            .groups = "drop") %>%
  arrange(median_trend)
cat("\n== 2. Route-level trends of non-native species (%/yr, converged routes) ==\n")
print(as.data.frame(nn_trend_summary), digits = 3)

cat("\nAll converged route x species rows, native vs non-native:\n")
print(trends %>%
        mutate(non_native = species_code %in% nn_codes) %>%
        group_by(non_native) %>%
        summarise(n = n(), median = median(trend), mean = mean(trend), pct_neg = 100 * mean(trend < 0)) %>%
        as.data.frame(), digits = 3)

# 3. Non-native species without trend estimates --------------------------------
missing_sp <- non_native %>% filter(!Code %in% trends$species_code)
cat("\n== 3. Non-native species with no trend estimates:", paste(missing_sp$Common.Name, collapse = ", "), "==\n")
if (file.exists(fit_log)) {
  log_lines <- readLines(fit_log, warn = FALSE)
  for (sp in missing_sp$Common.Name) {
    err <- grep(paste0("Data prep failed for ", sp, ":"), log_lines, value = TRUE, fixed = TRUE)
    cat(" ", sp, ":", if (length(err)) trimws(err[1]) else "no error line found in log", "\n")
  }
} else {
  cat("  fitting log not found:", fit_log, "\n")
}
cat("  (functions/covariate_model_fitting.R stops when a species has < 50 route-year observations,\n",
    "  before or after dropping rows without the Anthro covariate.)\n")

# 4. Non-native species in the tracking analysis --------------------------------
s1 <- read_csv(file.path(track_dir, paste0("sdm_tracking_stage1_species_", suffix, ".csv")), show_col_types = FALSE) %>%
  filter(status == "ok") %>%
  left_join(traits %>% distinct(Code, .keep_all = TRUE) %>% select(Code, invasive, redlistCategory, migration, diet),
            by = c(species_code = "Code"))
s1_nn <- s1 %>% filter(species_code %in% nn_codes)
cat("\n== 4. Non-native species in the tracking analysis ==\n")
cat(sprintf("  %d of %d species (%.1f%%); %d of %d route-by-species pairs (%.1f%%)\n",
            nrow(s1_nn), nrow(s1), 100 * nrow(s1_nn) / nrow(s1),
            sum(s1_nn$n_routes), sum(s1$n_routes), 100 * sum(s1_nn$n_routes) / sum(s1$n_routes)))
print(s1_nn %>% select(species, species_code, invasive, n_routes, slope_mean, slope_sd, slope_p_pos) %>%
        arrange(slope_mean) %>% as.data.frame(), digits = 3)

# 5. Category means (4c_sdm_tracking_meta.R) ---------------------------------------------
# All-species files are compared only when present (run 4a_sdm_tracking_routes.R /
# 4c_sdm_tracking_meta.R with native_only <- FALSE to produce them).
cat_nat <- read_if(file.path(track_dir, paste0("sdm_tracking_category_levels_", suffix, "_native.csv")))
cat_all <- read_if(file.path(track_dir, paste0("sdm_tracking_category_levels_", suffix, ".csv")))
# 4a labelled the main species set "all" before it was renamed "main"
if (!is.null(cat_all)) cat_all <- cat_all %>% mutate(species_set = recode(species_set, all = "main"))
if (!is.null(cat_nat)) cat_nat <- cat_nat %>% mutate(species_set = recode(species_set, all = "main"))
cat("\n== 5. Category means (%/yr), native only ==\n")
if (!is.null(cat_nat)) print(as.data.frame(cat_nat), digits = 3) else cat("  native category-level file not found\n")
if (!is.null(cat_all) && !is.null(cat_nat)) {
  cat_cmp <- cat_all %>% select(species_set, quantity, n_all = n_species, mu_all = mu) %>%
    inner_join(cat_nat %>% select(species_set, quantity, n_nat = n_species, mu_nat = mu),
               by = c("species_set", "quantity")) %>%
    mutate(change = mu_nat - mu_all)
  cat("\nAll species vs native only:\n")
  print(as.data.frame(cat_cmp), digits = 3)
} else cat("  (no all-species category-level file: comparison skipped)\n")

# 6. Tracking-slope models (4c_sdm_tracking_meta.R) ---------------------------
keep_cols <- c("model", "response", "n_species", "mu_mean", "mu_lo", "mu_hi", "p_mu_pos", "tau_mean",
               "pct_species_negative", "elpd_diff", "se_diff")
m_nat <- read_if(file.path(track_dir, paste0("sdm_tracking_stage2_models_", suffix, "_native.csv")))
m_all <- read_if(file.path(track_dir, paste0("sdm_tracking_stage2_models_", suffix, ".csv")))
cat("\n== 6. Stage-2 models, native only ==\n")
if (!is.null(m_nat)) print(as.data.frame(m_nat %>% select(any_of(keep_cols))), digits = 3) else cat("  native stage-2 file not found\n")
if (!is.null(m_all)) {
  cat("\nStage-2 models, all species:\n")
  print(as.data.frame(m_all %>% select(any_of(keep_cols))), digits = 3)
} else cat("  (no all-species stage-2 file: comparison skipped)\n")

# Trait models (4d_sdm_tracking_traits.R); the all-species run includes invasive status
t_nat <- read_if(file.path(track_dir, paste0("sdm_tracking_traits_models_", suffix, "_native.csv")))
t_all <- read_if(file.path(track_dir, paste0("sdm_tracking_traits_models_", suffix, ".csv")))
cat("\nTrait models, native only:\n")
if (!is.null(t_nat)) print(as.data.frame(t_nat %>% select(any_of(keep_cols))), digits = 3) else cat("  native trait-model file not found\n")
if (!is.null(t_all)) {
  cat("\nTrait models, all species (incl. invasive status):\n")
  print(as.data.frame(t_all %>% select(any_of(keep_cols))), digits = 3)
}

g_nat  <- read_if(file.path(track_dir, paste0("sdm_tracking_traits_groups_", suffix, "_native.csv")))
sp_nat <- read_if(file.path(track_dir, paste0("sdm_tracking_stage2_species_", suffix, "_native.csv")))
if (!is.null(g_nat) && !is.null(sp_nat)) {
  cat("\nHabitat-group slopes, native only (main.tex Table 2):\n")
  hab_tab <- g_nat %>%
    filter(model == "Habitat group", grouping == "g1") %>%
    select(level, n_species, mean, lo, hi, p_pos, class) %>%
    left_join(sp_nat %>% count(group, class_studentt) %>%
                pivot_wider(names_from = class_studentt, values_from = n, values_fill = 0),
              by = c(level = "group")) %>%
    arrange(desc(mean))
  print(as.data.frame(hab_tab), digits = 3)
  cat("\nSpecies classes (Student-t species-level model), native only:\n")
  print(table(factor(sp_nat$class_studentt, levels = c("Tracks", "Inconclusive", "Contrary"))))
}

# 7. SDM-matched routes, all vs native ---------------------------------------
# Same filters as 4a_sdm_tracking_routes.R: converged routes, classes 1-7,
# species with routes in at least min_classes classes.
sdm_files <- list.files(sdm_dir, pattern = paste0("_", model_tag, "_route_trends_sdm\\.csv$"), full.names = TRUE)
sdm_routes <- map_dfr(sdm_files, ~ read.csv(.x, stringsAsFactors = FALSE) %>%
                        transmute(species_code, route_converged, trend, cls = .data[[scenario]])) %>%
  filter(route_converged %in% TRUE)
route_counts <- function(d, label) {
  m <- d %>% filter(cls %in% 1:7)
  few <- m %>% group_by(species_code) %>% summarise(k = n_distinct(cls), .groups = "drop") %>%
    filter(k < min_classes) %>% pull(species_code)
  m2 <- m %>% filter(!species_code %in% few) %>%
    mutate(category = factor(case_when(cls <= 3 ~ "Decreasing (1-3)", cls == 4 ~ "Neutral (4)", TRUE ~ "Increasing (5-7)"),
                             levels = c("Decreasing (1-3)", "Neutral (4)", "Increasing (5-7)")))
  cat(sprintf("  %-7s converged %d | class 0 / no SDM %d | species with < %d classes: %d (%d routes) | kept %d routes from %d species\n",
              label, nrow(d), nrow(d) - nrow(m), min_classes, length(few), nrow(m) - nrow(m2),
              nrow(m2), n_distinct(m2$species_code)))
  print(m2 %>% group_by(category) %>%
          summarise(routes = n(), species = n_distinct(species_code), pooled_mean = round(mean(trend), 2)) %>%
          as.data.frame())
}
cat("\n== 7. SDM-matched route counts (", scenario, ") ==\n")
route_counts(sdm_routes, "all")
route_counts(sdm_routes %>% filter(!species_code %in% nn_codes), "native")

# 8. Trait counts, native species in the trait analysis -------------------------
nat_ok <- s1 %>% filter(!species_code %in% nn_codes)
cat("\n== 8. Trait counts for the", nrow(nat_ok), "native species in the trait analysis ==\n")
print(table(case_when(nat_ok$redlistCategory %in% c("Vulnerable", "Endangered", "Critically Endangered") ~ "Threatened",
                      TRUE ~ nat_ok$redlistCategory)))
print(table(nat_ok$migration))
print(sort(table(nat_ok$diet), decreasing = TRUE))
