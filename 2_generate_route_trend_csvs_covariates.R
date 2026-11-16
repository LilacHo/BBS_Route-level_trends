## =============================================================================
## Combine per-route trend output across ALL fitted models into one CSV, for
## cross-model comparison.
## Time period: 2010-2025
##
## Covers EVERY species in data/spp_names_codes_group_aou.csv (in_bbs == TRUE,
## all 12 groups), two models each -- "base" (no covariate) and "anthro" --
## as fit by 1a_species_iCAR_covariates.R. Each row carries its Group.
##
## NOTE: "base" is fit on the anthro-reduced dataset (only route-years where
## Anthro is non-NA), so it is directly comparable to "anthro".
##
## Reads, per species/model:
##   output/rds/<species>_iCAR_<tag>_<firstYear>_<lastYear>_summ_fit.rds
##   data/route_info/<species>_<tag>_<firstYear>_<lastYear>_route_info.rds
## and derives per-route trend, trend_lci, trend_uci, rel_abundance from
## beta/alpha.
##
## Writes:
##
##   1) Per-route file -- one row per route/species/model:
##      output/species_routes_covariates/all_route_trends_all_species_anthro_<firstYear>_<lastYear>.csv
##      columns: species, species_code, group, model, route, latitude,
##               longitude, alpha, beta, route_converged, trend, trend_lci,
##               trend_uci, rel_abundance
##
##      route_converged is TRUE only if this route's alpha_raw[r] AND
##      beta_raw_space[r] (the directly-sampled non-centered parameters) both
##      meet Rhat < 1.01 and bulk ESS > 400; FALSE if either fails or is NA.
##      It is independent of whole-model convergence
##      (helper/model_convergence.R), and it is this per-route flag that
##      drives which rows get excluded downstream. See the comment above
##      beta_raw_summ for why the transformed alpha[r]/beta[r] aren't used.
##
##   2) Model-level file -- one row per species/model:
##      output/species_routes_covariates/all_model_level_summary_all_species_anthro_<firstYear>_<lastYear>.csv
##      columns: species, species_code, group, model, n_routes, gamma1,
##               gamma1_lci, gamma1_uci, gamma1_excludes_zero
##
##      gamma1 is the covariate effect on log(lambda), one value per
##      species/model (NA for "base"). gamma1_lci/gamma1_uci are q5/q95 (90%
##      CI) on the log scale (no exp(), unlike trend_lci/trend_uci).
##      gamma1_excludes_zero = 90% CI does not span zero.
##
##   3) One per-route CSV per species/model (same columns as file 1), the
##      input format 3_add_SDM_covariates.R needs (one species_code per file):
##      output/species_routes_covariates/per_species/<species>_<model_tag>_route_trends.csv
##
## Post-processing only; safe to re-run at any point. Species/tags 1a hasn't
## fit yet are skipped and counted, not treated as errors.
## =============================================================================

library(tidyverse)
library(here)

here::i_am("2_generate_route_trend_csvs_covariates.R")

# Settings (match 1a_species_iCAR_covariates.R) ----------------------------
firstYear  <- 2010
lastYear   <- 2025
model_tags <- c("base", "anthro")

output_dir       <- here::here("output")
rds_dir          <- here::here("output", "rds")
route_info_dir   <- here::here("data", "route_info")
combined_out_dir <- here::here("output", "species_routes_covariates")
if (!dir.exists(combined_out_dir)) dir.create(combined_out_dir, recursive = TRUE)
per_species_out_dir <- here::here("output", "species_routes_covariates", "per_species")
if (!dir.exists(per_species_out_dir)) dir.create(per_species_out_dir, recursive = TRUE)

run_label <- "all_species_anthro"   # used only to name output files

# Full species list ---------------------------------------------------------
spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"),
                   stringsAsFactors = FALSE)

target_spp <- spp_df %>%
  filter(in_bbs == TRUE) %>%
  distinct(Common.Name, Code, Group, .keep_all = TRUE) %>%
  arrange(Common.Name)

cat("=== Combine route-trend CSVs across all species and both models ===\n")
cat("Period:", firstYear, "-", lastYear, "\n")
cat("Species (n =", nrow(target_spp), ") across", length(unique(target_spp$Group)), "groups\n")
cat("Models:", paste(model_tags, collapse = ", "), "\n")

# Helper: convert species name to file-safe format ------------------------
species_to_f <- function(sp) {
  gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
}

# Column order every row is coerced to before combining --------------------
route_cols <- c("species", "species_code", "group", "model", "route",
                "latitude", "longitude", "alpha", "beta", "route_converged",
                "trend", "trend_lci", "trend_uci", "rel_abundance")
model_cols <- c("species", "species_code", "group", "model", "n_routes",
                "gamma1", "gamma1_lci", "gamma1_uci", "gamma1_excludes_zero")

all_route_rows <- list()
all_model_rows <- list()
n_missing <- 0

for (i in seq_len(nrow(target_spp))) {
  sp       <- target_spp$Common.Name[i]
  sp_f     <- species_to_f(sp)
  sp_code  <- target_spp$Code[i]
  sp_group <- target_spp$Group[i]

  if (i %% 25 == 0 || i == 1) cat("\n[", i, "/", nrow(target_spp), "]", sp, "(", sp_group, ")\n")

  for (tag in model_tags) {

    out_base        <- paste0(sp_f, "_iCAR_", tag, "_", firstYear, "_", lastYear)
    summ_file       <- file.path(rds_dir, paste0(out_base, "_summ_fit.rds"))
    route_info_file <- file.path(route_info_dir,
                                 paste0(sp_f, "_", tag, "_", firstYear, "_", lastYear, "_route_info.rds"))

    if (!file.exists(summ_file) || !file.exists(route_info_file)) {
      n_missing <- n_missing + 1
      next   # not yet fit by 1a; counted in the summary at the end
    }

    summ       <- readRDS(summ_file)
    route_info <- readRDS(route_info_file)   # route, routeF, latitude, longitude

    # beta (slope) per route -> annual % trend + 90% CI ----------------------
    beta_summ <- summ %>%
      filter(str_detect(variable, "^beta\\[")) %>%
      transmute(routeF    = as.integer(str_extract(variable, "\\d+")),
                beta      = mean,
                trend     = 100 * (exp(mean) - 1),
                trend_lci = 100 * (exp(q5)   - 1),
                trend_uci = 100 * (exp(q95)  - 1))

    # alpha (intercept) per route -> relative abundance ---------------------
    alpha_summ <- summ %>%
      filter(str_detect(variable, "^alpha\\[")) %>%
      transmute(routeF        = as.integer(str_extract(variable, "\\d+")),
                alpha         = mean,
                rel_abundance = exp(mean))

    # route_converged uses the directly-sampled beta_raw_space[r]/alpha_raw[r],
    # NOT the transformed beta[r]/alpha[r]. beta[r] = sdbeta_space *
    # beta_raw_space[r] + BETA: when sdbeta_space is small, a poorly mixing
    # beta_raw_space[r] barely moves beta[r], so beta[r]'s rhat/ess can look
    # fine and mask the ICAR/funnel non-convergence this flag exists to catch.
    beta_raw_summ <- summ %>%
      filter(str_detect(variable, "^beta_raw_space\\[")) %>%
      transmute(routeF        = as.integer(str_extract(variable, "\\d+")),
                rhat_beta     = rhat,
                ess_bulk_beta = ess_bulk)

    alpha_raw_summ <- summ %>%
      filter(str_detect(variable, "^alpha_raw\\[")) %>%
      transmute(routeF         = as.integer(str_extract(variable, "\\d+")),
                rhat_alpha     = rhat,
                ess_bulk_alpha = ess_bulk)

    # Species/model-level covariate effect with 90% CI; NA for "base".
    if ("gamma1" %in% summ$variable) {
      gamma1_row        <- summ[summ$variable == "gamma1", ]
      gamma1_val        <- gamma1_row$mean[1]
      gamma1_lci_val    <- gamma1_row$q5[1]
      gamma1_uci_val    <- gamma1_row$q95[1]
      gamma1_excludes_zero_val <- !is.na(gamma1_lci_val) && !is.na(gamma1_uci_val) &&
        (gamma1_lci_val > 0 || gamma1_uci_val < 0)
    } else {
      gamma1_val               <- NA_real_
      gamma1_lci_val            <- NA_real_
      gamma1_uci_val            <- NA_real_
      gamma1_excludes_zero_val <- NA
    }

    route_trends <- beta_summ %>%
      left_join(alpha_summ, by = "routeF") %>%
      left_join(beta_raw_summ, by = "routeF") %>%
      left_join(alpha_raw_summ, by = "routeF") %>%
      left_join(route_info, by = "routeF") %>%
      mutate(species      = sp,
             species_code = sp_code,
             group        = sp_group,
             model        = tag,
             # NA rhat/ess counts as NOT converged
             route_converged = !is.na(rhat_alpha) & !is.na(ess_bulk_alpha) &
               !is.na(rhat_beta) & !is.na(ess_bulk_beta) &
               rhat_alpha < 1.01 & ess_bulk_alpha > 400 &
               rhat_beta  < 1.01 & ess_bulk_beta  > 400) %>%
      arrange(routeF) %>%
      select(all_of(route_cols))

    all_route_rows[[paste(sp, tag, sep = " | ")]] <- route_trends

    all_model_rows[[paste(sp, tag, sep = " | ")]] <- data.frame(
      species              = sp,
      species_code         = sp_code,
      group                = sp_group,
      model                = tag,
      n_routes             = nrow(route_trends),
      gamma1               = gamma1_val,
      gamma1_lci           = gamma1_lci_val,
      gamma1_uci           = gamma1_uci_val,
      gamma1_excludes_zero = gamma1_excludes_zero_val
    )

    # Per-species-per-model CSV (input for 3_add_SDM_covariates.R)
    per_species_csv <- file.path(per_species_out_dir,
                                 paste0(sp_f, "_", tag, "_route_trends.csv"))
    write.csv(route_trends, per_species_csv, row.names = FALSE)
  }
}

if (length(all_route_rows) == 0) {
  stop("Nothing to combine — run 1a_species_iCAR_covariates.R first.")
}

all_route_trends <- bind_rows(all_route_rows)
all_model_summary <- bind_rows(all_model_rows)

route_csv <- file.path(combined_out_dir,
                       paste0("all_route_trends_", run_label, "_",
                              firstYear, "_", lastYear, ".csv"))
model_csv <- file.path(combined_out_dir,
                       paste0("all_model_level_summary_", run_label, "_",
                              firstYear, "_", lastYear, ".csv"))

write.csv(all_route_trends, route_csv, row.names = FALSE)
write.csv(all_model_summary, model_csv, row.names = FALSE)

cat("\n=== Summary ===\n")
cat("Per-route rows:", nrow(all_route_trends), "\n")
cat("Species x model combinations written:", length(all_route_rows), "\n")
cat("Species x model combinations not yet fit by 1a (skipped):", n_missing, "\n")
print(table(all_route_trends$model))
cat("\nBy group:\n")
print(all_model_summary %>% distinct(species, group) %>% count(group, sort = TRUE))
cat("\nPer-route CSV written to:", route_csv, "\n")
cat("Model-level CSV written to:", model_csv, "\n")
cat("Per-species-per-model CSVs written to:", per_species_out_dir,
    "(", length(all_route_rows), "files )\n")

cat("\nSpecies with anthro's gamma1 90% CI excluding zero (credible effect):\n")
print(all_model_summary %>%
        filter(model == "anthro") %>%
        summarise(n_species = n(), n_excludes_zero = sum(gamma1_excludes_zero, na.rm = TRUE)))

cat("\nSame, broken down by group:\n")
print(all_model_summary %>%
        filter(model == "anthro") %>%
        group_by(group) %>%
        summarise(n_species = n(), n_excludes_zero = sum(gamma1_excludes_zero, na.rm = TRUE),
                  .groups = "drop") %>%
        arrange(desc(n_species)))

cat("\nNext: run 3_add_SDM_covariates.R to add climate-suitability columns.\n")
