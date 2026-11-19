## =============================================================================
## 1b_refit_nonconverged_species.R
##
## Targeted re-fit of ONLY the species/tag combos that fail whole-model
## convergence (Rhat < 1.01 and bulk ESS > 400 across EVERY parameter, every
## model tag including "base"; see helper/model_convergence.R).
##
##   Step 1. Find failing combos by reading every *_summ_fit.rds off disk
##           (not the diagnostics CSV, which only covers one run session).
##   Step 2. Refit those with more iterations and a higher adapt_delta, on the
##           same anthro-non-NA dataset 1a uses. The current fit is first
##           COPIED to *_under-converged/ backup folders, then overwritten at
##           the production filenames, so step 2 picks up the refit unchanged.
##   Step 3. Re-check convergence.
##   Step 4. Extra push (8000/8000, adapt_delta = 0.95) for a named set of
##           known "stubborn" species/tags (stubborn_species_tags below)
##           still failing after Step 3.
##   Step 5. (OPTIONAL, commented out) route-level drill-down for whatever
##           still fails.
##
## For gamma1's posterior summary on the post-refit fits, run
## s1_gamma.R.
##
## Resumable: each species/tag is re-checked right before refitting and
## skipped if an earlier interrupted run already fixed it.
##
## Note: a chain crashing outright (nonzero return code, e.g. a "differing
## number of rows: 4, 3" error) is a different failure mode from slow mixing,
## and more iterations alone won't fix it. fit_one_covariate_model() reports
## which chain failed.
##
## IMPORTANT: afterwards, re-run 2_generate_route_trend_csvs_covariates.R so
## the refits propagate into all_route_trends_*.csv /
## all_model_level_summary_*.csv.
## =============================================================================

library(bbsBayes2)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(sf)
library(spdep)
library(concaveman)
library(here)

here::i_am("1b_refit_nonconverged_species.R")

# species_to_f(), load_covariate(), prepare_species_data(),
# fit_one_covariate_model() -- shared with 1a_species_iCAR_covariates.R.
source(here::here("functions", "covariate_model_fitting.R"))

# Source without triggering its auto-run block; called below with this
# script's settings ------------------------------------------------------------
model_convergence_skip_autorun <- TRUE
source(here::here("helper", "model_convergence.R"))

cmdstanr_output_dir <- file.path(tempdir(), "cmdstan_output_refit")
if (!dir.exists(cmdstanr_output_dir)) dir.create(cmdstanr_output_dir, recursive = TRUE)

# Settings -- match 1a_species_iCAR_covariates.R -----------------------------
firstYear  <- 2010
lastYear   <- 2025
strat      <- "bcr"
model_tags <- c("base", "anthro")

# Refit settings (1a uses 2000/2000, adapt_delta 0.8, show_exceptions FALSE)
refit_iter_warmup    <- 4000
refit_iter_sampling  <- 4000
refit_adapt_delta    <- 0.9
refit_max_treedepth  <- 10     # raise if treedepth warnings show up
refit_show_exceptions <- TRUE  # surface per-chain crash text for these known-problem fits

# Step 4 extra-push settings, for species/tags that stay non-convergent after
# the standard refit despite being well sampled. All 12 combos converged with
# these settings in a 2026-09-11 test run.
stubborn_iter_warmup   <- 8000
stubborn_iter_sampling <- 8000
stubborn_adapt_delta   <- 0.95
stubborn_species_tags <- list(
  "Sharp-shinned Hawk" = c("base", "anthro"),
  "Hooded Merganser"   = c("base", "anthro"),
  "Broad-winged Hawk"  = c("base", "anthro"),
  "Belted Kingfisher"  = c("base", "anthro"),
  "Western Meadowlark" = c("anthro"),   # base already converges
  "Gray Catbird"       = c("anthro"),
  "Red Crossbill"      = c("anthro"),
  "Wood Duck"          = c("anthro")
)

# Output directories -- same as 1a; fits are overwritten in place ------------
output_dir     <- here::here("output")
rds_dir        <- here::here("output", "rds")
route_info_dir <- here::here("data", "route_info")
stan_data_dir  <- here::here("data", "stan_data")
if (!dir.exists(stan_data_dir)) dir.create(stan_data_dir, recursive = TRUE)
if (!dir.exists(here::here("data", "maps"))) dir.create(here::here("data", "maps"), recursive = TRUE)

# Backup directories for the pre-refit (under-converged) fits ---------------
rds_backup_dir        <- here::here("output", "rds_under-converged")
route_info_backup_dir <- here::here("data", "route_info_under-converged")
stan_data_backup_dir  <- here::here("data", "stan_data_under-converged")
if (!dir.exists(rds_backup_dir))        dir.create(rds_backup_dir, recursive = TRUE)
if (!dir.exists(route_info_backup_dir)) dir.create(route_info_backup_dir, recursive = TRUE)
if (!dir.exists(stan_data_backup_dir))  dir.create(stan_data_backup_dir, recursive = TRUE)

#' Archive the CURRENT (about-to-be-overwritten) fit files for one
#' species/tag into the *_under-converged backup directories, before
#' fit_one_covariate_model() overwrites them with a refit.
#'
#' COPIES rather than moves, so a refit that errors partway leaves the
#' original fit in place. Keeps only the most recent pre-refit snapshot.
archive_existing_fit <- function(species_f, model_tag, firstYear, lastYear) {
  out_base <- paste0(species_f, "_iCAR_", model_tag, "_", firstYear, "_", lastYear)

  copy_if_exists <- function(from, to_dir) {
    if (file.exists(from)) {
      to <- file.path(to_dir, basename(from))
      file.copy(from, to, overwrite = TRUE)
      cat("    Backed up pre-refit file:", basename(from), "->", to_dir, "\n")
    }
  }

  copy_if_exists(file.path(rds_dir, paste0(out_base, "_stanfit.rds")), rds_backup_dir)
  copy_if_exists(file.path(rds_dir, paste0(out_base, "_summ_fit.rds")), rds_backup_dir)
  copy_if_exists(file.path(route_info_dir,
                           paste0(species_f, "_", model_tag, "_", firstYear, "_", lastYear, "_route_info.rds")),
                 route_info_backup_dir)
  copy_if_exists(file.path(stan_data_dir,
                           paste0(species_f, "_", model_tag, "_", firstYear, "_", lastYear, "_stan_data.RData")),
                 stan_data_backup_dir)
}

# Covariate spec -- same as 1a ----------------------------------------------
covariate_specs <- list(
  anthro = list(file = "Anthro.csv", value_col = "Anthro", rescale = FALSE)
)

# Loaded even when only "base" needs refitting: base is fit on the same
# anthro-non-NA dataset.
cat("Loading covariates...\n")
covariate_lookups <- list(
  anthro = load_covariate(covariate_specs$anthro$file, covariate_specs$anthro$value_col, "anthro")
)
cat("  anthro covariate:", nrow(covariate_lookups$anthro), "route-year rows (",
    length(unique(covariate_lookups$anthro$route_key)), "routes)\n")

# Species list -- same as 1a -------------------------------------------------
spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"),
                   stringsAsFactors = FALSE)
target_spp <- spp_df %>%
  filter(in_bbs == TRUE) %>%
  distinct(Common.Name, Code, .keep_all = TRUE) %>%
  arrange(Common.Name)

# ==========================================================================
# Step 1: find who currently fails whole-model convergence ------------------
# ==========================================================================
cat("\n=== Checking current convergence status across all species/tags ===\n")
convergence_now <- model_convergence_table(target_spp = target_spp, model_tags = model_tags,
                                           firstYear = firstYear, lastYear = lastYear,
                                           rds_dir = rds_dir, bird_group = "all_species_anthro",
                                           write_csv = FALSE)

if (is.null(convergence_now)) {
  stop("No existing fits found at all -- run 1a_species_iCAR_covariates.R first.")
}

failing <- convergence_now %>% filter(!model_converged)

cat("\n", nrow(failing), "species/tag combo(s) currently fail the convergence criterion ",
    "(Rhat < 1.01 and bulk ESS > 400) and will be refit with increased iterations.\n", sep = "")

if (nrow(failing) == 0) {
  cat("Nothing to refit -- every existing fit already meets the criterion. Exiting.\n")
  quit(save = "no", status = 0)
}

# Reset na.print temporarily: an invalid session value makes
# print.data.frame() error with "invalid 'na.print' specification".
old_na_print <- getOption("na.print")
options(na.print = "NA")
print(as.data.frame(failing %>% select(species, species_code, model, model_max_rhat, model_min_ess_bulk)))
options(na.print = old_na_print)

# Compile only the Stan model(s) actually needed for the failing tags -------
tags_needed  <- unique(failing$model)
model_base   <- if ("base" %in% tags_needed) {
  cmdstan_model(here::here("models", "slope_iCAR_route_NB_New.stan"), stanc_options = list("O1"))
} else NULL
model_single <- if (any(tags_needed != "base")) {
  cmdstan_model(here::here("models", "slope_iCAR_route_NB_New_covariate.stan"), stanc_options = list("O1"))
} else NULL

# ==========================================================================
# Step 2: refit each failing species/tag -------------------------------------
# ==========================================================================
results_list <- list()
refit_diagnostics_list <- list()

# One prepare_species_data() call (BBS pull + Voronoi) per species, even if
# both tags failed ------------------------------------------------------------
failing_species <- failing %>% distinct(species, species_code) %>% arrange(species)

for (i in seq_len(nrow(failing_species))) {
  sp      <- failing_species$species[i]
  sp_code <- failing_species$species_code[i]
  sp_f    <- species_to_f(sp)
  sp_row  <- target_spp %>% filter(Common.Name == sp, Code == sp_code)
  sp_bbs   <- sp_row$bbs_english[1]
  sp_group <- sp_row$Group[1]

  tags_for_this_species <- failing %>%
    filter(species == sp, species_code == sp_code) %>%
    pull(model)

  cat("\n================================================================\n")
  cat("  [", i, "/", nrow(failing_species), "]", sp, " (", sp_group, ") -- tags:",
      paste(tags_for_this_species, collapse = ", "), "\n")
  cat("================================================================\n")

  prepped <- tryCatch(
    prepare_species_data(species = sp, species_bbs = sp_bbs, strat = strat,
                         firstYear = firstYear, lastYear = lastYear,
                         covariate_lookups = covariate_lookups),
    error = function(e) {
      message("  [ERROR] Data prep failed for ", sp, ": ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(prepped)) next

  for (tag in tags_for_this_species) {

    # Skip if an earlier (interrupted) run of this script already fixed it
    current <- model_convergence_check(sp_f, tag, firstYear, lastYear, rds_dir)
    if (!is.null(current) && isTRUE(current$model_converged)) {
      cat("  Skipping", tag, "-- already converged (fixed by a previous run of this script).\n")
      next
    }

    # Preserve the current (under-converged) fit before it gets overwritten
    archive_existing_fit(sp_f, tag, firstYear, lastYear)

    diagnostics <- tryCatch(
      fit_one_covariate_model(species = sp, species_f = sp_f, model_tag = tag,
                              prepped = prepped, firstYear = firstYear, lastYear = lastYear,
                              model_base = model_base, model_single = model_single,
                              rds_dir = rds_dir, route_info_dir = route_info_dir,
                              cmdstanr_output_dir = cmdstanr_output_dir,
                              stan_data_dir = stan_data_dir,
                              chains = 4, iter_warmup = refit_iter_warmup,
                              iter_sampling = refit_iter_sampling,
                              adapt_delta = refit_adapt_delta,
                              max_treedepth = refit_max_treedepth,
                              show_exceptions = refit_show_exceptions,
                              overwrite_route_info = TRUE),
      error = function(e) {
        message("  [ERROR] Refit '", tag, "' failed for ", sp, ": ", conditionMessage(e))
        return(NULL)
      }
    )
    if (is.null(diagnostics)) next

    diagnostics$group <- sp_group
    if (!is.null(current)) {
      diagnostics$before_max_rhat     <- current$model_max_rhat
      diagnostics$before_min_ess_bulk <- current$model_min_ess_bulk
    }

    results_list[[paste(sp, tag, sep = " | ")]] <- TRUE
    refit_diagnostics_list[[paste(sp, tag, sep = " | ")]] <- diagnostics

    cat("  Done — refit", tag, "for", sp, "\n")
  }

  if (length(refit_diagnostics_list) > 0 && i %% 10 == 0) {
    diag_csv <- file.path(output_dir,
                          paste0("diagnostics_refit_nonconverged_", firstYear, "_", lastYear, ".csv"))
    write.csv(bind_rows(refit_diagnostics_list), diag_csv, row.names = FALSE)
    cat("  [checkpoint] refit diagnostics written to:", diag_csv,
        "(", i, "/", nrow(failing_species), "species processed )\n")
  }
}

# Final refit diagnostics CSV -------------------------------------------------
if (length(refit_diagnostics_list) > 0) {
  diag_csv <- file.path(output_dir,
                        paste0("diagnostics_refit_nonconverged_", firstYear, "_", lastYear, ".csv"))
  write.csv(bind_rows(refit_diagnostics_list), diag_csv, row.names = FALSE)
  cat("\nRefit diagnostics written to:", diag_csv, "\n")
}

# ==========================================================================
# Step 3: re-check convergence after refitting --------------------------------
# ==========================================================================
cat("\n=== Re-checking convergence after refit ===\n")
convergence_after <- model_convergence_table(target_spp = target_spp, model_tags = model_tags,
                                             firstYear = firstYear, lastYear = lastYear,
                                             rds_dir = rds_dir,
                                             bird_group = "all_species_anthro_postrefit")

still_failing <- convergence_after %>% filter(!model_converged)
cat("\n", nrow(still_failing), "species/tag combo(s) still fail the convergence criterion after refit.\n", sep = "")
if (nrow(still_failing) > 0) {
  old_na_print <- getOption("na.print")
  options(na.print = "NA")
  print(as.data.frame(still_failing %>% select(species, species_code, model, model_max_rhat, model_min_ess_bulk)))
  options(na.print = old_na_print)
  cat("\nThese may need an even longer refit, a different adapt_delta, or individual investigation",
      "-- e.g. a chain that's outright crashing is a different",
      "failure mode than slow mixing, and more iterations alone won't necessarily fix it.\n")
}

# ==========================================================================
# Step 4: extra-push refit (stubborn_iter_* / stubborn_adapt_delta) for the
# species/tags in stubborn_species_tags that still fail after Step 3.
# ==========================================================================
stubborn_lookup <- bind_rows(lapply(names(stubborn_species_tags), function(sp) {
  data.frame(species = sp, model = stubborn_species_tags[[sp]], stringsAsFactors = FALSE)
}))

stubborn_still_failing <- still_failing %>% semi_join(stubborn_lookup, by = c("species", "model"))

if (nrow(stubborn_still_failing) > 0) {
  stubborn_species <- stubborn_still_failing %>% distinct(species, species_code) %>% arrange(species)

  cat("\n=== Step 4: extra-push refit for stubborn species (still failing:",
      paste(stubborn_species$species, collapse = ", "), ") ===\n")
  cat("iter_warmup =", stubborn_iter_warmup, "| iter_sampling =", stubborn_iter_sampling,
      "| adapt_delta =", stubborn_adapt_delta, "\n")

  # Reuse already-compiled models; compile only what's missing.
  if (is.null(model_single) && any(stubborn_still_failing$model != "base")) {
    model_single <- cmdstan_model(here::here("models", "slope_iCAR_route_NB_New_covariate.stan"), stanc_options = list("O1"))
  }
  if (is.null(model_base) && any(stubborn_still_failing$model == "base")) {
    model_base <- cmdstan_model(here::here("models", "slope_iCAR_route_NB_New.stan"), stanc_options = list("O1"))
  }

  for (i in seq_len(nrow(stubborn_species))) {
    sp       <- stubborn_species$species[i]
    sp_code  <- stubborn_species$species_code[i]
    sp_f     <- species_to_f(sp)
    sp_row   <- target_spp %>% filter(Common.Name == sp, Code == sp_code)
    sp_bbs   <- sp_row$bbs_english[1]
    sp_group <- sp_row$Group[1]
    sp_fails <- stubborn_still_failing %>% filter(species == sp, species_code == sp_code)

    cat("\n================================================================\n")
    cat("  [Extra-push refit]", sp, "-- tags:", paste(sp_fails$model, collapse = ", "), "\n")
    cat("================================================================\n")

    prepped <- tryCatch(
      prepare_species_data(species = sp, species_bbs = sp_bbs, strat = strat,
                           firstYear = firstYear, lastYear = lastYear,
                           covariate_lookups = covariate_lookups),
      error = function(e) {
        message("  [ERROR] Data prep failed for ", sp, ": ", conditionMessage(e))
        return(NULL)
      }
    )
    if (is.null(prepped)) next

    for (j in seq_len(nrow(sp_fails))) {
      tag <- sp_fails$model[j]

      diagnostics <- tryCatch(
        fit_one_covariate_model(species = sp, species_f = sp_f, model_tag = tag,
                                prepped = prepped, firstYear = firstYear, lastYear = lastYear,
                                model_base = model_base, model_single = model_single,
                                rds_dir = rds_dir, route_info_dir = route_info_dir,
                                cmdstanr_output_dir = cmdstanr_output_dir,
                                stan_data_dir = stan_data_dir,
                                chains = 4, iter_warmup = stubborn_iter_warmup,
                                iter_sampling = stubborn_iter_sampling,
                                adapt_delta = stubborn_adapt_delta,
                                max_treedepth = refit_max_treedepth,
                                show_exceptions = refit_show_exceptions,
                                overwrite_route_info = TRUE),
        error = function(e) {
          message("  [ERROR] Extra-push refit '", tag, "' failed for ", sp, ": ", conditionMessage(e))
          return(NULL)
        }
      )
      if (is.null(diagnostics)) next

      diagnostics$group               <- sp_group
      diagnostics$before_max_rhat     <- sp_fails$model_max_rhat[j]
      diagnostics$before_min_ess_bulk <- sp_fails$model_min_ess_bulk[j]
      results_list[[paste(sp, tag, "extra-push", sep = " | ")]] <- TRUE
      refit_diagnostics_list[[paste(sp, tag, "extra-push", sep = " | ")]] <- diagnostics

      cat("  Done — extra-push refit", tag, "for", sp, "\n")
    }
  }

  if (length(refit_diagnostics_list) > 0) {
    diag_csv <- file.path(output_dir,
                          paste0("diagnostics_refit_nonconverged_", firstYear, "_", lastYear, ".csv"))
    write.csv(bind_rows(refit_diagnostics_list), diag_csv, row.names = FALSE)
    cat("\nRefit diagnostics (incl. extra push) written to:", diag_csv, "\n")
  }

  cat("\n=== Re-checking convergence after stubborn-species extra-push refit ===\n")
  convergence_after <- model_convergence_table(target_spp = target_spp, model_tags = model_tags,
                                               firstYear = firstYear, lastYear = lastYear,
                                               rds_dir = rds_dir,
                                               bird_group = "all_species_anthro_postrefit")
  still_failing <- convergence_after %>% filter(!model_converged)
  cat("\n", nrow(still_failing), "species/tag combo(s) still fail the convergence criterion after all refit steps.\n", sep = "")
  if (nrow(still_failing) > 0) {
    old_na_print <- getOption("na.print")
    options(na.print = "NA")
    print(as.data.frame(still_failing %>% select(species, species_code, model, model_max_rhat, model_min_ess_bulk)))
    options(na.print = old_na_print)
  }
} else {
  cat("\nNone of the stubborn species/tags (", paste(stubborn_lookup$species, "|", stubborn_lookup$model, collapse = ", "),
      ") are in the failing set -- skipping the extra-push step.\n", sep = "")
}

# ==========================================================================
# Step 5 (OPTIONAL): for whatever still fails after Steps 2-4, find WHICH
# routes are offending and why, so specific routes can be flagged rather
# than dropping the whole species.
#
# Chain (each needs the previous step's output CSV):
#   helper/diagnose_nonconvergence.R        -- classifies EVERY offending
#     parameter (route_alpha / route_beta / gamma1 / hyperparameter_sd /
#     other) per still-failing species/tag; NOT sourced here, run it first
#     if output/files/nonconvergence_diagnosis_<firstYear>_<lastYear>.csv
#     doesn't exist yet.
#   helper/translate_offending_routes.R     -- translates the route_alpha/
#     route_beta rows' internal routeF index into real BBS route IDs +
#     coordinates, via each species/tag's own route_info.rds.
#   helper/check_flagged_route_sparsity.R   -- for those flagged routes,
#     tests count-volatility and spatial neighbor-mismatch against each
#     species' other (converged) routes (Wilcoxon rank-sum).
#
# Uncomment to run the last two here:
# translate_offending_routes_skip_autorun <- TRUE
# source(here::here("helper", "translate_offending_routes.R"))
# translate_offending_routes(firstYear = firstYear, lastYear = lastYear)
#
# check_flagged_route_sparsity_skip_autorun <- TRUE
# source(here::here("helper", "check_flagged_route_sparsity.R"))
# check_flagged_route_sparsity_table(firstYear = firstYear, lastYear = lastYear)
# ==========================================================================

cat("\n\nIMPORTANT: re-run 2_generate_route_trend_csvs_covariates.R next so these improved fits",
    "propagate into all_route_trends_*.csv / all_model_level_summary_*.csv.\n")

cat("\n=== DONE ===\n")
cat("Species/tag combos refit this run:", length(results_list), "\n")
