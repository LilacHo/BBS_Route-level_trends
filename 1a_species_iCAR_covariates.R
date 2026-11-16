## =============================================================================
## Full production run: iCAR route-level slope model, base vs. anthro
## Time period: 2010-2025
##
## Fits TWO models for EVERY BBS species in this project's species list
## (data/spp_names_codes_group_aou.csv, in_bbs == TRUE -- all groups:
## aridlands, boreal_forests, eastern_forests, western_forests,
## subtropical_forests, grasslands, marshlands, coastal, waterbirds,
## generalists, arctic, urban_suburban -- currently 607 species):
##
##   base   -- no covariate    (models/slope_iCAR_route_NB_New.stan)
##   anthro -- data/Anthro.csv  (models/slope_iCAR_route_NB_New_covariate.stan)
##     log(lambda[r,t]) = alpha[r] + beta[r]*t + gamma1*Anthro[r,t]
##
## Why only anthro: in pilot runs the "*_to_anthro" transition covariates
## never produced a credible gamma1 (too sparse at route/year resolution),
## group-specific habitat covariates aren't available consistently across
## all 12 groups, and anthro was the one covariate with a consistently
## credible, well-powered effect.
##
## Both models are fit on the SAME reduced dataset per species: counts with
## a missing (NA) Anthro value are dropped once up front, so the comparison
## is fair.
##
## Covariate data: data/Anthro.csv, one row per BBS route per year. Not
## rescaled (standing-stock proportion already spans most of [0,1]).
##
## strat = "bcr": bbsBayes2::prepare_data() only formats new_data$route as
## "<StateNum>-<Route>" (the key the covariate join needs) when stratified
## by "bcr".
##
## Species are sorted alphabetically for a stable, resumable run order.
## force_refit = FALSE skips any species/tag that already has a saved
## summ_fit.rds, so an interrupted run resumes where it left off.
##
## This is a LONG run: ~600 species x 2 models x (BBS data pull + Voronoi
## neighbours + 4-chain/4000-iteration Stan fit). Expect hours to days.
##
## This script fits the models and saves, per species and per model tag:
##   output/rds/<species>_iCAR_<tag>_<firstYear>_<lastYear>_stanfit.rds
##   output/rds/<species>_iCAR_<tag>_<firstYear>_<lastYear>_summ_fit.rds
##   data/route_info/<species>_<tag>_<firstYear>_<lastYear>_route_info.rds
##   data/stan_data/<species>_<tag>_<firstYear>_<lastYear>_stan_data.RData
## Plus one combined diagnostics CSV and (via helper/model_convergence.R)
## one whole-model convergence table across all species and tags, including
## "base" (Rhat < 1.01 and bulk ESS > 400 required across every parameter).
## For gamma1's posterior summary and credibility (90% CI excludes zero)
## across species, run 2a_anthro_cover_gamma_supplement.R.
##
## Requires (sourced below):
##   functions/covariate_model_fitting.R
##     species_to_f()            file-safe species name for output filenames
##     load_covariate()          reads data/Anthro.csv into a route-year lookup
##     prepare_species_data()    BBS data pull, Anthro join, NA-drop, and
##                               spatial neighbours for one species
##     fit_one_covariate_model() fits "base" or "anthro" for one species and
##                               saves the output files listed above
##     That file itself sources:
##       functions/neighbours_define_voronoi.R
##         neighbours_define_voronoi()  Voronoi route neighbourhood for the
##                                      iCAR prior (used by prepare_species_data())
##       functions/posterior_summary_functions.R  (loaded, not currently used)
##   helper/model_convergence.R
##     model_convergence_table() whole-model convergence table, run at the end
## =============================================================================

library(bbsBayes2)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(sf)
library(spdep)
library(concaveman)
library(here)

here::i_am("1a_species_iCAR_covariates.R")

# species_to_f(), load_covariate(), prepare_species_data(),
# fit_one_covariate_model() -- shared with 1b_refit_nonconverged_species.R.
source(here::here("functions", "covariate_model_fitting.R"))

# Each species/model_tag gets its own subfolder here (created/deleted inside
# fit_one_covariate_model()) so per-chain CSVs don't pile up in %TEMP%.
cmdstanr_output_dir <- file.path(tempdir(), "cmdstan_output")
if (!dir.exists(cmdstanr_output_dir)) dir.create(cmdstanr_output_dir, recursive = TRUE)

# Settings ----------------------------------------------------------------
firstYear  <- 2010
lastYear   <- 2025

strat <- "bcr"   # NOT "bbs_usgs" -- see header

force_refit <- FALSE   # TRUE = re-fit everything, ignoring saved fits

model_tags <- c("base", "anthro")

# Output directories --------------------------------------------------------
output_dir <- here::here("output")
if (!dir.exists(output_dir)) dir.create(output_dir)
rds_dir <- here::here("output", "rds")
if (!dir.exists(rds_dir)) dir.create(rds_dir, recursive = TRUE)
if (!dir.exists(here::here("data"))) dir.create(here::here("data"))
if (!dir.exists(here::here("data", "maps"))) dir.create(here::here("data", "maps"), recursive = TRUE)
if (!dir.exists(here::here("data", "stan_data"))) dir.create(here::here("data", "stan_data"), recursive = TRUE)
route_info_dir <- here::here("data", "route_info")
if (!dir.exists(route_info_dir)) dir.create(route_info_dir, recursive = TRUE)

# Thin wrapper around load_covariate(): optional rescale by observed max,
# plus a summary printout. out_name becomes both the joined column name in
# new_data AND the model_tag used to select it in fit_one_covariate_model().
load_group_covariate <- function(file, value_col, out_name, rescale) {
  lk <- load_covariate(file, value_col, out_name)
  if (rescale) {
    cscale <- max(lk[[out_name]], na.rm = TRUE)
    lk[[out_name]] <- lk[[out_name]] / cscale
    cat("  ", out_name, "covariate rescaled by its observed max (", round(cscale, 6), ") to span [0,1]\n")
  }
  cat("  ", out_name, "covariate:", nrow(lk), "route-year rows (",
      length(unique(lk$route_key)), "routes)\n")
  lk
}

# Covariate spec -- the list name is both the model_tag and the joined
# column name in new_data.
covariate_specs <- list(
  anthro = list(file = "Anthro.csv", value_col = "Anthro", rescale = FALSE)
)

# Species list ------------------------------------------------------------
spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"),
                   stringsAsFactors = FALSE)

target_spp <- spp_df %>%
  filter(in_bbs == TRUE) %>%
  distinct(Common.Name, Code, .keep_all = TRUE) %>%
  arrange(Common.Name)

cat("=== Full run: ", nrow(target_spp), " species x ", length(model_tags),
    " models (", paste(model_tags, collapse = ", "), ") ===\n", sep = "")
cat("Groups represented:", paste(sort(unique(target_spp$Group)), collapse = ", "), "\n")

# Compile the Stan models once (reused across every species) ---------------
model_base   <- cmdstan_model("models/slope_iCAR_route_NB_New.stan",
                              stanc_options = list("O1"))
model_single <- cmdstan_model("models/slope_iCAR_route_NB_New_covariate.stan",
                              stanc_options = list("O1"))

# Load covariates once, up front, for every species ------------------------
cat("Loading covariates...\n")
covariate_lookups <- list()
for (cn in setdiff(model_tags, "base")) {
  spec <- covariate_specs[[cn]]
  covariate_lookups[[cn]] <- load_group_covariate(spec$file, spec$value_col, cn, spec$rescale)
}

# ==========================================================================
# Main loop: every species, both model tags.
# ==========================================================================
# model_convergence_table() is called after the loop; skip its autorun on
# source() so it uses the settings above.
model_convergence_skip_autorun <- TRUE
source(here::here("helper", "model_convergence.R"))

results_list <- list()
diagnostics_list <- list()

for (i in seq_len(nrow(target_spp))) {
  sp       <- target_spp$Common.Name[i]
  sp_f     <- species_to_f(sp)
  sp_code  <- target_spp$Code[i]
  sp_bbs   <- target_spp$bbs_english[i]
  sp_group <- target_spp$Group[i]

  cat("\n================================================================\n")
  cat("  [", i, "/", nrow(target_spp), "]", sp, " (", sp_group, ")\n")
  cat("================================================================\n")

  all_exist <- all(vapply(model_tags, function(tag) {
    out_base  <- paste0(sp_f, "_iCAR_", tag, "_", firstYear, "_", lastYear)
    summ_file <- file.path(rds_dir, paste0(out_base, "_summ_fit.rds"))
    file.exists(summ_file)
  }, logical(1)))

  if (all_exist && !force_refit) {
    cat("  Skipping (both", paste(model_tags, collapse = " and "), "already fitted for", sp, ")\n")
    next
  }

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

  for (tag in model_tags) {
    out_base       <- paste0(sp_f, "_iCAR_", tag, "_", firstYear, "_", lastYear)
    summ_file      <- file.path(rds_dir, paste0(out_base, "_summ_fit.rds"))
    stan_data_file <- here::here("data", "stan_data",
                                 paste0(sp_f, "_", tag, "_", firstYear, "_", lastYear, "_stan_data.RData"))
    if (file.exists(summ_file) && file.exists(stan_data_file) && !force_refit) {
      cat("  Skipping (already fitted):", basename(summ_file), "\n")
      next
    }

    diagnostics <- tryCatch(
      fit_one_covariate_model(species = sp, species_f = sp_f, model_tag = tag,
                              prepped = prepped, firstYear = firstYear, lastYear = lastYear,
                              model_base = model_base, model_single = model_single,
                              rds_dir = rds_dir, route_info_dir = route_info_dir,
                              cmdstanr_output_dir = cmdstanr_output_dir,
                              chains = 4, iter_warmup = 2000, iter_sampling = 2000,
                              adapt_delta = 0.8, max_treedepth = 10,
                              show_exceptions = FALSE),
      error = function(e) {
        message("  [ERROR] Model '", tag, "' failed for ", sp, ": ", conditionMessage(e))
        return(NULL)
      }
    )
    if (is.null(diagnostics)) next

    diagnostics$group <- sp_group
    results_list[[paste(sp, tag, sep = " | ")]] <- summ_file
    diagnostics_list[[paste(sp, tag, sep = " | ")]] <- diagnostics

    cat("  Done —", tag, "fit + diagnostics saved for", sp, "\n")
  }

  # Checkpoint the diagnostics CSV every 25 species so an interruption
  # doesn't lose all progress info.
  if (length(diagnostics_list) > 0 && i %% 25 == 0) {
    diag_csv <- file.path(output_dir,
                          paste0("diagnostics_covariates_all_species_anthro_",
                                 firstYear, "_", lastYear, ".csv"))
    write.csv(bind_rows(diagnostics_list), diag_csv, row.names = FALSE)
    cat("  [checkpoint] diagnostics written to:", diag_csv, "(", i, "/", nrow(target_spp), "species processed )\n")
  }
}

# Final combined diagnostics CSV ---------------------------------------------
if (length(diagnostics_list) > 0) {
  diagnostics_all <- bind_rows(diagnostics_list)
  diag_csv <- file.path(output_dir,
                        paste0("diagnostics_covariates_all_species_anthro_",
                               firstYear, "_", lastYear, ".csv"))
  write.csv(diagnostics_all, diag_csv, row.names = FALSE)
  cat("\nDiagnostics written to:", diag_csv, "\n")
}

# Whole-model convergence table across every species AND every tag,
# including "base" -- Rhat < 1.01 and bulk ESS > 400 required across every
# parameter, not just gamma1 (see helper/model_convergence.R).
model_convergence_table(target_spp = target_spp, model_tags = model_tags,
                        firstYear = firstYear, lastYear = lastYear,
                        rds_dir = rds_dir, bird_group = "all_species_anthro")

cat("\n\n=== DONE ===\n")
cat("Species x model fits this run:", length(results_list), "\n")
cat("Fits saved in:", output_dir, "\n")
