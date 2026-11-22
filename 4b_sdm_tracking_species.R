# 4b_sdm_tracking_species.R
#
# Step 2 of the SDM-tracking analysis (Methods, Species tracking slope):
# species-level quantities from the posterior draws of the route trends.
# Run this where the Stan fits (*_stanfit.rds from 1a/1b) are stored; step 3
# (4c_sdm_tracking_meta.R) only needs this script's small CSV output.
#
# For each species, in EVERY one of n_draws posterior draws of the route-level
# slopes beta[r] (anthro iCAR model; a random subset, the same for every quantity):
#   1. each route's slope is converted to a trend, 100 * (exp(beta) - 1) %/yr;
#   2. TRACKING SLOPE: route trend is regressed on the route's SDM class score by
#      ordinary least squares (routes weighted equally), in %/yr per class step.
#      Positive = trends are better where the SDM projects improving suitability;
#   3. the same slope on the worsening half (classes 1-4) and the improving half
#      (classes 4-7) of the gradient, for species with >= 2 classes on that half;
#   4. CATEGORY MEANS: the mean trend of the species' routes in each range-change
#      category (Decreasing suitability = 1-3, Neutral = 4, Increasing suitability = 5-7),
#      for the categories in which the species has routes.
# Each quantity is then summarized across draws by its mean (y) and SD (se), the
# inputs to the meta-analysis in 4c. Because the draws are joint across routes,
# se carries both each route's trend uncertainty and the spatial (iCAR)
# correlation between neighbouring routes.
#
# Class score (class_score): ordinal steps centred on Neutral,
#   1 extirpation -3 | 2 worsening -2 | 3 slightly worsening -1 | 4 neutral 0 |
#   5 slightly improving +1 | 6 improving +2 | 7 colonization +3.
# The percent-change bands behind the classes are open-ended at the top
# (improving = 50%+, colonization = newly suitable), so equal steps are used.
#
# Route selection as in 4a: converged routes in classes 1-7; species need routes
# in at least min_classes classes. All species are processed, native or not;
# 4c drops the non-native species (native_only), so the same file also serves
# the all-species sensitivity analysis (s2_non_native.R).
#
# Reads:  <rds_dir>/<species>_iCAR_<tag>_<firstYear>_<lastYear>_stanfit.rds   (cmdstanr save_object)
#         <route_info_dir>/<species>_<tag>_<firstYear>_<lastYear>_route_info.rds (route <-> routeF)
#         output/species_routes_covariates/per_species_sdm/<species>_<tag>_route_trends_sdm.csv  (step 3)
# Writes: output/sdm_tracking/sdm_tracking_stage1_species_<suffix>.csv
#           one row per species: status, route counts per class and category, and for
#           each quantity its posterior mean, sd, 2.5% / 97.5% quantiles and P(> 0):
#           slope_*, slope_low_*, slope_up_*, mean_decreasing_*, mean_neutral_*, mean_increasing_*
#         output/sdm_tracking/sdm_tracking_stage1_draws_<suffix>.rds   (optional, save_draws)
#         output/sdm_tracking/sdm_tracking_stage1_log_<suffix>.txt
#
# Resumable: species already finished (status "ok") in an existing output CSV
# are skipped, so the script can be stopped and restarted. An output CSV from
# before the category means were added here is ignored (all species redone).
#
# Needs: tidyverse, posterior, cmdstanr (to read the saved CmdStanMCMC objects).

library(here)
library(tidyverse)
library(posterior)
library(cmdstanr)

here::i_am("4b_sdm_tracking_species.R")

# Settings -------------------------------------------------------------------
model_tag  <- "anthro"
scenario   <- "rcp45"
firstYear  <- 2010
lastYear   <- 2025
rds_dir        <- here::here("output", "rds")          # matches 1a_species_iCAR_covariates.R (*_stanfit.rds)
route_info_dir <- here::here("data", "route_info")     # matches 1a_species_iCAR_covariates.R (*_route_info.rds)
sdm_dir        <- here::here("output", "species_routes_covariates", "per_species_sdm")
out_dir        <- here::here("output", "sdm_tracking")

n_draws    <- 2000      # posterior draws used per species (random subset of the 8000 or more)
save_draws <- TRUE      # also save a thinned set of per-species slope draws (n_draws_saved)
n_draws_saved <- 500
seed       <- 20261014
min_classes <- 2        # species need routes in at least this many SDM classes (as in 4a)
class_score <- c(`1` = -3, `2` = -2, `3` = -1, `4` = 0, `5` = 1, `6` = 2, `7` = 3)
categories  <- list(decreasing = 1:3, neutral = 4L, increasing = 5:7)
only_species <- NULL    # e.g. c("Song Sparrow", "Wild Turkey") for a quick test run
# SDM classes to drop (sensitivity analysis). integer(0) = main analysis (all classes 1-7);
# c(7) = drop colonization, adding "_no7" to every output file name. Must match 4a and 4c.
exclude_classes <- integer(0)

run_tag  <- if (length(exclude_classes) > 0) paste0("_no", paste(sort(exclude_classes), collapse = "")) else ""
suffix   <- paste0(paste(model_tag, scenario, firstYear, lastYear, sep = "_"), run_tag)
out_csv  <- file.path(out_dir, paste0("sdm_tracking_stage1_species_", suffix, ".csv"))
out_rds  <- file.path(out_dir, paste0("sdm_tracking_stage1_draws_", suffix, ".rds"))
log_file <- file.path(out_dir, paste0("sdm_tracking_stage1_log_", suffix, ".txt"))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

species_to_f <- function(sp) gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)

log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), "  ", paste0(...))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

# OLS slope of each draw's route trends on x: one value per draw
slope_draws <- function(Y, x) {
  xc <- x - mean(x)
  as.vector(Y %*% xc) / sum(xc^2)
}
# Posterior summary of one quantity's draws (NA when the quantity is not defined)
summ <- function(v, nm) {
  if (is.null(v) || all(is.na(v))) {
    return(setNames(as.list(rep(NA_real_, 5)), paste0(nm, c("_mean", "_sd", "_q025", "_q975", "_p_pos"))))
  }
  setNames(list(mean(v), sd(v), unname(quantile(v, 0.025)), unname(quantile(v, 0.975)), mean(v > 0)),
           paste0(nm, c("_mean", "_sd", "_q025", "_q975", "_p_pos")))
}

# One species ----------------------------------------------------------------
process_species <- function(sdm_file) {
  sdm_all <- read.csv(sdm_file, stringsAsFactors = FALSE)
  sp   <- sdm_all$species[1]
  code <- sdm_all$species_code[1]

  # Routes: converged, in classes 1-7 (minus any excluded classes)
  sdm <- sdm_all %>%
    filter(route_converged %in% TRUE, .data[[scenario]] %in% setdiff(1:7, exclude_classes)) %>%
    transmute(route, cls = as.integer(.data[[scenario]]))
  class_counts <- setNames(as.list(tabulate(sdm$cls, nbins = 7)), paste0("n_class", 1:7))
  cat_counts   <- setNames(as.list(map_int(categories, ~ sum(sdm$cls %in% .x))), paste0("n_", names(categories)))
  base <- bind_cols(tibble(species = sp, species_code = code), as_tibble(class_counts), as_tibble(cat_counts),
                    tibble(n_routes = nrow(sdm), n_classes = n_distinct(sdm$cls)))
  too_few <- paste("fewer than", min_classes, "SDM classes")
  if (nrow(sdm) == 0 || n_distinct(sdm$cls) < min_classes) return(list(row = mutate(base, status = too_few), draws = NULL))

  # Posterior draws of the route slopes
  sp_f <- species_to_f(sp)
  fit_file <- file.path(rds_dir, paste0(sp_f, "_iCAR_", model_tag, "_", firstYear, "_", lastYear, "_stanfit.rds"))
  ri_file  <- file.path(route_info_dir, paste0(sp_f, "_", model_tag, "_", firstYear, "_", lastYear, "_route_info.rds"))
  if (!file.exists(fit_file)) return(list(row = mutate(base, status = "stanfit file not found"), draws = NULL))
  if (!file.exists(ri_file))  return(list(row = mutate(base, status = "route_info file not found"), draws = NULL))

  ri  <- readRDS(ri_file)
  sdm <- sdm %>% inner_join(ri %>% distinct(route, routeF), by = "route")
  if (nrow(sdm) == 0 || n_distinct(sdm$cls) < min_classes) {
    return(list(row = mutate(base, status = paste(too_few, "after route match")), draws = NULL))
  }

  fit <- readRDS(fit_file)
  dm  <- fit$draws(variables = "beta", format = "draws_matrix")
  rm(fit); invisible(gc())
  set.seed(seed)
  keep <- if (nrow(dm) > n_draws) sort(sample(nrow(dm), n_draws)) else seq_len(nrow(dm))
  cols <- paste0("beta[", sdm$routeF, "]")
  miss <- setdiff(cols, colnames(dm))
  if (length(miss) > 0) return(list(row = mutate(base, status = paste("beta columns missing:", length(miss))), draws = NULL))
  Y <- 100 * (exp(unclass(dm)[keep, cols, drop = FALSE]) - 1)   # draws x routes, %/yr
  rm(dm); invisible(gc())

  # Tracking slopes, in every draw
  x <- unname(class_score[as.character(sdm$cls)])
  s_all <- slope_draws(Y, x)
  low <- sdm$cls <= 4; up <- sdm$cls >= 4
  s_low <- if (n_distinct(x[low]) >= 2) slope_draws(Y[, low, drop = FALSE], x[low]) else NULL
  s_up  <- if (n_distinct(x[up])  >= 2) slope_draws(Y[, up,  drop = FALSE], x[up])  else NULL

  # Category means, in every draw
  m_cat <- map(categories, function(cls) {
    idx <- sdm$cls %in% cls
    if (any(idx)) rowMeans(Y[, idx, drop = FALSE]) else NULL
  })

  row <- bind_cols(base, as_tibble(c(
    summ(s_all, "slope"), summ(s_low, "slope_low"), summ(s_up, "slope_up"),
    summ(m_cat$decreasing, "mean_decreasing"), summ(m_cat$neutral, "mean_neutral"),
    summ(m_cat$increasing, "mean_increasing"),
    list(cor_low_up = if (!is.null(s_low) && !is.null(s_up)) cor(s_low, s_up) else NA_real_,
         n_draws_used = length(keep))
  ))) %>% mutate(status = "ok")

  dr <- NULL
  if (save_draws) {
    k <- seq_len(min(n_draws_saved, length(s_all)))
    dr <- list(slope = s_all[k], slope_low = s_low[k], slope_up = s_up[k])
  }
  list(row = row, draws = dr)
}

# Run ------------------------------------------------------------------------
sdm_files <- list.files(sdm_dir, pattern = paste0("_", model_tag, "_route_trends_sdm\\.csv$"), full.names = TRUE)
if (length(sdm_files) == 0) stop("No *_", model_tag, "_route_trends_sdm.csv files in ", sdm_dir)
if (!is.null(only_species)) {
  sdm_files <- sdm_files[basename(sdm_files) %in% paste0(species_to_f(only_species), "_", model_tag, "_route_trends_sdm.csv")]
}

done <- character()
results <- list(); draws_out <- list()
if (file.exists(out_csv)) {
  prev <- read.csv(out_csv, stringsAsFactors = FALSE)
  if (!"mean_neutral_mean" %in% names(prev)) {
    log_msg("Existing ", basename(out_csv), " has no category means (older version): recomputing all species")
  } else {
    results <- split(prev, seq_len(nrow(prev))) %>% map(as_tibble)
    done <- prev$species[prev$status == "ok"]
    if (save_draws && file.exists(out_rds)) draws_out <- readRDS(out_rds)
    log_msg("Resuming: ", length(done), " species already done")
  }
}

log_msg("Species-level quantities: ", length(sdm_files), " species files; draws from ", rds_dir)
for (i in seq_along(sdm_files)) {
  f <- sdm_files[i]
  sp <- read.csv(f, nrows = 1, stringsAsFactors = FALSE)$species
  if (sp %in% done) next
  res <- tryCatch(process_species(f), error = function(e) {
    list(row = tibble(species = sp, status = paste("error:", conditionMessage(e))), draws = NULL)
  })
  results <- results[!vapply(results, function(r) identical(r$species[1], sp), logical(1))]
  results[[length(results) + 1]] <- res$row
  if (!is.null(res$draws)) draws_out[[res$row$species_code[1]]] <- res$draws
  log_msg(sprintf("[%d/%d] %s: %s%s", i, length(sdm_files), sp, res$row$status[1],
                  if (res$row$status[1] == "ok") sprintf(" (slope %.2f, sd %.2f)", res$row$slope_mean, res$row$slope_sd) else ""))
  if (i %% 25 == 0 || i == length(sdm_files)) {          # checkpoint
    write.csv(bind_rows(results), out_csv, row.names = FALSE)
    if (save_draws) saveRDS(draws_out, out_rds)
  }
}
out <- bind_rows(results)
write.csv(out, out_csv, row.names = FALSE)
if (save_draws) saveRDS(draws_out, out_rds)

# Short summary for the log ------------------------------------------------------
ok <- out %>% filter(status == "ok")
log_msg("Done. Species ok: ", nrow(ok), " of ", nrow(out))
log_msg("Status counts: ", paste(names(table(out$status)), table(out$status), sep = " = ", collapse = "; "))
if (nrow(ok) > 0) {
  log_msg(sprintf("Tracking slope across species (raw, unpooled): median %.3f; %d of %d species with P(slope > 0) > 0.9; %d with < 0.1",
                  median(ok$slope_mean), sum(ok$slope_p_pos > 0.9), nrow(ok), sum(ok$slope_p_pos < 0.1)))
  log_msg(sprintf("Species with routes in Decreasing / Neutral / Increasing suitability: %d / %d / %d (all three: %d)",
                  sum(is.finite(ok$mean_decreasing_mean)), sum(is.finite(ok$mean_neutral_mean)),
                  sum(is.finite(ok$mean_increasing_mean)),
                  sum(is.finite(ok$mean_decreasing_mean) & is.finite(ok$mean_neutral_mean) & is.finite(ok$mean_increasing_mean))))
}
cat("\nWrote:\n ", out_csv, "\n ", if (save_draws) out_rds else "", "\n ", log_file, "\n")
