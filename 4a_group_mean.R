# 4a_group_mean.R
#
# Average-species trend in each SDM range-change group:
#   "The average species trends at X %/yr in its Decreasing-suitability areas,
#    Y %/yr in its Neutral areas and Z %/yr in its Increasing-suitability areas."
#
# Groups (SDM classified-change classes):
#   Decreasing suitability = classes 1-3 | Neutral = class 4 | Increasing suitability = classes 5-7
#   Class 0 (never suitable) and routes without a class are excluded.
#
# Part 0 -- pooled route means (runs anywhere; no Stan).
#   Route counts and the plain mean trend of all routes in each group, pooled
#   across species (posterior-mean route trends), as reported at the start of
#   the Results. Also the numbers of converged and excluded routes.
#
# Species filters (native_only, min_classes) are applied in all three parts, so
# Parts 0, 1 and 2 describe the same set of species.
#
# Part 1 -- per-species group means (needs the Stan fits, *_stanfit.rds from 1a/1b).
#   For each species, in every retained posterior draw of the route-level slopes
#   beta[r], each route's slope is converted to a trend, 100 * (exp(beta) - 1)
#   %/yr, and the trends of the species' routes in each group are averaged.
#   Repeating this over the draws gives a posterior for each group mean, which
#   carries each route's trend uncertainty and the spatial (iCAR) correlation
#   between neighbouring routes. Uses the same route selection, draw subset
#   (n_draws, seed) and min_classes rule as 4b_sdm_tracking_species.R.
#   Skipped when rds_dir does not exist (the fits are stored on another
#   machine); Part 2 then uses the saved means_csv.
#
# Part 2 -- across species (runs anywhere).
#   Each species' group mean (y) and its posterior SD (se) enter the same
#   random-effects meta-analysis as 4c (models/sdm_tracking_meta.stan), with a
#   Student-t distribution of species values (theta integrated out):
#     y[s] ~ Student-t(nu, mu, sqrt(se[s]^2 + tau^2))
#   mu = average species' trend in that group (a robust centre: a few species
#   with extreme trends pull it less than they would a Normal mean); tau = scale
#   of between-species variation; nu = tail weight, estimated. The Student-t is
#   used because species' group means are heavy-tailed: in leave-one-species-out
#   cross-validation it predicted held-out species clearly better than a Normal
#   for the Neutral and Increasing-suitability means and no worse elsewhere.
#   Species sets: "all" = every species with routes in that group;
#                 "common" = species with routes in all three groups.
#   Differences vs Neutral (Decreasing - Neutral, Increasing - Neutral) are
#   fitted on the "common" set with se = sqrt(se_a^2 + se_b^2). The two groups
#   use disjoint routes, so their draws are close to independent; this is an
#   approximation.
#
# Reads:  <rds_dir>/<species>_iCAR_<tag>_<firstYear>_<lastYear>_stanfit.rds          (Part 1)
#         <route_info_dir>/<species>_<tag>_<firstYear>_<lastYear>_route_info.rds    (Part 1)
#         output/species_routes_covariates/per_species_sdm/<species>_<tag>_route_trends_sdm.csv  (step 3; Parts 0, 1)
#         data/spp_names_codes_group_aou_taxa_redlist_diet_migration_invasive.csv  (native_only; Parts 0, 2)
#         models/sdm_tracking_meta.stan                                             (Part 2)
# Writes: output/sdm_tracking/sdm_tracking_pooled_means_<suffix>[_native].csv  (Part 0)
#           one row per group (+ all classified routes): routes, species, pooled mean trend;
#           plus rows with the numbers of converged and excluded routes
#         output/sdm_tracking/sdm_tracking_group_means_<suffix>.csv           (Part 1)
#           one row per species: status, route counts per group, and for each
#           group mean its posterior mean, sd, 2.5% / 97.5% quantiles and P(> 0)
#         output/sdm_tracking/sdm_tracking_group_means_log_<suffix>.txt       (Part 1)
#         output/sdm_tracking/sdm_tracking_category_levels_<suffix>[_native].csv  (Part 2)
#           one row per quantity x species set: mu [95% CrI], P(mu < 0), tau, nu,
#           convergence checks
#
# Part 1 is resumable: species already finished (status "ok") in an existing
# means CSV are skipped, so it can be stopped and restarted.

library(here)
library(tidyverse)
library(posterior)
library(cmdstanr)

here::i_am("4a_group_mean.R")

# Settings -------------------------------------------------------------------
model_tag  <- "anthro"
scenario   <- "rcp45"
firstYear  <- 2010
lastYear   <- 2025
rds_dir        <- here::here("output", "rds")
route_info_dir <- here::here("data", "route_info")
sdm_dir        <- here::here("output", "species_routes_covariates", "per_species_sdm")
out_dir        <- here::here("output", "sdm_tracking")
traits_csv     <- here::here("data", "spp_names_codes_group_aou_taxa_redlist_diet_migration_invasive.csv")

parts <- c(0, 1, 2)     # which parts to run; c(0) runs only the pooled means (no Stan needed)

# Species filters, applied in Parts 0, 1 and 2
native_only <- TRUE     # TRUE = drop US-RIIS listed non-native species (Established / Widespread invasive)
# Species need routes in at least this many distinct SDM classes. 2 matches
# 4b_sdm_tracking_species.R; 1 also includes species found in a single class.
min_classes <- 2

# Part 1
n_draws     <- 2000     # posterior draws kept per species (random subset)
seed        <- 20261014 # same as 4b_sdm_tracking_species.R
only_species <- NULL    # e.g. c("Song Sparrow", "Wild Turkey") for a quick test run
# SDM classes to drop (sensitivity analysis). integer(0) = main analysis;
# c(7) = drop colonization, adding "_no7" to the output file names.
exclude_classes <- integer(0)

# Part 2
prior_scale <- 3        # %/yr
chains <- 4; iter_warmup <- 1000; iter_sampling <- 2000; adapt_delta <- 0.99
meta_seed <- 20261112

groups <- list(decreasing = 1:3, neutral = 4L, increasing = 5:7)

run_tag   <- if (length(exclude_classes) > 0) paste0("_no", paste(sort(exclude_classes), collapse = "")) else ""
suffix    <- paste0(paste(model_tag, scenario, firstYear, lastYear, sep = "_"), run_tag)
means_csv <- file.path(out_dir, paste0("sdm_tracking_group_means_", suffix, ".csv"))
log_file  <- file.path(out_dir, paste0("sdm_tracking_group_means_log_", suffix, ".txt"))
native_tag <- if (native_only) "_native" else ""
pooled_csv <- file.path(out_dir, paste0("sdm_tracking_pooled_means_", suffix, native_tag, ".csv"))
levels_csv <- file.path(out_dir, paste0("sdm_tracking_category_levels_", suffix, native_tag, ".csv"))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

species_to_f <- function(sp) gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)

log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), "  ", paste0(...))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

summ <- function(v, nm) {
  if (is.null(v) || all(is.na(v))) {
    return(setNames(as.list(rep(NA_real_, 5)), paste0(nm, c("_mean", "_sd", "_q025", "_q975", "_p_pos"))))
  }
  setNames(list(mean(v), sd(v), unname(quantile(v, 0.025)), unname(quantile(v, 0.975)), mean(v > 0)),
           paste0(nm, c("_mean", "_sd", "_q025", "_q975", "_p_pos")))
}

non_native <- if (native_only) {
  read.csv(traits_csv, stringsAsFactors = FALSE) %>%
    filter(in_bbs %in% TRUE, invasive != "Not listed") %>% pull(Code)
} else character()

# ============================================================================
# Part 0: pooled route means (no Stan)
# ============================================================================
if (0 %in% parts) {
  route_files <- list.files(sdm_dir, pattern = paste0("_", model_tag, "_route_trends_sdm\\.csv$"), full.names = TRUE)
  routes <- map_dfr(route_files, function(f) {
    read.csv(f, stringsAsFactors = FALSE) %>%
      transmute(species_code, route_converged, trend, cls = .data[[scenario]])
  })
  dropped <- intersect(unique(routes$species_code), non_native)
  routes <- routes %>% filter(!species_code %in% non_native)

  conv <- routes %>% filter(route_converged %in% TRUE)
  classified <- conv %>%
    filter(cls %in% setdiff(1:7, exclude_classes)) %>%
    mutate(group = case_when(cls <= 3 ~ "Decreasing suitability",
                             cls == 4 ~ "Neutral",
                             TRUE     ~ "Increasing suitability"),
           group = factor(group, levels = c("Decreasing suitability", "Neutral", "Increasing suitability")))
  few_classes <- classified %>% group_by(species_code) %>%
    summarise(k = n_distinct(cls), .groups = "drop") %>% filter(k < min_classes) %>% pull(species_code)
  n_classified_all <- nrow(classified)
  classified <- classified %>% filter(!species_code %in% few_classes)

  pooled <- bind_rows(
    classified %>% group_by(group) %>%
      summarise(n_routes = n(), n_species = n_distinct(species_code), pooled_mean_trend = mean(trend), .groups = "drop") %>%
      mutate(group = as.character(group)),
    tibble(group = "All classified routes", n_routes = nrow(classified),
           n_species = n_distinct(classified$species_code), pooled_mean_trend = mean(classified$trend)),
    tibble(group = "Converged routes (any class)", n_routes = nrow(conv),
           n_species = n_distinct(conv$species_code), pooled_mean_trend = NA_real_),
    tibble(group = "Excluded: class 0, no class or excluded class", n_routes = nrow(conv) - n_classified_all,
           n_species = NA_integer_, pooled_mean_trend = NA_real_),
    tibble(group = paste("Excluded: species with fewer than", min_classes, "SDM classes"),
           n_routes = n_classified_all - nrow(classified), n_species = length(few_classes), pooled_mean_trend = NA_real_)
  )
  write.csv(pooled, pooled_csv, row.names = FALSE)
  cat("\n=== Part 0: pooled route means", if (native_only) "(non-native species excluded:" else "",
      if (native_only) paste0(length(dropped), " species: ", paste(dropped, collapse = ", "), ")") else "", "===\n")
  print(as.data.frame(pooled), digits = 4, row.names = FALSE)
  cat("Wrote:", pooled_csv, "\n")
}

# ============================================================================
# Part 1: per-species group means from the posterior draws
# ============================================================================
process_species <- function(sdm_file) {
  sdm_all <- read.csv(sdm_file, stringsAsFactors = FALSE)
  sp   <- sdm_all$species[1]
  code <- sdm_all$species_code[1]

  sdm <- sdm_all %>%
    filter(route_converged %in% TRUE, .data[[scenario]] %in% setdiff(1:7, exclude_classes)) %>%
    transmute(route, cls = as.integer(.data[[scenario]]))
  if (code %in% non_native) return(tibble(species = sp, species_code = code, status = "non-native species (excluded)"))
  counts <- map_int(groups, ~ sum(sdm$cls %in% .x))
  base <- tibble(species = sp, species_code = code,
                 n_decreasing = counts[["decreasing"]], n_neutral = counts[["neutral"]],
                 n_increasing = counts[["increasing"]], n_routes = nrow(sdm),
                 n_classes = n_distinct(sdm$cls))
  if (n_distinct(sdm$cls) < min_classes || nrow(sdm) == 0) {
    return(mutate(base, status = paste("fewer than", min_classes, "SDM classes")))
  }

  sp_f <- species_to_f(sp)
  fit_file <- file.path(rds_dir, paste0(sp_f, "_iCAR_", model_tag, "_", firstYear, "_", lastYear, "_stanfit.rds"))
  ri_file  <- file.path(route_info_dir, paste0(sp_f, "_", model_tag, "_", firstYear, "_", lastYear, "_route_info.rds"))
  if (!file.exists(fit_file)) return(mutate(base, status = "stanfit file not found"))
  if (!file.exists(ri_file))  return(mutate(base, status = "route_info file not found"))

  ri  <- readRDS(ri_file)
  sdm <- sdm %>% inner_join(ri %>% distinct(route, routeF), by = "route")
  if (n_distinct(sdm$cls) < min_classes || nrow(sdm) == 0) {
    return(mutate(base, status = paste("fewer than", min_classes, "SDM classes after route match")))
  }

  fit <- readRDS(fit_file)
  dm  <- fit$draws(variables = "beta", format = "draws_matrix")
  rm(fit); invisible(gc())
  set.seed(seed)
  keep <- if (nrow(dm) > n_draws) sort(sample(nrow(dm), n_draws)) else seq_len(nrow(dm))
  cols <- paste0("beta[", sdm$routeF, "]")
  miss <- setdiff(cols, colnames(dm))
  if (length(miss) > 0) return(mutate(base, status = paste("beta columns missing:", length(miss))))
  Y <- 100 * (exp(unclass(dm)[keep, cols, drop = FALSE]) - 1)   # draws x routes, %/yr
  rm(dm); invisible(gc())

  # mean trend of the species' routes in each group, in every draw
  group_draws <- map(groups, function(cls) {
    idx <- sdm$cls %in% cls
    if (any(idx)) rowMeans(Y[, idx, drop = FALSE]) else NULL
  })

  bind_cols(base, as_tibble(c(
    summ(group_draws$decreasing, "mean_decreasing"),
    summ(group_draws$neutral,    "mean_neutral"),
    summ(group_draws$increasing, "mean_increasing"),
    list(n_draws_used = length(keep))
  ))) %>% mutate(status = "ok")
}

if (1 %in% parts && dir.exists(rds_dir)) {
  sdm_files <- list.files(sdm_dir, pattern = paste0("_", model_tag, "_route_trends_sdm\\.csv$"), full.names = TRUE)
  if (length(sdm_files) == 0) stop("No *_", model_tag, "_route_trends_sdm.csv files in ", sdm_dir)
  if (!is.null(only_species)) {
    sdm_files <- sdm_files[basename(sdm_files) %in% paste0(species_to_f(only_species), "_", model_tag, "_route_trends_sdm.csv")]
  }

  results <- list()
  done <- character()
  if (file.exists(means_csv)) {
    prev <- read.csv(means_csv, stringsAsFactors = FALSE)
    results <- split(prev, seq_len(nrow(prev))) %>% map(as_tibble)
    done <- prev$species[prev$status == "ok"]
    log_msg("Resuming: ", length(done), " species already done")
  }

  log_msg("Group means start: ", length(sdm_files), " species files; draws from ", rds_dir)
  for (i in seq_along(sdm_files)) {
    f <- sdm_files[i]
    first <- read.csv(f, nrows = 1, stringsAsFactors = FALSE)
    sp <- first$species
    if (sp %in% done && !first$species_code %in% non_native) next
    row <- tryCatch(process_species(f), error = function(e) {
      tibble(species = sp, status = paste("error:", conditionMessage(e)))
    })
    results <- results[!vapply(results, function(r) identical(r$species[1], sp), logical(1))]
    results[[length(results) + 1]] <- row
    log_msg(sprintf("[%d/%d] %s: %s", i, length(sdm_files), sp, row$status[1]))
    if (i %% 25 == 0 || i == length(sdm_files)) write.csv(bind_rows(results), means_csv, row.names = FALSE)   # checkpoint
  }
  write.csv(bind_rows(results), means_csv, row.names = FALSE)
} else if (1 %in% parts) {
  cat("No fits in", rds_dir, "-- Part 1 skipped; Part 2 will use the saved", means_csv, "\n")
}

# ============================================================================
# Part 2: average species in each group (meta-analysis across species)
# ============================================================================
if (2 %in% parts) {
  if (!file.exists(means_csv)) {
    stop("No saved ", means_csv, " -- run Part 1 where the model fits are stored.")
  }
  means <- read.csv(means_csv, stringsAsFactors = FALSE)
  ok <- means %>% filter(status == "ok")
  cat("Species with group means:", nrow(ok), "of", nrow(means),
      "(with Decreasing / Neutral / Increasing routes:",
      sum(is.finite(ok$mean_decreasing_mean)), "/", sum(is.finite(ok$mean_neutral_mean)), "/",
      sum(is.finite(ok$mean_increasing_mean)), ")\n")

  s1 <- ok
  if (native_only) {
    cat("Dropping non-native species:", paste(intersect(s1$species_code, non_native), collapse = ", "), "\n")
    s1 <- s1 %>% filter(!species_code %in% non_native)
  }

  mod <- cmdstan_model(here::here("models", "sdm_tracking_meta.stan"))

  fit_mu <- function(y, se) {
    ok <- is.finite(y) & is.finite(se) & se > 0
    run <- function() mod$sample(data = list(N = sum(ok), y = y[ok], se = se[ok], K1 = 0, g1 = rep(0L, sum(ok)),
                                  K2 = 0, g2 = rep(0L, sum(ok)), student = 1, explicit_t = 0,
                                  fixed1 = 0, prior_scale = prior_scale),
                      chains = chains, parallel_chains = chains, iter_warmup = iter_warmup,
                      iter_sampling = iter_sampling, adapt_delta = adapt_delta, seed = meta_seed,
                      refresh = 0, show_messages = FALSE)
    fit <- run()
    # Windows can block a freshly compiled model executable on its first launch; retry once
    if (all(fit$return_codes() != 0)) { Sys.sleep(10); fit <- run() }
    dr <- fit$draws(c("mu", "tau", "nu"), format = "draws_df")
    tibble(n_species = sum(ok), nu = median(dr$nu),
           mu = mean(dr$mu), mu_lo = quantile(dr$mu, 0.025), mu_hi = quantile(dr$mu, 0.975),
           p_mu_neg = mean(dr$mu < 0), tau = mean(dr$tau),
           max_rhat = max(fit$summary(c("mu", "tau"))$rhat), divergences = sum(fit$diagnostic_summary()$num_divergent))
  }

  common <- s1 %>% filter(is.finite(mean_decreasing_mean), is.finite(mean_neutral_mean), is.finite(mean_increasing_mean))
  resp <- list(
    list("all", "Decreasing suitability (classes 1-3)", s1, "mean_decreasing_mean", "mean_decreasing_sd"),
    list("all", "Neutral (class 4)",                    s1, "mean_neutral_mean",    "mean_neutral_sd"),
    list("all", "Increasing suitability (classes 5-7)", s1, "mean_increasing_mean", "mean_increasing_sd"),
    list("common", "Decreasing suitability (classes 1-3)", common, "mean_decreasing_mean", "mean_decreasing_sd"),
    list("common", "Neutral (class 4)",                    common, "mean_neutral_mean",    "mean_neutral_sd"),
    list("common", "Increasing suitability (classes 5-7)", common, "mean_increasing_mean", "mean_increasing_sd")
  )
  common <- common %>% mutate(cs_mean = mean_decreasing_mean - mean_neutral_mean, cs_sd = sqrt(mean_decreasing_sd^2 + mean_neutral_sd^2),
                              es_mean = mean_increasing_mean - mean_neutral_mean, es_sd = sqrt(mean_increasing_sd^2 + mean_neutral_sd^2))
  resp <- c(resp, list(
    list("common", "Decreasing suitability - Neutral", common, "cs_mean", "cs_sd"),
    list("common", "Increasing suitability - Neutral", common, "es_mean", "es_sd")
  ))

  res <- map_dfr(resp, function(r) {
    cat("Fitting:", r[[1]], r[[2]], "\n")
    fit_mu(r[[3]][[r[[4]]]], r[[3]][[r[[5]]]]) %>%
      mutate(species_set = r[[1]], quantity = r[[2]], .before = 1)
  })
  res <- res %>% mutate(across(where(is.double), unname))
  write.csv(as.data.frame(res), levels_csv, row.names = FALSE)
  print(res %>% mutate(across(where(is.double), ~ round(.x, 3))), n = Inf, width = Inf)
  cat("\nWrote:", levels_csv, "\n")
}
