# 4a_sdm_tracking_routes.R
#
# Step 1 of the SDM-tracking analysis (Methods, Statistical analysis): the data.
# Selects the route-by-species trend estimates that enter the analysis and
# reports how many routes and species remain, plus the plain (pooled) mean
# trend of all routes in each range-change category, as quoted at the start of
# the Results. No Stan: runs anywhere.
#
# Selection (the same rules as 4b and 4c):
#   - route-level trends from the anthropogenic-cover model (model_tag), RCP 4.5 classes;
#   - routes whose own parameters converged (route_converged);
#   - routes in SDM classes 1-7 (class 0 = never suitable, and routes without a class, are excluded);
#   - native species only (native_only; US-RIIS non-native species dropped);
#   - species with routes in at least min_classes distinct classes.
#
# Range-change categories: Decreasing suitability = classes 1-3 | Neutral = class 4 |
# Increasing suitability = classes 5-7.
#
# Reads:  output/species_routes_covariates/per_species_sdm/<species>_<tag>_route_trends_sdm.csv  (step 3)
#         data/spp_names_codes_group_aou_taxa_redlist_diet_migration_invasive.csv  (non-native list)
# Writes: output/sdm_tracking/sdm_tracking_pooled_means_<suffix>[_native].csv
#           one row per category (+ all classified routes): routes, species, pooled mean trend;
#           plus rows with the numbers of converged and excluded routes and species
#
# Next: 4b_sdm_tracking_species.R (species-level quantities from the posterior draws).

library(here)
library(tidyverse)

here::i_am("4a_sdm_tracking_routes.R")

# Settings -------------------------------------------------------------------
model_tag <- "anthro"
scenario  <- "rcp45"
firstYear <- 2010
lastYear  <- 2025
sdm_dir    <- here::here("output", "species_routes_covariates", "per_species_sdm")
out_dir    <- here::here("output", "sdm_tracking")
traits_csv <- here::here("data", "spp_names_codes_group_aou_taxa_redlist_diet_migration_invasive.csv")

native_only <- TRUE   # TRUE = drop US-RIIS listed non-native species (Established / Widespread invasive)
min_classes <- 2      # species need routes in at least this many SDM classes (as in 4b)
# SDM classes to drop (sensitivity analysis). integer(0) = main analysis;
# c(7) = drop colonization, adding "_no7" to the output file name. Must match 4b and 4c.
exclude_classes <- integer(0)

run_tag    <- if (length(exclude_classes) > 0) paste0("_no", paste(sort(exclude_classes), collapse = "")) else ""
suffix     <- paste0(paste(model_tag, scenario, firstYear, lastYear, sep = "_"), run_tag)
pooled_csv <- file.path(out_dir, paste0("sdm_tracking_pooled_means_", suffix, if (native_only) "_native" else "", ".csv"))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

non_native <- if (native_only) {
  read.csv(traits_csv, stringsAsFactors = FALSE) %>%
    filter(in_bbs %in% TRUE, invasive != "Not listed") %>% pull(Code)
} else character()

# Routes -----------------------------------------------------------------------
route_files <- list.files(sdm_dir, pattern = paste0("_", model_tag, "_route_trends_sdm\\.csv$"), full.names = TRUE)
if (length(route_files) == 0) stop("No *_", model_tag, "_route_trends_sdm.csv files in ", sdm_dir)
routes <- map_dfr(route_files, function(f) {
  read.csv(f, stringsAsFactors = FALSE) %>%
    transmute(species_code, route_converged, trend, cls = .data[[scenario]])
})
dropped <- intersect(unique(routes$species_code), non_native)
routes  <- routes %>% filter(!species_code %in% non_native)

conv <- routes %>% filter(route_converged %in% TRUE)
classified <- conv %>%
  filter(cls %in% setdiff(1:7, exclude_classes)) %>%
  mutate(category = case_when(cls <= 3 ~ "Decreasing suitability",
                              cls == 4 ~ "Neutral",
                              TRUE     ~ "Increasing suitability"),
         category = factor(category, levels = c("Decreasing suitability", "Neutral", "Increasing suitability")))
n_classified_all <- nrow(classified)
n_sp_classified  <- n_distinct(classified$species_code)
few_classes <- classified %>% group_by(species_code) %>%
  summarise(k = n_distinct(cls), .groups = "drop") %>% filter(k < min_classes) %>% pull(species_code)
classified <- classified %>% filter(!species_code %in% few_classes)

# Counts and pooled means ---------------------------------------------------------
pooled <- bind_rows(
  classified %>% group_by(category) %>%
    summarise(n_routes = n(), n_species = n_distinct(species_code), pooled_mean_trend = mean(trend), .groups = "drop") %>%
    transmute(group = as.character(category), n_routes, n_species, pooled_mean_trend),
  tibble(group = "All classified routes", n_routes = nrow(classified),
         n_species = n_distinct(classified$species_code), pooled_mean_trend = mean(classified$trend)),
  tibble(group = "Converged routes (any class)", n_routes = nrow(conv),
         n_species = n_distinct(conv$species_code), pooled_mean_trend = NA_real_),
  tibble(group = "Excluded: class 0, no class or excluded class", n_routes = nrow(conv) - n_classified_all,
         n_species = n_distinct(conv$species_code) - n_sp_classified, pooled_mean_trend = NA_real_),
  tibble(group = paste("Excluded: species with fewer than", min_classes, "SDM classes"),
         n_routes = n_classified_all - nrow(classified), n_species = length(few_classes), pooled_mean_trend = NA_real_)
)
write.csv(pooled, pooled_csv, row.names = FALSE)

cat("\n=== Routes and species in the SDM-tracking analysis ===\n")
if (native_only) cat("Non-native species excluded (", length(dropped), "): ", paste(dropped, collapse = ", "), "\n", sep = "")
print(as.data.frame(pooled), digits = 4, row.names = FALSE)
cat("Wrote:", pooled_csv, "\n")
