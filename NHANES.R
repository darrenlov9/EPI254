########################################################################
# Class:       2026 Winter Nutritional Epidemiology Lab 1 (R Version)
# Converted from SAS by: Claude
# Original SAS by: Cheng-Tzu Hsieh (partially by Xiang Li)
# Instructor:  Dr. Liwei Chen
########################################################################

library(haven)
library(dplyr)
library(purrr)
setwd("C:/Users/Darren/Downloads/NHANES Data")

# -----------------------------------------------------------------------
# Load all local NHANES XPT/SAS files from wave subdirectories
# -----------------------------------------------------------------------
wave_dirs <- list.files(path = ".", full.names = TRUE)
wave_dirs <- wave_dirs[dir.exists(wave_dirs)]

nhanes_data <- list()
for (wave in wave_dirs) {
  files <- list.files(
    wave,
    pattern = "\\.(xpt|sas7bdat)$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(files) == 0) next
  wave_name <- basename(wave)
  nhanes_data[[wave_name]] <- lapply(files, function(f) {
    if (grepl("\\.xpt$", f, ignore.case = TRUE)) read_xpt(f) else read_sas(f)
  })
  names(nhanes_data[[wave_name]]) <- tools::file_path_sans_ext(basename(files))
}

waves <- c("2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2020",'2021-2023')

# -----------------------------------------------------------------------
# Helper: pool a component across waves
#   - Looks for the dataset by name, accounting for CDC naming conventions
#   - Adds a `cycle` column
# -----------------------------------------------------------------------

# Wave name → file suffix (matches CDC naming convention)
wave_suffix <- c(
  "2009-2010" = "_F",
  "2011-2012" = "_G",
  "2013-2014" = "_H",
  "2015-2016" = "_I",
  "2017-2020" = "P_",   # pre-pandemic cycle uses P_ prefix instead of suffix
  "2021-2023" = '_L'
)

get_component <- function(wave_name, prefix) {
  wave_files <- nhanes_data[[wave_name]]
  if (is.null(wave_files)) return(NULL)
  
  # 2017-2020: P_PREFIX naming; earlier cycles: PREFIX_X
  if (wave_name == "2017-2020") {
    pattern <- paste0("^P_", prefix, "$")
  } else {
    sfx     <- wave_suffix[[wave_name]]
    pattern <- paste0("^", prefix, sfx, "$")
  }
  
  match_name <- names(wave_files)[grepl(pattern, names(wave_files), ignore.case = TRUE)]
  
  if (length(match_name) == 0) {
    message("  [WARNING] '", prefix, "' not found in wave '", wave_name, "'")
    return(NULL)
  }
  wave_files[[match_name[1]]] |> mutate(cycle = wave_name)
}

pool_component <- function(prefix, wave_list = waves) {
  message("Pooling: ", prefix)
  map(wave_list, \(w) get_component(w, prefix)) |>
    compact() |>
    bind_rows()
}

# -----------------------------------------------------------------------
# Pool all components
# -----------------------------------------------------------------------
DEMO_pool   <- pool_component("DEMO")
DR1IFF_pool <- pool_component("DR1IFF")
DR2IFF_pool <- pool_component("DR2IFF")
DR1TOT_pool <- pool_component("DR1TOT")
DR2TOT_pool <- pool_component("DR2TOT")
BMX_pool    <- pool_component("BMX")
HDL_pool    <- pool_component("HDL")
TRIGLY_pool <- pool_component("TRIGLY")
TCHOL_pool  <- pool_component("TCHOL")
BPQ_pool    <- pool_component("BPQ")
DIQ_pool    <- pool_component("DIQ")
MCQ_pool    <- pool_component("MCQ")
PAQ_pool    <- pool_component("PAQ")
RHQ_pool    <- pool_component("RHQ")
SMQ_pool    <- pool_component("SMQ")
GHB_pool    <- pool_component("GHB")
FSQ_pool    <- pool_component("FSQ")

# Special case: Blood Pressure (file names changed between cycles)
pool_bpx <- function(wave_name) {
  wave_files <- nhanes_data[[wave_name]]
  if (is.null(wave_files)) return(NULL)
  
  # 2009-2016: BPX_F/G/H/I  |  2017-2020: P_BPXO
  if (wave_name == "2017-2020") {
    pattern <- "^P_BPXO$"
  } else {
    sfx     <- wave_suffix[[wave_name]]
    pattern <- paste0("^BPX", sfx, "$")
  }
  
  match_name <- names(wave_files)[grepl(pattern, names(wave_files), ignore.case = TRUE)]
  if (length(match_name) == 0) {
    message("  [WARNING] BPX not found in wave '", wave_name, "'")
    return(NULL)
  }
  wave_files[[match_name[1]]] |> mutate(cycle = wave_name)
}

BPX_pool <- map(waves, pool_bpx) |> compact() |> bind_rows()


# -----------------------------------------------------------------------
# Dietary: Targeted food consumption grams (Day 1 & Day 2)
# NOTE: Update food code vectors below for your own food of interest.
#       See the NHANES Food Code Book for details.
# -----------------------------------------------------------------------

# Helper: extract the first N digits of an 8-digit zero-padded food code
first_n_digits <- function(foodcode, n) {
  as.integer(substr(sprintf("%08d", as.integer(foodcode)), 1, n))
}

# Example 1: Regular carbonated soft drinks (page 100, food code book)
softdrink_codes <- c(
  9241011, 9241031, 9241033, 9241034, 9241036,
  9241039, 9241041, 9241051, 9241055, 9241061,
  9241071, 9241081, 9241151, 9241152, 9241601,
  9241701, 9243100, 9243200, 9243300
)

# Example 2: Low-calorie carbonated soft drinks (page 101, food code book)
lcsoftdrink_codes <- c(
  9240010, 9241021, 9241025, 9241030, 9241032,
  9241035, 9241037, 9241040, 9241042, 9241052,
  9241056, 9241062, 9241072, 9241082, 9241761,
  9241162
)

# -- Day 1 --
DR1IFF_sumgrams <- DR1IFF_pool |>
  dplyr::select(cycle, SEQN, WTDRD1, WTDR2D, DR1IFDCD, DR1IGRMS) |>
  mutate(
    foodcode_7dig      = first_n_digits(DR1IFDCD, 7),
    softdrink_day1_g   = if_else(foodcode_7dig %in% softdrink_codes,   DR1IGRMS, 0),
    lcsoftdrink_day1_g = if_else(foodcode_7dig %in% lcsoftdrink_codes, DR1IGRMS, 0)
  ) |>
  group_by(SEQN) |>
  summarise(
    softdrink_day1_g   = sum(softdrink_day1_g,   na.rm = TRUE),
    lcsoftdrink_day1_g = sum(lcsoftdrink_day1_g, na.rm = TRUE),
    .groups = "drop"
  )

# -- Day 2 --
DR2IFF_sumgrams <- DR2IFF_pool |>
  dplyr::select(cycle, SEQN, WTDRD1, WTDR2D, DR2IFDCD, DR2IGRMS) |>
  mutate(
    foodcode_7dig      = first_n_digits(DR2IFDCD, 7),
    softdrink_day2_g   = if_else(foodcode_7dig %in% softdrink_codes,   DR2IGRMS, 0),
    lcsoftdrink_day2_g = if_else(foodcode_7dig %in% lcsoftdrink_codes, DR2IGRMS, 0)
  ) |>
  group_by(SEQN) |>
  summarise(
    softdrink_day2_g   = sum(softdrink_day2_g,   na.rm = TRUE),
    lcsoftdrink_day2_g = sum(lcsoftdrink_day2_g, na.rm = TRUE),
    .groups = "drop"
  )

# -- Merge Day 1 + Day 2, compute 2-day means --
DR_IFF_2day <- full_join(DR1IFF_sumgrams, DR2IFF_sumgrams, by = "SEQN") |>
  mutate(
    D12mean_g_softdrink = case_when(
      !is.na(softdrink_day1_g) & !is.na(softdrink_day2_g) ~
        (softdrink_day1_g + softdrink_day2_g) / 2,
      !is.na(softdrink_day1_g) ~ softdrink_day1_g,
      !is.na(softdrink_day2_g) ~ softdrink_day2_g
    ),
    D12mean_g_lcsoftdrink = case_when(
      !is.na(lcsoftdrink_day1_g) & !is.na(lcsoftdrink_day2_g) ~
        (lcsoftdrink_day1_g + lcsoftdrink_day2_g) / 2,
      !is.na(lcsoftdrink_day1_g) ~ lcsoftdrink_day1_g,
      !is.na(lcsoftdrink_day2_g) ~ lcsoftdrink_day2_g
    )
  )


# -----------------------------------------------------------------------
# Total Nutrients: Day 1 & Day 2, then 2-day means
# -----------------------------------------------------------------------
nutrient_stems <- c("TKCAL","TPROT","TCARB","TSUGR","TFIBE",
                    "TTFAT","TSFAT","TMFAT","TPFAT","TCHOL",
                    "TFA","TCALC","TIRON","TSODI","TCAFF","TALCO")

nutrient_labels <- c(
  TKCAL_2dmean = "2-day mean total energy (kcal)",
  TPROT_2dmean = "2-day mean total protein (g)",
  TCARB_2dmean = "2-day mean total carbohydrate (g)",
  TSUGR_2dmean = "2-day mean total sugars (g)",
  TFIBE_2dmean = "2-day mean total fiber (g)",
  TTFAT_2dmean = "2-day mean total fat (g)",
  TSFAT_2dmean = "2-day mean saturated fat (g)",
  TMFAT_2dmean = "2-day mean monounsaturated fat (g)",
  TPFAT_2dmean = "2-day mean polyunsaturated fat (g)",
  TCHOL_2dmean = "2-day mean cholesterol (mg)",
  TFA_2dmean   = "2-day mean Folic acid (mcg)",
  TCALC_2dmean = "2-day mean calcium (mg)",
  TIRON_2dmean = "2-day mean iron (mg)",
  TSODI_2dmean = "2-day mean sodium (mg)",
  TCAFF_2dmean = "2-day mean caffeine (mg)",
  TALCO_2dmean = "2-day mean alcohol (g)"
)

nutrient_vars_d1 <- paste0("DR1", nutrient_stems)
nutrient_vars_d2 <- paste0("DR2", nutrient_stems)

DR1TOT_pool_select <- DR1TOT_pool |>
  filter(DR1DRSTZ %in% c(1, 4)) |>   # 1 = reliable; 4 = breast-milk
  dplyr::select(cycle, SEQN, DR1DRSTZ, WTDRD1, WTDR2D, any_of(nutrient_vars_d1))

DR2TOT_pool_select <- DR2TOT_pool |>
  filter(DR2DRSTZ %in% c(1, 4)) |>
  dplyr::select(cycle, SEQN, DR2DRSTZ, WTDRD1, WTDR2D, any_of(nutrient_vars_d2))

DR_TOT_2day <- full_join(
  DR1TOT_pool_select,
  DR2TOT_pool_select,
  by = c("cycle", "SEQN")
) |>
  filter(!is.na(DR1DRSTZ) | !is.na(DR2DRSTZ))

# Compute 2-day mean for each nutrient
for (stem in nutrient_stems) {
  v1  <- paste0("DR1", stem)
  v2  <- paste0("DR2", stem)
  avg <- paste0(stem, "_2dmean")
  DR_TOT_2day[[avg]] <- case_when(
    !is.na(DR_TOT_2day[[v1]]) & !is.na(DR_TOT_2day[[v2]]) ~
      (DR_TOT_2day[[v1]] + DR_TOT_2day[[v2]]) / 2,
    !is.na(DR_TOT_2day[[v1]]) ~ DR_TOT_2day[[v1]],
    !is.na(DR_TOT_2day[[v2]]) ~ DR_TOT_2day[[v2]],
    TRUE ~ NA_real_
  )
  attr(DR_TOT_2day[[avg]], "label") <- nutrient_labels[[avg]]
}

DR2TOT_pool <- DR2TOT_pool |>
  mutate(
    WTDR2D_unified = case_when(
      cycle == "2017-2020" ~ WTDR2DPP / 2,  # convert 4-year to 2-year equivalent
      TRUE ~ WTDR2D
    )
  )

# -----------------------------------------------------------------------
# Merge all pooled components into one analytic dataset
# DEMO_pool is the backbone; all others are left-joined by SEQN
# -----------------------------------------------------------------------

drop_cycle <- function(df) dplyr::select(df, -any_of("cycle"))

nhanes_link <- DEMO_pool |>
  left_join(DR_IFF_2day,                by = "SEQN") |>
  left_join(DR_TOT_2day  |> drop_cycle(), by = "SEQN") |>
  left_join(DR2TOT_pool |> drop_cycle(), by = "SEQN") |>
  left_join(BMX_pool     |> drop_cycle(), by = "SEQN") |>
  # left_join(HDL_pool     |> drop_cycle(), by = "SEQN") |>
  left_join(TRIGLY_pool  |> drop_cycle(), by = "SEQN") |>
  left_join(TCHOL_pool   |> drop_cycle(), by = "SEQN") |>
  left_join(GHB_pool     |> drop_cycle(), by = "SEQN") |>
  left_join(FSQ_pool     |> drop_cycle(), by = "SEQN") |>
  # left_join(GLU_pool     |> drop_cycle(), by = "SEQN") |>
  left_join(BPQ_pool     |> drop_cycle(), by = "SEQN") |>
  left_join(DIQ_pool     |> drop_cycle(), by = "SEQN") |>
  # left_join(MCQ_pool     |> drop_cycle(), by = "SEQN") |>
  left_join(PAQ_pool     |> drop_cycle(), by = "SEQN") |>
  # left_join(RHQ_pool     |> drop_cycle(), by = "SEQN") |>
  left_join(SMQ_pool     |> drop_cycle(), by = "SEQN") |>
  left_join(BPX_pool     |> drop_cycle(), by = "SEQN")

message("Final dataset: ", nrow(nhanes_link), " rows x ", ncol(nhanes_link), " cols")

nhanes_link <- nhanes_link %>% filter(RIDAGEYR >= 18)

library(dplyr)

nhanes_link <- nhanes_link |>
  mutate(
    # Keep only plausible values (optional but helps prevent garbage driving tertiles)
    TFIBE_2dmean = if_else(TFIBE_2dmean < 0 | TFIBE_2dmean > 100, NA_real_, TFIBE_2dmean),
    TSFAT_2dmean = if_else(TSFAT_2dmean < 0 | TSFAT_2dmean > 200, NA_real_, TSFAT_2dmean),
    TCHOL_2dmean = if_else(TCHOL_2dmean < 0 | TCHOL_2dmean > 3000, NA_real_, TCHOL_2dmean)
  ) |>
  # Create tertiles (computed on non-missing)
  mutate(
    fiber_tertile = ntile(TFIBE_2dmean, 3),     # 1=low, 3=high
    sfat_tertile  = ntile(TSFAT_2dmean, 3),     # 1=low, 3=high
    chol_tertile  = ntile(TCHOL_2dmean, 3)      # 1=low, 3=high
  ) |>
  mutate(
    # Score: high fiber = good; low sat fat = good; low cholesterol = good
    plant_forward_score =
      case_when(
        is.na(fiber_tertile) | is.na(sfat_tertile) | is.na(chol_tertile) ~ NA_real_,
        TRUE ~ (fiber_tertile - 1) + (3 - sfat_tertile) + (3 - chol_tertile)
      ),
    # Range check: 0 (least plant-forward) to 6 (most plant-forward)
    plant_forward_score = as.integer(plant_forward_score),
    
    # Optional: make groups for easier interpretation
    plant_forward_tertile = ntile(plant_forward_score, 3)
  )

nhanes_link <- nhanes_link |>
  mutate(
    animal_forward_score =
      case_when(
        is.na(fiber_tertile) | is.na(sfat_tertile) | is.na(chol_tertile) ~ NA_real_,
        TRUE ~ (3 - fiber_tertile) + (sfat_tertile - 1) + (chol_tertile - 1)
      ),
    # Range: 0 (least animal-forward) to 6 (most animal-forward)
    animal_forward_score = as.integer(animal_forward_score),
    
    # Tertiles for modeling
    animal_forward_tertile = ntile(animal_forward_score, 3)
  )

# convert typical NHANES dummy missing codes to NA
nhanes_na <- function(x) {
  dplyr::na_if(x, 7777) %>%
    dplyr::na_if(9999) %>%
    dplyr::na_if(77) %>%
    dplyr::na_if(99)
}

# compute weekly minutes for legacy style (typical week * typical day)
derive_lt_minutes_legacy <- function(df) {
  df %>%
    dplyr::mutate(
      mod_yes = nhanes_na(PAQ665),
      mod_days = nhanes_na(PAQ670),
      mod_min_day = nhanes_na(PAD675),
      
      vig_yes = nhanes_na(PAQ650),
      vig_days = nhanes_na(PAQ655),
      vig_min_day = nhanes_na(PAD660),
      
      lt_mod_min_wk = dplyr::case_when(
        mod_yes == 2 ~ 0,
        mod_yes == 1 & !is.na(mod_days) & !is.na(mod_min_day) ~ mod_days * mod_min_day,
        TRUE ~ NA_real_
      ),
      
      lt_vig_min_wk = dplyr::case_when(
        vig_yes == 2 ~ 0,
        vig_yes == 1 & !is.na(vig_days) & !is.na(vig_min_day) ~ vig_days * vig_min_day,
        TRUE ~ NA_real_
      )
    ) %>%
    dplyr::select(
      -mod_yes, -mod_days, -mod_min_day,
      -vig_yes, -vig_days, -vig_min_day
    )
}
# convert frequency per unit to per week
frequency_per_week <- function(freq, unit) {
  dplyr::case_when(
    is.na(freq) | is.na(unit) ~ NA_real_,
    unit == "W" ~ freq,
    unit == "D" ~ freq * 7,
    unit == "M" ~ freq * 12 / 52,
    unit == "Y" ~ freq / 52,
    TRUE ~ NA_real_
  )
}

# compute weekly minutes for PAQ_L
derive_lt_minutes_2021 <- function(df) {
  df %>%
    mutate(
      mod_freq_wk = frequency_per_week(nhanes_na(PAD790Q), PAD790U),
      vig_freq_wk = frequency_per_week(nhanes_na(PAD810Q), PAD810U),
      lt_mod_min_wk = mod_freq_wk * nhanes_na(PAD800),
      lt_vig_min_wk = vig_freq_wk * nhanes_na(PAD820)
    ) %>%
    dplyr::select(-mod_freq_wk, -vig_freq_wk)
}

# derive MVPA + guideline
derive_mvpa <- function(df) {
  df %>% mutate(
    mvpa_min_wk = lt_mod_min_wk + 2 * lt_vig_min_wk,
    meets_guideline = ifelse(!is.na(mvpa_min_wk), mvpa_min_wk >= 150, NA)
  )
}

derive_pa_all <- function(df) {
  
  df %>%
    mutate(
      lt_mod_min_wk = NA_real_,
      lt_vig_min_wk = NA_real_
    ) %>%
    
    # ---- legacy PAQ (2009–2020) ----
  {
    legacy <- . %>% filter(era %in% c("2009-2016","2017-2020"))
    legacy <- derive_lt_minutes_legacy(legacy)
    
    modern <- . %>% filter(era == "2021-2023")
    modern <- derive_lt_minutes_2021(modern)
    
    bind_rows(legacy, modern)
  } %>%
    
    derive_mvpa()
}

derive_lt_minutes_2021 <- function(df) {
  df %>%
    mutate(
      # clean inputs
      mod_freq = nhanes_na(PAD790Q),
      mod_unit = PAD790U,
      mod_min_each = nhanes_na(PAD800),
      
      vig_freq = nhanes_na(PAD810Q),
      vig_unit = PAD810U,
      vig_min_each = nhanes_na(PAD820),
      
      mod_freq_wk = frequency_per_week(mod_freq, mod_unit),
      vig_freq_wk = frequency_per_week(vig_freq, vig_unit),
      
      # moderate leisure PA
      lt_mod_min_wk = case_when(
        mod_freq == 0 ~ 0,
        mod_freq > 0 & !is.na(mod_freq_wk) & !is.na(mod_min_each) ~
          mod_freq_wk * mod_min_each,
        TRUE ~ NA_real_
      ),
      
      # vigorous leisure PA
      lt_vig_min_wk = case_when(
        vig_freq == 0 ~ 0,
        vig_freq > 0 & !is.na(vig_freq_wk) & !is.na(vig_min_each) ~
          vig_freq_wk * vig_min_each,
        TRUE ~ NA_real_
      )
    ) %>%
    dplyr::select(
      -mod_freq, -mod_unit, -mod_min_each, -mod_freq_wk,
      -vig_freq, -vig_unit, -vig_min_each, -vig_freq_wk
    )
}

derive_pa_all <- function(df) {
  
  # initialize PA variables
  df <- df %>%
    mutate(
      lt_mod_min_wk = NA_real_,
      lt_vig_min_wk = NA_real_
    )
  
  # ---- legacy PAQ (2009–2020) ----
  legacy <- df %>%
    filter(era %in% c("2009-2016", "2017-2020")) %>%
    derive_lt_minutes_legacy()
  
  # ---- modern PAQ_L (2021–2023) ----
  modern <- df %>%
    filter(era == "2021-2023") %>%
    derive_lt_minutes_2021()
  
  # ---- recombine + derive MVPA ----
  bind_rows(legacy, modern) %>%
    derive_mvpa()
}

nhanes_link <- nhanes_link %>%
  mutate(
    era = case_when(
      cycle %in% c("2009-2010","2011-2012","2013-2014","2015-2016") ~ "2009-2016",
      cycle == "2017-2020" ~ "2017-2020",
      cycle == "2021-2023" ~ "2021-2023"
    )
  )

df <- derive_pa_all(nhanes_link)

df$meets_ltpa <- df$meets_guideline


df$BMI <- df %>% 
  mutate(bmi = case_when(
    BMXBMI < 18 ~ 'Underweight',
    BMXBMI >=18 & BMXBMI <= 25  ~ 'Normal',
    BMXBMI > 25 & BMXBMI <30 ~ 'Overweight',
    BMXBMI >= 30 & BMXBMI < 35 ~ 'Obese',
    BMXBMI >= 35 ~ 'Very Obese',
    T ~ NA
  )) %>% pull(bmi)

df |>
  summarise(
    n_nonmissing_ltpa = sum(!is.na(meets_ltpa)),
    pct_meets = mean(meets_ltpa, na.rm = TRUE)
  )


# df <- df %>% 
#   mutate(WTDR2D = coalesce(WTDR2D.x, WTDR2D.y)) |>
#   dplyr::select(-WTDR2D.x, -WTDR2D.y)

df <- df |>
  mutate(
    WTDR2D_unified = case_when(
      cycle == "2017-2020" ~ WTDR2DPP,
      TRUE ~ WTDR2D
    ),
    WTDR2D_pooled = WTDR2D_unified / 6
  )

# df |>
#   group_by(cycle) |>
#   summarise(
#     mean_wt = mean(WTDR2D_pooled, na.rm = TRUE),
#     n = n()
#   )
# 
# df |>
#   group_by(cycle) |>
#   summarise(
#     total_pop = sum(WTDR2D_pooled, na.rm = TRUE)
#   )


df <- df |>
  filter(
    # !is.na(plant_forward_score),
    # !is.na(meets_ltpa),
    !is.na(WTDR2D_pooled),
    !is.na(SDMVPSU),
    !is.na(SDMVSTRA)
  )



make_quantile_safe <- function(x, n = 4) {
  
  # If all NA or <2 unique values, return NA
  if (all(is.na(x)) || length(unique(x[!is.na(x)])) < 2) {
    return(rep(NA_integer_, length(x)))
  }
  
  qs <- quantile(
    x,
    probs = seq(0, 1, length.out = n + 1),
    na.rm = TRUE,
    type = 7
  )
  
  # Remove duplicated cutpoints
  qs <- unique(qs)
  
  # If we can’t form intervals, return NA
  if (length(qs) <= 2) {
    return(rep(NA_integer_, length(x)))
  }
  
  as.integer(cut(
    x,
    breaks = qs,
    include.lowest = TRUE
  ))
}


make_quantile_by_group <- function(df, var, group, n = 4, new_name = NULL) {
  
  var   <- rlang::ensym(var)
  group <- rlang::ensym(group)
  
  if (is.null(new_name)) {
    new_name <- paste0(as.character(var), "_q", n)
  }
  
  df %>%
    group_by(!!group) %>%
    mutate(
      !!new_name := make_quantile_safe(!!var, n = n)
    ) %>%
    ungroup()
}


df <- make_quantile_by_group(
  df,
  var   = INDFMPIR,
  group = cycle,
  n     = 3,
  new_name = "pir_q"
)
df$pir_q <- as.factor(df$pir_q)

df <- make_quantile_by_group(
  df,
  var   = TKCAL_2dmean ,
  group = cycle,
  n     = 3,
  new_name = "calorie_q"
)

df$calorie_q <- as.factor(df$calorie_q)


df <- df %>%
  mutate(
    diabetes = case_when(
      DIQ010 == 1 ~ 1,
      DIQ010 == 2 ~ 0,
      TRUE ~ NA_real_
    )
  )

# df <- df %>%
#   mutate(
#     diabetes_c = case_when(
#       DIQ010 == 1 ~ 1,
#       DIQ010 == 2 & !is.na(hba1c) & hba1c >= 6.5 ~ 1,
#       DIQ010 == 2 ~ 0,
#       TRUE ~ NA_real_
#     )
#   )


df <- df %>%
  mutate(
    food_insecure = case_when(
      FSDAD %in% c(3, 4) ~ 1,
      FSDAD %in% c(1, 2) ~ 0,
      TRUE ~ NA_real_
    )
  )

df <- df %>%
  mutate(
    food_security_cat = case_when(
      FSDAD == 1 ~ "High",
      FSDAD == 2 ~ "Marginal",
      FSDAD == 3 ~ "Low",
      FSDAD == 4 ~ "Very low",
      TRUE ~ NA_character_
    )
  )

df$food_security_cat <- factor(
  df$food_security_cat,
  levels = c("High", "Marginal", "Low", "Very low")
)

df$age_cat <- df %>% 
  mutate(age = case_when(
    RIDAGEYR >= 18 & RIDAGEYR <= 29 ~ '18-29',
    RIDAGEYR >=30 & RIDAGEYR <= 49 ~ '40-49',
    RIDAGEYR >= 50 & RIDAGEYR <= 65 ~ '50-65',
    RIDAGEYR > 65 ~ '65+',
    T ~ NA
  )) %>% pull(age)

df <- df %>%
  mutate(
    meets_ltpa = factor(
      meets_ltpa,
      levels = c(FALSE, TRUE),
      labels = c("Does not meet", "Meets")
    )
  )

df$race <- df %>% 
  mutate(race = case_when(
    RIDRETH1 == 1 | RIDRETH1== 2 ~ 'Mexican/Hispanic',
    RIDRETH1 == 3 ~ 'Non-Hispanic White',
    RIDRETH1 == 4 ~ 'Non-Hispanic Black',
    RIDRETH1 == 5 ~ 'Other Race/Multi',
    T ~ NA
  )) %>% pull(race)


df <- df %>%
  mutate(
    # Avoid division by zero
    sodium_density = if_else(
      TKCAL_2dmean > 0,
      TSODI_2dmean / (TKCAL_2dmean / 1000),
      NA_real_
    ),
    
    sugar_pct_kcal = if_else(
      TKCAL_2dmean > 0,
      (TSUGR_2dmean * 4) / TKCAL_2dmean * 100,
      NA_real_
    )
  )


df <- make_quantile_by_group(
  df,
  var = sodium_density,
  group = cycle,
  n = 3,
  new_name = "sodium_t"
)

df <- make_quantile_by_group(
  df,
  var = sugar_pct_kcal,
  group = cycle,
  n = 3,
  new_name = "sugar_t"
)

df <- make_quantile_by_group(
  df,
  var = TSFAT_2dmean,
  group = cycle,
  n = 3,
  new_name = "sfat_t"
)

df <- make_quantile_by_group(
  df,
  var = TFIBE_2dmean,
  group = cycle,
  n = 3,
  new_name = "fiber_t"
)

df <- df %>%
  mutate(
    risk_sodium = if_else(sodium_t == 3, 1, 0, missing = NA_real_),
    risk_sugar  = if_else(sugar_t  == 3, 1, 0, missing = NA_real_),
    risk_sfat   = if_else(sfat_t   == 3, 1, 0, missing = NA_real_),
    risk_fiber  = if_else(fiber_t  == 1, 1, 0, missing = NA_real_)
  )

df <- df %>%
  mutate(
    cd_risk_score = risk_sodium +
      risk_sugar  +
      risk_sfat   +
      risk_fiber
  )

df <- df %>%
  mutate(
    cd_risk_cat = case_when(
      cd_risk_score %in% c(0,1) ~ "Low",
      cd_risk_score == 2        ~ "Moderate",
      cd_risk_score %in% c(3,4) ~ "High",
      TRUE ~ NA_character_
    )
  )

df$cd_risk_cat <- factor(df$cd_risk_cat,
                         levels = c("Low","Moderate","High"))

library(survey)
# 
# design <- svydesign(
#   id = ~SDMVPSU,
#   strata = ~SDMVSTRA,
#   weights = ~WTDR2D_pooled,
#   nest = TRUE,
#   data = df
# )
# 
# svytable(~cd_risk_score, design)
# svyby(~I(meets_ltpa=="Meets"),
#       ~cd_risk_score,
#       design,
#       svymean,
#       na.rm=TRUE)
# 
# svyby(~I(meets_ltpa=="Meets"),
#       ~interaction(cd_risk_score, era),
#       design,
#       svymean,
#       na.rm=TRUE)
# 
# model_int <- svyglm(
#   meets_ltpa ~ cd_risk_score * era +
#     age_cat + race + pir_q +
#     food_security_cat + calorie_q,
#   design = design,
#   family = quasibinomial()
# )
# summary(model_int)
#survey ----
library(survey)
df <- df |>
  mutate(
    strata_new = interaction(cycle, SDMVSTRA, drop = TRUE),
    psu_new    = interaction(cycle, SDMVPSU, drop = TRUE)
  )

des <- svydesign(
  ids = ~psu_new,
  strata = ~strata_new,
  weights = ~WTDR2D_pooled,
  nest = TRUE,
  data = df
)

# Table ----
library(gtsummary)
tbl_summary(
  df,
  by = pir_q,
  include = c(
    meets_ltpa,
    plant_forward_tertile,
    age_cat,
    RIAGENDR,
    race,
    cycle,
    food_insecure,
    BMXBMI,
    BMI
  ),
  # type = list(
  #   plant_forward_score ~ "continuous"
  # ),
  statistic = list(
    all_continuous() ~ "{mean} ({sd})",
    all_categorical() ~ "{n} ({p}%)"
  )
) %>% 
  add_overall() %>% 
  bold_labels()

# models ----


# des <- update(des,
#               pct_prot = (TPROT_2dmean*4)/TKCAL_2dmean,
#               pct_carb = (TCARB_2dmean*4)/TKCAL_2dmean,
#               pct_fat  = (TTFAT_2dmean*9)/TKCAL_2dmean,
#               fiber_density = TFIBE_2dmean / TKCAL_2dmean * 1000
# )

summary(svyglm(
  BMXBMI ~ pct_prot + pct_fat + fiber_density  + cycle + age_cat + race,
  design = des
))


m1 <- svyglm(
  meets_ltpa ~ factor(plant_forward_tertile)*
    age_cat + RIAGENDR + RIDRETH1 +
    pir_q + SMQ020 + calorie_q,
  design = des,
  family = quasibinomial()
)

parameters::parameters(m1,exponentiate = T, ci_method = 'wald',include_info = T)

m2 <- svyglm(
  meets_ltpa ~ plant_forward_score +
    RIDAGEYR + RIAGENDR + RIDRETH1 +
    INDFMPIR + SMQ020 + TKCAL_2dmean + BMI,
  design = des,
  family = quasibinomial()
)

parameters::parameters(m2, exponentiate = TRUE)





# Graphs ----



library(ggplot2)
library(marginaleffects)
library(scales)

plot_svy_trend <- function(model, design, group_var, time_var ,
                           x_lab = "NHANES Era",
                           y_lab = "Adjusted Probability",
                           legend_title = NULL,
                           title = NULL,
                           subtitle = NULL) {
  
  # mf <- model.frame(model) 
  # mf$wt_final <- weights(model$survey.design)
  model_data <- model$survey.design$variables
  group_vars <- c(time_var, group_var)
  
  model_data <- model_data[complete.cases(model_data[, group_vars]), ]
  
  pred <- predictions(
    model,
    newdata = model_data,
    by = group_vars,
    type = "response",
    wts = "WTDR2D_pooled"
  )
  
  # group_vars <- c(time_var, group_var)
  # 
  # mf2 <- mf[complete.cases(mf[, group_vars]), ]
  # 
  # pred <- predictions(
  #   model,
  #   newdata = mf2,
  #   by = group_vars,
  #   type = "response",
  #   wts = "wt_final"
  # )
  pred <- pred[!is.na(pred[[group_var]]), ]
  
  pred[[group_var]] <- as.factor(pred[[group_var]])
  pred[[time_var]] <- as.factor(pred[[time_var]])
  
  if (is.null(legend_title)) {
    legend_title <- group_var
  }
  
  ggplot(pred, aes(
    x = .data[[time_var]],
    y = estimate,
    group = .data[[group_var]],
    color = .data[[group_var]]
  )) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.5) +
    geom_errorbar(
      aes(ymin = conf.low, ymax = conf.high),
      width = 0.12,
      linewidth = 0.6
    ) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      x = x_lab,
      y = y_lab,
      color = legend_title,
      title = title,
      subtitle = subtitle
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 14)
    )
}

m1 <- svyglm(
  meets_ltpa ~ era+factor(plant_forward_tertile)*pir_q +
    age_cat + RIAGENDR + race+ food_insecure+
     calorie_q ,
  design = des,
  family = quasibinomial()
)
m2 <- svyglm(
  meets_ltpa ~ cycle + factor(animal_forward_tertile)*
    age_cat + RIAGENDR + race+ food_insecure+
    pir_q + calorie_q,
  design = des,
  family = quasibinomial()
)
parameters::parameters(m1,exponentiate = T, ci_method='wald',include_info= T)
regTermTest(m1, ~ cycle:factor(plant_forward_tertile))
emmeans(m1,
        ~ plant_forward_tertile | cycle,
        type = "response")
library(ggeffects)
ggaverage(m1, terms = c("plant_forward_tertile", "pir_q"))

m3 <- svyglm(
  BMXBMI ~ cycle * meets_ltpa+factor(plant_forward_tertile)+
    age_cat + RIAGENDR + race+ food_insecure+
    pir_q + calorie_q,
  design = des
)


pred <- ggpredict(
  m3,
  terms = c("plant_forward_tertile", "meets_ltpa", "cycle")
)

ggplot(pred, aes(x = x, y = predicted,
                 color = group, group = group)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  facet_wrap(~facet) +
  labs(
    x = "Plant-Forward Diet Tertile",
    y = "Predicted BMI",
    color = "Meets LTPA",
    title = "Predicted BMI by LTPA and Plant-Forward Diet"
  ) +
  theme_minimal()
parameters::parameters(m3,include_info= T)

plot_svy_trend(
  model = m1,
  design = des,
  group_var = "plant_forward_tertile",
  time_var = 'pir_q',
  legend_title = "Plant-Based Tertile",
  title = "Meeting Leisure-Time Physical Activity Standards by Plant-forward diet and Poverty",
  subtitle = "Survey-weighted, adjusted for food security, education, income, age, and race"
)
library(emmeans)
emmeans(m3, ~ meets_ltpa + cycle,rg.limit =20000)
emmeans(m3, ~ meets_ltpa + cycle,
        vcov. = vcov(m3),
        ,rg.limit =20000)

emm <- emmeans(m3,
               ~ meets_ltpa + cycle,
               options = list(survey.adjust = TRUE),
               rg.limit =20000)

# -----------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------
saveRDS(nhanes_link, "nhanes_link.rds")           # fast R-native format
haven::write_xpt(nhanes_link, "nhanes_link.xpt")  # SAS-compatible if needed

# Quick check
cat("Variables in dataset:\n")
print(names(nhanes_link))

# Optional descriptive summary (equivalent to PROC MEANS)
# summary(nhanes_link[, paste0(nutrient_stems, "_2dmean")])