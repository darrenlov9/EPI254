library(marginaleffects)
library(gtsummary)
library(survey)
library(tidyverse)

df <- data.table::fread('nhanes_data.csv')

#survey ----
library(survey)
df <- df |>
  mutate(
    strata_new = interaction(cycle, SDMVSTRA, drop = TRUE),
    psu_new    = interaction(cycle, SDMVPSU, drop = TRUE)
  )

df$meets_ltpa <- as.factor(df$meets_ltpa)

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
  include = c(
    meets_ltpa,
    plant_forward_tertile,
    pir_q,
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
  # add_overall() %>% 
  bold_labels()

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

# Visuals ----

# Weighted proportions
pf_ltpa <- svyby(
  ~I(meets_ltpa == "Meets"),
  ~plant_forward_tertile,
  des,
  svymean,
  na.rm = TRUE
)

pf_ltpa$percent <- pf_ltpa$`I(meets_ltpa == "Meets")TRUE` * 100

ggplot(pf_ltpa, aes(x = factor(plant_forward_tertile), y = percent)) +
  geom_col() +
  labs(
    x = "Plant-Forward Diet Tertile",
    y = "% Meeting LTPA Guidelines",
    title = "Weighted % Meeting LTPA by Plant-Forward Tertile"
  ) +
  theme_minimal()

pf_stack <- svytable(~plant_forward_tertile + meets_ltpa, des)

pf_prop <- prop.table(pf_stack, 1)  # row percentages
pf_df <- as.data.frame(pf_prop)

ggplot(pf_df, aes(x = factor(plant_forward_tertile), 
                  y = Freq,
                  fill = meets_ltpa)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Plant-Forward Diet Tertile",
    y = "Weighted %",
    fill = "LTPA Status"
  ) +
  theme_minimal()

pir_stack <- svytable(~pir_q + meets_ltpa, des)
pir_prop <- prop.table(pir_stack, 1)
pir_df <- as.data.frame(pir_prop)

ggplot(pir_df, aes(x = factor(pir_q), 
                   y = Freq,
                   fill = meets_ltpa)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Income-to-Poverty Ratio Tertile",
    y = "Weighted %",
    fill = "LTPA Status"
  ) +
  theme_minimal()


# Get % meeting LTPA for each combo
heat <- svyby(
  ~I(meets_ltpa == "Meets"),
  ~plant_forward_tertile + pir_q,
  des,
  svymean,
  na.rm = TRUE
)

heat$percent <- heat$`I(meets_ltpa == "Meets")TRUE` * 100

ggplot(heat, aes(
  x = factor(pir_q),
  y = factor(plant_forward_tertile),
  fill = percent
)) +
  geom_tile() +
  geom_text(aes(label = round(percent, 1))) +
  scale_fill_gradient(low = "white", high = "darkgreen") +
  labs(
    x = "PIR Tertile",
    y = "Plant-Forward Tertile",
    fill = "% Meets LTPA",
    title = "Weighted % Meeting LTPA by Diet and Income"
  ) +
  theme_minimal()

heat <- svyby(
  ~I(meets_ltpa == "Meets"),
  ~plant_forward_tertile + pir_q,
  des,
  svymean,
  vartype = "ci",
  na.rm = TRUE
)



# Weighted % meeting LTPA by plant tertile
pf_trend <- svyby(
  ~I(meets_ltpa == "Meets"),
  ~plant_forward_tertile,
  des,
  svymean,
  vartype = "ci",
  na.rm = TRUE
)

pf_trend <- pf_trend %>%
  mutate(
    percent = `I(meets_ltpa == "Meets")TRUE` * 100,
    ci_low = `ci_l.I(meets_ltpa == "Meets")TRUE` * 100,
    ci_high = `ci_u.I(meets_ltpa == "Meets")TRUE` * 100
  )

ggplot(pf_trend, aes(x = as.numeric(plant_forward_tertile),
                     y = percent)) +
  geom_line() +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.1) +
  scale_x_continuous(breaks = 1:3) +
  labs(
    x = "Plant-Forward Diet Tertile",
    y = "% Meeting LTPA (Weighted)",
    title = "Trend of LTPA Across Plant-Forward Diet Tertiles"
  ) +
  theme_minimal()


pir_trend <- svyby(
  ~I(meets_ltpa == "Meets"),
  ~pir_q,
  des,
  svymean,
  vartype = "ci",
  na.rm = TRUE
)

pir_trend <- pir_trend %>%
  mutate(
    percent = `I(meets_ltpa == "Meets")TRUE` * 100,
    ci_low = `ci_l.I(meets_ltpa == "Meets")TRUE` * 100,
    ci_high = `ci_u.I(meets_ltpa == "Meets")TRUE` * 100
  )

ggplot(pir_trend, aes(x = as.numeric(pir_q),
                      y = percent)) +
  geom_line() +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.1) +
  scale_x_continuous(breaks = 1:3) +
  labs(
    x = "Income-to-Poverty Ratio Tertile",
    y = "% Meeting LTPA (Weighted)",
    title = "Trend of LTPA Across Income Tertiles"
  ) +
  theme_minimal()


pf_trend$exposure <- "Plant-Forward"
pf_trend$tertile <- as.numeric(pf_trend$plant_forward_tertile)

pir_trend$exposure <- "PIR"
pir_trend$tertile <- as.numeric(pir_trend$pir_q)

combined <- bind_rows(
  pf_trend %>% select(tertile, percent, ci_low, ci_high, exposure),
  pir_trend %>% select(tertile, percent, ci_low, ci_high, exposure)
)

ggplot(combined, aes(x = tertile,
                     y = percent,
                     group = exposure,
                     color = exposure)) +
  geom_line() +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.1) +
  scale_x_continuous(breaks = 1:3) +
  labs(
    x = "Tertile Level",
    y = "% Meeting LTPA (Weighted)",
    color = "Exposure",
    title = "Trends in LTPA Across Diet and Income Tertiles"
  ) +
  theme_minimal()



heat_trend <- svyby(
  ~I(meets_ltpa == "Meets"),
  ~plant_forward_tertile + pir_q,
  des,
  svymean,
  na.rm = TRUE
)

heat_trend$percent <- heat_trend$`I(meets_ltpa == "Meets")TRUE` * 100

ggplot(heat_trend,
       aes(x = as.numeric(plant_forward_tertile),
           y = percent,
           group = pir_q,
           color = factor(pir_q))) +
  geom_line() +
  geom_point() +
  labs(
    x = "Plant-Forward Tertile",
    y = "% Meeting LTPA (Weighted)",
    color = "PIR Tertile",
    title = "Plant-Forward Diet and LTPA Stratified by Income"
  ) +
  theme_minimal()



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