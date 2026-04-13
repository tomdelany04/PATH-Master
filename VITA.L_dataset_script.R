#Vital Dataset

library(haven)
library(dplyr)
library(tidyverse)
library(broom)

#Data load
vital <- read_sas("VITAL_trial_NEJM_2022.sas7bdat")
names(vital)
view(vital)
dim(vital)

#Data preparation
vital.hte <- vital %>%
  select(cvdeath,fishoilactive, ageyr, bmi, sex, htnmed, diabetes, currsmk, parhxmi) %>% 
  drop_na() %>% 
  mutate( 
    age = as.numeric(ageyr),
    bmi = as.numeric(bmi),
    death = cvdeath,
    treatment = factor(fishoilactive,
                       levels = c(0,1),
                       labels = c("Placebo", "Active Omega-3")),
    sex = factor(sex,
                 levels = c(1,2),
                 labels = c("Male", "Female")),
    htnmed = factor(htnmed,
                    levels = c(0,1),
                    labels = c("No", "Yes")),
    diabetes = factor(diabetes,
                      levels = c(0,1),
                      labels = c("No", "Yes")),
    currsmk = factor(currsmk, 
                     levels = c(0,1),
                     labels = c("Former", "Current")),
    familyHist = factor(parhxmi,
                         levels = c(0,1),
                         labels = c("Former", "Current")))
  
View(vital.hte)


# Descriptive statistics 
# crude 2x2 table
vital_tab <- table(vital.hte$death, vital.hte$treatment)
vital_tab

# Descriptive table
library(table1)
vital_tab2 <- table1(~age + bmi + sex + htnmed + currsmk + familyHist, data =vital.hte,  digits =2)
vital_tab2
          
#Fit logistic regression without the treatment
glm.baseline.risk.vital <- glm(death ~ age + bmi + sex + htnmed + currsmk + familyHist,
                  data = vital.hte,
                  family = binomial)
summary(glm.baseline.risk.vital)

# Linear predictor; predict baseline risk; Individual
vital.hte$risk <- predict(glm.baseline.risk.vital, type = "link")

# Lake Wobegon Effect
hist(vital.hte$risk, main="Baseline Risk Distribution", xlab="Linear Predictor")
# It is skewed to the left,

# Centering risk for interpretability
vital.hte$risk.center <- vital.hte$risk - mean(vital.hte$risk, na.rm = T)

# Testing Treatment interaction with risk
vital.hte.inter <- glm(death ~ treatment * risk.center, 
                        data = vital.hte, family = binomial())
summary(vital.hte.inter) # Interaction term is not significant(P = 0.249)

# Baseline stratification 
# Dividing patients into four equal sized risk groups
risk.quarter <- quantile(vital.hte$risk, probs = c(0,0.25, 0.50, 0.75,1), na.rm = T) 

# Check if new column is created for risk quarter
vital.hte <- vital.hte %>% 
  mutate(risk.quarter = ntile(risk,4)) %>% 
    mutate(risk.quarter = factor(risk.quarter,
                          levels= c(1,2,3,4),
                          labels = c("Q(low)", "Q2","Q3", "Q4(High)")))
glimpse(vital.hte)


# Absolute risk Difference(ARD)
library(marginaleffects)
avg_comparisons(vital.hte.inter, variables = "treatment", by = "risk.quarter", newdata = vital.hte)


#overall risk: with interaction
overall_risk <- avg_comparisons(
  vital.hte.inter,
  variables = "treatment",
  comparison = "lnoravg",
  transform = "exp") %>% 
  as_tibble() %>% 
  mutate(subgroup = "overall") %>% 
  select(subgroup , estimate, conf.low, conf.high)

# Treatment effect per risk quarter
risk_per_quarter <- avg_comparisons(
  vital.hte.inter,
  variables = "treatment",
  by = "risk.quarter",
  comparison = "lnoravg",
  transform = "exp") %>% 
  as_tibble() %>% 
  select(subgroup = risk.quarter, estimate, conf.low, conf.high)


# Individual variables
vital.subgroup <- vital.hte %>% 
  mutate(
    age.categ = ifelse(age >= 70, "Age >= 70", "Age < 70"),
    bmi.categ = ifelse(bmi >= 30, "Obese(BMI>=30)", "Non-obese(BMI < 30)")
  )
subgroup.variables <- c( "age.categ", "bmi.categ")

# for age
risk.by.age <- avg_comparisons(
  vital.hte.inter,
  variables = "treatment",
  newdata = subset(vital.subgroup,age.categ == "Age >= 70"),
  comparison = "lnoravg",
  transform = "exp") %>% 
  as_tibble() %>% 
  mutate(subgroup = "Age >= 70") %>% 
  select(subgroup , estimate, conf.low, conf.high)

# for BMI
risk.by.bmi <- avg_comparisons(
  vital.hte.inter,
  variables = "treatment",
  newdata = subset(vital.subgroup,bmi.categ == "Obese(BMI>=30)"),
  comparison = "lnoravg",
  transform = "exp") %>% 
  as_tibble() %>% 
  mutate(subgroup = "Obese(BMI>=30)") %>% 
  select(subgroup , estimate, conf.low, conf.high)

# Forest Plot
library(ggplot2)
library(plotly)

vital.plot.HTE <- bind_rows(overall_risk, risk_per_quarter,risk.by.age, risk.by.bmi, ) 

p <- ggplot(vital.plot.HTE, aes(x = estimate, y = subgroup)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_point(size = 3, color = "darkred") +
  scale_x_log10() + 
  labs(title = "Forest Plot of Treatment Heterogeneity (VITAL)",
       x = "Odds Ratio (95% CI)",
       y = "Risk") +
  theme_minimal()

ggplotly(p)
