#Vital Dataset

library(haven)
library(dplyr)
library(tidyverse)
library(broom)
library(marginaleffects)
library(table1)
library(ggplot2)
library(plotly)
library(survival)
library(survminer)



#Data load
vital <- read_sas("VITAL_trial_NEJM_2022.sas7bdat")
names(vital)
# view(vital)
dim(vital)

# Data preparation
vital.hte <- vital %>%
  select(randyrs,cvdeath,vitdactive,fishoilactive, ageyr, bmi, sex, htnmed, diabetes, currsmk, parhxmi) %>% 
  drop_na() %>% 
  mutate( 
    time = as.numeric(randyrs),
    age = as.numeric(ageyr),
    bmi = as.numeric(bmi),
    death = cvdeath,
    vitd = factor(vitdactive,
                  levels = c(0,1),
                  labels = c("Placebo", "Active Vit-D")),
    fishoil = factor(fishoilactive,
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
vital_tab <- table(vital.hte$death, vital.hte$fishoilonly)
vital_tab

# Descriptive table
vital_tab2 <- table1(~age + bmi + sex + htnmed + currsmk + diabetes + familyHist, data =vital.hte,  digits =2)
vital_tab2

# Event rates by treatment group
vital.hte %>%
  group_by(vitd, fishoil) %>%
  summarise(
    n = n(),
    events = sum(death),
    event_rate = mean(death),
    .groups = "drop"
  ) %>% 
  print()

# Proportional hazards assumption 
ph.test <- cox.zph(cox.vital.hte)
print(ph.test)

# create age and BMI strata
vital.hte <- vital.hte %>%
  mutate(
    age.categ = ifelse(age >= 70, "Age >= 70", "Age < 70"),
    bmi.categ = ifelse(bmi >= 30, "Obese(BMI>=30)", "Non-obese")
  )

# Unadjusted factorial analysis with cox proportional hazards model
cox.vital.hte <- coxph(Surv(time,death) ~ vitd * fishoil, data = vital.hte)
summary(cox.vital.hte)
tidy(cox.vital.hte, conf.int = T, exponentiate=T)

# Adjusted factorial analysis with cox proportional hazards model
cox.vitd.fact <- coxph(Surv(time,death) ~ vitd + strata(age) + strata(bmi) + sex + htnmed + diabetes + currsmk + familyHist, data = vital.hte)
summary(cox.vitd.fact)
tidy(cox.vitd.fact, conf.int = T, exponentiate=T)

# Vitamin D and Omega-3 group compared to Double Placebo
vital.both1 <- vital.hte %>% 
  filter(
    (vitd == "Placebo" & fishoil == "Placebo")|
      (vitd == "Active Vit-D" & fishoil == "Active Omega-3")
  ) %>% 
  mutate(
    treatment1 = factor(
      ifelse(vitd == "Active Vit-D" & fishoil == "Active Omega-3", "vitd_with_omega3", "Double_placebo"),
      levels = c("Double_placebo", "vitd_with_omega3")
    )
  )

table(vital.both1$treatment1)

# Fit cox model for vitD and Omega-3 group compared to Double Placebo group
cox.both1 <- coxph(Surv(time,death) ~ treatment1 + strata(age) + strata(bmi) + sex + htnmed + diabetes + currsmk + familyHist,data = vital.both1)
summary(cox.both1)
tidy(cox.both1, conf.int = T, exponentiate=T)


# Vitamin D and Omega-3 group compared to Omega-3 and Placebo group
vital.both2 <- vital.hte %>% 
  filter(
    (vitd == "Placebo" & fishoil == "Active Omega-3")|
      (vitd == "Active Vit-D" & fishoil == "Active Omega-3")
  ) %>% 
  mutate(
    treatment2 = factor(
      ifelse(vitd == "Active Vit-D" & fishoil == "Active Omega-3", "vitd_with_omega3", "Placebo_vitd_active_omega3"),
      levels = c("Placebo_vitd_active_omega3", "vitd_with_omega3")
    )
  )

table(vital.both2$treatment2)

# Fit cox model for vitD and Omega-3 group compared to Omega-3 and Placebo group
cox.both2 <- coxph(Surv(time,death) ~ treatment2 + strata(age) + strata(bmi) + sex + htnmed + diabetes + currsmk + familyHist, data = vital.both2)
summary(cox.both2)
tidy(cox.both2, conf.int = T, exponentiate=T)




#****PATH approach********

# Fit a Cox proportional hazards model with out treatment 
# doi: 10.7326/M18-3668
cox.vital.baseline <- coxph(Surv(time,death) ~ strata(age) + strata(bmi)+ sex + htnmed + currsmk + diabetes + familyHist, 
                            data = (vital.hte))
summary(cox.vital.baseline)

#linear predictor; predict baseline risk; Individual with cox model
vital.hte$lp <- predict(cox.vital.baseline, newdata= vital.hte, type = "lp")

# # Centering risk for interpretability
vital.hte$risk.center.cox <- vital.hte$lp - mean(vital.hte$lp, na.rm = T)

# Baseline risk stratification: Dividing patients into four equal sized risk groups 
risk.quarter <- quantile(vital.hte$risk.center.cox, probs = c(0,0.25, 0.50, 0.75,1), na.rm = T) 
vital.hte <- vital.hte %>%
  mutate(risk = ntile(risk.center.cox,4),
         risk = factor(risk,
                       levels= c(1,2,3,4),
                       labels = c("Q(low)", "Q2","Q3", "Q4(High)")))

# fit Cox model with interaction between vitamin D and risk quarter
cox.vital.vitd <- coxph(Surv(time,death) ~ vitd * risk, data = vital.hte)
summary(cox.vital.vitd)


# fit cox model with interaction between Omega-3 and risk quarter
cox.vital.omega3 <- coxph(Surv(time,death) ~ fishoil * risk, data = vital.hte)
summary(cox.vital.omega3)

# fit cox model with interaction between vitamin D, Omega-3 and risk quarter
cox.vital.interaction <- coxph(Surv(time, death)~ vitd* fishoil * risk, data = vital.hte)
# summary(cox.interaction.hte)


library(marginaleffects)
library(ggplot2)
library(plotly)

# Overall Hazard lnratioavg
overall.risk <- avg_comparisons(
  cox.both1.risk,
  variables = "treatment1",
  comparison = "lnratioavg",   
  transform = "exp"
) %>% 
  as_tibble() %>% 
  mutate(subgroup = "Overall") %>% 
  select(subgroup, estimate, conf.low, conf.high)

# Hazard lnratioavg by Risk Quarter
hr.by.risk <- avg_comparisons(
  cox.both1.risk,
  variables = "treatment1",
  by = "risk",
  comparison = "lnratioavg",   
  transform = "exp"
) %>% 
  as_tibble() %>% 
  select(subgroup = risk, estimate, conf.low, conf.high)

#Subgroup Hazard lnratioavg by Age and BMI
# Age >= 70
hr.by.age <- avg_comparisons(
  cox.both1.risk,
  variables = "treatment1",
  newdata = subset(vital.subgroup, age.categ == "Age >= 70"),
  comparison = "lnratioavg",   
  transform = "exp"
) %>% 
  as_tibble() %>% 
  mutate(subgroup = "Age >= 70") %>% 
  select(subgroup, estimate, conf.low, conf.high)

# Age < 70
hr.by.age2 <- avg_comparisons(
  cox.both1.risk,
  variables = "treatment1",
  newdata = subset(vital.subgroup, age.categ == "Age < 70"),
  comparison = "lnratioavg",   
  transform = "exp"
) %>% 
  as_tibble() %>% 
  mutate(subgroup = "Age < 70") %>% 
  select(subgroup, estimate, conf.low, conf.high)

# Bmi >= 30 : Obese 
hr.by.bmi <- avg_comparisons(
  cox.both1.risk,
  variables = "treatment1",
  newdata = subset(vital.subgroup, bmi.categ == "Obese(BMI>=30)"),
  comparison = "lnratioavg",   
  transform = "exp"
) %>% 
  as_tibble() %>% 
  mutate(subgroup = "Obese(BMI>=30)") %>% 
  select(subgroup, estimate, conf.low, conf.high)

# BMI < 30
hr.by.bmi2 <- avg_comparisons(
  cox.both1.risk,
  variables = "treatment1",
  newdata = subset(vital.subgroup, age == "Non-obese"),
  comparison = "lnratioavg",   
  transform = "exp"
) %>% 
  as_tibble() %>% 
  mutate(subgroup = "Non-obese") %>% 
  select(subgroup, estimate, conf.low, conf.high)





# Forest plot

vital.plot.HTE <- bind_rows(overall.risk,
                            hr.by.risk, 
                            hr.by.age, 
                            hr.by.age2, 
                            hr.by.bmi,
                            hr.by.bmi2) %>% 
  mutate(subgroup = factor(subgroup, 
                           levels = c("Overall", "Q(low)", "Q2", "Q3", "Q4(High)", "Age >= 70", "Age < 70", "Obese(BMI>=30)", "Non-obese")))
print(vital.plot.HTE)

p1 <- ggplot(vital.plot.HTE,
             aes(x = estimate, y = subgroup)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, color = "grey30") +
  geom_point(size = 3, color = "maroon") +
  scale_x_log10() + 
  labs(title = "Forest Plot of Treatment Heterogeneity ",
       x = "Hazard ratio (95% CI)",
       y = "Risk") +
  theme_minimal() 

ggplotly(p1)