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

#Data preparation
vital.hte <- vital %>%
  select(randyrs, cvdeath,fishoilactive, ageyr, bmi, sex, htnmed, diabetes, currsmk, parhxmi) %>% 
  drop_na() %>% 
  mutate( 
    time = as.numeric(randyrs),
    age = as.numeric(ageyr),
    bmi = as.numeric(bmi),
    death = cvdeath,
    fishoilonly = factor(fishoilactive,
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
# glimpse(vital.hte)


# Absolute risk Difference(ARD)
avg_comparisons(vital.hte.inter, variables = "fishoilonly", by = "risk.quarter", newdata = vital.hte)


#overall risk: with interaction
overall_risk <- avg_comparisons(
  vital.hte.inter,
  variables = "fishoilonly",
  comparison = "lnoravg",
  transform = "exp") %>% 
  as_tibble() %>% 
  mutate(subgroup = "overall") %>% 
  select(subgroup , estimate, conf.low, conf.high)

# Treatment effect per risk quarter
risk_per_quarter <- avg_comparisons(
  vital.hte.inter,
  variables = "fishoilonly",
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
  variables = "fishoilonly",
  newdata = subset(vital.subgroup,age.categ == "Age >= 70"),
  comparison = "lnoravg",
  transform = "exp") %>% 
  as_tibble() %>% 
  mutate(subgroup = "Age >= 70") %>% 
  select(subgroup , estimate, conf.low, conf.high)

# for BMI
risk.by.bmi <- avg_comparisons(
  vital.hte.inter,
  variables = "fishoilonly",
  newdata = subset(vital.subgroup,bmi.categ == "Obese(BMI>=30)"),
  comparison = "lnoravg",
  transform = "exp") %>% 
  as_tibble() %>% 
  mutate(subgroup = "Obese(BMI>=30)") %>% 
  select(subgroup , estimate, conf.low, conf.high)

# Forest Plot
vital.plot.HTE <- bind_rows(overall_risk, risk_per_quarter,risk.by.age, risk.by.bmi ) 

p <- ggplot(vital.plot.HTE, aes(x = estimate, y = subgroup)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_point(size = 3, color = "darkred") +
  scale_x_log10() + 
  labs(title = "Forest Plot of Treatment Heterogeneity (Omega-3)",
       x = "Odds Ratio (95% CI)",
       y = "Risk") +
  theme_minimal()

ggplotly(p)

#Time to event for vitD only
# Events by treatment(vitD only)
vital.hte %>% 
  group_by(fishoilonly) %>% 
  summarise(
    total_n = n(),
    Deaths = sum(death == 1),
    Censored = sum(death == 0),
    Death_rate = mean(death==1)
  )
#Kaplan-Meier not accounting for risk quarter
surv.vital <- survfit(Surv(time,death)~ fishoilonly, data = vital.hte)

# Survival distribution
ggsurvplot(surv.vital,
           data = vital.hte,
           legend.title="vit-D_only",
           risk.table = TRUE, 
           break.x.by = 1, 
           surv.median.line = "hv",
           conf.int = TRUE, pval = T,
           ylab = "Probability \n no death", xlab = "Time (year)")

# Significant test 
survdiff(Surv(time,death)~ fishoilonly, data = vital.hte)

# Cox proportional hazard 
cox.vital <- coxph(Surv(time,death) ~ fishoilonly, data = vital.hte)
cox.results <- tidy(cox.vital, conf.int = T, exponentiate=T)

#################*Adjusting for Vitamin D treatment*##################
# Data preparation
vital.hte2 <- vital %>%
  select(randyrs,cvdeath,vitdactive, ageyr, bmi, sex, htnmed, diabetes, currsmk, parhxmi) %>% 
  drop_na() %>% 
  mutate( 
    time = as.numeric(randyrs),
    age = as.numeric(ageyr),
    bmi = as.numeric(bmi),
    death = cvdeath,
    vitdonly = factor(vitdactive,
                      levels = c(0,1),
                      labels = c("Placebo", "Active Vit-D")),
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
head(vital.hte2)

# Descriptive statistics 
# crude 2x2 table
vital_tab <- table(vital.hte2$death, vital.hte2$vitdonly)
vital_tab

# Descriptive table
vital_tab2 <- table1(~age + bmi + sex + htnmed + currsmk + familyHist, data =vital.hte2,  digits =2)
vital_tab2

#Fit logistic regression without the treatment
glm.baseline.risk.vital2 <- glm(death ~ age + bmi + sex + htnmed + currsmk + familyHist,
                                data = vital.hte2,
                                family = binomial)
summary(glm.baseline.risk.vital2)

# Linear predictor; predict baseline risk; Individual
vital.hte2$risk2 <- predict(glm.baseline.risk.vital2, type = "link")

# Lake Wobegon Effect
hist(vital.hte2$risk2, main="Baseline Risk Distribution", xlab="Linear Predictor")
# It is skewed to the left,


# Centering risk for interpretability
vital.hte2$risk.center2 <- vital.hte2$risk2 - mean(vital.hte2$risk2, na.rm = T)

# Testing Treatment interaction with risk
vital.hte.inter2 <- glm(death ~ vitdonly * risk.center2, 
                        data = vital.hte2, family = binomial())
summary(vital.hte.inter2) # Interaction term is not significant(P = 0.249)

# Baseline stratification 
# Dividing patients into four equal sized risk groups
risk.quarter2 <- quantile(vital.hte2$risk2, probs = c(0,0.25, 0.50, 0.75,1), na.rm = T) 

# Check if new column is created for risk quarter
vital.hte2 <- vital.hte2 %>%
  mutate(risk.quarter = ntile(risk2,4),
         risk.quarter = factor(risk.quarter,
                               levels= c(1,2,3,4),
                               labels = c("Q(low)", "Q2","Q3", "Q4(High)")))
# glimpse(vital.hte2)


# Absolute risk Difference(ARD)

avg_comparisons(vital.hte.inter2, variables = "vitdonly", by = "risk.quarter", newdata = vital.hte2)


#overall risk: with interaction
overall_risk2 <- avg_comparisons(
  vital.hte.inter2,
  variables = "vitdonly",
  comparison = "lnoravg",
  transform = "exp") %>% 
  as_tibble() %>% 
  mutate(subgroup = "overall") %>% 
  select(subgroup , estimate, conf.low, conf.high)

# Treatment effect per risk quarter
risk_per_quarter2 <- avg_comparisons(
  vital.hte.inter2,
  variables = "vitdonly",
  by = "risk.quarter",
  newdata = vital.hte2,
  comparison = "lnoravg",
  transform = "exp") %>% 
  as_tibble() %>% 
  select(subgroup = risk.quarter, estimate, conf.low, conf.high)


# Individual variables
vital.subgroup2 <- vital.hte2 %>% 
  mutate(
    age.categ = ifelse(age >= 70, "Age >= 70", "Age < 70"),
    bmi.categ = ifelse(bmi >= 30, "Obese(BMI>=30)", "Non-obese(BMI < 30)")
  )
subgroup.variables <- c( "age.categ", "bmi.categ")

# for age
risk.by.age <- avg_comparisons(
  vital.hte.inter2,
  variables = "vitdonly",
  newdata = subset(vital.subgroup2,age.categ == "Age >= 70"),
  comparison = "lnoravg",
  transform = "exp") %>% 
  as_tibble() %>% 
  mutate(subgroup = "Age >= 70") %>% 
  select(subgroup , estimate, conf.low, conf.high)

# for BMI
risk.by.bmi <- avg_comparisons(
  vital.hte.inter2,
  variables = "vitdonly",
  newdata = subset(vital.subgroup2,bmi.categ == "Obese(BMI>=30)"),
  comparison = "lnoravg",
  transform = "exp") %>% 
  as_tibble() %>% 
  mutate(subgroup = "Obese(BMI>=30)") %>% 
  select(subgroup , estimate, conf.low, conf.high)

# Forest Plot
library(ggplot2)
library(plotly)

vital.plot.HTE2 <- bind_rows(overall_risk2, risk_per_quarter2,risk.by.age, risk.by.bmi) 

p2 <- ggplot(vital.plot.HTE2, aes(x = estimate, y = subgroup)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_point(size = 3, color = "darkred") +
  scale_x_log10() + 
  labs(title = "Forest Plot of Treatment Heterogeneity (Vitamin D)",
       x = "Odds Ratio (95% CI)",
       y = "Risk") +
  theme_minimal()

ggplotly(p2)

#Time to event for vitD only
# Events by treatment(vitD only)
vital.hte2 %>% 
  group_by(vitdonly) %>% 
  summarise(
    total_n = n(),
    Deaths = sum(death == 1),
    Censored = sum(death == 0),
    Death_rate = mean(death==1)
  )
#Kaplan-Meier not accounting for risk quarter
surv.vital <- survfit(Surv(time,death)~ vitdonly, data = vital.hte2)

# Survival distribution
ggsurvplot(surv.vital,
           data = vital.hte2,
           legend.title="vit-D_only",
           risk.table = TRUE, 
           break.x.by = 1, 
           pval = TRUE,
           censor.size = 1, 
           ylab = "Probability \n no death", xlab = "Time (year)")

# Survival distribution
ggsurvplot(surv.vital,
           data = vital.hte2,
           legend.title="vit-D_only",
           risk.table = TRUE, 
           break.x.by = 1, 
           pval = TRUE,
           censor.size = 1,
           ylim = c(0.985, 1.00), 
           ylab = "Probability \n no death", xlab = "Time (year)")

# Significant test 
survdiff(Surv(time,death)~ vitdonly, data = vital.hte2)

# Cox proportional hazard 
cox.vital <- coxph(Surv(time,death) ~ vitdonly, data = vital.hte2)
cox.results <- tidy(cox.vital, conf.int = T, exponentiate=T)



#################*Adjusting for both treatments(Vitamin-D and Omega-3 fish oil)*#################

vital.hte.both <- vital %>%
  select(randyrs,cvdeath,vitdactive, fishoilactive, ageyr, bmi, sex, htnmed, diabetes, currsmk, parhxmi) %>% 
  drop_na() %>% 
  mutate( 
    time = as.numeric(randyrs),
    age = as.numeric(ageyr),
    bmi = as.numeric(bmi),
    death = cvdeath,
    vitdandomega3= ifelse(vitdactive == 1 & fishoilactive == 1, "Active_vitd_omega3", "Double_placebo"),
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
                        labels = c("Former", "Current"))) %>% 
  select(-c(vitdactive, fishoilactive))
View(vital.hte.both)

#Fit logistic regression without the treatment
glm.baseline.risk.vital3 <- glm(death ~ age + bmi + sex + htnmed + currsmk + familyHist + diabetes,
                                data = vital.hte.both,
                                family = binomial)
summary(glm.baseline.risk.vital3)

# Linear predictor; predict baseline risk; Individual
vital.hte.both$risk3 <- predict(glm.baseline.risk.vital3, type = "link")

# Centering risk for interpretability
vital.hte.both$risk.center3 <- vital.hte.both$risk3 - mean(vital.hte.both$risk3, na.rm = T)

# Check if new column is created for risk quarter
vital.hte.both <- vital.hte.both %>%
  mutate(risk.quarter = ntile(risk3,4),
         risk.quarter = factor(risk.quarter,
                               levels= c(1,2,3,4),
                               labels = c("Q(low)", "Q2","Q3", "Q4(High)")))
# glimpse(vital.hte.both)

# Testing Treatment interaction with risk
vital.hte.inter3 <- glm(death ~ vitdandomega3 * risk.center3, 
                        data = vital.hte.both, family = binomial())
summary(vital.hte.inter3) # Interaction term is not significant(P = 0.249)

# Baseline stratification 
# Dividing patients into four equal sized risk groups
risk.quarter3 <- quantile(vital.hte.both$risk3, probs = c(0,0.25, 0.50, 0.75,1), na.rm = T) 

# Check if new column is created for risk quarter
vital.hte.both <- vital.hte.both %>%
  mutate(risk.quarter = ntile(risk3,4),
         risk.quarter = factor(risk.quarter,
                               levels= c(1,2,3,4),
                               labels = c("Q(low)", "Q2","Q3", "Q4(High)")))
# glimpse(vital.hte.both)


# Absolute risk Difference(ARD)
library(marginaleffects)
avg_comparisons(vital.hte.inter3, variables = "vitdandomega3", by = "risk.quarter", newdata = vital.hte.both)


#overall risk: with interaction
overall_risk3 <- avg_comparisons( 
  vital.hte.inter3,
  variables = "vitdandomega3",
  comparison = "lnoravg",
  transform = "exp") %>% 
  as_tibble() %>% 
  mutate(subgroup = "overall") %>% 
  select(subgroup , estimate, conf.low, conf.high)

# Treatment effect per risk quarter
risk_per_quarter3 <- avg_comparisons(
  vital.hte.inter3,
  variables = "vitdandomega3",
  by = "risk.quarter",
  newdata = vital.hte.both,
  comparison = "lnoravg",
  transform = "exp") %>% 
  as_tibble() %>% 
  select(subgroup = risk.quarter, estimate, conf.low, conf.high)


# Individual variables
vital.subgroup3 <- vital.hte.both %>% 
  mutate(
    age.categ = ifelse(age >= 70, "Age >= 70", "Age < 70"),
    bmi.categ = ifelse(bmi >= 30, "Obese(BMI>=30)", "Non-obese(BMI < 30)")
  )
subgroup.variables <- c( "age.categ", "bmi.categ")

# for age
risk.by.age <- avg_comparisons(
  vital.hte.inter3,
  variables = "vitdandomega3",
  newdata = subset(vital.subgroup3,age.categ == "Age >= 70"),
  comparison = "lnoravg",
  transform = "exp") %>% 
  as_tibble() %>% 
  mutate(subgroup = "Age >= 70") %>% 
  select(subgroup , estimate, conf.low, conf.high)

# for BMI
risk.by.bmi <- avg_comparisons(
  vital.hte.inter3,
  variables = "vitdandomega3",
  newdata = subset(vital.subgroup3,bmi.categ == "Obese(BMI>=30)"),
  comparison = "lnoravg",
  transform = "exp") %>% 
  as_tibble() %>% 
  mutate(subgroup = "Obese(BMI>=30)") %>% 
  select(subgroup , estimate, conf.low, conf.high)

# Forest Plot
library(ggplot2)
library(plotly)

vital.plot.HTE3 <- bind_rows(overall_risk3, risk_per_quarter3,risk.by.age, risk.by.bmi) 

p3 <- ggplot(vital.plot.HTE3, aes(x = estimate, y = subgroup)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_point(size = 3, color = "darkred") +
  scale_x_log10() + 
  labs(title = "Forest Plot of Treatment Heterogeneity (Vitamin-D and Omega-3)",
       x = "Odds Ratio (95% CI)",
       y = "Risk") +
  theme_minimal()

ggplotly(p3)

#Time to event for vitD Plus Omega3
# Events by treatment
vital.hte.both %>% 
  group_by(vitdandomega3) %>% 
  summarise(
    total_n = n(),
    Deaths = sum(death == 1),
    Censored = sum(death == 0),
    Death_rate = mean(death==1)
  )
#Kaplan-Meier not accounting for risk quarter
surv.vital <- survfit(Surv(time,death)~ vitdandomega3, data = vital.hte.both)

# Survival distribution
ggsurvplot(surv.vital,
           data = vital.hte.both,
           legend.title="vitdandomega3",
           risk.table = TRUE, 
           break.x.by = 1, 
           surv.median.line = "hv",
           conf.int = TRUE, pval = T,
           ylab = "Probability \n of no death", xlab = "Time (year)")

# Significant test 
survdiff(Surv(time,death)~ vitdandomega3, data = vital.hte.both)

# Cox proportional hazard 
cox.vital <- coxph(Surv(time,death) ~ vitdandomega3, data = vital.hte.both)
cox.results <- tidy(cox.vital, conf.int = T, exponentiate=T)

ggplot(cox.results,
       aes(x = term, y =estimate, ymin = conf.low, ymax = conf.high))+
  geom_pointrange()+
  geom_hline(yintercept = 1, linetype = "dashed")+
  coord_flip()+
  labs(
    x = "",
    y = "Hazard Ratio",
    title = "Cox proportional hazards model"
  )+
  theme_minimal()

# # Coxproportional hazard_ Interaction
cox.vital2 <- coxph(Surv(time,death) ~ vitdandomega3 * risk.quarter, data = vital.hte.both)
cox.results2 <- tidy(cox.vital, conf.int = T, exponentiate=T)

ggplot(cox.results2,
       aes(x = term, y =estimate, ymin = conf.low, ymax = conf.high))+
  geom_pointrange()+
  geom_hline(yintercept = 1, linetype = "dashed")+
  coord_flip()+
  labs(
    x = "",
    y = "Hazard Ratio",
    title = "Cox proportional hazards model"
  )+
  theme_minimal()






##### Random Forest
vital.hte.rf <- vital.hte$dea %>%
  select(-time)

library(randomForest)

## reproducibility 
set.seed(1) 
n <- nrow(vital.hte)
## 
train <- runif(n) < 0.7 # train the model on 70% random sample of the row numbers 
test <- (!train)

train <- vital.hte[train,]
test <- vital.hte[test, ]

## create x and matrices for glmnet (this could be done later)
x.train <- model.matrix(death ~ ., data = train)[, -1]
x.test <- model.matrix(death ~ ., data = test)[, -1]
y.train <- train$death
y.test <- test$death

set.seed(1)
oob_values <- vector("double", length = (ncol(train) - 1))
for(i in 1:(ncol(train)-1)) {
  rf0 <- randomForest(death ~ ., data = train, mtry = i, ntree = 500)
  oob_values[i] <- rf0$err.rate[rf0$ntree, "OOB"] # the OOB error rate
}

##
best.mtry <- which.min(oob_values)
best.mtry


# Fitting Random Forest with optimized mtry
rf.vital <- randomForest(death ~ vitdand, data = train, 
                         mtry = best.mtry, importance = TRUE)

# Check variable importance
importance(rf.vital)

varImpPlot(rf.vital, main = "Variable Importance")


# Random forest Performance 
rf.mal.pred <- predict(rf.vital, newdata = test, type = "prob")[ , "yes"]

