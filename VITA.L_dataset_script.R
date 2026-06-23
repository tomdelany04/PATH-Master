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
library(Hmisc)
library(splines)
library(scales)



#Data load
vital <- read_sas("VITAL_trial_NEJM_2022.sas7bdat")
names(vital)
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
                     labels = c("No", "Yes")),
    familyHist = factor(parhxmi,
                        levels = c(0,1),
                        labels = c("No", "Yes")))


# Descriptive statistics
table1(~ as.factor(death) + age + bmi + sex + htnmed + diabetes + currsmk + familyHist| vitd, data =vital.hte,  digits =2)

# Descriptive statistics by treatment group
table(vital.hte$vitd, vital.hte$fishoil)

#***kaplan Meier curve***#
# kaplan Meier curve for vitd
vital.kp <- vital.hte %>% 
  select(time, death, vitd, fishoil) %>% 
  drop_na()

#Kaplan Meier for vit-D
kp_vitd <- survfit(Surv(time, death)~ vitd, data = vital.kp)

ggsurvplot(
  kp_vitd, 
  data = vital.kp,
  legend.title= "Vitamin D", 
  legend.labs = c("Placebo", "Active Vit-D"),
  risk.table = TRUE, 
  break.x.by = 1,
  xlim = c(0,6),
  ylim = c(0.986, 1.00),
  conf.int = FALSE, 
  pval = TRUE,
  ylab = "Probability \n no CVD death", xlab = "Time (Years)")

#Kaplan Meier for Omega-3
kp_omega3 <- survfit(Surv(time, death)~ fishoil, data = vital.kp)

ggsurvplot(
  kp_vitd, 
  data = vital.kp,
  legend.title= "Omega-3", 
  legend.labs = c("Placebo", "Active Omega-3"),
  risk.table = TRUE, 
  break.x.by = 1,
  xlim = c(0,6),
  ylim = c(0.986, 1.00),
  pval = TRUE,
  ylab = "Probability \n no CVD death", xlab = "Time (Years)")





#***Factorial analysis: overall treatment effect**#  
# Coxph model for overall treatment effect 
cox.factorial <- coxph(Surv(time,death) ~ vitd + fishoil + age + bmi + sex + htnmed + diabetes + currsmk + familyHist ,data = vital.hte)

# Fit a Cox proportional hazards model with-out treatment 
cox.vital.baseline <- coxph(Surv(time,death) ~ age + bmi + sex + htnmed + currsmk + diabetes + familyHist, data = vital.hte)

#linear predictor; predict baseline individual risk 
vital.hte$lp <- predict(cox.vital.baseline, newdata= vital.hte, type = "lp")

#Centering risk for interpretability
vital.hte$risk.center.cox <- vital.hte$lp - mean(vital.hte$lp, na.rm = T)

# Baseline risk stratification 
# Dividing patients into four equal sized risk groups 
vital.hte <- vital.hte %>%
  mutate(risk = ntile(risk.center.cox,4),
         risk = factor(risk,
                       levels= c(1,2,3,4),
                       labels = c("Q(low)", "Q2","Q3", "Q4(High)")))


# Fit cox with treatment and risk quarter
cox.tx.risk <- coxph(Surv(time,death) ~ vitd + fishoil + risk, data = vital.hte, ties = "breslow")

# Spline interaction effect model
cox.spline <- coxph(Surv(time,death) ~ (vitd + fishoil) * ns(risk.center.cox, df =3) , data = vital.hte, ties = "breslow")

# Individual variables
vital.subgroup <- vital.hte %>% 
  mutate(
    age.categ = ifelse(age >= 70, "Age >= 70", "Age < 70"),
    bmi.categ = ifelse(bmi >= 30, "Obese(BMI>=30)", "Non-obese(BMI < 30)"),
    sex
  )
subgroup.variables <- c( "age.categ", "bmi.categ")



#****marginal effects for treatment heterogeneity**# 
#**Vitd**

# Overall risk difference
overall.ard <- avg_comparisons(
  cox.tx.risk,
  variables = "vitd",
  type = "survival",
  comparison = "differenceavg",
  conf_level = 0.95,
  time = 5,
  vcov = TRUE
) %>%
  as_tibble() %>%
  mutate(subgroup = "Overall") %>%
  select(subgroup, estimate, conf.low, conf.high, p.value)

# Risk difference by risk quarter
ard.by.risk <- avg_comparisons(
  cox.tx.risk,
  variables = "vitd",
  by = "risk",
  type = "survival",
  comparison = "differenceavg",
  conf_level = 0.95,
  time = 5,
  vcov = TRUE
) %>%
  as_tibble() %>%
  select(subgroup = risk, estimate, conf.low, conf.high, p.value)

# ARD at 5 years among age groups
# Age >= 70
ard.by.age1 <- avg_comparisons(
  cox.tx.risk,
  variables = "vitd",
  newdata = subset(vital.subgroup, age.categ == "Age >= 70"),
  type = "survival",
  comparison = "differenceavg",
  conf_level = 0.95,
  time = 5,
  vcov = TRUE
) %>%
  as_tibble() %>%
  mutate(subgroup = "Age >= 70") %>%
  select(subgroup, estimate, conf.low, conf.high, p.value)

# Age <= 70
ard.by.age2 <- avg_comparisons(
  cox.tx.risk,
  variables = "vitd",
  newdata = subset(vital.subgroup, age.categ == "Age < 70"),
  type = "survival",
  comparison = "differenceavg",
  conf_level = 0.95,
  time = 5,
  vcov = TRUE
) %>%
  as_tibble() %>%
  mutate(subgroup = "Age < 70") %>%
  select(subgroup, estimate, conf.low, conf.high, p.value)

# ARD among at 5 years BMI group
#BMI >= 30
ard.by.bmi1 <- avg_comparisons(
  cox.tx.risk,
  variables = "vitd",
  newdata = subset(vital.subgroup, bmi.categ == "Obese(BMI>=30)"),
  type = "survival",
  comparison = "differenceavg",
  conf_level = 0.95,
  time = 5,
  vcov = TRUE
) %>%
  as_tibble() %>%
  mutate(subgroup = "Obese(BMI>=30)") %>%
  select(subgroup, estimate, conf.low, conf.high, p.value)

# BMI < 30
ard.by.bmi2 <- avg_comparisons(
  cox.tx.risk,
  variables = "vitd",
  newdata = subset(vital.subgroup, bmi.categ == "Non-obese(BMI < 30)"),
  type = "survival",
  comparison = "differenceavg",
  conf_level = 0.95,
  time = 5,
  vcov = TRUE
) %>%
  as_tibble() %>%
  mutate(subgroup = "Non-obese(BMI < 30)") %>%
  select(subgroup, estimate, conf.low, conf.high, p.value)

#***Forest plot: Risk modelling***#
vital.plot.HTE <- bind_rows(overall.ard,
                            ard.by.risk,
                            ard.by.age1,
                            ard.by.age2,
                            ard.by.bmi1,
                            ard.by.bmi2
) %>%
  mutate(subgroup = factor(subgroup,
                           levels = c("Overall", "Q(low)", "Q2", "Q3",
                                      "Q4(High)", "Age >= 70", "Age < 70",
                                      "Obese(BMI>=30)", "Non-obese(BMI < 30)")))

ggplot(vital.plot.HTE,
       aes(x = estimate, y = subgroup)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.2, color = "grey30") +
  geom_point(size = 3, color = "maroon") +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.01)) +
  labs(title = "Forest Plot: Absolute Risk Difference for Vitamin D",
       x = "Absolute Risk Difference at 5 years (99% CI)",
       y = NULL) +
  theme_minimal()



#***marginal effects for treatment heterogeneity fo Omega-3**#

# Overall risk difference
overall.ard.omega <- avg_comparisons(
  cox.tx.risk,
  variables = "fishoil",
  type = "survival",
  comparison = "differenceavg",
  conf_level = 0.95,
  time = 5,
  vcov = TRUE
) %>%
  as_tibble() %>%
  mutate(subgroup = "Overall") %>%
  select(subgroup, estimate, conf.low, conf.high, p.value)

# Risk difference by risk quarter
ard.by.risk.omega <- avg_comparisons(
  cox.tx.risk,
  variables = "fishoil",
  by = "risk",
  type = "survival",
  comparison = "differenceavg",
  conf_level = 0.95,
  time = 5,
  vcov = TRUE
) %>%
  as_tibble() %>%
  select(subgroup = risk, estimate, conf.low, conf.high, p.value)

# ARD at 5 years among age groups
# Age >= 70
ard.by.age1.omega <- avg_comparisons(
  cox.tx.risk,
  variables = "fishoil",
  newdata = subset(vital.subgroup, age.categ == "Age >= 70"),
  type = "survival",
  comparison = "differenceavg",
  conf_level = 0.95,
  time = 5,
  vcov = TRUE
) %>%
  as_tibble() %>%
  mutate(subgroup = "Age >= 70") %>%
  select(subgroup, estimate, conf.low, conf.high, p.value)

# Age <= 70
ard.by.age2.omega <- avg_comparisons(
  cox.tx.risk,
  variables = "fishoil",
  newdata = subset(vital.subgroup, age.categ == "Age < 70"),
  type = "survival",
  comparison = "differenceavg",
  conf_level = 0.95,
  time = 5,
  vcov = TRUE
) %>%
  as_tibble() %>%
  mutate(subgroup = "Age < 70") %>%
  select(subgroup, estimate, conf.low, conf.high, p.value)

# ARD among at 5 years BMI group
#BMI >= 30
ard.by.bmi1.omega <- avg_comparisons(
  cox.tx.risk,
  variables = "fishoil",
  newdata = subset(vital.subgroup, bmi.categ == "Obese(BMI>=30)"),
  type = "survival",
  comparison = "differenceavg",
  conf_level = 0.95,
  time = 5,
  vcov = TRUE
) %>%
  as_tibble() %>%
  mutate(subgroup = "Obese(BMI>=30)") %>%
  select(subgroup, estimate, conf.low, conf.high, p.value)

# BMI < 30
ard.by.bmi2.oemga <- avg_comparisons(
  cox.tx.risk,
  variables = "fishoil",
  newdata = subset(vital.subgroup, bmi.categ == "Non-obese(BMI < 30)"),
  type = "survival",
  comparison = "differenceavg",
  conf_level = 0.95,
  time = 5,
  vcov = TRUE
) %>%
  as_tibble() %>%
  mutate(subgroup = "Non-obese(BMI < 30)") %>%
  select(subgroup, estimate, conf.low, conf.high, p.value)

# Forest plot
vital.plot.HTE.omega <- bind_rows(overall.ard.omega,
                                  ard.by.risk.omega,
                                  ard.by.age1.omega,
                                  ard.by.age2.omega,
                                  ard.by.bmi1.omega,
                                  ard.by.bmi2.oemga,
) %>%
  mutate(subgroup = factor(subgroup,
                           levels = c("Overall", "Q(low)", "Q2", "Q3",
                                      "Q4(High)", "Age >= 70", "Age < 70",
                                      "Obese(BMI>=30)", "Non-obese(BMI < 30)")))
# ARD forest plot for Omega-3 fish oil
ggplot(vital.plot.HTE.omega,
       aes(x = estimate, y = subgroup)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.2, color = "grey30") +
  geom_point(size = 3, color = "darkblue") +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.01)) +
  labs(title = "Forest Plot: Absolute Risk Difference for Omega-3",
       x = "Absolute Risk Difference at 5 years (99% CI)",
       y = NULL) +
  theme_minimal()


#**Effect modelling***#

# survival at t=5 from cox.vital.baseline
surv_baseline <- survfit(cox.vital.baseline)

surv_t5 <- summary(
  surv_baseline,
  times  = 5,
  extend = TRUE
)$surv

surv_baseline_time5 <- surv_t5^exp(vital.hte$lp)

# Baseline risk per patient
vital.hte$baseline_risk <- 1 - surv_baseline_time5

hist(vital.hte$baseline_risk,
     main = "Baseline risk distribution",
     xlab = "Baseline risk at 5 years")

# Proportional effect of VitD 
cox.prop.vitd <- coxph(Surv(time, death)~ vitd + fishoil + risk.center.cox,
                       data = vital.hte,
                       ties = "breslow")

# Hazard ratio vitd from cox.prop.vitd
hr_vitd <- exp(coef(cox.prop.vitd)["vitdActive Vit-D"])

proportional_line <- data.frame(
  baseline_risk = vital.hte$baseline_risk,
  benefit       = surv_baseline_time5^hr_vitd - surv_baseline_time5
)

proportional_line <- proportional_line %>%
  arrange(baseline_risk)

# Grouped ARD by baseline risk quarter
grouped.ard.by.risk <- avg_comparisons(
  cox.tx.risk,
  newdata    = vital.hte,
  variables  = "vitd",
  by         = "risk",
  type       = "survival",
  comparison = "differenceavg",
  conf_level = 0.95,
  time       = 5,
  vcov       = TRUE
) %>%
  as_tibble() %>%
  mutate(
    ard      = -estimate,
    ard.low  = -conf.high,
    ard.high = -conf.low
  ) %>%
  left_join(
    vital.hte %>%
      group_by(risk) %>%
      summarise(med.baseline = median(baseline_risk)),
    by = "risk"
  )

#**Effect modelling for omega-3 fishoil**
cox.prop.omega <- cox.prop.vitd
hr_omega <- exp(coef(cox.prop.omega)["fishoilActive Omega-3"])

# Proportional effect of omega-3 
proportional_line_omega <- data.frame(
  baseline_risk = vital.hte$baseline_risk,
  benefit       = surv_baseline_time5^hr_omega - surv_baseline_time5
)

proportional_line_omega <- proportional_line_omega %>%
  arrange(baseline_risk)

# Grouped ARD by baseline risk quarter
grouped.ard.by.risk.omega <- avg_comparisons(
  cox.tx.risk,
  newdata    = vital.hte,
  variables  = "fishoil",
  by         = "risk",
  type       = "survival",
  comparison = "differenceavg",
  conf_level = 0.95,
  time       = 5,
  vcov       = TRUE
) %>%
  as_tibble() %>%
  mutate(
    ard      = -estimate,
    ard.low  = -conf.high,
    ard.high = -conf.low
  ) %>%
  left_join(
    vital.hte %>%
      group_by(risk) %>%
      summarise(med.baseline = median(baseline_risk)),
    by = "risk"
  )



#*** Absolute benefit versus baseline risk plot***
#**VitD**

xmax <- quantile(vital.hte$baseline_risk, 0.95, na.rm =T)
# X-axis 
xp <- seq(
  min(vital.hte$baseline_risk ,na.rm = T), xmax,
  length.out = 300
)

# baseline risk probability scale
p1exp <- (1 - xp) - (1 - xp)^hr_vitd 

# y-axis
# yax <- range(
#   c(p1exp, grouped.ard.by.risk$ard.low, grouped.ard.by.risk$ard.high)
#   , na.rm =T
# )


png("VITAL_plot_Vit.png", width=1800, height=1200, res=220)
# Proportional effect line
VITAL_plot_Vit <- plot(
  x    = xp,
  y    = p1exp,
  type = "l",
  lty  = 2,
  lwd  = 3,
  xlim = range(0, xmax),
  ylim = range(-0.006,0.009),
  col  = "maroon",
  xlab = "Predicted baseline risk of CVD death at 5 years",
  ylab = "Absolute risk benefit at 5 years: VitD",
  cex.lab = 1.2,
  las  = 1,
  bty  = "l"
)

# reference line (ARD = 0)
abline(h = 0, lty = 1, col ="black")

# Baseline risk distribution 
histSpike(
  vital.hte$baseline_risk,
  add  = TRUE,
  side = 1,
  nint = 300,
  frac = 0.15
)

# Grouped ARD 
points(
  x   = grouped.ard.by.risk$med.baseline,
  y   = grouped.ard.by.risk$ard,
  pch = 1,
  cex = 2,
  lwd = 2,
  col = "darkblue"
)

# 99% CI 
arrows(
  x0    = grouped.ard.by.risk$med.baseline,
  x1    = grouped.ard.by.risk$med.baseline,
  y0    = grouped.ard.by.risk$ard.low,
  y1    = grouped.ard.by.risk$ard.high,
  angle = 90,
  code  = 3,
  length = 0.1,
  col   = "darkblue"
)


legend(
  "topleft",
  lty    = c(2, NA),
  pch    = c(NA, 1),
  lwd    = c(3, 2),
  bty    = "n",
  col    = c("maroon", "darkblue"),
  cex    = 1.1,
  legend = c("Expected with proportional effect", "Grouped patients")
)
VITAL_plot_Vit
dev.off()

#*** Absolute benefit versus baseline risk plot***
#**fishoil/omega-3**

xmax <- quantile(vital.hte$baseline_risk, 0.95, na.rm =T)

# X-axis 
xp <- seq(
  min(vital.hte$baseline_risk ,na.rm = T), xmax,
  length.out = 300
)

# baseline risk probability scale
p1exp <- (1 - xp) - (1 - xp)^hr_omega 

# y-axis
yax <- range(
  c(p1exp, grouped.ard.by.risk.omega$ard.low, grouped.ard.by.risk.omega$ard.high)
  , na.rm =T
)


png("VITAL_plot_omg.png", width=1800, height=1200, res=220)
# Proportional effect line
VITAL_plot_Omg <- plot(
  x    = xp,
  y    = p1exp,
  type = "l",
  lty  = 2,
  lwd  = 3,
  xlim = range(0, xmax),
  ylim = yax,
  col  = "maroon",
  xlab = "Predicted baseline risk of CVD death at 5 years",
  ylab = "Absolute risk benefit at 5yrs:fishoil",
  cex.lab = 1.2,
  las  = 1,
  bty  = "l"
)

# reference line (ARD = 0)
abline(h = 0, lty = 1, col ="black")

# Baseline risk distribution 
histSpike(
  vital.hte$baseline_risk,
  add  = TRUE,
  side = 1,
  nint = 300,
  frac = 0.15
)

# Grouped ARD 
points(
  x   = grouped.ard.by.risk.omega$med.baseline,
  y   = grouped.ard.by.risk.omega$ard,
  pch = 1,
  cex = 2,
  lwd = 2,
  col = "darkblue"
)

# 99% CI 
arrows(
  x0    = grouped.ard.by.risk.omega$med.baseline,
  x1    = grouped.ard.by.risk.omega$med.baseline,
  y0    = grouped.ard.by.risk.omega$ard.low,
  y1    = grouped.ard.by.risk.omega$ard.high,
  angle = 90,
  code  = 3,
  length = 0.1,
  col   = "darkblue"
)


legend(
  "topleft",
  lty    = c(2, NA),
  pch    = c(NA, 1),
  lwd    = c(3, 2),
  bty    = "n",
  col    = c("maroon", "darkblue"),
  cex    = 1.2,
  legend = c("Expected with proportional effect", "Grouped patients")
)
VITAL_plot_Omg
dev.off()


#**Standardized plot for presentation and report*
#**Forest plot: VitD*
library(metafor)
library(dplyr)

png(
  "VITAL_presn_forest_vitd.png",
  width  = 1800,
  height = 1200,
  res    = 220
)


# Category order
vitd.order <- c(
  "Overall",
  "Q4(High)",
  "Q3",
  "Q2",
  "Q1(low)",
  "Age >= 70",
  "Age < 70",
  "Obese (BMI >= 30)",
  "Non-obese (BMI < 30)"
)

# Labels on the left side 
vitd.labels <- c(
  "Overall (Trial Average)",
  "Risk Quarter 4",
  "Risk Quarter 3",
  "Risk Quarter 2",
  "Risk Quarter 1",
  "Age >= 70",
  "Age < 70",
  "Obese (BMI >= 30)",
  "Non-obese (BMI < 30)"
)

# Row positions 
vitd.rows <- c(13, 11, 10, 9, 8, 6, 5, 3, 2)

# data prep
vitd.fp <- vital.plot.HTE %>%
  mutate(
    subgroup = as.character(subgroup),
    estimate = 100 * estimate,
    low      = 100 * conf.low,
    high     = 100 * conf.high
  )

# CI text
ci.txt <- sprintf(
  "%.2f [%.2f, %.2f]",
  vitd.fp$estimate,
  vitd.fp$low,
  vitd.fp$high
)

count_vital <- function(data, tx_var, active_level) {
  tab <- table(data[[tx_var]], data$death)
  data.frame(
    n_active = sum(tab[active_level, ]),
    event_active = tab[active_level, "1"],
    n_placebo = sum(tab["Placebo", ]),
    event_placebo = tab["Placebo", "1"]
  )
}

vitd.counts <- bind_rows(
  count_vital(vital.subgroup, "vitd", "Active Vit-D"),
  count_vital(filter(vital.subgroup, risk == "Q4(High)"), "vitd", "Active Vit-D"),
  count_vital(filter(vital.subgroup, risk == "Q3"), "vitd", "Active Vit-D"),
  count_vital(filter(vital.subgroup, risk == "Q2"), "vitd", "Active Vit-D"),
  count_vital(filter(vital.subgroup, risk == "Q(low)"), "vitd", "Active Vit-D"),
  count_vital(filter(vital.subgroup, age.categ == "Age >= 70"), "vitd", "Active Vit-D"),
  count_vital(filter(vital.subgroup, age.categ == "Age < 70"), "vitd", "Active Vit-D"),
  count_vital(filter(vital.subgroup, bmi.categ == "Obese(BMI>=30)"), "vitd", "Active Vit-D"),
  count_vital(filter(vital.subgroup, bmi.categ == "Non-obese(BMI < 30)"), "vitd", "Active Vit-D")
)


# Plot limits
x.left  <- -8
x.right <-  3.2

par(mar = c(4, 4, 2, 2))

forest(
  x       = vitd.fp$estimate,
  ci.lb   = vitd.fp$low,
  ci.ub   = vitd.fp$high,
  slab    = vitd.labels,
  rows    = vitd.rows,
  xlim    = c(x.left, x.right),
  alim    = c(-1, 1),
  at      = seq(-1, 1, by = 0.5),
  refline = 0,
  ilab = cbind(
    vitd.counts$n_active,
    vitd.counts$event_active,
    vitd.counts$n_placebo,
    vitd.counts$event_placebo
  ),
  ilab.xpos= c(-5.8, -4.8, -3.5, -2.5),
  xlab    = "Vitamin D: Absolute Risk Difference at 5 years, percentage points",
  mlab    = "",
  psize   = 1.4,
  lwd     = 1.8,
  col     = "black",
  cex     = 0.9,
  ylim    = c(0, 16),
  header  = FALSE,
  annotate = FALSE
)


# Top headers
text(x.left, 15.5, "Categories", pos = 4, font = 2, cex = 1.0)
text(0, 15.5, "VITAL Trial: Vitamin D", font = 2, cex = 1.0)
text(x.right, 15.5, "ARD [95% CI]", pos = 2, font = 2, cex = 1.0)
text(c(-5.8, -4.8, -3.5, -2.5),
     15.5,
     c("Vit-D", "Events", "Placebo", "Events"),
     font = 2,
     cex = 0.85)

# Horizontal line under headers
segments(x.left, 15.0, x.right, 15.0, lwd = 1.5)

# Section headers
text(x.left, 12.2, "Baseline risk (quartiles)", pos = 4, font = 2, cex = 0.9)
text(x.left, 6.5, "Age", pos = 4, font = 2, cex = 0.9)
text(x.left, 3.5, "BMI", pos = 4, font = 2, cex = 0.9)
text(x.right, vitd.rows, ci.txt, pos = 2, cex = 0.85, col = "black")


# Right-side CI values
text(
  x.right,
  vitd.rows,
  ci.txt,
  pos = 2,
  cex = 0.85,
  col = "black"
)

dev.off()


#**Forest plot: Omega-3 fish oil*
png(
  "VITAL_presn_forest_omega.png",
  width  = 1800,
  height = 1200,
  res    = 220
)


# Category order
omega.order <- c(
  "Overall",
  "Q4(High)",
  "Q3",
  "Q2",
  "Q1(low)",
  "Age >= 70",
  "Age < 70",
  "Obese (BMI >= 30)",
  "Non-obese (BMI < 30)"
)

# Labels on the left side 
omega.labels <- c(
  "Overall (Trial Average)",
  "Risk Quarter 4",
  "Risk Quarter 3",
  "Risk Quarter 2",
  "Risk Quarter 1",
  "Age >= 70",
  "Age < 70",
  "Obese (BMI >= 30)",
  "Non-obese (BMI < 30)"
)

# Row positions
omega.rows <- c(13, 11, 10, 9, 8, 6, 5, 3, 2)

# data prep
omega.fp <- vital.plot.HTE %>%
  mutate(
    subgroup = as.character(subgroup),
    estimate = 100 * estimate,
    low      = 100 * conf.low,
    high     = 100 * conf.high
  )


# CI text
ci.txt.omega <- sprintf(
  "%.2f [%.2f, %.2f]",
  omega.fp$estimate,
  omega.fp$low,
  omega.fp$high
)

count_vital <- function(data, tx_var, active_level) {
  tab <- table(data[[tx_var]], data$death)
  data.frame(
    n_active = sum(tab[active_level, ]),
    event_active = tab[active_level, "1"],
    n_placebo = sum(tab["Placebo", ]),
    event_placebo = tab["Placebo", "1"]
  )
}

omega.counts <- bind_rows(
  count_vital(vital.subgroup, "fishoil", "Active Omega-3"),
  count_vital(filter(vital.subgroup, risk == "Q4(High)"), "fishoil", "Active Omega-3"),
  count_vital(filter(vital.subgroup, risk == "Q3"), "fishoil", "Active Omega-3"),
  count_vital(filter(vital.subgroup, risk == "Q2"), "fishoil", "Active Omega-3"),
  count_vital(filter(vital.subgroup, risk == "Q(low)"), "fishoil", "Active Omega-3"),
  count_vital(filter(vital.subgroup, age.categ == "Age >= 70"), "fishoil", "Active Omega-3"),
  count_vital(filter(vital.subgroup, age.categ == "Age < 70"), "fishoil", "Active Omega-3"),
  count_vital(filter(vital.subgroup, bmi.categ == "Obese(BMI>=30)"), "fishoil", "Active Omega-3"),
  count_vital(filter(vital.subgroup, bmi.categ == "Non-obese(BMI < 30)"), "fishoil", "Active Omega-3")
)

# Plot limits
x.left  <- -9.5
x.right <-  4.2
par(mar = c(4, 4, 2, 2))


forest(
  x       = omega.fp$estimate,
  ci.lb   = omega.fp$low,
  ci.ub   = omega.fp$high,
  slab    = omega.labels,
  rows    = omega.rows,
  xlim    = c(x.left, x.right),
  alim    = c(-1, 1),
  at      = seq(-1, 1, by = 0.5),
  refline = 0,
  ilab = cbind(
    omega.counts$n_active,
    omega.counts$event_active,
    omega.counts$n_placebo,
    omega.counts$event_placebo
  ),
  ilab.xpos= c(-5.6, -4.5, -3.3, -2.2),
  xlab    = "Omega-3: Absolute Risk Difference at 5 years, percentage points",
  mlab    = "",
  psize   = 1.4,
  lwd     = 1.8,
  col     = "black",
  cex     = 0.9,
  ylim    = c(0, 16),
  header  = FALSE,
  annotate = FALSE
)


# Top headers
text(x.left, 15.5, "Categories", pos = 4, font = 2, cex = 1.0)
text(0, 15.5, "VITAL Trial: Fish oil(Omega-3)", font = 2, cex = 1.0)
text(x.right, 15.5, "ARD [95% CI]", pos = 2, font = 2, cex = 1.0)
text(
  c(-5.8, -4.8, -3.5, -2.5),
  15.5,
  c("Omega-3", "Events", "Placebo", "Events"),
  font = 2,
  cex = 0.9
)

segments(x.left, 15.0, x.right, 15.0, lwd = 1.5)

# Section headers
text(x.left, 12.2, "Baseline risk (quartiles)", pos = 4, font = 2, cex = 0.9)
text(x.left, 6.5, "Age", pos = 4, font = 2, cex = 0.9)
text(x.left, 3.5, "BMI", pos = 4, font = 2, cex = 0.9)


# Right-side CI values
text(
  x.right,
  omega.rows,
  ci.txt,
  pos = 2,
  cex = 0.85,
  col = "black"
)


dev.off()

