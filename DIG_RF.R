library(dplyr)
library(Hmisc)
library(table1)
library(ranger)
library(meta)
library(DescTools)
library(gt)
library(splines)
library(readr)

dig_raw <- read_csv("DIG.csv", show_col_types = FALSE)

# Baseline variables are taken from the pre-outcome section of DIG.csv.
# Follow-up event indicators and event-day variables are excluded.
baseline_predictors <- c(
  "AGE", "RACE", "SEX", "EJF_PER", "EJFMETH", "CHESTX", "BMI",
  "KLEVEL", "CREAT", "CHFDUR", "RALES", "ELEVJVP", "PEDEMA",
  "RESTDYS", "EXERTDYS", "ACTLIMIT", "S3", "PULCONG", "NSYM",
  "HEARTRTE", "DIABP", "SYSBP", "FUNCTCLS", "CHFETIOL", "PREVMI",
  "ANGINA", "DIABETES", "HYPERTEN", "DIGUSE", "DIURETK", "DIURET",
  "KSUPP", "ACEINHIB", "NITRATES", "HYDRAL", "VASOD"
)

factor_predictors <- c(
  "RACE", "SEX", "EJFMETH", "RALES", "ELEVJVP", "PEDEMA", "RESTDYS",
  "EXERTDYS", "ACTLIMIT", "S3", "PULCONG", "FUNCTCLS", "CHFETIOL",
  "PREVMI", "ANGINA", "DIABETES", "HYPERTEN", "DIGUSE", "DIURETK",
  "DIURET", "KSUPP", "ACEINHIB", "NITRATES", "HYDRAL", "VASOD"
)

numeric_predictors <- setdiff(baseline_predictors, factor_predictors)

dig_rf <- dig_raw %>%
  select(ID, TRTMT, DEATH, all_of(baseline_predictors)) %>%
  mutate(
    digoxin = as.integer(TRTMT),
    death = factor(as.integer(DEATH), levels = c(0, 1), labels = c("No", "Yes")),
    deathn = as.integer(DEATH)
  )

dig_rf[factor_predictors] <- lapply(dig_rf[factor_predictors], factor)
dig_rf[numeric_predictors] <- lapply(dig_rf[numeric_predictors], as.numeric)

dig_rf <- dig_rf %>%
  filter(complete.cases(across(c(digoxin, deathn, all_of(baseline_predictors)))))

label(dig_rf$AGE)      <- "Age (years)"
label(dig_rf$RACE)     <- "Race"
label(dig_rf$SEX)      <- "Sex"
label(dig_rf$EJF_PER)  <- "Ejection fraction (%)"
label(dig_rf$EJFMETH)  <- "EF measurement method"
label(dig_rf$CHESTX)   <- "Cardiothoracic ratio"
label(dig_rf$BMI)      <- "Body mass index"
label(dig_rf$KLEVEL)   <- "Potassium"
label(dig_rf$CREAT)    <- "Creatinine"
label(dig_rf$CHFDUR)   <- "CHF duration"
label(dig_rf$HEARTRTE) <- "Heart rate"
label(dig_rf$DIABP)    <- "Diastolic BP (mmHg)"
label(dig_rf$SYSBP)    <- "Systolic BP (mmHg)"
label(dig_rf$FUNCTCLS) <- "NYHA functional class"
label(dig_rf$death)    <- "All-cause mortality"

table1(~ death | digoxin, data = dig_rf)

rf_formula <- as.formula(
  paste("death ~", paste(baseline_predictors, collapse = " + "))
)

rf_fit_dig <- ranger(
  rf_formula,
  data = dig_rf,
  probability = TRUE,
  num.trees = 1000,
  seed = 88
)

# Out-of-bag predicted baseline mortality risk.
dig_risk_base <- rf_fit_dig$predictions[, "Yes"]

dig_rf$DIG_rf_groups0 <- cut2(dig_risk_base, g = 4)
dig_rf$DIG_rf_groups0 <- factor(dig_rf$DIG_rf_groups0, ordered = TRUE)

groups_rf_DIG <- dig_rf$DIG_rf_groups0

group_placebo_rf <- groups_rf_DIG[dig_rf$digoxin == 0]
group_digoxin_rf <- groups_rf_DIG[dig_rf$digoxin == 1]

rate_placebo <- prop.table(
  table(group_placebo_rf, dig_rf$death[dig_rf$digoxin == 0]), 1
)[, "Yes"]

rate_digoxin <- prop.table(
  table(group_digoxin_rf, dig_rf$death[dig_rf$digoxin == 1]), 1
)[, "Yes"]

ratediff <- rate_placebo - rate_digoxin

placebo_tab <- table(group_placebo_rf, dig_rf$death[dig_rf$digoxin == 0])
placebo_events <- placebo_tab[, "Yes"]
placebo_nonevents <- placebo_tab[, "No"]
placebo_n <- placebo_events + placebo_nonevents

digoxin_tab <- table(group_digoxin_rf, dig_rf$death[dig_rf$digoxin == 1])
digoxin_events <- digoxin_tab[, "Yes"]
digoxin_nonevents <- digoxin_tab[, "No"]
digoxin_n <- digoxin_events + digoxin_nonevents

rf_df <- data.frame(
  subgroup = paste("RF Quartile", 1:4),
  event_digoxin = digoxin_events,
  n_digoxin = digoxin_n,
  event_placebo = placebo_events,
  n_placebo = placebo_n
)

rf_meta <- metabin(
  event.e = event_digoxin,
  n.e = n_digoxin,
  event.c = event_placebo,
  n.c = placebo_n,
  studlab = subgroup,
  data = rf_df,
  sm = "OR",
  method = "Inverse",
  incr = 0.5,
  random = FALSE
)

RF_DIG_forest <- forest(
  rf_meta,
  layout = "RevMan5",
  leftcols = c("studlab", "event.e", "n.e", "event.c", "n.c"),
  leftlabs = c("Subgroup", "Digoxin events", "Digoxin n", "Placebo events", "Placebo n"),
  rightlabs = c("OR", "95% CI"),
  xlab = "Odds Ratio (digoxin vs placebo)",
  col.diamond = "maroon",
  col.square = "maroon",
  overall.hetstat = FALSE
)

CI <- BinomDiffCI(
  x1 = placebo_events,
  n1 = placebo_n,
  x2 = digoxin_events,
  n2 = digoxin_n,
  method = "scorecc"
)


RF_DIG_forest


colnames(CI) <- c("Absolute difference", "Lower CI", "Upper CI")
rownames(CI) <- names(placebo_events)

result_df <- data.frame(
  Group = c("Lowest risk", "Low risk", "Moderate risk", "Highest risk"),
  round(CI, 3)
)

result_df %>% gt()

prop_formula <- as.formula(
  paste("death ~ digoxin +", paste(baseline_predictors, collapse = " + "))
)

model1 <- glm(
  prop_formula,
  data = dig_rf,
  family = binomial
)

dig_rf$rf_lp <- log(dig_risk_base / (1 - dig_risk_base))

xp <- seq(0.002, 0.7, by = 0.001)
logxp0 <- log(xp / (1 - xp))

# Expected absolute risk difference if proportional treatment effect holds.
p1exp <- plogis(logxp0) - plogis(logxp0 + coef(model1)["digoxin"])

plot(
  x = xp,
  y = p1exp,
  type = "l",
  lty = 2,
  lwd = 3,
  xlim = c(0, 0.7),
  ylim = c(min(-0.02, CI[, 2], na.rm = TRUE), max(0.08, CI[, 3], na.rm = TRUE)),
  col = "maroon",
  xlab = "Baseline risk",
  ylab = "Benefit by digoxin",
  cex.lab = 1.2,
  las = 1,
  bty = "l"
)

lines(x = c(0, 0.7), y = c(0, 0), col = "grey50")
histSpike(
  dig_risk_base,
  add = TRUE,
  side = 1,
  nint = 300,
  frac = 0.15,
  col = "grey80"
)

quartile_mean_risk <- tapply(
  dig_risk_base,
  dig_rf$DIG_rf_groups0,
  mean
)

points(
  x = quartile_mean_risk,
  y = ratediff,
  pch = 1,
  cex = 2,
  lwd = 2,
  col = "maroon"
)
arrows(
  x0 = quartile_mean_risk,
  x1 = quartile_mean_risk,
  y0 = CI[, 2],
  y1 = CI[, 3],
  angle = 90,
  code = 3,
  len = 0.1,
  col = "grey50"
)

legend(
  "topleft",
  lty = c(2, NA),
  pch = c(NA, 1),
  lwd = c(3, 2),
  bty = "n",
  col = c("maroon", "maroon"),
  cex = 1.1,
  legend = c("Expected with proportional effect", "Grouped patients")
)
A random forest model was also fitted to estimate baseline mortality risk. Consistent with the baseline characteristics collected in the DIG trial, all prognostic variables available prior to treatment allocation were included as candidate predictors [@InternationalRandomizedTrial1993]. These comprised demographic characteristics, cardiac function measures, clinical signs and symptoms of heart failure, laboratory measurements, cardiovascular history, comorbidities, and concomitant medication use. The model was implemented using the `ranger` package with 1,000 trees and probability estimation enabled.


