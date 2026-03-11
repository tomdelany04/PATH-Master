IST_df <- read.csv("IST_corrected.csv")

head(IST_df)

library(Hmisc)
library(dplyr)
library(table1)
library(ranger)
library(meta)

# Treatment and outcome coding for aspirin vs control
IST_df$aspirin <- ifelse(IST_df$RXASP == "Y", 1, 0)
IST_df$death14 <- factor(IST_df$ID14, levels = c(0, 1), labels = c("No", "Yes"))
IST_df$death14n <- as.integer(IST_df$death14 == "Yes")

# Ensure categorical predictors are handled as factors
IST_df$SEX <- as.factor(IST_df$SEX)
IST_df$RCONSC <- as.factor(IST_df$RCONSC)
IST_df$RATRIAL <- as.factor(IST_df$RATRIAL)
IST_df$RSLEEP <- as.factor(IST_df$RSLEEP)
IST_df$RCT <- as.factor(IST_df$RCT)
IST_df$RVISINF <- as.factor(IST_df$RVISINF)

label(IST_df$AGE)     <- "Age (years)"
label(IST_df$SEX)     <- "Sex"
label(IST_df$RCONSC)  <- "Conscious state"
label(IST_df$RSBP)    <- "Systolic BP (mmHg)"
label(IST_df$RATRIAL) <- "Atrial fibrillation"
label(IST_df$RDELAY)  <- "Delay to randomisation (hours)"
label(IST_df$RSLEEP)  <- "Symptoms on waking"
label(IST_df$RCT)     <- "CT before randomisation"
label(IST_df$RVISINF) <- "Infarct visible on CT"
label(IST_df$death14) <- "14-day mortality"

################################################################################################################
#Logistic Approach










###################################################################
############################################ Random Forest Approach

rf_fit_IST <- ranger(
  death14 ~ AGE + SEX + RCONSC + RSBP + RATRIAL + RDELAY +
    RSLEEP + RCT + RVISINF + RDEF1 + RDEF2 + RDEF3 + RDEF4 + RDEF5,
  data = IST_df,
  probability = TRUE,
  num.trees = 1000,
  seed = 88
)

# Out-of-bag predicted baseline risk of 14-day mortality
IST_risk_base <- rf_fit_IST$predictions[, "Yes"]

# Divide patients into quartiles of predicted baseline risk
IST_df$IST_rf_groups0 <- cut2(IST_risk_base, g = 4)
IST_df$IST_rf_groups0 <- factor(IST_df$IST_rf_groups0, ordered = TRUE)

groups_rf_IST <- IST_df$IST_rf_groups0

# Split risk groups by treatment arm: control vs aspirin
group_control_rf <- groups_rf_IST[IST_df$aspirin == 0]
group_aspirin_rf <- groups_rf_IST[IST_df$aspirin == 1]

#########################################################################################################################

# Mortality rates per RF risk group within each treatment arm
rate_control <- prop.table(table(group_control_rf, IST_df$death14[IST_df$aspirin == 0]), 1)[, "Yes"]
rate_aspirin <- prop.table(table(group_aspirin_rf, IST_df$death14[IST_df$aspirin == 1]), 1)[, "Yes"]

# Absolute benefit of aspirin by RF risk group
ratediff <- rate_control - rate_aspirin

# Absolute event counts by RF risk group
# Control arm
control_tab <- table(group_control_rf, IST_df$death14[IST_df$aspirin == 0])
control_events <- control_tab[, "Yes"]
control_nonevents <- control_tab[, "No"]
control_n <- control_events + control_nonevents

# Aspirin arm
aspirin_tab <- table(group_aspirin_rf, IST_df$death14[IST_df$aspirin == 1])
aspirin_events <- aspirin_tab[, "Yes"]
aspirin_nonevents <- aspirin_tab[, "No"]
aspirin_n <- aspirin_events + aspirin_nonevents

########################################################################################
# Forest plot inputs
rf_df <- data.frame(
  subgroup = paste("RF Quartile", 1:4),
  event_aspirin = aspirin_events,
  n_aspirin = aspirin_n,
  event_control = control_events,
  n_control = control_n
)

rf_df

rf_meta <- metabin(
  event.e = event_aspirin,
  n.e     = n_aspirin,
  event.c = event_control,
  n.c     = n_control,
  studlab = subgroup,
  data    = rf_df,
  sm      = "OR",
  method  = "Inverse",
  incr    = 0.5,
  random  = FALSE
)

rf_meta

forest(
  rf_meta,
  layout = "RevMan5",
  leftcols = c("studlab", "event.e", "n.e", "event.c", "n.c"),
  leftlabs = c("Subgroup", "Aspirin events", "Aspirin n", "Control events", "Control n"),
  rightlabs = c("OR", "95% CI"),
  xlab = "Odds Ratio (aspirin vs control)",
  col.diamond = "maroon",
  col.square = "maroon",
  overall.hetstat = FALSE
)

######################################################################################
# Absolute benefit table by RF quartile
library(DescTools)
library(gt)

CI <- BinomDiffCI(
  x1 = control_events,
  n1 = control_n,
  x2 = aspirin_events,
  n2 = aspirin_n,
  method = "scorecc"
)

colnames(CI) <- c("Absolute difference", "Lower CI", "Upper CI")
rownames(CI) <- names(control_events)

groups <- levels(IST_df$IST_rf_groups0)

result_df <- data.frame(
  Group = groups,
  round(CI, 3)
)

result_df$Group <- c(
  "Lowest risk",
  "Low risk",
  "Moderate risk",
  "Highest risk"
)

result_df %>% gt()

##########################################
# Plot under a proportional treatment effect assumption
library(splines)

model1 <- glm( ##Model not tuned at all.
  death14 ~
    aspirin +
    AGE +
    SEX +
    RCONSC +
    RSBP +
    RATRIAL +
    RDELAY +
    RSLEEP +
    RCT +
    RVISINF +
    RDEF1 +
    RDEF2 +
    RDEF3 +
    RDEF4 +
    RDEF5,
  data = IST_df,
  family = binomial
)

IST_df$rf_lp <- log(IST_risk_base / (1 - IST_risk_base))

xp <- seq(0.002, 0.5, by = 0.001)
logxp0 <- log(xp / (1 - xp))

# Expected absolute risk difference if the proportional treatment effect model holds
p1exp <- plogis(logxp0) - plogis(logxp0 + coef(model1)[2])

plot(
  x = xp,
  y = p1exp,
  type = "l",
  lty = 2,
  lwd = 3,
  xlim = c(0, 0.35),
  ylim = c(-0.01, 0.05),
  col = "red",
  xlab = "Baseline risk",
  ylab = "Benefit by aspirin",
  cex.lab = 1.2,
  las = 1,
  bty = "l"
)

# Add horizontal reference line
lines(x = c(0, 0.5), y = c(0, 0))

# Distribution of predicted baseline risk from the RF model
histSpike(IST_risk_base, add = TRUE, side = 1, nint = 300, frac = 0.15)

# Grouped absolute risk differences by RF quartile with confidence intervals
points(x = rate_control, y = ratediff, pch = 1, cex = 2, lwd = 2, col = "blue")
arrows(x0 = rate_control, x1 = rate_control, y0 = CI[, 2], y1 = CI[, 3], angle = 90, code = 3, len = 0.1, col = "blue")

legend(
  "topleft",
  lty = c(2, NA),
  pch = c(NA, 1),
  lwd = c(3, 2),
  bty = "n",
  col = c("red", "blue"),
  cex = 1.2,
  legend = c(
    "Expected with proportional effect",
    "Grouped patients"
  )
)

