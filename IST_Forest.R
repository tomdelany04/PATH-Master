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

#Base variable Selection - The list of all included variables is as follows: age, sex, level of consciousness at presentation,
#presence of wake-up stroke,
# underlying atrial fibrillation, visible infarction on computed tomography, heparinization within 24 h, aspirin administration within
# 3 days, systolic blood pressure, presence of deficits (including face, upper and lower extremities, dysphasia,
#hemianopsia, visuospatial disorder, and other neurological deficits), and aspirin or heparin administration at
#presentation. Age and systolic blood pressure were continuous variables; sex and level of consciousness were
#categorical; all other variables were binomial.

#########################

table1(~death14|aspirin, data= IST_df)


library(broom)

mod_IST <- glm(death14 ~ aspirin, data = IST_df, family = binomial) #Perform logistic regression



or_tbl <- tidy(mod_IST, exponentiate = TRUE, conf.int = TRUE) #odds ratio table



or_clean <- or_tbl %>% #or_tbl will include standard error, statistic,
  select(Odds_Ratio = estimate,
         Lower_CI  = conf.low,
         Upper_CI  = conf.high,
         p.value
  )

or_clean


###Log model recoding
motor_deficit = ifelse(IST_df$RDEF2 == "Y" | IST_df$RDEF3 == "Y", 1, 0)
consc_bin = ifelse(IST_df$RCONSC == "F", 0, 1)
afib = ifelse(IST_df$RATRIAL == "Y", 1, 0)
ct_infarct = ifelse(IST_df$RVISINF == "Y", 1, 0)

#Adjusted based on Thompson et al 2015

# Done in logistic R script

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

model1 <- glm(
  death14 ~ aspirin +
    AGE +
    consc_bin +
    ct_infarct +
    motor_deficit +
    afib,
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
quartile_mean_risk <- tapply(
  IST_risk_base,
  IST_df$IST_rf_groups0,
  mean
)

points(
  x = quartile_mean_risk,
  y = ratediff,
  pch = 1,
  cex = 2,
  lwd = 2,
  col = "blue"
)
arrows(x0 = quartile_mean_risk, x1 = quartile_mean_risk, y0 = CI[, 2], y1 = CI[, 3], angle = 90, code = 3, len = 0.1, col = "blue")

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


############################################
# Forest plot for PATH approach, aligns with Logistic R script (Full plot for report)
# Only change was rf_quartile dataframe from rf groups


IST_df$consc <- ifelse(IST_df$RCONSC == "F", 0, 1)

IST_df$afib <- ifelse(IST_df$RATRIAL == "Y", 1, 0)

############################################################
# Forest plot inputs (need dfs for quartiles, classics and overall _ then bind)
overall_df <- data.frame(
  subgroup = "Overall",

  event_aspirin =
    sum(IST_df$death14 == "Yes" & IST_df$aspirin == 1),

  n_aspirin =
    sum(IST_df$aspirin == 1),

  event_control =
    sum(IST_df$death14 == "Yes" & IST_df$aspirin == 0),

  n_control =
    sum(IST_df$aspirin == 0)
)

# RF quartile subgroup dataframe

rf_quartile_df <- data.frame(

  subgroup = c(
    "Lowest risk",
    "Low risk",
    "Moderate risk",
    "Highest risk"
  ),

  event_aspirin = aspirin_events,
  n_aspirin = aspirin_n,

  event_control = control_events,
  n_control = control_n
)

# CLASSICS

IST_df$age_group <- cut(
  IST_df$AGE,
  breaks = c(0, 65, 75, 120),
  labels = c("<65", "65–75", "75+")
)

age_df <- IST_df %>%
  group_by(age_group) %>%
  summarise(

    event_aspirin =
      sum(death14 == "Yes" & aspirin == 1),

    n_aspirin =
      sum(aspirin == 1),

    event_control =
      sum(death14 == "Yes" & aspirin == 0),

    n_control =
      sum(aspirin == 0)

  ) %>%
  rename(subgroup = age_group)

IST_df$consc_group <- factor(
  IST_df$consc,
  levels = c(0,1),
  labels = c("Alert", "Drowsy/coma")
)


consc_df <- IST_df %>%
  group_by(consc_group) %>%
  summarise(

    event_aspirin =
      sum(death14 == "Yes" & aspirin == 1),

    n_aspirin =
      sum(aspirin == 1),

    event_control =
      sum(death14 == "Yes" & aspirin == 0),

    n_control =
      sum(aspirin == 0)

  ) %>%
  rename(subgroup = consc_group)


IST_df$afib_group <- factor(IST_df$afib, labels = c("No AF", "AF"))

afib_df <- IST_df %>%
  group_by(afib_group) %>%
  summarise(

    event_aspirin =
      sum(death14 == "Yes" & aspirin == 1),

    n_aspirin =
      sum(aspirin == 1),

    event_control =
      sum(death14 == "Yes" & aspirin == 0),

    n_control =
      sum(aspirin == 0)

  ) %>%
  rename(subgroup = afib_group)


#Bind into one df

forest_df <- bind_rows(
  afib_df,
  age_df,
  consc_df,
  rf_quartile_df,
  overall_df
)


###################################
# Plot for forest
library(metafor)

res <- rma(
  ai = event_aspirin,
  bi = n_aspirin - event_aspirin,
  ci = event_control,
  di = n_control - event_control,
  data = forest_df,
  measure = "OR",
  slab = subgroup,
  method = "ML"
)

par(mar = c(4,4,1,2))

forest(

  res,

  xlim = c(-8, 2.5),

  at = log(c(0.5, 1, 1.5)),

  alim = c(log(0.2), log(2)),

  atransf = exp,

  ilab = cbind(
    forest_df$n_aspirin,
    forest_df$event_aspirin,
    forest_df$n_control,
    forest_df$event_control
  ),

  ilab.xpos = c(-5,-4,-3,-2),

  slab = forest_df$subgroup,

  rows = c(1:2, 4:6, 8:9, 11:14, 16),

  xlab = "Odds Ratio",

  mlab = "",

  psize = 1.2,

  lwd = 1.5,

  col = "maroon",

  cex = 0.9
)


############################################################
# Subgroup headers

text(-8, 15, "Baseline risk (quartiles)", pos = 4, font = 2)

text(-8, 10, "Age", pos = 4, font = 2)

text(-8, 7, "Consciousness", pos = 4, font = 2)

text(-8, 3, "Atrial fibrillation", pos = 4, font = 2)


############################################################
# Column headers

text(
  c(-5,-4,-3,-2,3),
  18,
  c(
    "Aspirin",
    "Events",
    "Control",
    "Events",
    "OR [95% CI]"
  ),
  font = 2
)

text(-8, 19, "IST trial", pos = 4, font = 2)
############################################################

#####################################################
# Clean GGPlot for absolute
library(ggplot2)
library(microshades)

#Data prep
curve_df <- data.frame(risk = xp, prop = p1exp)


group_df <- data.frame(
  risk = quartile_mean_risk,
  benefit = ratediff,
  lower = CI[,2],
  upper = CI[,3]
)

hist_df <- data.frame(risk = IST_risk_base)

# plot
aspirin_plot_rf <- ggplot() +

  # baseline risk distribution (scaled density)
  geom_density(
    data = hist_df,
    aes(x = risk, y = after_stat(scaled) * 0.015),
    fill = "grey80",
    color = NA,
    alpha = 0.4
  ) +

  # proportional model
  geom_line(
    data = curve_df,
    aes(x = risk, y = prop),
    linetype = "dashed",
    linewidth = 1.7,
    colour = "#A80050"
  ) +

  # grouped estimates
  geom_point(
    data = group_df,
    aes(x = risk, y = benefit),
    size = 3,
    color = "#A80050"
  ) +

  geom_errorbar(
    data = group_df,
    aes(x = risk, ymin = lower, ymax = upper),
    width = 0.005,
    alpha = 0.2
  ) +

  # zero line
  geom_hline(yintercept = 0, linetype = "dotted") +

  coord_cartesian(xlim = c(0, 0.4), ylim = c(-0.02, 0.035)) +

  labs(
    x = "Baseline risk",
    y = "Benefit by aspirin (absolute risk difference)"
  ) +

  theme_classic(base_size = 13) +
  theme(
    axis.ticks = element_blank(),
    axis.line = element_line(color = "grey60", linewidth = 0.5),
  ) +

  annotate(
    "text",
    x = 0.3,
    y = 0.016,
    label = "Proportional effect",
    color = "#a80050",
    hjust = 0,
    size = 6
  ) +

  labs(
    title = "Absolute Benefit of Aspirin Across RF-Predicted Baseline Risk"
  )

aspirin_plot_rf


saveRDS(aspirin_plot_rf, "aspirin_plot_rf.rds")
