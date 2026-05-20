
IST_df <- read.csv("IST_corrected.csv")

head(IST_df)

library(Hmisc)
library(dplyr)
library(table1)
library(ranger)
library(meta)
library(gt)
library(DescTools)
library(gt)



IST_df$aspirin <- ifelse(IST_df$RXASP == "Y", 1, 0)
IST_df$death14 <- factor(IST_df$ID14, levels = c(0, 1), labels = c("No", "Yes"))
IST_df$death14n <- as.integer(IST_df$death14 == "Yes")

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
IST_df$motor <- ifelse(
  IST_df$RDEF2 == "Y" | IST_df$RDEF3 == "Y", 1, 0
)

IST_df$consc <- ifelse(IST_df$RCONSC == "F", 0, 1)

IST_df$afib <- ifelse(IST_df$RATRIAL == "Y", 1, 0)

IST_df$ct_infarct <- ifelse(IST_df$RVISINF == "Y", 1, 0)



#Adjusted based on Thompson et al., 2015
model <- glm(
  death14 ~ AGE + consc + ct_infarct + motor + afib,
  family = binomial(),
  data = IST_df
)




summary(model)
#Linear predictor for later plots
IST_df$lp <- predict(model, type = "link")
# Predicted baseline risk
IST_risk_base <- predict(model, type = "response")

############################## INT test, not needed.
model_lp_int <- glm(
  death14 ~ lp*aspirin,
  family = binomial(),
  data = IST_df
)

anova(model_lp_int, test = "LRT")

#No neccessity for Intercation


#Predict from model

IST_df$lp <- predict(model, type = "link")


library(ggplot2)

ggplot(IST_df, aes(x = plogis(lp))) +
  geom_histogram(bins = 15) + xlim(0,0.4)


############################################################
# Create quartiles of baseline risk

IST_df$IST_logit_groups <- cut2(IST_risk_base, g = 4)
IST_df$IST_logit_groups <- factor(IST_df$IST_logit_groups, ordered = TRUE)

groups_logit <- IST_df$IST_logit_groups

# Split by treatment arm
group_control <- groups_logit[IST_df$aspirin == 0]
group_aspirin <- groups_logit[IST_df$aspirin == 1]





############################################################
# Event rates per group

rate_control <- prop.table(
  table(group_control, IST_df$death14[IST_df$aspirin == 0]), 1
)[, "Yes"]

rate_aspirin <- prop.table(
  table(group_aspirin, IST_df$death14[IST_df$aspirin == 1]), 1
)[, "Yes"]

ratediff <- rate_control - rate_aspirin
############################################################
# Counts for meta-analysis

control_tab <- table(group_control, IST_df$death14[IST_df$aspirin == 0])
control_events <- control_tab[, "Yes"]
control_n <- rowSums(control_tab)

aspirin_tab <- table(group_aspirin, IST_df$death14[IST_df$aspirin == 1])
aspirin_events <- aspirin_tab[, "Yes"]
aspirin_n <- rowSums(aspirin_tab)





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

quartile_df <- data.frame(
  subgroup = c("Lowest risk",
                   "Low risk",
                   "Moderate risk",
                   "Highest risk"),
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
  quartile_df,
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


png("IST_forest_report.png", width = 2000, height = 1400, res = 300)
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

text(-8, 10, "Consciousness", pos = 4, font = 2)

text(-8, 7, "Age", pos = 4, font = 2)

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


dev.off()

############################################################
# Absolute benefit with CI

CI <- BinomDiffCI(
  x1 = control_events,
  n1 = control_n,
  x2 = aspirin_events,
  n2 = aspirin_n,
  method = "scorecc"
)

colnames(CI) <- c("Absolute difference", "Lower CI", "Upper CI")
rownames(CI) <- names(control_events)

groups <- levels(IST_df$IST_logit_groups)

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

############################################################
# Proportional treatment effect plot

model1 <- glm(
  death14 ~ aspirin + lp,
  data = IST_df,
  family = binomial
)

xp <- seq(0.002, 0.5, by = 0.001)
logxp0 <- log(xp / (1 - xp))

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
  las = 1,
  bty = "l"
)

lines(x = c(0, 0.5), y = c(0, 0))

# Distribution of predicted baseline risk (logistic)
histSpike(IST_risk_base, add = TRUE, side = 1, nint = 300, frac = 0.15)

# Grouped results
points(x = rate_control, y = ratediff, pch = 1, cex = 2, lwd = 2, col = "blue")
arrows(
  x0 = rate_control,
  x1 = rate_control,
  y0 = CI[, 2],
  y1 = CI[, 3],
  angle = 90,
  code = 3,
  len = 0.1,
  col = "blue"
)

legend(
  "topleft",
  lty = c(2, NA),
  pch = c(NA, 1),
  lwd = c(3, 2),
  bty = "n",
  col = c("red", "blue"),
  legend = c(
    "Expected with proportional effect",
    "Grouped patients"
  )
)

############################################################
# Risk distribution plot

ggplot(IST_df, aes(x = IST_risk_base)) +
  geom_histogram(bins = 15) +
  xlim(0, 0.4) +
  labs(x = "Predicted baseline risk", y = "Count")

############################################################
# Spline regression plot:
library(splines)
model1 <- glm(
  death14 ~ aspirin + lp,
  data = IST_df,
  family = binomial
)

#spline approach
knots <- quantile(IST_df$lp, probs = c(0.25, 0.5, 0.75))

model2 <- glm(
  death14 ~ aspirin + ns(lp, df = 3) + aspirin:lp,
  data = IST_df,
  family = binomial
)

#predictions for proportional
xp <- seq(0.002, 0.5, by = 0.001)
logxp0 <- log(xp / (1 - xp))

p1exp <- plogis(logxp0) - plogis(logxp0 + coef(model1)[2])

#Predictions for spline
quantiles <- quantile(IST_df$lp, probs = c(0.01, 0.99))

xp2 <- seq(quantiles[1], quantiles[2], length.out = 500)

newdata_control <- data.frame(lp = xp2, aspirin = 0)
newdata_treated <- data.frame(lp = xp2, aspirin = 1)

p_control <- predict(model2, newdata = newdata_control, type = "response")
p_treated <- predict(model2, newdata = newdata_treated, type = "response")

p2exp <- p_control - p_treated




plot(
  x = xp,
  y = p1exp,
  type = "l",
  lty = 2,
  lwd = 3,
  xlim = c(0, 0.35),
  ylim = c(-0.01, 0.04),
  col = "red",
  xlab = "Baseline risk",
  ylab = "Benefit by aspirin",
  las = 1,
  bty = "l"
)

lines(
  x = plogis(xp2),
  y = p2exp,
  col = "darkgreen",
  lwd = 3
)

lines(x = c(0, 0.4), y = c(0, 0))

# Distribution of predicted baseline risk (logistic)
histSpike(IST_risk_base, add = TRUE, side = 1, nint = 300, frac = 0.15)

# Grouped results
points(x = rate_control, y = ratediff, pch = 1, cex = 2, lwd = 2, col = "blue")
arrows(
  x0 = rate_control,
  x1 = rate_control,
  y0 = CI[, 2],
  y1 = CI[, 3],
  angle = 90,
  code = 3,
  len = 0.1,
  col = "blue"
)

legend(
  "topleft",
  lty = c(2, 1, NA),
  pch = c(NA, NA, 1),
  lwd = c(3, 3, 2),
  bty = "n",
  col = c("red", "darkgreen", "blue"),
  legend = c(
    "Proportional effect",
    "Spline model",
    "Grouped patients"
  )
)


range(IST_risk_base)



#########################################
#cleaner plot
library(ggplot2)
library(microshades)

#Data prep
curve_df <- data.frame(risk = xp, prop = p1exp)
spline_df <- data.frame(risk = plogis(xp2), spline = p2exp)

group_df <- data.frame(
  risk = rate_control,
  benefit = ratediff,
  lower = CI[,2],
  upper = CI[,3]
)

hist_df <- data.frame(risk = IST_risk_base)

# plot
aspirin_plot <- ggplot() +

  # baseline risk distribution (scaled density)
  geom_density(
    data = hist_df,
    aes(x = risk, y = after_stat(scaled) * 0.015),
    fill = "grey80",
    color = NA,
    alpha = 0.2
  ) +

  # proportional model
  geom_line(
    data = curve_df,
    aes(x = risk, y = prop),
    linetype = "dashed",
    linewidth = 1.7,
    colour = "#F09163"
  ) +

  # spline model
  geom_line(
    data = spline_df,
    aes(x = risk, y = spline),
    linewidth = 1.7,
    color = "#4292C6"
  ) +

  # grouped estimates
  geom_point(
    data = group_df,
    aes(x = risk, y = benefit),
    size = 3,
    color = "#238B45"
  ) +

  geom_errorbar(
    data = group_df,
    aes(x = risk, ymin = lower, ymax = upper),
    width = 0.005,
    alpha = 0.2
  ) +

  # zero line
  geom_hline(yintercept = 0, linetype = "dotted") +

  coord_cartesian(xlim = c(0, 0.5), ylim = c(-0.02, 0.035)) +

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
    x = 0.4,
    y = 0.016,
    label = "Proportional effect",
    color = "#F09163",
    hjust = 0,
    size = 6
  ) +

  annotate(
    "text",
    x = 0.3,
    y = 0.0085,
    label = "Spline (df = 3)",
    color = "#4292C6",
    hjust = 0,
    size = 6
  ) +

  labs(
    title = "Absolute Benefit of Aspirin Across Baseline Risk"
  )

aspirin_plot


saveRDS(aspirin_plot, "aspirin_plot.rds")




############################################################
# Forest for poster

png(
  "IST_pres_forest.png",
  width = 1800,
  height = 1200,
  res = 220
)

par(mar = c(4,4,1,2))

gusto_pres_forest <- forest(
    res,
    xlim = c(-5, 2.5),
    at = log(c(0.5, 1, 1.5)),
    alim = c(log(0.2), log(2)),
    atransf = exp,
    slab = forest_df$subgroup,
    rows = c(1:2, 4:6, 8:9, 11:14, 16),
    xlab = "Odds Ratio",
    mlab = "",
    psize = 1.2,
    lwd = 1.5,
    col = "red",
    cex = 0.9
  )

############################################################
# Subgroup headers
############################################################

text(-5, 15, "Baseline risk (quartiles)", pos = 4, font = 2)
text(-5, 10, "Consciousness", pos = 4, font = 2)
text(-5, 7, "Age", pos = 4, font = 2)
text(-5, 3, "Atrial fibrillation", pos = 4, font = 2)



############################################################
# Title
############################################################

text(-2, 18, "IST-1 trial", pos = 4, font = 2)

dev.off()







