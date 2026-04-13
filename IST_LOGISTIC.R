
IST_df <- read.csv("IST_corrected.csv")

head(IST_df)

library(Hmisc)
library(dplyr)
library(table1)
library(ranger)
library(meta)
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



#Adjusted based on HWANGBO et al., 2022

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
# Forest plot inputs

logit_df <- data.frame(
  subgroup = paste("Logit Quartile", 1:4),
  event_aspirin = aspirin_events,
  n_aspirin = aspirin_n,
  event_control = control_events,
  n_control = control_n
)

logit_meta <- metabin(
  event.e = event_aspirin,
  n.e     = n_aspirin,
  event.c = control_events,
  n.c     = control_n,
  studlab = subgroup,
  data    = logit_df,
  sm      = "OR",
  method  = "Inverse",
  incr    = 0.5,
  random  = FALSE
)

forest(
  logit_meta,
  layout = "RevMan5",
  leftcols = c("studlab", "event.e", "n.e", "event.c", "n.c"),
  leftlabs = c("Subgroup", "Aspirin events", "Aspirin n", "Control events", "Control n"),
  rightlabs = c("OR", "95% CI"),
  xlab = "Odds Ratio (aspirin vs control)",
  col.diamond = "maroon",
  col.square = "maroon",
  overall.hetstat = FALSE
)

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
