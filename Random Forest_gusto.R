data(gusto)   # loads original GUSTO dataset

gusto_raw <- gusto   # always keep this untouched

head(gusto_raw)

library(Hmisc)
library(dplyr)
library(table1)
library(ranger)
library(meta)


###### Relabelling
gusto <- gusto_raw %>% filter(tx %in% c("SK","tPA")) %>% mutate(
    day30 = factor(day30, levels = c(0, 1), labels = c("No", "Yes")), across(c(tx, Killip, pmi, miloc, sex), as.factor))

gusto$tpa <- ifelse(gusto$tx == "tPA", 1, 0)
head(gusto)



label(gusto$age)    <- "Age (years)"
label(gusto$Killip) <- "Killip class"
label(gusto$sysbp)  <- "Systolic BP (mmHg)"
label(gusto$pulse)  <- "Pulse (bpm)"
label(gusto$pmi)    <- "Previous MI"
label(gusto$miloc)  <- "MI location"
label(gusto$sex)    <- "Sex"
label(gusto$day30)  <- "30-day mortality"
label(gusto$tx)     <- "Treatment"

gusto$day30 <- factor(gusto$day30, levels = c("No","Yes"))

gusto$day30n <- as.integer(gusto$day30 == "Yes")

###################################################################
############################################ Random Forest Approach



rf_fit0 <- ranger(
  day30 ~ age + Killip + sysbp + pulse + pmi + miloc + sex, #Random forests already model nonlinear relationships automatically. No need for Splines
  data = gusto,
  probability = TRUE,
  num.trees = 1000,
  seed = 88
)


########################################################################################################################
rf_risk_base <- rf_fit0$predictions[, "Yes"] #Vector containing baseline mortality of death for each patient, predictions are "out of bag"

gusto$rf_groups0 <- cut2(rf_risk_base, g = 4) #Cutting the predictions into quartiles

gusto$rf_groups0 <- factor(gusto$rf_groups0, ordered = TRUE) # Ensuring these are ordered


groups_rf <- gusto$rf_groups0 #Convenient storing of grouped variable

#splitting risk groups by teatment arm
group0_rf <- groups_rf[gusto$tpa == 0]
group1_rf <- groups_rf[gusto$tpa == 1]

#########################################################################################################################

#Data Extraction for building forest plot
#Mortality rates per RF risk group
rate_sk <- prop.table(table(group0_rf, gusto$day30[gusto$tpa==0]),1 )[,"Yes"] #P(death | SK , RF risk quartile)
rate_tpa <- prop.table(table(group1_rf, gusto$day30[gusto$tpa==1]),1 )[,"Yes"] #P(death | tPA , RF risk quartile)

ratediff <- rate_sk - rate_tpa # benefit of tPA by group

tab0 <- table(group0_rf, gusto$day30[gusto$tpa==0])
tab1 <- table(group1_rf, gusto$day30[gusto$tpa==1])

#Event counts for forest plot

events1 <- tab0[,"Yes"]#Deaths in SK
nevents1 <- tab0[,"No"]#Survivors is SK

events2 <- tab1[,"Yes"]#Deaths in tPA
nevents2 <- tab1[,"No"]#survivors in tPA

n1 <- events1 + nevents1 #Total patients for group SK
n2 <- events2 + nevents2 #Total Patients for group tPA

#######################################################################################


#Dataframe for event counts
rf_df <- data.frame(
  subgroup = paste("RF Quartile", 1:4),
  event_tpa = events2,
  n_tpa     = n2,
  event_sk  = events1,
  n_sk      = n1
)

rf_df


#Metabin object is good for calculation of common effect and random effects estimates (e.g OR)
rf_meta <- metabin(
  event.e = event_tpa,
  n.e     = n_tpa,
  event.c = event_sk,
  n.c     = n_sk,
  studlab = subgroup,
  data    = rf_df,
  sm      = "OR",
  method  = "Inverse",
  incr    = 0.5,
  random = FALSE
)

#Basic Forest plot for random forest quartiles, no heterogeneity stats included
forest(
  rf_meta,
  layout = "RevMan5",
  leftcols = c("studlab","event.e","n.e","event.c","n.c"),
  leftlabs = c("Subgroup","tPA events","tPA n","SK events","SK n"),
  rightlabs = c("OR","95% CI"),
  xlab = "Odds Ratio (tPA vs SK)",
  col.diamond = "maroon",
  col.square  = "maroon",
  overall.hetstat = FALSE
)

######################################################################################

#Absolute benefit Table - same as FHarrel, just using gt object, will spruce up for a final report
library(DescTools)
library(gt)
CI      <-
  BinomDiffCI(
    x1 = events1,
    n1 = n1,
    x2 = events2,
    n2 = n2,
    method = "scorecc")

colnames(CI) <- c("Absolute difference", "Lower CI", "Upper CI")
rownames(CI) <- names(events1)


groups <- levels(gusto$rf_groups0)

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
#Basic plot with proportional
# ------------------------------------------------------------------------------
# Under a proportional treatment effect assumption, this relationship does not
# hold as well with the random forest model, but the absolute benefit still
# appears greater among higher-risk groups.

# We had to model the proprtional effect using a standard log regression curve
# The histogram is taken from RF
# And of course, so are the point estimates for absolute benefit and CIs
# ------------------------------------------------------------------------------

library(splines)
sysbp_maxed = with(gusto,pmin(sysbp, 120))

model1 <- glm(
  day30 ~
    tpa +
    age +
    Killip +
    sysbp_maxed +
    ns(pulse, knots = 50) +
    pmi +
    miloc,
  data = gusto,
  family = binomial
)

gusto$rf_lp <- log(rf_risk_base/(1-rf_risk_base))

lp.no.tx <- gusto$rf_lp

xp <- seq(0.002,.5,by=0.001)
logxp0 <- log(xp/(1-xp))

# expected difference, if covariate adjusted model holds
p1exp <- plogis(logxp0) - plogis(logxp0+coef(model1)[2]) # proportional effect assumed

plot(x=xp, y=p1exp, type='l', lty=2, lwd=3, xlim=c(0,.35), ylim=c(-0.007,.05), col="red",
     xlab="Baseline risk", ylab="Benefit by tPA", cex.lab=1.2, las=1, bty='l' )

# add horizontal line
lines(x=c(0,.5), y=c(0,0))
# distribution of predicted
histSpike(rf_risk_base, add=T, side=1, nint=300, frac=.15)

points(x=rate_sk, y=ratediff, pch=1, cex=2, lwd=2, col="blue")
arrows(x0=rate_sk, x1=rate_sk, y0=CI[,2], y1=CI[,3], angle=90, code=3,len=.1, col="blue")

legend("topleft", lty=c(2,NA), pch=c(NA,1), lwd=c(3,2), bty='n',col=c("red", "blue"), cex=1.2,
       legend=c("Expected with proportional effect",
                "Grouped patients"))

##################################################################################################

#Relaxation of this effect
library(rms)
h  <- lrm(day30 ~ tpa + tpa * rf_lp, data=gusto, eps=0.005, maxit=30)
h2 <- lrm(day30 ~ tpa + rcs(rf_lp,3)*tpa, data=gusto, eps=0.005, maxit=99)
# loess smoothing
l0 <- loess(day30n ~ rf_risk_base, data=gusto, subset=tpa==0)
l1 <- loess(day30n ~ rf_risk_base, data=gusto, subset=tpa==1)

# subtract predicted risks with from without tx
p1 <- plogis(Predict(h,  tpa=0, rf_lp = logxp0)[,3]) -
  plogis(Predict(h,  tpa=1, rf_lp = logxp0)[,3])
p2 <- plogis(Predict(h2, tpa=0, rf_lp = logxp0)[,3]) -
  plogis(Predict(h2, tpa=1, rf_lp = logxp0)[,3])
l  <- predict(l0, data.frame(rf_risk_base = xp)) -
  predict(l1, data.frame(rf_risk_base = xp))

plot(x=xp, y=p1exp, type='l', lty=1, lwd=4, xlim=c(0,.35), ylim=c(-0.007,.05), col="red",
     xlab="Baseline risk", ylab="Benefit by tPA", cex.lab=1.2, las=1, bty='l' )
# benefit with interaction terms
lines(x=xp, y=p1, type='l', lty=2, lwd=3, col="darkblue")
lines(x=xp, y=p2, type='l', lty=3, lwd=2, col="purple")
lines(x=xp, y=l,  type='l', lty=1, lwd=3, col="black")

# horizontal line
lines(x=c(0,.5), y=c(0,0))
# distribution of predicted
histSpike(rf_risk_base, add=T, side=1, nint=300, frac=.15)

points(x=rate_sk, y=ratediff, pch=1, cex=2, lwd=2, col="blue")
arrows(x0=rate_sk, x1=rate_sk, y0=CI[,2], y1=CI[,3], angle=90, code=3,len=.1, col="blue")

legend("topleft", lty=c(1,2,3,1), pch=c(NA,NA,NA,NA), lwd=c(4,3,2,3), bty='n',
       col=c("red", "darkblue", "purple","black"), cex=1.2,
       legend=c("Expected with proportional effect",
                "Linear interaction", "Spline smoothing, 2 df",
                "Loess"))

############################################################################
##############################################################
# Clean ggplot version
############################################################

library(ggplot2)
library(microshades)

############################################################
# Data prep
############################################################

curve_df <- data.frame(
  risk = xp,
  prop = p1exp
)

spline_df <- data.frame(
  risk = plogis(logxp0),
  spline = p2
)

group_df <- data.frame(
  risk = rate_sk,
  benefit = ratediff,
  lower = CI[,2],
  upper = CI[,3]
)

hist_df <- data.frame(
  risk = rf_risk_base
)

############################################################
# Plot
############################################################

gusto_plot_rf <- ggplot() +

  ##########################################################
# Baseline risk distribution
##########################################################

geom_density(
  data = hist_df,
  aes(
    x = risk,
    y = after_stat(scaled) * 0.015
  ),
  fill = "grey80",
  color = NA,
  alpha = 0.4
) +

  ##########################################################
# Proportional model
##########################################################

  geom_line(
  data = curve_df,
  aes(
    x = risk,
    y = prop
  ),
  linetype = "dashed",
  linewidth = 1.7,
  colour = "darkblue"
) +

  ##########################################################
# Spline model
##########################################################

  geom_line(
  data = spline_df,
  aes(
    x = risk,
    y = spline
  ),
  linewidth = 1.7,
  color = "#a80050"
) +

  ##########################################################
# Grouped estimates
##########################################################

  geom_point(
  data = group_df,
  aes(
    x = risk,
    y = benefit
  ),
  size = 3,
  color = "#a80050"
) +

  geom_errorbar(
    data = group_df,
    aes(
      x = risk,
      ymin = lower,
      ymax = upper
    ),
    width = 0.005,
    alpha = 0.2
  ) +

  ##########################################################
# Zero line
##########################################################

geom_hline(
  yintercept = 0,
  linetype = "dotted"
) +

  ##########################################################
# Axes
##########################################################

coord_cartesian(
  xlim = c(0, 0.3),
  ylim = c(-0.01, 0.05)
) +

  ##########################################################
# Labels
##########################################################

labs(
  x = "Baseline risk",
  y = "Benefit by tPA (absolute risk difference)",
  title = "Absolute Benefit of tPA Across Baseline Risk"
) +

  ##########################################################
# Theme
##########################################################

theme_classic(base_size = 13) +

  theme(
    axis.ticks = element_blank(),
    axis.line = element_line(
      color = "grey60",
      linewidth = 0.5
    )
  ) +

  ##########################################################
# Annotations
##########################################################

  annotate(
  "text",
  x = 0.20,
  y = 0.04,
  label = "Proportional effect",
  color = "darkblue",
  hjust = 0,
  size = 6
) +

  annotate(
    "text",
    x = 0.20,
    y = 0.015,
    label = "Restricted cubic spline (3 knots)",
    color = "#a80050",
    hjust = 0,
    size = 6
  )

gusto_plot_rf

############################################################
# Save plot object
############################################################

saveRDS(
  gusto_plot_rf,
  "gusto_plot_rf.rds"
)

############################################################
# Metafor forest using random forest risk quartiles
############################################################

gusto$age_group <- factor(
  gusto$age >= 75,
  levels = c(FALSE, TRUE),
  labels = c("<75", "75+")
)

gusto$sex_group <- factor(gusto$sex)

gusto$mi_group <- factor(
  gusto$miloc == "Anterior",
  levels = c(FALSE, TRUE),
  labels = c("Other MI", "Anterior")
)

overall_df <- data.frame(
  subgroup = "Overall",
  event_tpa = sum(gusto$day30 == "Yes" & gusto$tpa == 1),
  n_tpa = sum(gusto$tpa == 1),
  event_sk = sum(gusto$day30 == "Yes" & gusto$tpa == 0),
  n_sk = sum(gusto$tpa == 0)
)

quartile_df <- data.frame(
  subgroup = c(
    "Lowest risk",
    "Low risk",
    "Moderate risk",
    "Highest risk"
  ),
  event_tpa = events2,
  n_tpa = n2,
  event_sk = events1,
  n_sk = n1
)

age_df <- gusto %>%
  group_by(age_group) %>%
  summarise(
    event_tpa = sum(day30 == "Yes" & tpa == 1),
    n_tpa = sum(tpa == 1),
    event_sk = sum(day30 == "Yes" & tpa == 0),
    n_sk = sum(tpa == 0),
    .groups = "drop"
  ) %>%
  rename(subgroup = age_group)

sex_df <- gusto %>%
  group_by(sex_group) %>%
  summarise(
    event_tpa = sum(day30 == "Yes" & tpa == 1),
    n_tpa = sum(tpa == 1),
    event_sk = sum(day30 == "Yes" & tpa == 0),
    n_sk = sum(tpa == 0),
    .groups = "drop"
  ) %>%
  rename(subgroup = sex_group)

mi_df <- gusto %>%
  group_by(mi_group) %>%
  summarise(
    event_tpa = sum(day30 == "Yes" & tpa == 1),
    n_tpa = sum(tpa == 1),
    event_sk = sum(day30 == "Yes" & tpa == 0),
    n_sk = sum(tpa == 0),
    .groups = "drop"
  ) %>%
  rename(subgroup = mi_group)

forest_df_rf <- bind_rows(
  mi_df,
  age_df,
  sex_df,
  quartile_df,
  overall_df
)

library(metafor)
par(fg = "black")

res_rf <- rma(
  ai = event_tpa,
  bi = n_tpa - event_tpa,
  ci = event_sk,
  di = n_sk - event_sk,
  data = forest_df_rf,
  measure = "OR",
  slab = subgroup,
  method = "ML"
)

png(
  "gusto_rf_report_forest.png",
  width = 1800,
  height = 1200,
  res = 220
)

par(mar = c(4, 4, 1, 2))

forest(
  res_rf,
  addfit = FALSE,
  xlim = c(-8, 2.5),
  at = log(c(0.5, 1)),
  alim = c(log(0.2), log(2)),
  atransf = exp,
  ilab = cbind(
    forest_df_rf$n_tpa,
    forest_df_rf$event_tpa,
    forest_df_rf$n_sk,
    forest_df_rf$event_sk
  ),
  ilab.xpos = c(-5, -4, -3, -2),
  slab = forest_df_rf$subgroup,
  rows = c(1:2, 4:5, 7:8, 10:13, 15),
  xlab = "Odds Ratio",
  mlab = "",
  psize = 1.2,
  lwd = 1.5,
  colout = c(rep("black", 10), "#a80050"),
  cex = 0.9
)

text(-8, 14, "Baseline risk (quartiles)", pos = 4, font = 2)
text(-8, 9, "Sex", pos = 4, font = 2)
text(-8, 6, "Age", pos = 4, font = 2)
text(-8, 3, "MI location", pos = 4, font = 2)

text(
  c(-5, -4, -3, -2, 3),
  17,
  c("tPA", "Events", "SK", "Events", "OR [95% CI]"),
  font = 2
)

text(-8, 18, "GUSTO-I trial", pos = 4, font = 2)

dev.off()

png(
  "gusto_rf_pres_forest.png",
  width = 1800,
  height = 1200,
  res = 220
)

par(mar = c(4, 4, 1, 2))

gusto_rf_pres_forest <- forest(
  res_rf,
  addfit = FALSE,
  xlim = c(-1.5, 1.5),
  at = log(c(0.5, 1, 2)),
  alim = c(log(0.2), log(2.5)),
  atransf = exp,
  slab = forest_df_rf$subgroup,
  rows = c(1:2, 4:5, 7:8, 10:13, 15),
  xlab = "Odds Ratio",
  mlab = "",
  psize = 1.4,
  lwd = 1.8,
  colout = c(rep("black", 10), "#a80050"),
  cex = 1
)

text(-1.5, 14, "Baseline risk (quartiles)", pos = 4, font = 2)
text(-1.5, 9, "Sex", pos = 4, font = 2)
text(-1.5, 6, "Age", pos = 4, font = 2)
text(-1.5, 3, "MI location", pos = 4, font = 2)
text(-0.4, 17, "GUSTO-I trial", pos = 4, font = 2)

dev.off()
