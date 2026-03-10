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
gusto$rf_groups1 <- cut2(rf_risk_base, g = 10) #For deciles at the end

gusto$rf_groups0 <- factor(gusto$rf_groups0, ordered = TRUE) # Ensuring these are ordered
gusto$rf_groups1 <- factor(gusto$rf_groups1, ordered = TRUE)


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
CI      <- BinomDiffCI(x1 = events1, n1 = n1, x2 = events2, n2 = n2, method = "scorecc")

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

# ------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------
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


