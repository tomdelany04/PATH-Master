IST_df <- read.csv("IST_corrected.csv")

head(IST_df)

library(Hmisc)
library(dplyr)
library(table1)
library(ranger)
library(meta)


IST_df$aspirin <- ifelse(IST_df$RXASP == "Y", 1, 0)
IST_df$death14 <- factor(IST_df$ID14, levels=c(0,1), labels=c("No","Yes"))
IST_df$death14n <- as.integer(IST_df$death14=="Yes")

label(IST_df$AGE)    <- "Age (years)"
label(IST_df$SEX)    <- "Sex"
label(IST_df$RCONSC) <- "Conscious state"
label(IST_df$RSBP)   <- "Systolic BP (mmHg)"
label(IST_df$RATRIAL)<- "Atrial fibrillation"
label(IST_df$RDELAY) <- "Delay to randomisation (hours)"
label(IST_df$RSLEEP) <- "Symptoms on waking"
label(IST_df$RCT)    <- "CT before randomisation"
label(IST_df$RVISINF)<- "Infarct visible on CT"
label(IST_df$death14)<- "14-day mortality"


rf_fit_IST <- ranger(
  death14 ~ AGE + SEX + RCONSC + RSBP + RATRIAL + RDELAY +
    RSLEEP + RCT + RVISINF +
    RDEF1 + RDEF2 + RDEF3 + RDEF4 + RDEF5, #Random forests already model nonlinear relationships automatically. No need for Splines
  data = IST_df,
  probability = TRUE,
  num.trees = 1000,
  seed = 88
)


IST_risk_base <- rf_fit_IST$predictions[, "Yes"]


IST_df$IST_rf_groups0 <- cut2(IST_risk_base, g = 4) #Cutting the predictions into quartiles

IST_df$IST_rf_groups0 <- factor(IST_df$IST_rf_groups0, ordered = TRUE) # Ensuring these are ordered


groups_rf_IST <- gusto$rf_groups0 #Convenient storing of grouped variable

#splitting risk groups by teatment arm
group0_rf_IST <- groups_rf_IST[IST_df$aspirin == 0]
group1_rf_IST <- groups_rf_IST[IST_df$aspirin == 1]




#########################################################################################################################

#Data Extraction for building forest plot
#Mortality rates per RF risk group
rate_asp <- prop.table(table(group0_rf_IST, IST_df$death14[IST_df$aspirin==0]),1 )[,"Yes"] #P(death | SK , RF risk quartile)
rate_plc <- prop.table(table(group1_rf_IST, IST_df$death14[IST_df$aspirin==1]),1 )[,"Yes"] #P(death | tPA , RF risk quartile)

ratediff <- rate_asp - rate_plc # benefit of asp by group

tab0 <- table(group0_rf, IST_df$death14[IST_df$aspirin==0])
tab1 <- table(group1_rf, IST_df$death14[IST_df$aspirin==1])

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













