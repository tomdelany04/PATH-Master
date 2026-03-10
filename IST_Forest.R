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
group0_rf_IST <- groups_rf_IST[IST_df$tpa == 0]
group1_rf_IST <- groups_rf_IST[IST_df$tpa == 1]


















