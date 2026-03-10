data(gusto)   # loads original GUSTO dataset

gusto_raw <- gusto   # always keep this untouched

library(Hmisc)

library(dplyr)
library(table1)
library(ranger)
library(metafor)


###### Relabelling
gusto <- gusto_raw %>% filter(tx %in% c("SK","tPA")) %>%
  select(day30, tx, age, Killip, sysbp, pulse, pmi, miloc, sex) %>% mutate(
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
