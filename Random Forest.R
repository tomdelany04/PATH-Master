data(gusto)   # loads original GUSTO dataset

gusto_raw <- gusto   # always keep this untouched

library(Hmisc)

library(dplyr)

gusto <- gusto_raw %>% filter(tx %in% c("SK","tPA")) %>%
  select(day30, tx, age, Killip, sysbp, pulse, pmi, miloc, sex) %>% mutate(
    day30 = factor(day30, levels = c(0, 1), labels = c("No", "Yes")), across(c(tx, Killip, pmi, miloc, sex), as.factor))

gusto$tpa <- ifelse(gusto$tx == "tPA", 1, 0)
head(gusto)



library(table1)

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





table(gusto$day30, useNA = "ifany")



library(ranger)

rf_fit0 <- ranger(
  day30 ~ age + Killip + sysbp + pulse + pmi + miloc + sex, #Random forests already model nonlinear relationships automatically. No need for Splines
  data = gusto,
  probability = TRUE,
  num.trees = 1000,
  seed = 88
)

rf_risk_base <- rf_fit0$predictions[, "Yes"]

gusto$rf_groups0 <- cut2(rf_risk_base, g = 4)
gusto$rf_groups0 <- factor(gusto$rf_groups0, ordered = TRUE)

groups_rf <- gusto$rf_groups0

group0_rf <- groups_rf[gusto$tpa == 0]
group1_rf <- groups_rf[gusto$tpa == 1]

gusto$day30n <- as.integer(gusto$day30 == "Yes")
#Mortality rates per RF risk group
rate0 <- prop.table(table(group0_rf, gusto$day30[gusto$tpa==0]),1 )[,"Yes"]
rate1 <- prop.table(table(group1_rf, gusto$day30[gusto$tpa==1]),1 )[,"Yes"]
ratediff <- rate0-rate1 # benefit of tPA by group


#Event counts for forest plot
events1   <- table(group0_rf, gusto$day30[gusto$tpa==0])[,2]
nevents1  <- table(group0_rf, gusto$day30[gusto$tpa==0])[,1]

events2   <- table(group1_rf, gusto$day30[gusto$tpa==1])[,2]
nevents2  <- table(group1_rf, gusto$day30[gusto$tpa==1])[,1]

n1 <- events1 + nevents1
n2 <- events2 + nevents2

gusto$day30n <- as.integer(gusto$day30 == "Yes")

#Build subgroup dataset
data.subgroups <- as.data.frame(matrix(nrow=(4+6+1), ncol=10))

colnames(data.subgroups) <- c(
  "tevent", "tnoevent",
  "cevent", "cnoevent",
  "name", "type",
  "tn", "pt",
  "cn", "pc"
)

#Overall results
data.subgroups[11,1:4] <- table(gusto$tpa, gusto$day30)[4:1]


#RF risk groups
data.subgroups[10:7,1:4] <- cbind(events2, nevents2, events1, nevents1)


#Classical subgroups
#sex
data.subgroups[5,1:4] <- table(1-gusto$day30n,1-gusto$tpa, gusto$sex)[1:4]
data.subgroups[6,1:4] <- table(1-gusto$day30n,1-gusto$tpa, gusto$sex)[5:8]

#age
data.subgroups[3,1:4] <- table(1-gusto$day30n,1-gusto$tpa, gusto$age>=75)[1:4]
data.subgroups[4,1:4] <- table(1-gusto$day30n,1-gusto$tpa, gusto$age>=75)[5:8]


#MI location
data.subgroups[1,1:4] <- table(1-gusto$day30n,1-gusto$tpa, gusto$miloc=="Anterior")[1:4]
data.subgroups[2,1:4] <- table(1-gusto$day30n,1-gusto$tpa, gusto$miloc=="Anterior")[5:8]

#Labelling for Forest plot
data.subgroups[11,5]   <- "Overall"
data.subgroups[10:7,5] <- paste("RF Risk Quartile",1:4)

data.subgroups[5:6,5]  <- c("Male sex","Female sex")
data.subgroups[3:4,5]  <- c("Age <75","Age>=75")
data.subgroups[1:2,5]  <- c("Other MI","Anterior")


#subgroup type
data.subgroups[11,6]   <- ""
data.subgroups[10:7,6] <- rep("RF Risk-based subgroups",4)

data.subgroups[1:6,6] <- c(
  rep("Location",2),
  rep("Age",2),
  rep("Sex",2)
)


#Mortality %

data.subgroups[,7] <- data.subgroups[,1] + data.subgroups[,2]

data.subgroups[,8] <- paste(
  round(100*data.subgroups[,1] / data.subgroups[,7] , 1),
  "%",
  sep=""
)

data.subgroups[,9] <- data.subgroups[,3] + data.subgroups[,4]

data.subgroups[,10] <- paste(
  round(100*data.subgroups[,3] / data.subgroups[,9] , 1),
  "%",
  sep=""
)


#FIT FOREST
library(metafor)

res <- rma(
  ai = tevent,
  bi = tnoevent,
  ci = cevent,
  di = cnoevent,
  data = data.subgroups,
  measure = "OR",
  slab = name,
  method = "ML"
)


#Forest Plot

par(fg="maroon")

forest(
  res,
  xlim=c(-8, 2.5),
  at=log(c(0.5, 1)),
  alim=c(log(0.2), log(2)),
  atransf=exp,
  ilab=cbind(data.subgroups$tn,
             data.subgroups$pt,
             data.subgroups$cn,
             data.subgroups$pc),
  ilab.xpos=c(-5,-4,-3,-2),
  adj=1,
  cex=.9,
  ylim=c(0, 24),
  rows=c(1:2, (4:5)-.5, 6:7, 10:13, 15),
  xlab="",
  mlab="",
  psize=1,
  lwd=1.5,
  addfit=FALSE,
  col="maroon",
  shade=FALSE
)

text(c(-5,-4,-3,-2, 2.2), 18,
     c("n", "%mort", "n", "%mort", "OR [95% CI]"),
     font=2, adj=1, cex=.9)

text(-8, 18, "GUSTO-I trial", font=2, adj=0, cex=.9)

text(c(-4.5,-2.5), 19,
     c("tPA", "SK"),
     font=2, adj=1)



######## Code for the (Loess, proprtional, spline curve)


