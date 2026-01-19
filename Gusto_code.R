data(gusto)   # loads original GUSTO dataset

gusto_raw <- gusto   # always keep this untouched

library(Hmisc)

library(dplyr)
gusto <- gusto_raw %>% filter(tx %in% c("SK","tPA")) %>%
  select(day30, tx, age, Killip, sysbp, pulse, pmi, miloc, sex) %>% mutate(
    day30 = factor(day30, levels = c(0, 1), labels = c("No", "Yes")), across(c(tx, Killip, pmi, miloc, sex), as.factor))


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



tab1 <- table1( ~ age + Killip + sysbp + pulse + pmi + miloc + sex + day30 | tx, data = gusto, test = TRUE)

tab1






#**PATH**
 # First we look at the primary outcome

table1(~day30|tx, data=gusto)



gusto$tx <- droplevels(gusto$tx) # Drop the SK + tPA category which has 0 obs
gusto$tpa <- as.numeric(gusto$tx) - 1 #code for logistic regression

label(gusto$tpa) <- "tPA"

levels(gusto$tpa) <- c("SK", "tPA")

tab2 <- table(gusto$day30, gusto$tpa)
tab2

library(broom)

mod <- glm(day30 ~ tpa, data = gusto, family = binomial) #Perform logistic regression



or_tbl <- tidy(mod, exponentiate = TRUE, conf.int = TRUE) #odds ratio table



or_clean <- or_tbl %>% #or_tbl will include standard error, statistic, drop these
  filter(term == "tpa") %>%
  select(Odds_Ratio = estimate,
         Lower_CI  = conf.low,
         Upper_CI  = conf.high,
         p.value
  )

or_clean



#Extract values from count table to compute absolute differences
library(DescTools)
library(kableExtra)
CI <- BinomDiffCI(
  x1 = tab2[2, 1],
  n1 = sum(tab2[, 1]),
  x2 = tab2[2, 2],
  n2 = sum(tab2[, 2]),
  method = "scorecc"
)

colnames(CI) <- c("Absolute difference", "Lower CI", "Upper CI") #rename columns

result <- round(CI, 3)

result %>%
  as.data.frame() %>%
  kable() %>%
  kable_styling(full_width = FALSE, position = "left")




library(splines)
#Recreate splines shown in harrel's model using glm and ns()

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

summary(model1)



##Check for interaction

model2 <- glm(
  day30 ~ tx*(
    tpa +
      age +
      Killip +
      sysbp_maxed +
      ns(pulse, knots = 50) +
      pmi +
      miloc),
  data = gusto,
  family = binomial
)


anova(model2, test = "LRT")



#The Treatment effect appears constant across all covariates.




#Check for interaction with linear predictor
model_no_tx <- glm(
  day30 ~ age + Killip + sysbp_maxed + ns(pulse, 50) + pmi + miloc,
  data = gusto,
  family = binomial
)

gusto$lp <- predict(model_no_tx, type = "link")

model_lp_int <- glm(
  day30 ~ tx * lp,
  data = gusto,
  family = binomial
)

anova(model_lp_int, test = "LRT")




#There is no significant interaction between treatment and the linear predictor.

#**No interaction is needed**
# We have no evidence against the assumption that the overall effect of treatment is applicable to all patients.


library(ggplot2)

ggplot(gusto, aes(x = plogis(lp))) +
  geom_histogram(bins = 15) + xlim(0,0.4)


#PATH Principle
#The Path Statement is concerned with Reporting RCT results
#stratified by a risk model when we have positive overall trial results.
#This allows us to understand the distribution of effects across the trial population


## The following is an expansion of the classical forest plot for subgroup effects with risk-based subgroups.
#4 risk-based subgroups and 3 classical (Sex, age, type of interaction)

groups <- cut2(gusto$lp, g=4)
group0 <- groups[gusto$tpa==0]  # SK gropup
group1 <- groups[gusto$tpa==1]  # tPA group

rate0 <- prop.table(table(group0, gusto$day30[gusto$tpa==0]),1 )[,"Yes"]
rate1 <- prop.table(table(group1, gusto$day30[gusto$tpa==1]),1 )[,"Yes"]
ratediff <- rate0-rate1 # benefit of tPA by group

# Make a data frame for the results
data.subgroups <- as.data.frame(matrix(nrow=(4+6+1), ncol=10))
colnames(data.subgroups) <- c("tevent", "tnoevent", "cevent", "cnoevent",
                              "name", "type", "tn", "pt", "cn", "pc")

data.subgroups[11,1:4] <- table(gusto$tpa,gusto$day30)[4:1] # overall results
# define event and non-event numbers
events1   <- table(group0, gusto$day30[gusto$tpa==0])[,2]
nevents1  <- table(group0, gusto$day30[gusto$tpa==0])[,1]
events2   <- table(group1, gusto$day30[gusto$tpa==1])[,2]
nevents2  <- table(group1, gusto$day30[gusto$tpa==1])[,1]
n1      <- events1 + nevents1
n2      <- events2 + nevents2

data.subgroups[10:7,1:4] <- cbind(events2,nevents2,events1,nevents1)

gusto$day30n <- as.integer(gusto$day30 == "Yes")


# Use `table`  to get the summary of cell numbers, by subgroup
# SEX
data.subgroups[5,1:4] <- table(1-gusto$day30n,1-gusto$tpa, gusto$sex)[1:4]
data.subgroups[6,1:4] <- table(1-gusto$day30n,1-gusto$tpa, gusto$sex)[5:8]
# AGE
data.subgroups[3,1:4] <- table(1-gusto$day30n,1-gusto$tpa, gusto$age>=75)[1:4]
data.subgroups[4,1:4] <- table(1-gusto$day30n,1-gusto$tpa, gusto$age>=75)[5:8]
# ANT
data.subgroups[1,1:4] <- table(1-gusto$day30n,1-gusto$tpa, gusto$miloc=="Anterior")[1:4]
data.subgroups[2,1:4] <- table(1-gusto$day30n,1-gusto$tpa, gusto$miloc=="Anterior")[5:8]

# Names
data.subgroups[11,5]   <- "Overall"
data.subgroups[10:7,5] <- paste("Quarter",1:4, sep=" ")
data.subgroups[5:6,5]  <- c("Male sex","Female sex")
data.subgroups[3:4,5]  <- c("Age <75","Age>=75")
data.subgroups[1:2,5]  <- c("Other MI","Anterior")

# Type of subgroup
data.subgroups[11,6]   <- ""
data.subgroups[10:7,6] <- c(rep("Risk-based subgroups", length(ratediff)))
data.subgroups[1:6,6] <- c(rep("Location",2), rep("Age",2), rep("Sex",2))

data.subgroups[,7] <- data.subgroups[,1] + data.subgroups[,2]
data.subgroups[,8] <- paste(round(100*data.subgroups[,1] / data.subgroups[,7] , 1),"%", sep="")
data.subgroups[,9] <- data.subgroups[,3] + data.subgroups[,4]
data.subgroups[,10] <- paste(round(100*data.subgroups[,3] / data.subgroups[,9] , 1),"%", sep="")

# Show the data
kable(as.data.frame((data.subgroups))) %>% kable_styling(full_width=F, position = "left")

library(metafor)
par(mar=c(4,4,1,2))
### fit random-effects model (use slab argument to define "study" labels)
res <- rma(ai=tevent, bi=tnoevent, ci=cevent, di=cnoevent, data=data.subgroups, measure="OR",
           slab=name, method="ML")

### set up forest plot (with 2x2 table counts added); rows argument is used
### to specify exactly in which rows the outcomes will be plotted)
par(fg="maroon")
forest(res, xlim=c(-8, 2.5), at=log(c(0.5, 1)), alim=c(log(0.2), log(2)), atransf=exp,
       ilab=cbind(data.subgroups$tn, data.subgroups$pt, data.subgroups$cn, data.subgroups$pc),
       ilab.xpos=c(-5,-4,-3,-2), adj=1,
       cex=.9, ylim=c(0, 24),
       rows=c(1:2, (4:5)-.5, 6:7, 10:13, 15),
       xlab="", mlab="", psize=1, lwd = 1.5, addfit=F,
       col = "maroon",
       shade = F)
# lines(x=c(-.15, -.15), y=c(0, 17)) ## could add a reference line of the overall treatment effect

text(c(-5,-4,-3,-2, 2.2), 18, c("n", "%mort", "n", "%mort", "OR    [95% CI]"),
     font=2, adj=1, cex=.9)
text(-8, 18, c("GUSTO-I trial"), font=2, adj=0, cex=.9)
text(c(-4.5,-2.5),  19, c("tPA", "SK"), font=2, adj=1)



# function for adjustment
library(rms)
options(digits=3)

subgroup.adj <- function(data=gusto, subgroup=gusto$sex) {
  subgroup <- factor(subgroup)

  coef.unadj <- sapply(split(data, subgroup), function(x) lrm.fit(y=x$day30n, x=x$tpa)$coef[2])
  var.unadj  <- sapply(split(data, subgroup), function(x) as.numeric(lrm.fit(y=x$day30n, x=x$tpa)$var[2,2]))
  coef.adj   <- sapply(split(data, subgroup), function(x) lrm.fit(y=x$day30n, x=x$tpa, offset=x$lp)$coef[2])
  var.adj    <- sapply(split(data, subgroup), function(x) as.numeric(lrm.fit(y=x$day30n, x=x$tpa, offset=x$lp)$var[2,2]))

  coef.unadj <- as.numeric(coef.unadj)
  coef.adj   <- as.numeric(coef.adj)
  var.unadj  <- as.numeric(var.unadj)
  var.adj    <- as.numeric(var.adj)

  result <- cbind(coef.unadj, coef.adj, coef.ratio=coef.adj/coef.unadj,
                  SEunadj=sqrt(var.unadj), SEadj=sqrt(var.adj),
                  SEratio=sqrt(var.adj)/ sqrt(var.unadj))
  rownames(result) <- levels(subgroup)
  result
}


kable(as.data.frame(subgroup.adj(gusto, gusto$tpa>-1))) %>%
  kable_styling(full_width=F, position = "left")

kable(as.data.frame(subgroup.adj(gusto, groups))) %>%
  kable_styling(full_width=F, position = "left")

kable(as.data.frame(subgroup.adj(gusto, gusto$age>=75))) %>%
  kable_styling(full_width=F, position = "left")

CI <- BinomDiffCI(x1 = events1, n1 = n1, x2 = events2, n2 = n2, method = "scorecc")
colnames(CI) <- c("Absolute difference", "Lower CI", "Upper CI")
rownames(CI) <- names(events1)
result <- round(CI, 3)
kable(as.data.frame(result)) %>% kable_styling(full_width=F, position = "left")

xp <- seq(0.002,.5,by=0.001)
logxp0 <- log(xp/(1-xp))

lp.no.tx <- gusto$lp
f <- model1

p1exp <- plogis(logxp0) - plogis(logxp0+coef(f)[2])

plot(x=xp, y=p1exp, type='l', lty=2, lwd=3, xlim=c(0,.35), ylim=c(-0.007,.05), col = "maroon",
     xlab="Baseline risk", ylab="Benefit by tPA", cex.lab=1.2, las=1, bty='l' )

lines(x=c(0,.5), y=c(0,0))
histSpike(plogis(lp.no.tx), add=T, side=1, nint=300, frac=.15)

points(x=rate0, y=ratediff, pch=1, cex=2, lwd=2, col = "blue")
arrows(x0=rate0, x1=rate0, y0=CI[,2], y1=CI[,3], angle=90, code=3, len=.1, col = "blue")

legend("topleft", lty=c(2,NA), pch=c(NA,1), lwd=c(3,2), bty='n', col = c("maroon", "blue"),cex=1.2,
       legend=c("Expected with proportional effect", "Grouped patients"))





#Relaxation of the proportional effect
#Inclusion of the linear predictor,
#Benefit was the differences between these 2 risk groups conditional on baseline risk.

h  <- lrm(day30 ~ tpa + tpa * lp, data=gusto, eps=0.005, maxit=30)
h2 <- lrm(day30 ~ tpa + rcs(lp,3)*tpa, data=gusto, eps=0.005, maxit=99)

l0 <- loess(day30n ~ lp, data=gusto, subset=tpa==0)
l1 <- loess(day30n ~ lp, data=gusto, subset=tpa==1)

p1 <- plogis(Predict(h,  tpa=0, lp = logxp0)[,3]) -
  plogis(Predict(h,  tpa=1, lp = logxp0)[,3])
p2 <- plogis(Predict(h2, tpa=0, lp = logxp0)[,3]) -
  plogis(Predict(h2, tpa=1, lp = logxp0)[,3])
l  <- predict(l0, data.frame(lp = logxp0)) -
  predict(l1, data.frame(lp = logxp0))

plot(x=xp, y=p1exp, type='l', lty=1, lwd=4, xlim=c(0,.35), ylim=c(-0.007,.05), col = "maroon",
     xlab="Baseline risk", ylab="Benefit by tPA", cex.lab=1.2, las=1, bty='l' )

lines(x=xp, y=p1, type='l', lty=2, lwd=3, col = "blue")
lines(x=xp, y=p2, type='l', lty=3, lwd=2, col = "green")
lines(x=xp, y=l,  type='l', lty=1, lwd=3, col = "black")

lines(x=c(0,.5), y=c(0,0))
histSpike(plogis(lp.no.tx), add=T, side=1, nint=300, frac=.1)

legend("topleft", lty=c(1,2,3,1), pch=c(NA,NA,NA,NA), lwd=c(4,3,2,3), bty='n',
       col=c("maroon", "blue", "green","black"),
       cex=1.2,
       legend=c("Expected with proportional effect",
                "Linear interaction", "Spline smoothing, 2 df",
                "Loess"))


g <- 4  # quarters

plot(x=xp, y=p1exp, type='l', lty=1, lwd=4,
     xlim=c(0,.35), ylim=c(-0.007,.05),
     col="maroon",
     xlab="Baseline risk", ylab="Benefit by tPA",
     cex.lab=1.2, las=1, bty='l' )

lines(x=xp, y=p1, type='l', lty=2, lwd=3, col="blue")
lines(x=xp, y=p2, type='l', lty=3, lwd=2, col="green")
lines(x=xp, y=l,  type='l', lty=1, lwd=3, col="black")

lines(x=c(0,.5), y=c(0,0))

histSpike(plogis(lp.no.tx), add=TRUE, side=1, nint=300, frac=.1)

groups <- cut2(lp.no.tx, g=g)
group0 <- groups[gusto$tpa==0]
group1 <- groups[gusto$tpa==1]

rate0 <- prop.table(table(group0, gusto$day30[gusto$tpa==0]),1)[,"Yes"]
rate1 <- prop.table(table(group1, gusto$day30[gusto$tpa==1]),1)[,"Yes"]
ratediff <- rate0 - rate1

events1  <- table(group0, gusto$day30[gusto$tpa==0])[,"Yes"]
nevents1 <- table(group0, gusto$day30[gusto$tpa==0])[,"No"]
events2  <- table(group1, gusto$day30[gusto$tpa==1])[,"Yes"]
nevents2 <- table(group1, gusto$day30[gusto$tpa==1])[,"No"]

n1 <- events1 + nevents1
n2 <- events2 + nevents2

CI <- BinomDiffCI(x1=events1, n1=n1, x2=events2, n2=n2, method="scorecc")

points(x=rate0, y=ratediff, pch=1, cex=1, lwd=1, col="black")
arrows(x0=rate0, x1=rate0, y0=CI[,2], y1=CI[,3],
       angle=90, code=3, len=.1, col="black", lwd=.5)

legend("topleft",
       lty=c(1,2,3,1,NA),
       pch=c(NA,NA,NA,NA,1),
       lwd=c(4,3,2,3,1),
       bty='n',
       col=c("maroon","blue","green","black","black"),
       cex=1.2,
       legend=c("Expected with proportional effect",
                "Linear interaction",
                "Spline smoothing, 2 df",
                "Loess",
                "Grouped patients"))




