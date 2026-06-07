
library(dplyr)
library(table1)
library(kableExtra)
library(Hmisc)      
library(DescTools)  
library(metafor)
library(readr)

#Loading and preparing dataset

dig <- read_csv("DIG.csv")   


dig.df <- dig %>%
  select(TRTMT, AGE, EJF_PER, FUNCTCLS, DEATH, SYSBP, HEARTRTE, PREVMI) %>%
  mutate(
    TRTMT = as.integer(TRTMT),
    DEATH = as.integer(DEATH),
    TRTMT_f =factor(TRTMT, levels = c(0,1), labels = c("Placebo", "Digoxin")),
    DEATH_f = factor(DEATH, levels = c(0,1), labels = c("Alive", "Dead")),
    DEATHn = DEATH, 
    Digoxin = TRTMT, 
    Classification = factor(
      FUNCTCLS, levels = c(1, 2, 3, 4),
      labels = c("NYHA_Class1", "NYHA_Class2", "NYHA_Class3", "NYHA_Class4")
    ),
    Ejection_Fraction = EJF_PER,  # EF(0% - 100%)
    sysbp_maxed = pmin(SYSBP, 120)
  )
dig.df <- na.omit(dig.df)



# simple cross-table

table1(~ factor(DEATH) | factor(TRTMT), data=dig, digits=2) 


#**Baseline Descriptive table(Group wise summary by treatment group)**
dig.discr <- table1(
  ~ DEATH + AGE + Ejection_Fraction + Classification,
  data = dig.df, digits = 2
)


#Drop unused labels
dig$TRTMT_f <- droplevels(factor(dig$TRTMT,
                                 levels= c(0,1),
                                 labels= c("Placebo", "Digoxin")))

dig$digoxin01 <- as.numeric(dig$TRTMT_f) - 1
label(dig$digoxin01) <- "Digoxin"
levels(dig$digoxin01) <- c("Placebo", "Digpxon")

# crude table
dig_tab2 <- table(dig$DEATH, dig$digoxin01)

# Odds ration and 95% CI
or_dig <- OddsRatio(dig_tab2, conf.level = 0.95)
names(or_dig) <- c("Odds Ratio", "Lower CI", "Upper CI")

kable(as.data.frame(t(or_dig)))%>%
  kable_styling(full_width = FALSE, position = "left")

# Absolute risk difference (Baseline risk - Digoxin Risk) and 95% CI
library(DescTools)

CI_dig <- BinomDiffCI(
  x1 = dig_tab2["1", "0"],
  n1 = sum(dig_tab2[, "0"]),
  x2 = dig_tab2["1", "1"],
  n2 = sum(dig_tab2[, "1"]),
  method = "scorecc"
)

colnames(CI_dig) <- c("Absolute difference", "Lower CI", "Upper CI")

round(CI_dig, 3) %>%
  as.data.frame() %>%
  kable() %>%
  kable_styling(full_width = FALSE, position = "left")
# Patients taking Digoxin have about 1.4% lower odds of DEATH compared with placebo.
# As the CI includes 1, no statistical evidence agains Ho(OR=1).
# Therefore, No evidence that Digoxin affects mortality(Digoxing does not reduce mortality).

# Adjusted outcome model including treatment (for btx)
model_with_tx <- glm(
  DEATHn ~ Digoxin + AGE + Classification + SYSBP + PREVMI + EJF_PER,
  data = dig.df,
  x = T, y= T,
  maxit = 99
)

print(model_with_tx)
coef(model_with_tx)["Digoxin"] #-0.00224  


#Check for interaction
# 1. Traditional approach
g <- glm(DEATHn ~ Digoxin * (AGE + Classification + SYSBP
                             + PREVMI + EJF_PER),
         data = dig.df,
         x=T , y=T,
         maxit = 100)

print(anova(g, test = 'LRT'))

# 2. PATH statement: linear interaction with linear predictor; baseline risk, tx=ref, SK
# Baseline risk model excluding treatment (PATH lp model)




model_no_tx <- glm(
  DEATHn ~ AGE + Classification + SYSBP +
    PREVMI + EJF_PER,
  data = dig.df,
  family = binomial()
)

# Predict lp for ALL rows 
dig.df$lp <- predict(model_no_tx, newdata = dig.df, type = "link")

# Use complete-case dataset for PATH interaction and risk grouping 
dig_clean <- dig.df %>%
  select(DEATHn, Digoxin, lp) %>%
  filter(complete.cases(.))
# PATH interaction test: Digoxin * lp  (compare main vs interaction)
m_lp_main <- glm(DEATHn ~ Digoxin + lp, data = dig_clean, family = binomial())
m_lp_int  <- glm(DEATHn ~ Digoxin * lp, data = dig_clean, family = binomial())

anova(m_lp_main, m_lp_int, test = "LRT")

# Risk quarters (4 rows): use events2/nevents2 for Digoxin, events1/nevents1 for Placebo
library(Hmisc)
groups <- cut2(dig_clean$lp, g = 4)

# Placebo and Digoxin groups (by quartile)
group0 <- groups[dig_clean$Digoxin == 0]  # Placebo
group1 <- groups[dig_clean$Digoxin == 1]  # Digoxin

# Event rates by quartile (Dead)
rate0 <- prop.table(table(group0, dig_clean$DEATHn[dig_clean$Digoxin == 0]), 1)[, "1"]
rate1 <- prop.table(table(group1, dig_clean$DEATHn[dig_clean$Digoxin == 1]), 1)[, "1"]
ratediff <- rate0 - rate1  # Benefit = Placebo risk - Digoxin risk

#Subgroup results coubtainer
data.subgroups <- as.data.frame(matrix(nrow = 11, ncol = 10))
colnames(data.subgroups) <- c("tevent","tnoevent","cevent","cnoevent",
                              "name","type","tn","pt","cn","pc")
data.subgroups <- na.omit(data.subgroups)
#Overall results
tab_all <- table(dig_clean$Digoxin, dig_clean$DEATHn)
data.subgroups[11, 1:4] <- c(
  tevent   = tab_all["1", "1"], #Digoxin DEATHs
  tnoevent = tab_all["1", "0"], # Digoxin Alive
  cevent   = tab_all["0", "1"], #Placebo DEATHs
  cnoevent = tab_all ["0", "0"] # Placebo Alive
)
# counts for CI
events1  <- table(group0, dig_clean$DEATHn[dig_clean$Digoxin == 0])[, "1"] #Placebo DEATHs
nevents1 <- table(group0, dig_clean$DEATHn[dig_clean$Digoxin == 0])[, "0"] # Placebo Alive

events2  <- table(group1, dig_clean$DEATHn[dig_clean$Digoxin == 1])[, "1"] # Digoxin DEATHs
nevents2 <- table(group1, dig_clean$DEATHn[dig_clean$Digoxin == 1])[, "0"] # Digoxin Alive

n1 <- events1 + nevents1
n2 <- events2 + nevents2

data.subgroups[10:7, 1:4] <- cbind(
  tevent   = as.numeric(events2), 
  tnoevent = as.numeric(nevents2), 
  cevent   = as.numeric(events1),
  cnoevent = as.numeric(nevents1)
)

# DIG Risk quarter forest plot
# NYHA Functional Class (Variable 26: Classes I/II vs III/IV)
# Note: FUNCTCLS > 2 captures NYHA Class III and IV
data.subgroups[5,1:4] <- table(dig.df$DEATHn, dig.df$TRTMT, dig.df$FUNCTCLS > 2)[1:4]
data.subgroups[6,1:4] <- table(dig.df$DEATHn, dig.df$TRTMT, dig.df$FUNCTCLS > 2)[5:8]

# AGE (Using 70 years as the primary DIG risk threshold)
data.subgroups[3,1:4] <- table(dig.df$DEATHn, dig.df$TRTMT, dig.df$AGE >= 70)[1:4]
data.subgroups[4,1:4] <- table(dig.df$DEATHn, dig.df$TRTMT, dig.df$AGE >= 70)[5:8]

# Ejection Fraction (Using 0.25 as the high-risk cut-off)
data.subgroups[1,1:4] <- table(dig.df$DEATHn, dig.df$TRTMT, dig.df$EJF_PER <= 25)[1:4]
data.subgroups[2,1:4] <- table(dig.df$DEATHn, dig.df$TRTMT, dig.df$EJF_PER <= 25)[5:8]

# 2. Subgroup Names
data.subgroups[8, 9]   <- "Overall"
data.subgroups[10:7,5] <- paste("Quarter", 1:4, sep=" ")
data.subgroups[5:6,5]  <- c("NYHA Class I or II", "NYHA Class III or IV")
data.subgroups[3:4,5]  <- c("Age < 70", "Age >= 70")
data.subgroups[1:2,5]  <- c("EF > 0.25", "EF <= 0.25")

# 3. Type of Subgroup
data.subgroups[8, 10]   <- ""
data.subgroups[10:7,6] <- c(rep("Risk-based subgroups", 4))
data.subgroups[1:6,6]  <- c(rep("Ejection Fraction", 2), 
                            rep("Age", 2), 
                            rep("Functional Class", 2))

# 4. Calculate Totals and Mortality Percentages (Placebo and Digoxin)
data.subgroups[,7] <- data.subgroups[,1] + data.subgroups[,2] # N Placebo
data.subgroups[,8] <- paste(round(100 * data.subgroups[,2] / data.subgroups[,7], 1), "%", sep="") # % Mort Placebo
data.subgroups[,9] <- data.subgroups[,3] + data.subgroups[,4] # N Digoxin
data.subgroups[,10] <- paste(round(100 * data.subgroups[,4] / data.subgroups[,9], 1), "%", sep="") # % Mort Digoxin

# Table
kable(as.data.frame(data.subgroups), 
      col.names = c("P_Alive", "P_Dead", "D_Alive", "D_Dead", "Subgroup", "Type", "N_Plac", "%_Plac", "N_Dig", "%_Dig")) %>% 
  kable_styling(full_width = F, position = "left")
colnames(data.subgroups)


par(mar=c(4,4,1,2))
data.subgroups[1:2, "name"]  <- c("EF > 0.25", "EF <= 0.25")
data.subgroups[3:4, "name"]  <- c("Age < 70", "Age >= 70")
data.subgroups[5:6, "name"]  <- c("NYHA Class I/II", "NYHA Class III/IV")
data.subgroups[7:10, "name"] <- paste("Risk Quarter", 1:4)
data.subgroups[11, "name"]   <- "Overall (Trial Average)"
res <- rma(ai=tevent, bi=tnoevent, ci=cevent, di=cnoevent, 
           data=data.subgroups, measure="OR",
           slab=name, method="ML")

forest(res, 
       xlim=c(-1.2,0.8), 
       at=log(c(0.5, 1, 2)), 
       alim=c(log(0.2), log(2)),
       atransf=exp,
       adj=1,
       cex=.9, ylim=c(0,22),
       rows=c(1:2, (4:5)-.5, 6:7, 10:13, 15),
       xlab="",
       mlab="", 
       psize=1, 
       lwd=1.5, 
       addfit=FALSE,
       colout = c(rep("black", 10), 
                  "#a80050"),
       shade=FALSE)


text(-1.2, 18, "DIG Trial Subgroups", font=2, adj=0, cex=.9)  # left header
text(0.8,  18, "OR [95% CI]",        font=2, adj=1, cex=.9)  # right header

# Save for poster use
saveRDS(res, file = "dig_forest_baseline_risk_model.rds")


#**Risk vs Benefit plot*

# Treatment effect model  
m_lp_main <- glm(DEATHn ~ Digoxin + lp, data = dig_clean, family = binomial())

beta_tx <- coef(m_lp_main)["Digoxin"]   # log-odd for Digoxin vs Placebo

# Expected benefit curve as a function of baseline risk xp
xp <- seq(0.002, 0.7, by = 0.001)      # baseline risk grid
logxp0 <- log(xp/(1-xp))               # baseline logit

# Expected ARD = P(Y=1|Placebo) - P(Y=1|Digoxin)
p1exp <- plogis(logxp0) - plogis(logxp0 + beta_tx)

# Distribution of baseline risk from the tx_no model 
lp.no.tx <- dig.df$lp   # from model_no_tx
baseline_risk_dist <- plogis(lp.no.tx)

# 95% CI for grouped patients
n1 <- as.numeric(events1 + nevents1)  
n2 <- as.numeric(events2 + nevents2)  

CI <- BinomDiffCI(
  x1 = as.numeric(events1), n1 = n1,   # Placebo worsening
  x2 = as.numeric(events2), n2 = n2,   # Digoxin Dead
  method = "scorecc"
)


# CI columns: est, lwr, upr 
ci_low <- CI[, 2]
ci_up  <- CI[, 3]

#Plot 
par(mar=c(5,5,2,2))
plot(x = xp,
     y = p1exp,
     type = "l",
     lty = 2, 
     lwd = 3, 
     colout = c(rep("black", 10), 
                "#a80050"),
     xlim = c(0.15, 0.5), 
     ylim = c(min(-0.01, min(ci_low, na.rm=TRUE)), 0.25),
     xlab = "Baseline risk",
     ylab = "Benefit by Digoxin",
     cex.lab = 1.2, 
     las = 1, 
     bty = "l")

# zero line
lines(x = c(0, 0.5), y = c(0, 0), col = "gray30")

# spike histogram along x-axis
histSpike(baseline_risk_dist, add = TRUE, side = 1, nint = 300, frac = 0.15)

# grouped patients (quartiles): x = placebo risk in quartile, y = ARR
points(x = as.numeric(rate0), y = as.numeric(ratediff),
       pch = 1, cex = 1.5, lwd = 2, col = "blue")

# CI bars
arrows(x0 = as.numeric(rate0), x1 = as.numeric(rate0),
       y0 = as.numeric(ci_low), y1 = as.numeric(ci_up),
       angle = 90, code = 3, len = 0.08, col = "blue")

legend("topleft",
       lty = c(2, NA), pch = c(NA, 1),
       lwd = c(3, 2), bty = "n",
       col = c("maroon", "blue"), cex = 1.1,
       legend = c("Expected with proportional effect", "Grouped patients"))




#***Plot for presentation**
#**Forest plot*

# For final Presentation (modified to only have Labels, OR and plot)
png(
  "DIG_presn_forest.png",
  width  = 1800,
  height = 1200,
  res    = 220
)

par(mar = c(4, 4, 1, 2))

DIG_presn_forest <- forest(
  res,
  xlim    = c(-1.5, 1.5),
  at      = log(c(0.5, 1, 2)),
  alim    = c(log(0.2), log(2.5)),
  atransf = exp,
  slab    = data.subgroups$name,         
  rows    = c(1:2, 4:5, 7:8, 10:13, 15), 
  xlab    = "Odds Ratio",
  mlab    = "",
  psize   = 1.4,
  lwd     = 1.8,
  colout = c(rep("black", 10), 
             "#a80050"),
  addfit = FALSE,
  cex     = 1
)

#Subgroup 
text(-1.5, 14, "Baseline risk (quartiles)", pos = 4, font = 2)
text(-1.5,  9, "NYHA Class",               pos = 4, font = 2)
text(-1.5,  6, "Age",                       pos = 4, font = 2)
text(-1.5,  3, "Ejection Fraction",         pos = 4, font = 2)

#Title 
text(-0.4, 17, "DIG Trial", pos = 4, font = 2)

DIG_presn_forest

dev.off()

png(
  "DIG_report_forest.png",
  width  = 1800,
  height = 1200,
  res    = 220
)

par(mar = c(4, 4, 1, 2))

DIG_report_forest <- forest(
  res,
  addfit  = FALSE,
  xlim    = c(-8, 2.5),
  at      = log(c(0.5, 1, 2)),
  alim    = c(log(0.2), log(2)),
  atransf = exp,
  ilab    = cbind(
    data.subgroups$tn,
    data.subgroups$tevent,
    data.subgroups$cn,
    data.subgroups$cevent
  ),
  ilab.xpos = c(-5, -4, -3, -2),
  slab      = data.subgroups$name,
  rows      = c(1:2, 4:5, 7:8, 10:13, 15),
  xlab      = "Odds Ratio",
  mlab      = "",
  psize     = 1.2,
  lwd       = 1.5,
  colout    = c(rep("black", 10), "#a80050"),
  cex       = 0.9
)

text(-8, 14, "Baseline risk (quartiles)", pos = 4, font = 2)
text(-8,  9, "NYHA Class",               pos = 4, font = 2)
text(-8,  6, "Age",                      pos = 4, font = 2)
text(-8,  3, "Ejection Fraction",        pos = 4, font = 2)

text(
  c(-5, -4, -3, -2, 3),
  17,
  c(
    "Digoxin",
    "Events",
    "Placebo",
    "Events",
    "OR [95% CI]"
  ),
  font = 2
)

text(-8, 18, "DIG Trial", pos = 4, font = 2)

DIG_report_forest

dev.off()

#**Absolute benefit vs baseline risk plot*
#**Preparation*
library(ggplot2)
library(microshades)
hist_df <- data.frame(
  risk = baseline_risk_dist
)

curve_df <- data.frame(
  risk = xp,
  prop = p1exp
)

spline_df <- data.frame(
  risk = xp,
  spline = plogis(logxp0) - plogis(logxp0 + beta_tx)
)

group_df <- data.frame(
  risk = as.numeric(rate0),
  benefit = as.numeric(ratediff),
  lower = as.numeric(ci_low),
  upper = as.numeric(ci_up)
)


DIG_plot <- ggplot() +
  
  # Baseline risk distribution 
  geom_density(
    data    = hist_df,
    mapping = aes(
      x = risk,
      y = ..scaled..* 0.015   # FIX 3: ..scaled.. → after_stat(scaled)
    ),
    fill   = "grey80",
    colour = NA,
    alpha  = 0.4
  ) +
  
  # Spline model
  geom_line(
    data      = spline_df,
    mapping   = aes(x = risk, y = spline),
    linewidth = 1.7,
    colour    = "darkblue"
  ) +
  
  # Risk-quartile point estimates
  geom_point(
    data    = group_df,
    mapping = aes(x = risk, y = benefit),
    size    = 3,
    colour  = "#a80050"
  ) +
  
  # 95% CI bars
  geom_errorbar(
    data    = group_df,
    mapping = aes(x = risk, ymin = lower, ymax = upper),
    width   = 0.005,
    alpha   = 0.6,
    colour  = "grey80"
  ) +
  
  # Zero-benefit reference line 
  geom_hline(yintercept = 0, linetype = "dotted") +
  
  # Axis limits
  coord_cartesian(
    xlim = c(0, 0.7),
    ylim = c(-0.07, 0.07)
  ) +
  
  # Labels
  labs(
    title = "Absolute Benefit of Digoxin Across Baseline Mortality Risk",
    x     = "Baseline risk",
    y     = "Benefit by Digoxin (absolute risk difference)"
  ) +
  
  # Theme
  theme_classic(base_size = 17) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.ticks = element_blank(),
    axis.line  = element_line(colour = "grey60", linewidth = 0.5)
  ) +
  
  # annotations without legend needed
  annotate(
    "text",

    x = 0.43, y = -0.012,
    label    = "Proportional Effect",
    colour   = "darkblue",
    hjust    = 0,
    size     = 6,
    fontface = "bold"
  ) 

DIG_plot

# saving
saveRDS(DIG_plot, "dig_benefit_baseline_risk_plot.rds")




