```{r}
data(gusto)   # loads original GUSTO dataset

gusto_raw <- gusto   # always keep this untouched



```

```{r}
library(dplyr)
gusto <- gusto_raw %>% filter(tx %in% c("SK","tPA")) %>%
  select(day30, tx, age, Killip, sysbp, pulse, pmi, miloc, sex) %>% mutate(
    day30 = factor(day30, levels = c(0, 1), labels = c("No", "Yes")), across(c(tx, Killip, pmi, miloc, sex), as.factor))


head(gusto)



```
```{r}
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



```

```{r}
tab1 <- table1( ~ age + Killip + sysbp + pulse + pmi + miloc + sex + day30 | tx, data = gusto, test = TRUE)

tab1




```


**PATH**
  First we look at the primary outcome
```{r}

table1(~day30|tx, data=gusto)

```

```{r}

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
```


```{r}
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
```




```{r}
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

```
```{r}
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

```

The Treatment effect appears constant across all covariates. 



```{r}
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