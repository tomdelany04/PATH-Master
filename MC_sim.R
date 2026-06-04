

#Propogate uncertainty in model using bootstraps

#################
#The simulation evaluates whether risk-based
#interaction testing can indirectly recover underlying HTE.
#########################



# grid of sample sizes and vector of interaction effects 

# For type I error do a grid parameters and sample 

# Expand_grid (n = c(1,1000) p = c(10,20))

# Get type one error rates <- geom_surface

#using ADEMP framework (Morris et al., 2019)
library(magick)
library(foreach)
library(doParallel)
library(ggplot2)

ADEMP <-image_read("ADEMP.webp")

print(ADEMP)


p = c(1:10)

n = n <- c(
  50, 75, 100, 150, 200,
  300, 500, 750,
  1000, 2000, 5000
)

(grid = expand.grid(p = p, n = n))


null_simulation <- function(n = 5000, p = 3) {

    ##########################
    # Covariates
    ##########################

    X <- matrix(
      rnorm(n * p),
      nrow = n,
      ncol = p
    )
    
    
    colnames(X) <- paste0("X", 1:p) #Names of X1, X2 ETC
    ##########################
    # Treatment assignment
    ##########################

    beta <- rep(0.5, p) # Coef vector
    
    
    
    Z <- rbinom(n, 1, 0.5)
    
    xb_null <- -2 -
      0.5 * Z +
      X %*% beta
    
    
  

    ##########################
    # NULL truth
    ##########################


    p_null <- plogis(xb_null)

    Y_null <- rbinom(n, 1, p_null)

    ##########################
    # Dataset
    ##########################

    df <- data.frame(
      X,
      Z,
      Y_null
    )

    
    
    
    
    ##########################
    # Train/test split
    ##########################

    train_id <- sample(
      seq_len(n),
      size = floor(n/2)
    )

    train_df <- df[train_id, ]
    test_df  <- df[-train_id, ]

    ##########################
    # Baseline risk model
    # TRAINING 
    ##########################

    
    
    covariates <- paste0("X", 1:p)
    paste(covariates, collapse = " + ")
    
    
    risk_formula <- as.formula( #R formula object
      paste(
        "Y_null ~",
        paste(covariates,
              collapse = " + ")
      )
    )
    
    
    risk_model <- glm(
      risk_formula,
      family = binomial(),
      data = train_df
    )

    ##########################
    # Predict LP in TEST 
    ##########################

    test_df$lp <- predict(
      risk_model,
      newdata = test_df,
      type = "link"
    )

    ##########################
    # Interaction models
    # TEST 
    ##########################

    model_main <- glm(
      Y_null ~ Z + lp,
      family = binomial(),
      data = test_df
    )

    model_interaction <- glm(
      Y_null ~ Z + lp + Z:lp,
      family = binomial(),
      data = test_df
    )

    ##########################
    # Likelihood ratio test
    ##########################
    lrt <- anova(
      model_main,
      model_interaction,
      test = "LRT"
    )
    pvalue <- lrt$`Pr(>Chi)`[2]
    return(pvalue)
}

cl <- makeCluster(4)

registerDoParallel(cl)

############################################################
# Run simulation study
############################################################

set.seed(2140)

results <- foreach(
  i = 1:nrow(grid),
  .combine = rbind
) %dopar% {
  
  ##########################
  # Current simulation condition
  ##########################
  
  current_n <- grid$n[i]
  
  current_p <- grid$p[i]
  
  ##########################
  # Monte Carlo repetitions
  ##########################
  
  pvalues <- replicate(
    1000,
    null_simulation(
      n = current_n,
      p = current_p
    )
  )
  
  ##########################
  # Empirical Type I error
  ##########################
  
  data.frame(
    n = current_n,
    p = current_p,
    pvalue = pvalues,
    reject = pvalues < 0.05
  )
}


summary_results <- aggregate(
  reject ~ n + p,
  data = results,
  mean
)

colnames(summary_results)[3] <- "type1_error"

saveRDS(summary_results, "summary_results.rds")

############################################################
# Stop cluster
############################################################

stopCluster(cl)





############################################################
# Type I Error Surface
############################################################

png(
  "Heatmap_MC.png",
  width = 1800,
  height = 1200,
  res = 220
)

null_plot <- ggplot(
  summary_results,
  aes(
    x = n,
    y = p,
    z = type1_error
  )
) +
  geom_contour_filled(aes(fill = after_stat(level_mid))) +
  geom_contour(color = "black") +
  geom_point(size = 2) +
  scale_fill_viridis_c(name = "Type I Error") +
  scale_x_log10() +
  theme_minimal(base_size = 15) +  
  theme(
    panel.grid = element_blank()
  ) +
  labs(
    title = "Type I Error Surface Under the Null",
    x = "Sample Size (log(n))",
    y = "Parameters (p)",
    fill = "Type I Error"
  ) 




print(null_plot)
dev.off()


saveRDS(null_plot, "type1_error_surface_plot.rds")


null_plot


# geom surface, 2d kernal density plot as alternative
#plot dots instead and interpolate between 

# discrete points are mapped to a grid
#Nearby values are connected
#Intermediate values are estimated locally
#Contour polygons are drawn for ranges of z


ggplot(summary_results,
       aes(x = n,
           y = p,
           z = type1_error)) +
  geom_contour_filled(aes(fill = after_stat(level_mid))) +
  scale_fill_viridis_c(name = "Type I Error") +
  theme_minimal()


