
source("sim_utilities.R")

# function for classifying the process
classify_m <- function(mhat, tol = 0.09) {
  if (abs(mhat - 1) <= tol) return("critical")
  if (mhat < 1) return("subcritical")
  "supercritical"
}

# mc function
mc_sim_gw <- function(R = 500,
                      ngen = 5,
                      offspring = "poisson",
                      true_param,
                      tol = 0.05,
                      support = "unknown",
                      dirichlet_prior = "neutral",
                      alpha = 1,
                      g0 = "po",
                      lambda_g0 = 0.63,
                      prob_g0 = 0.29,
                      sample.dp = FALSE,
                      C = 50,
                      seed = NULL) {

  #scenario <- match.arg(scenario)
  if (!is.null(seed)) set.seed(seed)

  ## true mean m
  #m_true <- switch(offspring,
                   # "poisson"   = true_param,
                   # "geometric" = (1 - true_param) / true_param,
                   # stop("offspring not implemented yet")
  
  
  true_m <- m_true(
  offspring   = offspring,
  true_param  = true_param,
  prob_finite = prob_finite,
  xm          = xm,
  alpha       = alpha
)


  true_class <- classify_m(true_m)

  ## storage
  est_names <- c("mle", "heyde_mean","heyde_median","dir_mean", "dir_median", "dp_mean", "dp_median")
  if(sample.dp){est_names <- c("mle", "heyde_mean","heyde_median","dir_mean", "dir_median", "dp_mean", "dp_median","dp_mean.mc","dp_median.mc")}
  est_mat   <- matrix(NA, nrow = R, ncol = length(est_names))
  colnames(est_mat) <- est_names

  correct_class <- matrix(FALSE, nrow = R, ncol = length(est_names))
  colnames(correct_class) <- est_names
  
  N_tot_vec <- numeric(R)

  # mc loop
  for (r in 1:R) {

    ## simulate GW
    sim <- gw.sim(
      ngen = ngen,
      offspring = offspring,
      lambda = if (offspring == "poisson") true_param else NULL,
      p      = if (offspring == "geometric") true_param else NULL,
        xm     = if (offspring == "pareto") true_param$xm else NULL,
  alpha  = if (offspring == "pareto") true_param$alpha else NULL,
    prob_finite = if (offspring == "finite") true_param else NULL 
    )

    Zc <- sim$Zc
    N_tot_vec[r] <- sim$N_tot

    ## estimate m
    est <- estimate_m(
      data = Zc,
      support = support,
      dirichlet_prior = dirichlet_prior,
      C = C,
      alpha = alpha,
      g0 = g0,
      lambda_g0 = lambda_g0,
      prob_g0 = prob_g0,
      sample.dp = sample.dp
    )

    mhat <- c(
      mle        = est$m_estimates["mle", "mean"],
      heyde_mean        = est$m_estimates["heyde", "mean"],
      heyde_median        = est$m_estimates["heyde", "median"],
      dir_mean   = est$m_estimates["Dirichlet", "mean"],
      dir_median = est$m_estimates["Dirichlet", "median"],
      dp_mean    = est$m_estimates["DP", "mean"],
      dp_median  = est$m_estimates["DP", "median"]
    )
    if(sample.dp){
    	dp_mean.mc = est$m_estimates["DP.mc", "mean"] 
    	dp_median.mc = est$m_estimates["DP.mc", "median"]
    	mhat <- c(mhat, dp_mean.mc, dp_median.mc)
    	}

    est_mat[r, ] <- mhat

    for (j in seq_along(mhat)) {
      correct_class[r, j] <-
        classify_m(mhat[j], tol) == true_class
    }
  }
  
  N_tot_mc_mean <- mean(N_tot_vec)
  N_tot_mc_sd   <- sd(N_tot_vec)


  list(
    N_tot_mc_mean <- N_tot_mc_mean,
    true_m = true_m,
    scenario = true_class,
    classification_accuracy = colMeans(correct_class),
    mc_mean = colMeans(est_mat),
    mc_sd   = apply(est_mat, 2, sd)
  )
}
# example with poisson offspring
lambda_sub  <- 0.9
lambda_crit <- 1
lambda_sup  <- 1.1

poi_res_sub <- mc_sim_gw(
  R = 500,
  ngen = 3,
  offspring = "poisson",
  true_param = lambda_sub,
  g0 = "po",
  prob_g0 = 0.29,
  lambda_g0 = 0.63,sample.dp=TRUE,
  seed = 123
)
poi_res_sub

## ===============================
## Geometric offspring
## ===============================

# target means (same as Poisson)
m_sub  <- 0.9
m_crit <- 1.0
m_sup  <- 1.1
m_sup2  <- 1.5

# corresponding geometric parameters
p_sub  <- 1 / (m_sub  + 1)
p_crit <- 1 / (m_crit + 1)
p_sup  <- 1 / (m_sup  + 1)
p_sup2  <- 1 / (m_sup2  + 1)

# Monte Carlo simulation 
res_geo_sup <- mc_sim_gw(
  R = 500,
  ngen = 3,
  offspring = "geometric",
  true_param = p_sup,
  g0 = "geo",
  prob_g0 = 0.29, sample.dp=TRUE,   
  seed = 123
)

res_geo_sup

## ===============================
## Pareto offspring
## ===============================


# Monte Carlo simulation (sub critico m = 0.9)
xm=0.6;alpha=3
xm=0.6;alpha=6
(alpha*xm)/(alpha-1)
res_par_sub <- mc_sim_gw(
  R = 500,
  ngen = 3,
  offspring = "pareto",
  true_param = list(xm = xm, alpha = alpha),
  g0 = "geo", sample.dp=FALSE,         
  lambda_g0 = 0.63,
  #seed = 123
)
res_par_sub



# target means
m_sup  <- 1.1

# Pareto parameters (xm fixed)
xm <- 1

# Pareto mean exists only for m > 1
alpha_sup <- m_sup / (m_sup - 1)

# Monte Carlo simulation (supercritico m =1.1)
pareto_res_par_sup <- mc_sim_gw(
  R = 500,
  ngen = 3,
  offspring = "pareto",
  true_param = list(xm = xm, alpha = alpha_sup),
  g0 = "geo",  sample.dp=FALSE,        # oppure "geo" se vuoi stressare ancora di più
  seed = 123
)

pareto_res_par_sup

# Monte Carlo simulation (supercritico m = Inf)
pareto_res_par_sup_minf <- mc_sim_gw(
  R = 500,
  ngen = 3,
  offspring = "pareto",
  true_param = list(xm = xm, alpha = 0.8),
  g0 = "po",  sample.dp= T,         
  #lambda_g0 = 1,
  seed = 123
)

pareto_res_par_sup_minf



## ===============================
## Finite offspring (unknown support)
## ===============================

# true offspring distribution (finite but with heavy upper tail)
prob_finite_1<-  c(  0.4,  0.3,  0.01, 0.01, 0.01, 0.01, 0.01, 0.01,  0.01, 0.03)# 1.15
prob_finite_1<-  c(  0.45,  0.4,  0.01, 0.01, 0.01, 0.01, 0.01, 0.01,  0.01, 0.01)# 0.9
prob_finite_1<-  c(  0.4,  0.5,  0.01, 0.01,0.01, 0.01,  0.01, 0.01)# 0.8
prob_finite_1 <- c(0.5,0.2,0.2,0.1)# sub1
prob_finite_1 <- c(0.30, 0.60, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01) #0.9 sub 2
prob_finite_1 <- c(0.36, 0.31, 0.2,0.13)# 1.1 sup 1
prob_finite_1 <- c(0.2,0.3,0.35,0.2) # 1. 5  sup 2
prob_finite_1 <- c(0.1,0.1,0.25,0.2, 0.30) # 2.25  sup 3
prob_finite <- prob_finite_1 / sum(prob_finite_1)



# true mean
m_true_finite <- sum((0:(length(prob_finite)-1)) * prob_finite)
m_true_finite

res_finite_sub1_k.ignoto_n3 <- mc_sim_gw(
  R = 500,
  ngen = 5,
  offspring = "finite",
  true_param = prob_finite,
  support = "unknown",
  dirichlet_prior = "neutral",
  g0 = "po",          
  lambda_g0 = 0.63, sample.dp=T,     # prior neutrale
  C = 40,
  seed = 123
)

res_finite_sub1_k.ignoto_n5

