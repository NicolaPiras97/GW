### load packages
#require(LaplacesDemon)# for dirichlet
require(extraDistr)

# Galton–Watson process with flexible offspring distribution
gw.sim <- function(ngen,
                   offspring = c("poisson","geometric","nbinom","pareto","finite"),
                   lambda = NULL,# par poisson
                   p = NULL,#  par delle geom
                   size = NULL,# par binom 1
                   prob_binom = NULL, # par binom 2
                   prob_finite = NULL,# par offspring finita
                   xm = NULL,
                   alpha = NULL,
                   q = 0.99,M_max=1000){

  offspring <- match.arg(offspring)

  Z <- integer(ngen)

  ## determine M (max offspring size for bookkeeping)
  if(offspring == "geometric"){
    if(is.null(p))
      stop("p must be provided for geometric offspring")

    M <- qgeom(q, p)
    r_offspring <- function(n) rgeom(n, p)

  } else if(offspring == "poisson"){
    if(is.null(lambda))
      stop("lambda must be provided for poisson offspring")

    M <- qpois(q, lambda)
    r_offspring <- function(n) rpois(n, lambda)

  } else if(offspring == "nbinom"){
    if(is.null(size) || is.null(prob_binom))
      stop("size and prob must be provided for negative binomial offspring")

    M <- qnbinom(q, size = size, prob = prob_binom)
    r_offspring <- function(n) rnbinom(n, size = size, prob = prob_binom)

  }  else if(offspring == "pareto"){
  # if (is.null(xm) || is.null(alpha))
    # stop("xm and alpha must be provided for Pareto offspring")

  # ## Pareto quantile (discretized)
  # M <- floor(xm / (1 - q)^(1 / alpha))

  # r_offspring <- function(n){
    # y <- xm / runif(n)^(1 / alpha)
    # floor(y)
        if(is.null(xm) || is.null(alpha))
      stop("xm and alpha must be provided for Pareto offspring")

    M_raw <- floor(xm / (1 - q)^(1 / alpha))
    M     <- min(M_raw, M_max)

    r_offspring <- function(n){
      y <- xm / runif(n)^(1 / alpha)
      pmin(floor(y), M)
  }
}
else if(offspring == "finite"){
    if(is.null(prob_finite))
      stop("prob must be provided for finite offspring distribution")

    if(abs(sum(prob_finite) - 1) > 1e-8)
      stop("prob must sum to 1")

    M <- length(prob_finite) - 1
    r_offspring <- function(n){
      sample(0:M, size = n, replace = TRUE, prob = prob_finite)
    }
  }

  ## complete data matrix
  Zc <- matrix(0, nrow = ngen, ncol = M + 1)

  ## initialization
  Z[1] <- 1
  zc1 <- r_offspring(1)
  zc1 <- min(zc1, M)
  Zc[1, zc1 + 1] <- 1

  ## GW recursion
  for(i in 2:ngen){
    Zi <- sum((0:M) * Zc[i - 1, ])
    offsprings <- r_offspring(Zi)
    offsprings <- pmin(offsprings, M)
    offsprings2 <- factor(offsprings, levels = 0:M)
    Zc[i, ] <- as.numeric(table(offsprings2))
    Z[i] <- Zi
  }
#  extinct_immediately <- (Z[2] == 0) 
N_tot <- sum(Zc)

  return(list(Z = Z, Zc = Zc, M = M, N_tot = N_tot))
 # return(list(Z = Z, Zc = Zc, M = M, extinct_immediately = extinct_immediately)) 
 
 # sim <- ...
 # if(sim$extinct_immediately) next
}


complete.data <- function(x, k, max_iter = 1e5) {
  # x: vector of incomplete data (Z0, Z1, ..., Z_n)
  # k: truncation level (maximum offspring count)
  # max_iter: safety cap for accept-reject
  
  n_gen <- length(x) - 1
  if (n_gen < 1) stop("x must contain at least two generations")

  # initialize output matrix: rows = generations, cols = offspring counts 0,...,k
  T_mat <- matrix(0, nrow = n_gen, ncol = k + 1)
  colnames(T_mat) <- 0:k
  rownames(T_mat) <- paste0("gen_", 0:(n_gen - 1))

  # default offspring probabilities (flat prior)
  pi <- rep(1 / (k + 1), k + 1)

  for (i in 1:n_gen) {
    Zi <- x[i]
    Zi1 <- x[i + 1]

    if (Zi == 0) {
      # no individuals, all counts remain zero
      next
    }

    accepted <- FALSE
    iter <- 0

    while (!accepted && iter < max_iter) {
      iter <- iter + 1

      # propose from multinomial
      T_prop <- as.vector(rmultinom(1, size = Zi, prob = pi))

      # check constraint sum_j j T_ij = Z_{i+1}
      if (sum((0:k) * T_prop) == Zi1) {
        accepted <- TRUE
        T_mat[i, ] <- T_prop
      }
    }

    if (!accepted) {
      stop(paste(
        "Accept-reject failed at generation", i - 1,
        "after", max_iter, "iterations. Consider increasing k."
      ))
    }
  }

  return(T_mat)
}

estimate_m_gibbs <- function(
  Z,                    
  k,                    
  niter = 5000,
  burnin = 1000,
  prior = c("dirichlet", "dp"),
  a,                    
  alpha = NULL,
  g0,
  lambda_g0 = NULL,
  prob_g0 = NULL,
  sample.dp = TRUE
) {

  prior <- match.arg(prior)

  ngen <- length(Z) - 1
  C <- k + 1

  ## ------------------------------------------------------------
  ## Storage
  ## ------------------------------------------------------------
  m_chain <- numeric(niter)

  ## ------------------------------------------------------------
  ## Initialization
  ## ------------------------------------------------------------
  # inizializza p^(0) (Dirichlet o DP mean)
  if (prior == "dirichlet") {
    p_curr <- rep(1 / C, C)
  } else {
    p_curr <- rep(1 / C, C)   # G^(0) discreto troncato
  }

  ## ------------------------------------------------------------
  ## Gibbs sampler
  ## ------------------------------------------------------------
  for (l in seq_len(niter)) {

    ## ================================
    ## STEP 1: Imputation of T_ij
    ## ================================
    Zc <- matrix(0, nrow = ngen, ncol = C)

    for (i in seq_len(ngen)) {

      Zi   <- Z[i]
      Zip1 <- Z[i + 1]

      ## accept–reject multinomial
      repeat {
        Ti <- as.numeric(rmultinom(1, size = Zi, prob = p_curr))
        if (sum((0:k) * Ti) == Zip1) break
      }

      Zc[i, ] <- Ti
    }

    ## ================================
    ## STEP 2: Posterior update
    ## ================================

    if (prior == "dirichlet") {

      ## --- usa ESATTAMENTE il tuo codice ---
      tot.nj <- colSums(Zc)
      apost  <- a + tot.nj

      p_dir_mean <- apost / sum(apost)
      m_curr <- sum((0:k) * p_dir_mean)

      p_curr <- p_dir_mean

    } else {

      ## --- DP posterior mean ---
      G0n <- funG0n(
        g0 = g0,
        alpha = alpha,
        data = Zc,
        lambda_g0 = lambda_g0,
        prob_g0 = prob_g0
      )

      m_curr <- sum((0:k) * G0n)
      p_curr <- G0n
    }

    m_chain[l] <- m_curr
  }

  ## ------------------------------------------------------------
  ## Output
  ## ------------------------------------------------------------
  list(
    m_chain = m_chain[-seq_len(burnin)],
    m_mean  = mean(m_chain[-seq_len(burnin)]),
    m_sd    = sd(m_chain[-seq_len(burnin)])
  )
}


  m_true <- function(offspring,
                   true_param = NULL,# par di pois o geom
                   prob_finite = NULL,
                   xm = NULL,
                   alpha = NULL) {
  switch(offspring,
    "poisson" = {
      true_param
    },
    "geometric" = {
      (1 - true_param) / true_param
    },
    "finite" = {
      if (is.null(prob_finite))
        stop("prob_finite must be provided for finite offspring")
      sum((0:(length(prob_finite) - 1)) * prob_finite)
    },
  "pareto" = {
      if (is.null(true_param$xm) || is.null(true_param$alpha))
        stop("true_param must contain xm and alpha for Pareto offspring")

      xm    <- true_param$xm
      alpha <- true_param$alpha

      if (alpha <= 1)
        Inf
      else
        alpha * xm / (alpha - 1)
    },
    stop("offspring not implemented yet")
  )
}



# the function below estimates the support size (with complete data only)
#nc.est<-function(x){max(which(apply(x,2,sum)>0))}
nc.est <- function(x){
  idx <- which(colSums(x) > 0)
  if(length(idx) == 0) return(0)
  max(idx)
}


# # OLD funG0n calculates the DP posterior measure
# funG0n<-function(g0, alpha, data){
	# N=sum(data)
	# if(g0=="po"){
# lambda <-1#0.6954 # median-neutral
 # G0 <- dpois(0:(ncol(data)-1),lambda)
# G0n<- 1/(alpha + N)*apply(data,2,sum) + alpha/(alpha + N)*G0}
 	# if(g0=="geo"){
# lambda <-1#0.2929 # median-neutral
 # G0 <- dgeom(0:(ncol(data)-1),lambda)
 # G0n <- 1/(alpha + N)*apply(data,2,sum) + alpha/(alpha + N)*G0}
 # return(G0n)
# }


#NEW funG0n calculates the DP posterior measure
funG0n <- function(g0 = c("po", "geo"),
                   alpha,
                   data,
                   lambda_g0 = 0.695, # poisson # 0.6954 # median-neutral for poisson
                   prob_g0 = 0.292    ){# par geom #0.2929 # median-neutral for geometric

  g0 <- match.arg(g0)
  N  <- sum(data)

  # offspring values j = 0, ..., J
  j_vals <- 0:(ncol(data) - 1)
  freq   <- colSums(data)

  # base measure G0
  if(g0 == "po"){
    G0 <- dpois(j_vals, lambda_g0)
  }

  if(g0 == "geo"){
    G0 <- dgeom(j_vals, prob_g0)
  }

  # normalize G0 on the implicit (possibly truncated) support
  G0 <- G0 / sum(G0)

  # posterior base measure G0^n
  G0n <- (freq + alpha * G0) / (N + alpha)

  return(G0n)
}

#mle.p <- function(x){phat<- apply(x, 2, sum)/sum(x);return(phat)}
# funzione per ottenere mle of m (da complete data x)
mle.m <- function(x){(sum(x)-1)/sum(x[-nrow(x),])} 

# funzione per calcolare mediana della beta
beta.median <- function(a, b){
  if(a <= 0 || b <= 0) return(0)
  qnsbeta(0.5, a, b)
}


# function rdp for sampling from DP(alpha, G0)
# function rdp(size, alpha, data) with G0= Poisson, Geom or discrete on 0...k
# if data<-matrix(0,C,nrow=1) it samples from G(alpha, g0) ; if data <- Zc its samples from posterior G | Zc
# rdp( size = 10, alpha = alpha, g0 = g0, lambda=lambda_g0,prob=prob_g0,data = Zc)
rdp<-function(size,alpha,g0, lambda, prob ,data){
	n = sum(data) 
	alpha.post = alpha + n
	ifelse(n==0, p.epdf<-1 , p.epdf<- colSums(data)/n )
# truncation index
trunc=5000
 # stick breaking for weigths  p_j          
P=vector()
B=rbeta(trunc,1,alpha.post)
P[1]<-B[1]
for (j in 2:trunc){
P[j]=B[j]*prod(1-B[1:(j-1)])  
              }
Pcum=c(0,cumsum(P))
# generate V = sample from G | x
# sample from epdf with prob n/(n+a) and from G0 with prob a/(n+a)
u <- rbinom(trunc,  1,prob = n/(n+alpha))
V <- vector()
for (i in 1:trunc){
	# Go Poisson:
	if(g0=="po") {
		V[i] = ifelse(   u[i] ==1, 
		sample(0:(ncol(data)-1),size=1, prob=p.epdf), 
		rpois(1,lambda) )}
# Go Geometric:
	if(g0=="geo") {
   V[i] = ifelse(   u[i] ==1, 
   sample(0:(ncol(data)-1),size=1, prob=p.epdf),
    rgeom(1,prob) )}
# # G0 flat on 0...k:
 # if(g0=="finite.flat"){
 # probs<-rep(1/(ncol(data)),ncol(data))
 # V[i] = ifelse(   u[i] ==1, sample(0:(ncol(data)-1), size=1, prob=p.epdf),
                                    # # sample(0:(ncol(data)-1), size=1,prob=probs) )    
 # if(k+1==1){V[i] = ifelse(   u[i] ==1, 0,0 )  } }
 }
# Doss' algorithm (1993)
# initialize x = sample vector
x<-vector()
for (t in 1:size){ # start cycle filling x
    U=runif(1,0,1)
    for (i in 1:trunc){
        if (   Pcum[i]<U & U<=Pcum[1+i]   )
        {  x[t]=V[i]   }
                   }
              } # fine ciclo
    return(x)	
	}



 ## NEW funzione di stima	    
 	    estimate_m <- function(data,
                       support = c("known", "unknown"),
                       dirichlet_prior = c("flat", "neutral"),
                       C = 100,
                       alpha = 1,
                       g0 = "po",
                       lambda_g0 = 0.635,
                       prob_g0 = 0.2929,
                       prob_finite = NULL,
                       sample.dp = FALSE) {

  support <- match.arg(support)
  dirichlet_prior <- match.arg(dirichlet_prior)

  Zc <- data

## ------------------------------------------------------------
## Support handling
## ------------------------------------------------------------
if (support == "known") {
  nc <- ncol(Zc)
  C  <- nc
} else {
  nc <- nc.est(Zc)            # supporto osservato
  buffer <- 10                # buffer per la coda
  if (C <= nc) {
    C <- nc + buffer
  }
  if (C > ncol(Zc)) {
    Zc0 <- matrix(0, nrow = nrow(Zc), ncol = C - ncol(Zc))
    Zc <- cbind(Zc, Zc0)
  }
}


  k <- nc - 1   # maximum offspring value

  ## ------------------------------------------------------------
  ## Dirichlet prior specification
  ## ------------------------------------------------------------
  if (dirichlet_prior == "flat") {

    a <- c(rep(1, nc), rep(0, C - nc))

  } else if (dirichlet_prior == "neutral") {

    if (nc > 2) {
      a <- c(1,
             rep(0, nc - 2),
             log(2) / log(k),
             rep(0, C - nc))
    } else {
      a <- c(rep(1, nc), rep(0, C - nc))
    }
  }
  ## ----------------------------------------
  ## Maximum likelihood estimator
  ## ---------------------------------------
  #est.p.mle <- mle.p(Zc)
  #m_mle <- sum(0:(ncol(Zc)-1)*est.p.mle)
  m_mle <- mle.m(Zc)
  
  ## -----------------------------------------
  ## Heyde estimator (assumes poisson offspring and a non informative prior)
  ## ------------------------------------------
  # bayes estimator with improper prior (C. Heyde)
 # m ~ m^-1 (only for Poisson GW)
    yn <- sum(Zc)
    y0 <-1
   yn_minus_1 <- sum(Zc)-sum(Zc[nrow(Zc),])
   m_post_heyde <- (yn-y0)/yn_minus_1
      # media a posteriori
    m_heyde1 <- m_post_heyde
    # P(m | x >1) (normal approximation; see heyde pag. 285)
    pm_greater_1 <- 1-pnorm((yn_minus_1 - (yn-y0))/sqrt(yn-y0)   )
    m_heyde2  <- ifelse(pm_greater_1 >0.5, 1.5, 0.5)
    
  ## ------------------------------------------------------------
  ## Dirichlet posterior
  ## ------------------------------------------------------------
  tot.nj <- colSums(Zc)
  apost  <- a + tot.nj

  p_dir_mean <- apost / sum(apost)

  p_dir_median_raw <- sapply(1:C, function(j) {
    beta.median(apost[j], sum(apost) - apost[j])  
  })
  if (sum(p_dir_median_raw) > 0) {
  p_dir_median <- p_dir_median_raw / sum(p_dir_median_raw)
} else {
  p_dir_median <- p_dir_mean   # fallback sensato
}

  m_dir_mean   <- sum((0:(C - 1)) * p_dir_mean)
  m_dir_median <- sum((0:(C - 1)) * p_dir_median)

  ## ------------------------------------------------------------
  ## DP posterior (G0^n)
  ## ------------------------------------------------------------
  G0n <- funG0n(g0 = g0,
                alpha = alpha,
                data = Zc,
                lambda_g0 = lambda_g0,
                prob_g0 = prob_g0)

  m_dp_mean <- sum((0:(C - 1)) * G0n)

  N <- sum(Zc)

  p_dp_median_raw <- sapply(1:C, function(j) {
    beta.median((alpha + N) * G0n[j],
                (alpha + N) * (1 - G0n[j]))
  })
  
  ## normalize + fallback
if (sum(p_dp_median_raw) > 0) {
  p_dp_median <- p_dp_median_raw / sum(p_dp_median_raw)
} else {
  p_dp_median <- G0n   # fallback naturale: DP posterior mean
}

  #m_dp_median <- (0:(C - 1))[min(which(cumsum(G0n) >= 0.5))] # mediana di g0
  m_dp_median <- sum((0:(C - 1)) * p_dp_median)
  
  # mean and median from sample of posterior DP
  if (sample.dp){
  	dp.sample <- rdp( size = 100, alpha = alpha, g0 = g0, lambda=lambda_g0,prob=prob_g0,data = Zc)
  	m_dp_mean.mc <- mean(dp.sample)
  	m_dp_median.mc <-median(dp.sample)
  }

  ## ------------------------------------------------------------
  ## Output
  ## ------------------------------------------------------------
  results <- matrix(
    c(m_mle, m_mle,
        m_heyde1, m_heyde2,
       m_dir_mean, m_dir_median,
      m_dp_mean,  m_dp_median),
    nrow = 4, byrow = TRUE
  )
  
if (sample.dp) {
	results <- rbind(results,c(m_dp_mean.mc, m_dp_median.mc))
  rownames(results) <- c("mle","heyde","Dirichlet", "DP", "DP.mc")
} else {
  rownames(results) <- c("mle","heyde","Dirichlet", "DP")
}


  colnames(results) <- c("mean", "median")
 # rownames(results) <- c("mle","heyde","Dirichlet", "DP")

  list(
    m_estimates = results,
    support_true = if (!is.null(prob_finite)) length(prob_finite) else NA,
    support_estimated = nc,
    dirichlet_prior = dirichlet_prior
  )
}

# example offsping poisson
# set.seed(321)
# data <- gw.sim(offspring = "poisson", lambda = 0.9, ngen = 3)
# Z <- data[[1]]
# Zc <- data[[2]]
eZc <- complete.data(Z,10)

