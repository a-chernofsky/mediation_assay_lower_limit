############################################################################################################
#
# MCEM
# created: 06/07/2020
#
#
############################################################################################################


# load libraries ----------------------------------------------------------

library(truncnorm)
library(tidyverse)

# functions used in the mcem function -------------------------------------

#machine precision
LOGEPS <- log(.Machine$double.eps / 2)

#evluates log(1+exp(x)) in a numerically stable format
log1pe <- function (x) {
  l <- ifelse(x > 0, x, 0)
  x <- ifelse(x > 0, -x, x)
  ifelse(x < LOGEPS, l, l + log1p(exp(x)))
}

#inverse logit or expit function
expit <- function(eta)1/(1+exp(-eta))

#function for evaluating the posterior distribution M | Y, C, theta^{(t)}
bayes_prob <- function(alpha, beta, y, m, c, limit, sigma) {
  y_design <- cbind(rep(1, length(m)), m, rep(c, length(m)))
  m_design <- c(1,c)
  dbinom(y, 1, expit(y_design %*% beta))
  dtruncnorm(m, b = limit, mean = m_design %*% alpha, sd = sigma)
}


# mcem function for calculating the indirect effect -----------------------

mcem <- function(Y, M, C, delta, lod, shift, initial){
  limit <- lod
  #collect inputs into a dataframe
  data <- data.frame(id = seq(1, length(Y)), Y, M, C, delta)
  
  #values with delta = 0 are below the assay limit
  cen <- which(data$delta == 0)
  
  #intialize imputations to half limit of detection
  data$M0 <- ifelse(data$delta == 1, M, initial)
  
  #initialze alpha, sigma^2, beta
  mlm <- lm(M0 ~ C, data = data)
  alpha_est <- mlm$coef
  
  #estimate of sigmam
  sigmam <- sigma(mlm)
  
  #initialize outcome model estimates
  yglm <- glm(Y ~ M0 + C, data = data, 
              family = binomial())
  beta_est <- yglm$coef
  
  log_odds_init <- yglm$linear.predictors
  
  #initial log likelihood for convergence criteria
  loglik <- sum(dnorm(data$M0, mlm$fitted.values, sigmam, log = T)) +
    sum(data$Y * log_odds_init - log1pe(log_odds_init))
  
  #starting value for m samples
  J <- 500
  
  #split obs and cen data
  obsdat <- data[which(data$delta == 1), ] 
  
  #create a new mediator column to store imputed values
  cendat <- data[rep(cen, J), ]
  #initialize M0 column to fill with sampled M values
  cendat$M0 = NA
  
  #save vectors of censored variables
  t <- 1
  
  mvals <- seq(-10, limit, length.out = 1000)
  
  repeat{
    grid_sample <- rep(NA, length(cen)*J)
    
    for(i in 1:length(cen)){
      
      num <- bayes_prob(alpha = alpha_est, 
                        beta = beta_est, 
                        y = data$Y[cen][i],
                        m = mvals, 
                        c = data$C[cen][i], 
                        limit = limit,
                        sigma = sigmam)
      den <- sum(num)
      
      lik <- num/den
      grid_sample[(1+J*(i - 1)):(i*J)] <- sample(mvals, size = J, 
                                                 replace = T, prob = lik)
    }
    
    cendat$M0 <- grid_sample
    
    #update alpha, sigma^2, beta
    newdat <- rbind(obsdat, cendat)
    #create weights
    newdat$wts <- ifelse(newdat$delta == 1, 1, 1/J)
    #update M linear model
    mlm_new <- lm(M0 ~ C, data = newdat, weights = wts)
    #update sigmam
    sigmam_new <- sqrt(1/nrow(data) * 
                         sum(newdat$wts*(newdat$M0 - mlm_new$fitted.values)^2))
    
    #calculate log likelihood for convergence criteria
    yglm_new <- glm(Y ~ M0 + C, data = newdat,
                    weights = wts, family = binomial())
    
    #update log likelihood estimates
    loglik_new <- sum(newdat$wts*(dnorm(newdat$M0, mlm_new$fitted.values, 
                                        sigmam_new, log = T) + 
                                    newdat$Y * yglm_new$linear.predictors - 
                                    log1pe(yglm_new$linear.predictors)))
    
    #check for convergence
    if((abs((loglik_new - loglik)/loglik) < 1e-6)) break
    
    #if convergence not reached update values
    else{
      alpha_est <- mlm_new$coef
      sigmam <- sigmam_new
      loglik <- loglik_new
      beta_est <- yglm_new$coef
      t <- t + 1
    }
  }
  
  #repeat M0 the number of shifts and subtract each column by a shift
  M1 <- sweep(matrix(rep(newdat$M0, length(shift)), ncol = length(shift)), 2, shift)
  
  #log(odds) = b1 + b2 * M1 + b3 * C
  log_odds_final <- beta_est[1] + M1 * beta_est[2] + newdat$C * beta_est[3] 
  
  #phat estimate
  phat_final <- 1/(1+exp(-log_odds_final))
  
  #E[Y_0]
  ey0 <- mean(data$Y)
  
  #E[Y_{0, I = 1}]
  ey0i <- apply(phat_final, 2, 
                function(x) weighted.mean(x, w = newdat$wts))
  
  #indirect effect
  ey0i - ey0
  
}

