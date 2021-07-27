############################################################################################################
#
# Numerical Optimization
# created: 06/07/2020
#
#
############################################################################################################

library(tidyverse)
library(truncnorm)


noptim <- function(Y, M, C, delta, lod, shift, initial){
  
  #function to integrate for Y = 1 
  p1 <- function(mm, C){
    #linear predictor for mediator model
    meta <- a0 + a2*C
    #linear predictor for outcome model
    yeta <- b0 + b2* mm + b3 * C
    #f(y | a, m, c)*f(m | a, c)
    (1/(1+exp(-yeta))) * 
      dnorm(mm, meta, s)
  }
  # function to integrate for Y = 0
  p0 <- function(mm, C){
    #linear predictor for mediator model
    meta <- a0 + a2*C
    #linear predictor for outcome model
    yeta <- b0 + b2* mm + b3 * C
    (1/(1+exp(yeta))) * 
      dnorm(mm, meta, s)
  }
  # Y = 1 integrate function
  int1 <- function(C) integrate(p1, lower = -Inf, 
                                upper = lod,  C = C, stop.on.error = F)$value
  #vectorize function
  vint1 <- Vectorize(int1)
  
  # Y = 0 integrate function
  int0 <- function(C) integrate(p0, lower = -Inf, 
                                upper = lod, C = C)$value
  #vectorize integrate function
  vint0 <- Vectorize(int0)
  
  #loglik function
  ll <- function(par){
    a0 <<- par[1]
    a2 <<- par[2]
    s <<- par[3]
    b0 <<- par[4]
    b2 <<- par[5]
    b3 <<- par[6]
    
    meta_obs <- a0 + a2 * C
    yeta_obs <- b0 + b2 * M + b3 * C
    p <- exp(yeta_obs)/(1+exp(yeta_obs))
    
    l.obs <- sum(delta*(Y*log(p) + (1-Y)*log(1-p) + 
                          dnorm(M, meta_obs, s, log = T)), na.rm = T) 
    l.cen <- sum((1-delta)*(Y*log(vint1(C)) +
                              (1-Y)*log(vint0(C))), na.rm = T)
    -(l.obs + l.cen)
  }
  
  data <- data.frame(id = seq(1, length(Y)), Y, M, C, delta)
  
  cen <- which(delta == 0)
  obs <- which(delta == 1)
  
  #impute missing values with lod/2
  M.init <- M
  M.init[cen] <- initial
  
  
  #initial estimates for parameters
  mlm <- lm(M.init ~ C)
  sigmam <- sigma(mlm)
  yglm <- glm(Y ~ M.init + C, family = binomial())
  
  #initial values for parameters
  pars.init <- c(mlm$coefficients, sigmam, yglm$coefficients)
  
  #numerical optimization - constrain sigmam^2 >0
  fit <- optim(par=pars.init, fn=ll, method = "L-BFGS-B",
               lower = c(rep(-Inf, 2), 0.0001 , rep(-Inf,3)))
  
  #extract estimated mediator model coefficients and sigmam
  alpha.est <- fit$par[c(1,2)]
  sigma.est <- fit$par[3]
  beta.est <- fit$par[4:6]
  
  #simulate mediator for censored values
  J <- 1000
  
  obsdat <- data[obs,]
  obsdat$M0 <- obsdat$M
  
  cendat <- data[rep(cen, J), ]
  cendat$M0 <- NA
  
  mhat <- cbind(rep(1, length(cen)), C[cen]) %*% 
    alpha.est
  
  cendat$M0 <- rtruncnorm(length(cen)*J,
                          b = lod,
                          mean = mhat,
                          sd = sigma.est)
  
  newdat <- rbind(obsdat, cendat)
  newdat$wts <- ifelse(newdat$delta == 1, 1, 1/J) 
  
  M1 <- sweep(matrix(rep(newdat$M0, length(shift)), ncol = length(shift)), 2, shift)
  
  eta0i <- beta.est["(Intercept)"] + beta.est["M.init"]*M1 + beta.est["C"]*newdat$C
  phat0i <- 1/(1+exp(-eta0i))
  
  ey0i <- apply(phat0i, 2, function(x) weighted.mean(x, w = newdat$wts))
  ey0 <- mean(Y) 
  
  ey0i - ey0
}
