##############################################################
#
#
# Extrapolation
# created: 06/19/2020
#
#
#
###############################################################


library(tidyverse)
library(truncnorm)

extrapolation <- function(Y, M, C, delta, lod, shift, initial){
  
  data <- data.frame(id = seq(1, length(Y)), Y, M, C, delta)
  
  #save indices for censored and observed 
  cen <- which(data$delta == F)
  obs <- which(data$delta == T)
  
  n <- nrow(data)
  
  #standardization function
  z <- function(lod, mu, sigma){(lod - mu)/sigma}
  
  data.init <- data
  data.init$M[cen] <- initial 
  
  #initial values for coefficients for observed
  lm.cur <- lm(M ~ C, data = data.init)
  
  #initial value for sigma.sq
  sigma.cur <- sigma(lm.cur)
  
  #fitted mediator values for observed model
  mhat <- lm.cur$fitted.values
  
  
  #current censored standardized values
  z.cur <- z(lod, mhat[cen], sigma.cur)
  
  
  #calculate log likelihood 
  log.lik.cur <- sum(dnorm(data$M[obs], 
                           mean = mhat[obs], 
                           sd = sigma.cur, log = T)) + 
    sum(pnorm(z.cur, log.p = T))
  
  #initialize counter for iteration
  i = 1
  
  #iterative least squares
  repeat{
    
    #calculate update for sigma square
    num <- sum((data$M[obs] - mhat[obs])^2)
    den <- length(mhat[obs]) + sum(z.cur * (dnorm(z.cur)/pnorm(z.cur))) 
    
    sigma.new <- sqrt(num/den)
    
    #impute censored values
    z.new <- z(lod, mhat[cen], sigma.new)
    
    mimp <- mhat[cen] - (sigma.new * (dnorm(z.new))/pnorm(z.new))
    
    #create a new dataset to store imputed values
    data.new <- data 
    data.new$M[cen] <- mimp
    
    #fit lm on "new" complete data with imputed values
    lm.new <- lm(M ~ C, data = data.new)
    
    #fitted values based on new data
    mhat.new <- lm.new$fitted.values
    
    #calc new loglik to test convergence
    log.lik.new <- sum(dnorm(data$M[obs], 
                             mean = mhat.new[obs], 
                             sd = sigma.new, log = T)) + 
      sum(pnorm(z(lod, mhat.new[cen], sigma.new), log.p = T))
    
    #convergence criteria
    if(abs((log.lik.new - log.lik.cur)/log.lik.cur) < 1e-10){
      break
    } 
    #if convergence criteria not met set the new values to current values
    else{
      sigma.cur <- sigma.new
      mhat <- mhat.new
      log.lik.cur <- log.lik.new
      z.cur <- z.new
      i = i + 1
    }
  }
  
  Mest <- list(alpha = lm.cur$coeff, mhat = mhat, 
               sigma = sigma.cur)
  
  mu <- Mest$mhat
  sigma <- Mest$sigma
  alpha <- Mest$alpha
  J <- 1000
  
  obsdat <- data[obs,]
  obsdat$M0 <- obsdat$M
  
  cendat <- data[rep(cen, J), ]
  cendat$M0 <- NA
  
  cendat$M0 <- rtruncnorm(length(cen)*J, mean = mu[cen], sd = sigma, b = lod)
  
  newdat <- rbind(obsdat, cendat)
  
  beta.est <- glm(Y ~ M + C, family = binomial(), data = data, 
                  subset = which(delta == 1))$coeff
  newdat$wts <- ifelse(newdat$delta == 1, 1, 1/J)
  
  #estimate E[Y_{0, I=1}]
  M1 <- sweep(matrix(rep(newdat$M0, length(shift)), ncol = length(shift)), 2, shift)
  
  eta.final <- beta.est[1] + M1 * beta.est[2] + newdat$C * beta.est[3] 
  
  phat.final <- 1/(1+exp(-eta.final))
  
  ey0 <- mean(data$Y)
  
  ey0i <- apply(phat.final, 2, function(x) weighted.mean(x, w = newdat$wts))
  
  #indirect effect
  ey0i - ey0
}


