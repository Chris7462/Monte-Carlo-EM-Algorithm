rm(list=ls())
library("MASS")

set.seed(312)

# (a)
mwanza <- read.csv('mwanza.csv', header=TRUE)
mwanza$community <- as.factor(mwanza$community)
glm.fit <- glmmPQL(hiv~arm,random=~1|community,family='binomial',data=mwanza)
summary(glm.fit)

alpha.glm <- as.numeric(glm.fit$coefficients$fixed[1])
beta.glm <- as.numeric(glm.fit$coefficients$fixed[2])
sigma2.glm <- as.numeric(VarCorr(glm.fit)[1,1])
u.glm <- as.numeric(ranef(glm.fit)$"(Intercept)")

# (b)
# organize data
n <- length(unique(mwanza$community))
y <- vector("list",n)
x <- vector("list",n)

for ( i in 1:n ){
  index <- which(mwanza$community==i)
  y[[i]] <- mwanza$hiv[index]
  x[[i]] <- mwanza$arm[index]
}

## log-likelihood function ##
log.likelihood <- function(y,x,alpha,beta,u){
  pred <- alpha + beta*x + u 
  log.like <- sum(y*pred-log(1+exp(pred)))
  return(log.like)
}

## log-prior distribution ##
log.prior <- function(u,sigma2){
  u.log.prior <- dnorm(u, mean=0, sd=sqrt(sigma2), log=TRUE)
  return(u.log.prior)
}

## log-posterior distribution ##
log.posterior <- function(y,x,alpha,beta,u,sigma2){
  u.log.post <- log.likelihood(y,x,alpha,beta,u) + log.prior(u,sigma2) 
  return(u.log.post)
}

## Hasting-Metroplis algorithm ##
H_M <- function(y,x,alpha,beta,u,sigma2,size,burnin=500,intval=20){
  iter <- burnin+intval*size
  u.chain <- rep(0,iter)
  u.chain[1] <- u
  for (i in 1:(iter-1) ){
    # proposal function
    u.prop <- rnorm(1,mean=u.chain[i],sd=sqrt(sigma2))
    prob <- min(exp(log.posterior(y,x,alpha,beta,u.prop,sigma2)-log.posterior(y,x,alpha,beta,u.chain[i],sigma2)),1)
    if( runif(1) < prob ){
      u.chain[i+1] <- u.prop
    } else {
      u.chain[i+1] <- u.chain[i]
    }
  }
  # remove burnin sample
  u.chain <- u.chain[-(1:burnin)]
  # subsample
  u.sample <- u.chain[seq(1,iter-burnin,by=intval)]
  return(u.sample)
}

## EM function ##
Qm <- function(theta,y,x,u,w,size){
  alpha <- theta[1]
  beta <- theta[2]
  sigma2 <- theta[3]
  n <- length(y)
   
  Q <- 0
  for ( t in 1:size ){ 
    tmp.log <- 0
    for ( i in 1:n ){
      tmp.log <- tmp.log + log.posterior(y[[i]],x[[i]],alpha,beta,u[t,i],sigma2)
    }
    Q <- Q + w[t]*tmp.log
  }
  Q <- Q/sum(w)
  return(Q)
}

## partial derivative for log-likelihood
der.log.alpha <- function(y,x,alpha,beta,u){
  n <- length(y)
  tmp <- 0
  for ( i in 1:n ){
    pred_i <- alpha + beta*x[[i]]+u[i]
    tmp <- tmp + sum(y[[i]]-exp(pred_i)/(1+exp(pred_i)))
  }
  return(tmp)
}

der.log.beta <- function(y,x,alpha,beta,u){
  n <- length(y)
  tmp <- 0
  for ( i in 1:n ){
    pred_i <- alpha + beta*x[[i]]+u[i]
    tmp <- tmp + sum(y[[i]]*x[[i]]-x[[i]]*exp(pred_i)/(1+exp(pred_i)))
  }
  return(tmp)
}

der.log.sigma2 <- function(sigma2, u){
  n <- length(u)
  tmp <- crossprod(u)/(2*(sigma2^2))-n/(2*sigma2)
  return(as.numeric(tmp))
}

## start MCEM ##
## constant setting
nu <- 1
d <- 0.5
c <- 3
a <- 0.25
niter <- 100

# step 1. Initialize 
# number of random sample generated for u
mc.size <- 100
alpha.em <- alpha.glm
beta.em <- beta.glm
sigma2.em <- sigma2.glm
#alpha.em <- -0.5
#beta.em <- -0.5
#sigma2.em <- 0.5

# step 2. Generate u_1,...,u_m via Hasting-Metropolis Algorithm
u <- matrix(0, nrow = mc.size, ncol = n)
for ( i in 1:n ){
  u[,i] <- H_M(y[[i]],x[[i]],alpha.em[1],beta.em[1],u.glm[i],sigma2.em[1],mc.size)
}

# start EM iter step 3 to step 9
r <- 1
repeat{
  print(r)
  # step 3. compute importance weights
  w <- numeric()
  for ( t in 1:mc.size ){
    tmp <- 0
    for ( i in 1:n ){
      tmp <- tmp + (log.posterior(y[[i]],x[[i]],alpha.em[r],beta.em[r],u[t,i],sigma2.em[r])
                   - log.posterior(y[[i]],x[[i]],alpha.em[1],beta.em[1],u[t,i],sigma2.em[1]))
    }
    w[t] <- exp(tmp)
  }

  # step 4 & 5 E-M step
  par <- optim(c(alpha.glm,beta.glm,sigma2.glm), Qm, NULL, y = y, x = x, u = u, w = w, size=mc.size,
               lower = c(-6,-2,0.05), upper = c(-2,1,1), method = "L-BFGS-B",
               control = list(fnscale = -1 ), hessian = FALSE)$par
  alpha.em[r+1] <- par[1]
  beta.em[r+1] <- par[2]
  sigma2.em[r+1] <- par[3]
  print(par)

  # step 6 MC error estimation
  # (a)
  mu.alpha <- 0
  mu.beta <- 0
  mu.sigma2 <- 0
  for ( t in 1:mc.size ){
    mu.alpha <- mu.alpha + w[t]*der.log.alpha(y,x,alpha.em[r+1],beta.em[r+1],u[t,])
    mu.beta <- mu.beta + w[t]*der.log.beta(y,x,alpha.em[r+1],beta.em[r+1],u[t,])
    mu.sigma2 <- mu.sigma2 + w[t]*der.log.sigma2(sigma2.em[r+1],u[t,])
  }
  mu.alpha <- mu.alpha/sum(w)
  mu.beta <- mu.beta/sum(w)
  mu.sigma2 <- mu.sigma2/sum(w)
 
  # (b)
  var.alpha <- 0
  var.beta <- 0
  var.sigma2 <- 0
  for ( t in 1:mc.size ){
    var.alpha <- var.alpha + w[t]*(der.log.alpha(y,x,alpha.em[r+1],beta.em[r+1],u[t,])^2)
    var.beta <- var.beta + w[t]*(der.log.beta(y,x,alpha.em[r+1],beta.em[r+1],u[t,])^2)
    var.sigma2 <- var.sigma2 + w[t]*(der.log.sigma2(sigma2.em[r+1],u[t,])^2)
  }
  var.alpha <- var.alpha/sum(w) - mu.alpha^2
  var.beta <- var.beta/sum(w) - mu.beta^2
  var.simga2 <- var.sigma2/sum(w) - mu.sigma2^2

  # (c)
  alpha.CI <- c(mu.alpha - qnorm(1-a/2)*sqrt(var.alpha),mu.alpha + qnorm(1-a/2)*sqrt(var.alpha) )
  beta.CI <- c(mu.beta - qnorm(1-a/2)*sqrt(var.beta),mu.beta + qnorm(1-a/2)*sqrt(var.beta) )
  sigma2.CI <- c(mu.sigma2 - qnorm(1-a/2)*sqrt(var.sigma2),mu.sigma2 + qnorm(1-a/2)*sqrt(var.sigma2) )

  # step 7   
  tk <- numeric()
  k <- 1
  repeat{
    xk <- rpois(1, lambda=nu*(k^d)) + 1
    tk[k] <- xk + sum(tk)
    if ( tk[k] > mc.size ){
      break
    }
    k <- k+1
  }
  N <- k-1
  tk <- tk[1:N]

  # step 8
  if ( r > 1 ){
    if ( ((alpha.CI[1] <= Q.alpha) && (Q.alpha <= alpha.CI[2])) && 
        ((beta.CI[1] <= Q.beta) && (Q.beta <= beta.CI[2])) &&
        ((sigma2.CI[1] <= Q.sigma2) && (Q.sigma2 <= sigma2.CI[2])) ){
      m0.size <- mc.size
      mc.size <- m0.size + floor(m0.size/c)
      ## generate new sample via MCMC algorithm
      u.new <- matrix(0, nrow = (mc.size - m0.size), ncol = n)
      # run H-M Algorithm for each u_i
      for ( i in 1:n ){
        u.new[,i] <- H_M(y[[i]],x[[i]],alpha.em[1],beta.em[1],u.glm[i],sigma2.em[1],(mc.size - m0.size))
      }
      u <- rbind(u,u.new)
      print('MC sample size')
      print(mc.size)
    }
  }

  # step 9
  Q.alpha <- 0
  Q.beta <- 0
  Q.sigma2 <- 0

  for ( k in 1:N ){
    Q.alpha <- Q.alpha + w[tk[k]]*der.log.alpha(y,x,alpha.em[r+1],beta.em[r+1],u[tk[k],])
    Q.beta <- Q.beta + w[tk[k]]*der.log.beta(y,x,alpha.em[r+1],beta.em[r+1],u[tk[k],])
    Q.sigma2 <- Q.sigma2 + w[tk[k]]*der.log.sigma2(sigma2.em[r+1],u[tk[k],])
  }
  Q.alpha <- Q.alpha/sum(w[tk])
  Q.beta <- Q.beta/sum(w[tk])
  Q.sigma2 <- Q.sigma2/sum(w[tk])

  r <- r+1
  if ( r > niter ){
    break
  }

}

## fit model by R ##
library(lme4)
lmer.fit <- glmer(hiv~arm+(1|community),family=binomial(link = logit),data=mwanza)
summary(lmer.fit)
