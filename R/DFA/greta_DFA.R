# Greta Dynamic Factor Analysis:
#

#### Source:
# https://mdscheuerell.github.io/gretaDFA/



####  Packages  ####
library(gmRi)
library(targets)
library(tidyverse)
library(greta)
library(R6)
library(MASS)
library(viridis)
library(corrplot)

## R pkgs
installed.packages()[c("greta", "tensorflow"), "Version"]

## tensorflow
# tensorflow::install_tensorflow()
tensorflow::tf_version()


####  Process  ####



####  Functions  ####

"
We need to define a helper function to fill in the diagonal of the prior for 
the loadings matrix Z with the appropriate distribution (see below).
"

zd <- function(M, sigma) {
  zd <- zeros(M)
  for(i in 1:M) {
    zd[i] <- ld(i, M = M, sigma = sigma, dim = 1)
  }
  return(zd)
}

"
We use a unique distribution for the diagonal of Z, as defined by Leung & Drton (2016), 
that allows for parameter identification that is invariant to shifts in the rows of the data 
matrix y. At present, this distribution is not defined internally in greta, so we define it 
explicitly here via an R6 class.
"

distrib <- .internals$nodes$constructors$distrib
distribution_node <- .internals$nodes$node_classes$distribution_node
check_dims <- .internals$utils$checks$check_dims
as.greta_array <- .internals$greta_arrays$as.greta_array
is_scalar <- .internals$utils$misc$is_scalar
fl <- .internals$utils$misc$fl

# function to call 
ld <- function (i, M, sigma, dim = 1) {
  distrib("ld", i, M, sigma, dim)
}

# Leung & Drton distn for prior of diag(Z)
ld_distribution <- R6Class (
  "ld_distribution",
  inherit = distribution_node,
  public = list(
    
    initialize = function (i, M, sigma, dim) {
      
      # check if (i, M, dim) are in counting set
      # (i)
      if (length(i) > 1 ||
          i <= 0 ||
          !is.finite(i) ||
          i != floor(i)) {
        
        stop ("i must be a scalar positive integer, but was: ",
              capture.output(dput(i)),
              call. = FALSE)
        
      }
      
      # (M)
      if (length(M) > 1 ||
          M <= 0 ||
          !is.finite(M) ||
          M != floor(M)) {
        
        stop ("M must be a scalar positive integer, but was: ",
              capture.output(dput(M)),
              call. = FALSE)
        
      }
      
      # (dim)
      if (length(dim) > 1 ||
          dim <= 0 ||
          !is.finite(dim) ||
          dim != floor(dim)) {
        
        stop ("dim must be a scalar positive integer, but was: ",
              capture.output(dput(dim)),
              call. = FALSE)
        
      }
      
      # check if i > M
      if (i > M) {
        
        stop ("i can be no larger than M",
              call. = FALSE)
        
      }
      
      # check if sigma is scalar
      # (sigma)
      if (length(sigma) > 1) {
        
        stop ("sigma must be a scalar positive integer, but was: ",
              capture.output(dput(sigma)),
              call. = FALSE)
        
      }
      
      i <- as.greta_array(i)
      M <- as.greta_array(M)
      sigma <- as.greta_array(sigma)
      
      self$bounds <- c(0, Inf)
      super$initialize("ld", dim, truncation = c(0, Inf))
      self$add_parameter(i, "i")
      self$add_parameter(M, "M")
      self$add_parameter(sigma, "sigma")
      
    },
    
    tf_distrib = function (parameters, dag) {
      
      i <- parameters$i
      M <- parameters$M
      sigma <- parameters$sigma
      
      # log pdf(x | i, M, sigma)
      log_prob = function (x) {
        (M - i) * tf$log(x) - x ^ fl(2) / (fl(2) * sigma)
      }
      
      list(log_prob = log_prob, cdf = NULL, log_cdf = NULL)
      
    },
    
    # no CDF for discrete distributions
    tf_cdf_function = NULL,
    tf_log_cdf_function = NULL
    
  )
)


####  Simulate data
"Our general approach is to create a large number of time series, each of which to a 
greater or lesser degree share some temporal patterns with one another. 
We’ll use 30 time series that are 30 units long, each of which is a linear 
combination of 3 different latent trends.
"


NN <- 30
TT <- 30
MM <- 3


#### Latent factors
"In a DFA model, the rows in the matrix of latent factors x are generally 
assumed to be independent random walks, each of which is a cumulative sum 
of a sequence of independent process errors."



set.seed(123)
## MM x TT matrix of innovations
ww <- matrix(rnorm(MM*TT, 0, 1), MM, TT)
ww[,1] <- rnorm(MM, 0, sqrt(5))
## MM x TT matrix of scaled latent trends
xx <- t(scale(apply(ww,1,cumsum)))


matplot(t(xx), type="b", lty="solid", cex=0.7,
        xlab="Time", ylab=expression(italic(x)[italic(t)]),
        col=plasma(MM, end=0.8))

#### Loadings matrix
"The matrix Z maps the factors x onto the observations y. We draw each of the 
sub-diagonal elements from a Uniform(-1,1); the diagonal elements are drawn from a 
Uniform(0,1). We also sort the diagonal elements from largest to smallest.
"

ZZ <- matrix(runif(NN*MM, -1, 1), NN, MM)
diag(ZZ) <- rev(sort(abs(diag(ZZ))))
ZZ[upper.tri(ZZ)] <- 0
ZZ <- round(ZZ, 2)



####  Observed time series
"Now we can use the loadings and some observation errors to create the observed 
time series. Here we assume that the errors are IID, and their standard deviation 
is 0.2. We can ignore the additive effect of the offset vector a because the 
expectations of x and e are both zero."



## obs var
obs_var <- 0.2^2
## obs errors
ee <- t(mvrnorm(TT, matrix(0,NN,1), diag(obs_var,NN,NN)))
## NN x TT matrix of observed data
yy <- ZZ %*% xx + ee

matplot(t(yy), type="l", lty="solid",
        xlab="Time", ylab=expression(italic(y)[italic(i)]),
        col=plasma(NN, alpha=0.7, end=0.8))


"It's hard to tell from this plot how many of the `r NN` time series are correlated. 
Here is a plot of the correlation coefficients for all of the pairwise 
comparisons with them clustered by similarity. 
"

rho <- cor(t(yy))
par(mai=c(0.8,0.8,0.2,0.2))
corrplot(rho, 
         method="ellipse", 
         type="lower", 
         order = "hclust",
         tl.col = "black", 
         tl.srt = 0, 
         tl.cex = 0.6, 
         tl.offset = 0.7,
         cl.cex = 0.8, 
         cl.offset = 0.9, 
         cl.ratio = 0.1)




# Specify priors

Before we can estimate the model parameters we need to specify distributions for the priors and likelihoods.

## Loadings matrix

"To make the model identifiable when $M$ > 1, we need to impose several constraints 
on the model. First, the upper right triangle of $\mathbf{Z}$ must be set 
equal to zero (i.e., $z_{ij} = 0 ~ \forall ~ i < j$). For example, if $M$ = 3 then 

Second, when estimating a DFA model in a Bayesian framework, as with __greta__, 
we need to constrain the diagonal of $\mathbf{Z}$ to be positive to ensure 
convergence (i.e., $z_{ij} > 0 ~ \forall ~ i = j$). In particular, we implement 
the prior derived by Leung & Drton [-@leung2016], wherein each of the diagonal 
elements $z_{ii}$ has density proportional to
"
  

  
"Thus, we need to tell __greta__ which of the parameters in 
$\mathbf{Z}$ are to be estimated, and which should be fixed at 0. 
The order of operations is important here because __greta__ does 
not allow you to reassign elements in an array to fixed values after 
they have been declared as random. (Note that we follow the __greta__ 
convention of using the `=` assignment for stochastic nodes.)
"




## empty loadings matrix
ZZ_est <- zeros(NN,MM)
## define sigma
sigma_est = normal(0, 2, truncation = c(0, Inf))
## diagonal
idx_d <- row(ZZ_est) == col(ZZ_est)
ZZ_est_raw_d = zd(MM, sigma_est)
ZZ_est[idx_d] <- ZZ_est_raw_d
## sub-diagonal
idx_s <- lower.tri(ZZ_est)
ZZ_est_raw_s = normal(0, sigma_est, dim = sum(idx_s), truncation = c(-1,1))
ZZ_est[idx_s] <- ZZ_est_raw_s

