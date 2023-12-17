## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----cache=FALSE--------------------------------------------------------------
library(rstan)
library(bridgesampling)
library(lorad)

## ----cache=FALSE--------------------------------------------------------------
set.seed(12345)

mu <- 0
tau2 <- 0.5
sigma2 <- 1

n <- 20
theta <- rnorm(n, mu, sqrt(tau2))
y <- rnorm(n, theta, sqrt(sigma2))

## ----cache=FALSE--------------------------------------------------------------
mu0 <- 0
tau20 <- 1
alpha <- 1
beta <- 1

## ----cache=FALSE--------------------------------------------------------------
stancodeH0 <- 'data {
  int<lower=1> n; // number of observations
  vector[n] y; // observations
  real<lower=0> alpha;
  real<lower=0> beta;
  real<lower=0> sigma2;
}
parameters {
  real<lower=0> tau2; // group-level variance
  vector[n] theta; // participant effects
}
model {
  target += inv_gamma_lpdf(tau2 | alpha, beta);
  target += normal_lpdf(theta | 0, sqrt(tau2));
  target += normal_lpdf(y | theta, sqrt(sigma2));
}
'
stancodeH1 <- 'data {
  int<lower=1> n; // number of observations
  vector[n] y; // observations
  real mu0;
  real<lower=0> tau20;
  real<lower=0> alpha;
  real<lower=0> beta;
  real<lower=0> sigma2;
}
parameters {
  real mu;
  real<lower=0> tau2; // group-level variance
  vector[n] theta; // participant effects
}
model {
  target += normal_lpdf(mu | mu0, sqrt(tau20));
  target += inv_gamma_lpdf(tau2 | alpha, beta);
  target += normal_lpdf(theta | mu, sqrt(tau2));
  target += normal_lpdf(y | theta, sqrt(sigma2));
}
'
# compile models
stanmodelH0 <- stan_model(model_code = stancodeH0, model_name="stanmodel")
stanmodelH1 <- stan_model(model_code = stancodeH1, model_name="stanmodel")

## ----cache=FALSE--------------------------------------------------------------
# fit models
stanfitH0 <- sampling(stanmodelH0, data = list(y = y, n = n,
                                               alpha = alpha,
                                               beta = beta,
                                               sigma2 = sigma2),
                      iter = 50000, warmup = 1000, chains = 3, cores = 1, seed = 12345)
stanfitH1 <- sampling(stanmodelH1, data = list(y = y, n = n,
                                               mu0 = mu0,
                                               tau20 = tau20,
                                               alpha = alpha,
                                               beta = beta,
                                               sigma2 = sigma2),
                      iter = 50000, warmup = 1000, chains = 3, cores = 1, seed = 12345)

## ----cache=FALSE--------------------------------------------------------------
# compute log marginal likelihood via bridge sampling for H0
H0.bridge <- bridge_sampler(stanfitH0, silent = TRUE)

# compute log marginal likelihood via bridge sampling for H1
H1.bridge <- bridge_sampler(stanfitH1, silent = TRUE)

print(H0.bridge)
print(H1.bridge)

## ----cache=FALSE--------------------------------------------------------------
paramsH0 <- data.frame(extract(stanfitH0, permuted=TRUE))

paramsH1 <- data.frame(extract(stanfitH1, permuted=TRUE))

## ----cache=FALSE--------------------------------------------------------------
colspecH0 <- c("tau2"="unconstrained", "theta.1"="unconstrained", "theta.2"="unconstrained", "theta.3"="unconstrained", "theta.4"="unconstrained", "theta.5"="unconstrained", "theta.6"="unconstrained", "theta.7"="unconstrained", "theta.8"="unconstrained", "theta.9"="unconstrained", "theta.10"="unconstrained", "theta.11"="unconstrained", "theta.12"="unconstrained", "theta.13"="unconstrained", "theta.14"="unconstrained", "theta.15"="unconstrained", "theta.16"="unconstrained", "theta.17"="unconstrained", "theta.18"="unconstrained", "theta.19"="unconstrained", "theta.20"="unconstrained", "lp__"="posterior")
results0 <- lorad_estimate(paramsH0, colspecH0, 0.5, "random", 0.1)
lorad_summary(results0)

## ----cache=FALSE--------------------------------------------------------------
colspecH1 <- c("mu"="unconstrained", "tau2"="unconstrained", "theta.1"="unconstrained", "theta.2"="unconstrained", "theta.3"="unconstrained", "theta.4"="unconstrained", "theta.5"="unconstrained", "theta.6"="unconstrained", "theta.7"="unconstrained", "theta.8"="unconstrained", "theta.9"="unconstrained", "theta.10"="unconstrained", "theta.11"="unconstrained", "theta.12"="unconstrained", "theta.13"="unconstrained", "theta.14"="unconstrained", "theta.15"="unconstrained", "theta.16"="unconstrained", "theta.17"="unconstrained", "theta.18"="unconstrained", "theta.19"="unconstrained", "theta.20"="unconstrained", "lp__"="posterior")
results1 <- lorad_estimate(paramsH1, colspecH1, 0.5, "random", 0.1)
lorad_summary(results1)

