# 1_bayes_donation_model_gender_ordinal.R
# version with ordinal outcome
# September 2021
library(dplyr)
library(R2WinBUGS)

### Part 0: set up the data

# get the data
load('data/donate_analysis_ready.RData') # from 0_read_data_donation.R

# country numbers
countries = unique(data$country)
n.country = n_distinct(data$country) # number of countries

## set up the data as a matrix
small = select(data, donate_ordinal, country, countrynum) # other variables to not add to contrast matrix
# independent variable
formula = 'donate_ordinal ~ -1 + factor(gender)'
for.model = data.frame(model.matrix(as.formula(formula), data=data))
names(for.model) = c('gender_non_male','gender_male')
for.model = bind_cols(for.model, small)
# labels used in plots/tables
glabels = c('Non-Male','Male')


### Part 1: fixed effects model

# create external text file with bugs model
model.file = 'bugs_donation_gender_ordinal_fixed.txt'
bugs = file(model.file, 'w')
cat('model{
  ## Likelihood
  for (i in 1:N){
    donate[i, 1:ncat] ~ dmulti(pi[i, 1:ncat], 1)
    # Cumulative probability of > category k given cutpoint
    for (k in 1:(ncat-1)){
      logit(Q[i,k]) <- beta[1]*X[i,1] + 
      beta[2]*X[i,2] - C[k];
    }
    # Calculate probabilities
    pi[i,1] <- 1 - Q[i,1]; # Pr(cat=1) = 1 - Pr(cat>1);
    for (k in 2:ncat-1) {
      pi[i,k] <- Q[i,k-1] - Q[i,k]; # Pr(cat>k-1) - Pr(cat>k);
    }
    pi[i,ncat] <- Q[i,ncat-1]; # Pr(cat=k) = Pr(cat>k-1);
  }
  ## Priors
  for (k in 1:2){beta[k] ~ dnorm(0, 0.00001)}
  # ordered cut-offs
  C[1] <- 0; # for identifiability
  C[2] ~ dunif(C[1], C[3])
  C[3] ~ dunif(C[2], C[4])
  C[4] ~ dunif(C[3], 10)
}\n', file=bugs)
close(bugs)

# MCMC parameters
source('1_MCMC_parameters.R')

# prepare the random data
N = nrow(for.model) # 
X = as.matrix(select(for.model, starts_with('gender_'))) # extract variables
X[,1] = 1 # reference group is non-male
P = ncol(X)
# ordinal donation as a matrix
formula = 'countrynum ~ -1 + factor(donate_ordinal)'
donate = model.matrix(as.formula(formula), data=for.model) # one answer per row
ncat = ncol(donate)
bdata = list(N = N, 
             X = X, 
             donate = donate, 
             ncat = ncat)
inits = list(beta=rep(0,P), C=c(NA,0.1,0.2,0.3))
inits = rep(list(inits), n_chains)

# run BUGS
parms = c('beta','C')
bugs.results.fixed = bugs(data=bdata, inits=inits, parameters=parms, model.file=model.file, bugs.seed=189,
                     n.chains=n_chains, n.iter=n_samples*2*n_thin, n.thin=n_thin, DIC=dic, debug=debug,
                     bugs.directory="c:/Program Files/WinBUGS14")
bugs.results.fixed$summary

### part 2: model with between country heterogeneity

# create external text file with bugs model
model.file = 'bugs_donation_gender_ordinal_variable.txt'
bugs = file(model.file, 'w')
cat('model{
  ## Likelihood
  for (i in 1:N){
    donate[i, 1:ncat] ~ dmulti(pi[i, 1:ncat], 1)
    # Cumulative probability of > category k given cutpoint
    for (k in 1:(ncat-1)){
      logit(Q[i,k]) <- beta[countrynum[i],1]*X[i,1] + 
      beta[countrynum[i],2]*X[i,2] - C[countrynum[i],k];
    }
    # Calculate probabilities
    pi[i,1] <- 1 - Q[i,1]; # Pr(cat=1) = 1 - Pr(cat>1);
    for (k in 2:ncat-1) {
      pi[i,k] <- Q[i,k-1] - Q[i,k]; # Pr(cat>k-1) - Pr(cat>k);
    }
    pi[i,ncat] <- Q[i,ncat-1]; # Pr(cat=k) = Pr(cat>k-1);
  }
  ## Priors
  for (i in 1:M){ # loop through countries
    beta[i, 1:P] ~ dmnorm(mu.beta[1:P], omega[1:P, 1:P])
  }
  omega[1:P, 1:P] ~ dwish(R[1:P, 1:P], P)
  for (i in 1:P){ # 
    mu.beta[i] ~ dnorm(0, tau.beta[i])
    tau.beta[i] ~ dgamma(0.1, 0.1)
  }
  # ordered cut-offs
  for (i in 1:M){ # loop through countries
    C[i,1] <- 0; # for identifiability
    C[i,2] ~ dunif(C[i,1], C[i,3])
    C[i,3] ~ dunif(C[i,2], C[i,4])
    C[i,4] ~ dunif(C[i,3], 10)
  }
}\n', file=bugs)
close(bugs)

# prepare the additional data needed data for the random model
R = toeplitz(c(1,0))
P = ncol(R)
bdata = list(N = N, 
             M = n.country,
             countrynum = for.model$countrynum,
             X = X, 
             R = R,
             P = P,
             donate = donate, 
             ncat = ncat)
C = matrix(rep(c(NA,0.1,0.2,0.3), n.country), byrow=TRUE, ncol=ncat-1) # needs ordered start
inits = list(beta = matrix(rep(0, P*n.country), ncol=P), 
             C = C,
             tau.beta = rep(1,P), 
             mu.beta = rep(0,P))
inits = rep(list(inits), n_chains)

# run BUGS
parms = c('beta','C','mu.beta','tau.beta')
bugs.results =  bugs(data=bdata, inits=inits, parameters=parms, model.file=model.file, bugs.seed=89,
                           n.chains=n_chains, n.iter=n_samples*2*n_thin, n.thin=n_thin, DIC=dic, debug=debug,
                           bugs.directory="c:/Program Files/WinBUGS14")
bugs.results$summary

###
# save results
save(countries, n.country, glabels, bugs.results, bugs.results.fixed, file='results/bugs_ordinal_model_gender.RData')

