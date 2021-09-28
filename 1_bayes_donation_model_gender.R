# 1_bayes_donation_model_gender.R
# model examining gender
# March 2021
# question is does Income, political leaning or education predict the binary outcome any willingness to donate vaccine versus none
library(dplyr)
library(R2WinBUGS)

### Part 0: set up the data

# get the data
load('data/donate_analysis_ready.RData') # from 0_read_data_donation.R

# manage data
data = mutate(data, 
              male = gender %in% c('Homme','Male','Maschio','Masculino'), 
              countrynum = as.numeric(as.factor(country)))

# labels used in plots/tables
countries = unique(data$country)
n.country = n_distinct(data$country) # number of countries
glabels = c('Non-Male','Male')

## set up the data as a matrix
small = select(data, donateany, country, countrynum) # other variables to not add to contrast matrix
formula = 'donateany ~ -1 + factor(male)'
for.model = data.frame(model.matrix(as.formula(formula), data=data))
names(for.model) = c('gender_non_male','gender_male')
for.model = bind_cols(for.model, small)

### Part 1: model with between country heterogeneity

# create external text file with bugs model (version with varying country effect)
model.file = 'bugs_donation_gender.txt'
bugs = file(model.file, 'w')
cat('model{
for (j in 1:N) {
  donateany[j] ~ dbern(mu[j])
  logit(mu[j]) <- alpha[countrynum[j],1]*X[j,1] +
                  alpha[countrynum[j],2]*X[j,2] 
}
for (i in 1:M){ # loop through countries
  alpha[i, 1:C] ~ dmnorm(mu.alpha[1:C], omega[1:C, 1:C])
}
for (i in 1:C){ # 
  mu.alpha[i] ~ dnorm(0, tau.alpha[i])
  tau.alpha[i] ~ dgamma(0.1, 0.1)
}
omega[1:C, 1:C] ~ dwish(R[1:C, 1:C], C)
}\n', file=bugs)
close(bugs)

# MCMC parameters
source('1_MCMC_parameters.R')

# prepare the random data
N = nrow(data) # 
X = as.matrix(select(for.model, starts_with('gender_'))) # extract variables
X[,1] = 1 # reference group is left
C = ncol(X)
R = toeplitz(c(1,0))
bdata = list(N = N, M = n.country, C = C, X = X, R = R, donateany = for.model$donateany, countrynum=for.model$countrynum)
inits = list(alpha=matrix(rep(0, C*n.country), ncol=C), tau.alpha=rep(1,C), mu.alpha=rep(0,C))
inits = rep(list(inits), n_chains)

# run BUGS
parms = c('alpha','mu.alpha','tau.alpha')
bugs.results =  bugs(data=bdata, inits=inits, parameters=parms, model.file=model.file, bugs.seed=89,
                     n.chains=n_chains, n.iter=n_samples*2*n_thin, n.thin=n_thin, DIC=dic, debug=debug,
                     bugs.directory="c:/Program Files/WinBUGS14")
bugs.results$summary

### part 2: fixed effects model

# create external text file with bugs model (version with non-varying country effect, but with random intercept per country)
model.file = 'bugs_donation_gender_fixed.txt'
bugs = file(model.file, 'w')
cat('model{
for (j in 1:N) {
  donateany[j] ~ dbern(mu[j])
  logit(mu[j]) <- alpha[countrynum[j]] +
                  beta*X[j,2] 
}
for (i in 1:M){ # loop through countries
  alpha[i] ~ dnorm(0, tau.country)
}
tau.country ~ dgamma(0.1, 0.1)
beta ~ dnorm(0, 0.00001)
}\n', file=bugs)
close(bugs)

# data and initial values
bdata = list(N = N, M = n.country, X = X, donateany = for.model$donateany, countrynum=for.model$countrynum)
inits = list(alpha=rep(0, n.country), tau.country=1, beta=0)
inits = rep(list(inits), n_chains)

# run BUGS
parms = c('alpha','beta','tau.country')
bugs.results.fixed =  bugs(data=bdata, inits=inits, parameters=parms, model.file=model.file, bugs.seed=89,
                     n.chains=n_chains, n.iter=n_samples*2*n_thin, n.thin=n_thin, DIC=dic, debug=debug,
                     bugs.directory="c:/Program Files/WinBUGS14")
bugs.results.fixed$summary

###
# save results
save(countries, n.country, glabels, bugs.results, bugs.results.fixed, file='results/bugs_model_gender.RData')

