# 1_bayes_donation_model_age.R
# Feb 2021
# question is does Income, political leaning or education predict the binary outcome any willingness to donate vaccine versus none
library(ggplot2)
library(dplyr)
library(R2WinBUGS)

# get the data
load('data/donate_analysis_ready.RData') # from 0_read_data_donation.R

# make age categories
data = mutate(data, 
  countrynum = as.numeric(as.factor(country)),
#  age = ifelse(age< 18, NA, age), # missing, only 4 missing, so just leave as is
#  age = ifelse(age>105, NA, age), # missing
  age_group = 1 + as.numeric(age>=40) + as.numeric(age>=60))
n_missing = sum(is.na(data$age))

## set up the data as a matrix
small = select(data, donateany, country, countrynum) # other variables to not add to contrast matrix
formula = 'donateany ~ -1 + factor(age_group)'
for.model = data.frame(model.matrix(as.formula(formula), data=data))
names(for.model) = c('age_young','age_mid','age_old')
for.model = bind_cols(for.model, small)
n.country = n_distinct(data$country) # number of countries
# labels used in plots/tables
countries = unique(data$country)
age_labels = c('Up to 39','40 to 59','60+')

# create external text file with bugs model
model.file = 'bugs_donation_age.txt'
bugs = file(model.file, 'w')
cat('model{
for (j in 1:N) {
  donateany[j] ~ dbern(mu[j])
  logit(mu[j]) <- alpha[countrynum[j],1]*X[j,1] +
                  alpha[countrynum[j],2]*X[j,2] +
                  alpha[countrynum[j],3]*X[j,3] 
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
n_chains = 2
n_samples = 4000
n_thin = 3

# prepare the random data
N = nrow(for.model) # 
X = as.matrix(select(for.model, starts_with('age_'))) # extract age variables
X[,1] = 1 # reference group is young
C = ncol(X)
R = toeplitz(c(1,0,0))
bdata = list(N = N, M = n.country, C = C, X = X, R = R, donateany = for.model$donateany, countrynum=for.model$countrynum)
inits = list(alpha=matrix(rep(0, C*n.country), ncol=C), tau.alpha=rep(1,C), mu.alpha=rep(0,C))
inits = rep(list(inits), n_chains)

# run BUGS
parms = c('alpha','omega', 'mu.alpha')
bugs.results =  bugs(data=bdata, inits=inits, parameters=parms, model.file=model.file, bugs.seed=89,
                     n.chains=n_chains, n.iter=n_samples*2*n_thin, n.thin=n_thin, DIC=FALSE, debug=TRUE,
                     bugs.directory="c:/Program Files/WinBUGS14")
bugs.results$summary

# save results
save(countries, n.country, age_labels, bugs.results, file='results/bugs_model_age.RData')

