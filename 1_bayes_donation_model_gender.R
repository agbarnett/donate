# 1_bayes_donation_model_gender.R
# Feb 2021
# question is does Income, political leaning or education predict the binary outcome any willingness to donate vaccine versus none
library(ggplot2)
library(dplyr)
library(R2WinBUGS)

# get the data
load('data/donate_analysis_ready.RData') # from 0_read_data_donation.R

# manage data
data = mutate(data, 
              male = gender %in% c('Homme','Male','Maschio','Masculino'), 
              countrynum = as.numeric(as.factor(country)))

# labels used in plots/tables
n.country = n_distinct(data$country) # number of countries

countries = unique(data$country)
elabels = c('Less than primary','Primary','Secondary','University','Unknown')

# create external text file with bugs model
model.file = 'bugs_donation_gender.txt'
bugs = file(model.file, 'w')
cat('model{
for (j in 1:N) {
  donateany[j] ~ dbern(mu[j])
  logit(mu[j]) <- alpha[countrynum[j]] + beta[countrynum[j]]*male[j] # centre on zero
}
for (i in 1:M){ # loop through countries
  alpha[i] ~ dnorm(0, tau.country)
  beta[i] ~ dnorm(mu.beta, tau.beta)
}
  tau.beta ~ dgamma(0.1, 0.1)
  mu.beta ~ dnorm(0, 0.0001) # overall average
  tau.country ~ dgamma(0.1, 0.1)
}\n', file=bugs)
close(bugs)

# MCMC parameters
n_chains = 2
n_samples = 4000
n_thin = 3

# prepare the random data
N = nrow(data) # 
I = length(unique(data$ideology))
bdata = list(N = N, 
             M = n.country, 
             male = as.numeric(data$male), # plus one so that the lowest number is 1
             donateany = data$donateany, 
             countrynum = data$countrynum)
inits = list(alpha=rep(0,n.country), beta=rep(0,n.country), tau.beta=1, mu.beta=0)
inits = rep(list(inits), n_chains)

# run BUGS
parms = c('alpha','mu.beta','beta')
bugs.results =  bugs(data=bdata, inits=inits, parameters=parms, model.file=model.file, bugs.seed=89,
                     n.chains=n_chains, n.iter=n_samples*2*n_thin, n.thin=n_thin, DIC=FALSE, debug=TRUE,
                     bugs.directory="c:/Program Files/WinBUGS14")
bugs.results$summary

# save results
save(countries, n.country, elabels, bugs.results, file='results/bugs_model_gender.RData')

