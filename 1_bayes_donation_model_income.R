# 1_bayes_donation_model_income.R
# June 2021
# question is does Income, political leaning or education predict the binary outcome any willingness to donate vaccine versus none
library(ggplot2)
library(dplyr)
library(R2WinBUGS)

### Part 0: set up the data

# get the data
load('data/donate_analysis_ready.RData') # from 0_read_data_donation.R

# plot raw data as quick check
summary = group_by(data, incomecat) %>%
  summarise(mean = mean(donateany)) 
eplot = ggplot(data=summary, aes(x=incomecat, y=mean))+
  geom_point()+
  theme_bw()+
  coord_flip()
eplot

## set up the data as a matrix
small = select(data, donateany, country) # other variables to not add to contrast matrix
formula = 'donateany ~ -1 + incomecat'
for.model = data.frame(model.matrix(as.formula(formula), data=data))
names(for.model) = c('i_low','i_high','i_missing')
for.model = bind_cols(for.model, small) %>%
  mutate(countrynum = as.numeric(as.factor(country)))
n.country = n_distinct(data$country) # number of countries
# labels used in plots/tables
countries = unique(data$country)
ilabels = c('Low','High','Missing')

### Part 1: model with between country heterogeneity

# create external text file with bugs model
model.file = 'bugs_donation_income.txt'
bugs = file(model.file, 'w')
cat('model{
for (j in 1:N) {
  donateany[j] ~ dbern(mu[j])
  logit(mu[j]) <- alpha[countrynum[j],1]*X[j,1] +
                  alpha[countrynum[j],2]*X[j,2] +
                  alpha[countrynum[j],3]*X[j,3] 
}
for (i in 1:M){ # loop through countries
  alpha[i,1:C] ~ dmnorm(mu.alpha[1:C], omega[1:C,1:C])
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
N = nrow(for.model) # 
X = as.matrix(select(for.model, starts_with('i_'))) # extract income variables
X[,1] = 1 # reference group is less than primary
C = ncol(X)
R = toeplitz(c(1,0,0))
bdata = list(N = N, M = n.country, C = C, X = X, R = R, donateany = for.model$donateany, countrynum=for.model$countrynum)
inits = list(alpha=matrix(rep(0, C*n.country), ncol=C), tau.alpha=rep(1,C), mu.alpha=rep(0,C))
inits = rep(list(inits), n_chains)

# run BUGS
parms = c('alpha','omega', 'mu.alpha')
bugs.results =  bugs(data=bdata, inits=inits, parameters=parms, model.file=model.file, bugs.seed=89,
                     n.chains=n_chains, n.iter=n_samples*2*n_thin, n.thin=n_thin, DIC=dic, debug=debug,
                     bugs.directory="c:/Program Files/WinBUGS14")
bugs.results$summary

### Part 2: model with fixed effect 

# create external text file with bugs model
model.file = 'bugs_donation_income_fixed.txt'
bugs = file(model.file, 'w')
cat('model{
for (j in 1:N) {
  donateany[j] ~ dbern(mu[j])
  logit(mu[j]) <- alpha[countrynum[j]] +
                  beta[1]*X[j,2] +
                  beta[2]*X[j,3] 
}
for (i in 1:M){ # loop through countries
  alpha[i] ~ dnorm(0, tau.country)
}
tau.country ~ dgamma(0.1, 0.1)
for (k in 1:2){
   beta[k] ~ dnorm(0, 0.0001)
}
}\n', file=bugs)
close(bugs)

# data and initial values
bdata = list(N = N, M = n.country, X = X, donateany = for.model$donateany, countrynum=for.model$countrynum)
inits = list(alpha=rep(0, n.country), tau.country=1, beta=c(0,0))
inits = rep(list(inits), n_chains)

# run BUGS
parms = c('alpha','beta','tau.country')
bugs.results.fixed =  bugs(data=bdata, inits=inits, parameters=parms, model.file=model.file, bugs.seed=89,
                           n.chains=n_chains, n.iter=n_samples*2*n_thin, n.thin=n_thin, DIC=dic, debug=debug,
                           bugs.directory="c:/Program Files/WinBUGS14")
bugs.results.fixed$summary

###
# save results
save(countries, n.country, ilabels, bugs.results, bugs.results.fixed, file='results/bugs_model_income.RData')

