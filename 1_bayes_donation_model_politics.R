# 1_bayes_donation_model_politics.R
# using tertiles within country
# Feb 2021
library(ggplot2)
library(dplyr)
library(R2WinBUGS)

# get the data
load('data/donate_analysis_ready.RData') # from 0_read_data_donation.R

# plot raw data
summary = group_by(data, ideology) %>%
  summarise(mean = mean(donateany)) 
eplot = ggplot(data=summary, aes(x=ideology, y=mean))+
  geom_point()+
  theme_bw()+
  coord_flip()
eplot
# plot ideology distribution by country
bplot = ggplot(data=data, aes(y=ideology, x=country))+
  geom_boxplot()+
  theme_bw()
bplot

# make into tertiles per country
data = mutate(data, 
              countrynum = as.numeric(as.factor(country))) %>%
  group_by(country, countrynum) %>%
  mutate(lower = quantile(ideology, 0.33, na.rm=TRUE),
         upper = quantile(ideology, 0.67, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(group = 1 + as.numeric(ideology>=lower) + as.numeric(ideology>upper),
         group = ifelse(is.na(group)==TRUE, 4, group)) # separate group for missing
with(data, table(group, country)) # check distribution

# country numbers
n.country = n_distinct(data$country) # number of countries
countries = unique(data$country)

## set up the data as a matrix
small = select(data, donateany, country, countrynum) # other variables to not add to contrast matrix
formula = 'donateany ~ -1 + factor(group)'
for.model = data.frame(model.matrix(as.formula(formula), data=data))
names(for.model) = c('poll_left','poll_central','poll_right','poll_missing')
for.model = bind_cols(for.model, small)
# labels used in plots/tables
poll_labels = c('Left','Central','Right','Missing')

# create external text file with bugs model
model.file = 'bugs_donation_ideology.txt'
bugs = file(model.file, 'w')
cat('model{
for (j in 1:N) {
  donateany[j] ~ dbern(mu[j])
  logit(mu[j]) <- alpha[countrynum[j],1]*X[j,1] +
                  alpha[countrynum[j],2]*X[j,2] +
                  alpha[countrynum[j],3]*X[j,3] +
                  alpha[countrynum[j],3]*X[j,4]
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
N = nrow(data) # 
X = as.matrix(select(for.model, starts_with('poll_'))) # extract age variables
X[,1] = 1 # reference group is left
C = ncol(X)
R = toeplitz(c(1,0,0,0))
bdata = list(N = N, M = n.country, C = C, X = X, R = R, donateany = for.model$donateany, countrynum=for.model$countrynum)
inits = list(alpha=matrix(rep(0, C*n.country), ncol=C), tau.alpha=rep(1,C), mu.alpha=rep(0,C))
inits = rep(list(inits), n_chains)

# run BUGS
parms = c('alpha','mu.alpha')
bugs.results =  bugs(data=bdata, inits=inits, parameters=parms, model.file=model.file, bugs.seed=189,
                     n.chains=n_chains, n.iter=n_samples*2*n_thin, n.thin=n_thin, DIC=FALSE, debug=TRUE,
                     bugs.directory="c:/Program Files/WinBUGS14")
bugs.results$summary

# save results
save(countries, n.country, poll_labels, bugs.results, file='results/bugs_model_politics.RData')

