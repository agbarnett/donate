# 1_bayes_donation_model_altruism_ordinal.R
# using tertiles within country for altruism
# version with ordinal outcome
# September 2021
library(dplyr)
library(R2WinBUGS)

# select response categories, five with missing data, four without
categories = 'four' 
# results are saved
outfile = 'results/bugs_ordinal_model_altruism.RData'
if(categories == 'four'){
  outfile = 'results/bugs_ordinal_model_altruism_four.RData'
}

### Part 0: set up the data

# get the data
load('data/donate_analysis_ready.RData') # from 0_read_data_donation.R
# remove do not know / prefer not to say
if(categories == 'four'){
  data = filter(data, donate_ordinal != 2) 
}
# country numbers
data = mutate(data, 
              altruism = relevel(altruism, ref='Meagre'), # use Meagre as reference
              countrynum = as.numeric(as.factor(country)))
countries = unique(data$country)
n.country = n_distinct(data$country) # number of countries

## set up the data as a matrix
small = select(data, donate_ordinal, country, countrynum) # other variables to not add to contrast matrix
# independent variable
formula = 'donate_ordinal ~ -1 + altruism'
for.model = data.frame(model.matrix(as.formula(formula), data=data))
names(for.model) = c('altruism_meagre','altruism_noanswer','altruism_generous')
for.model = bind_cols(for.model, small)
# labels used in plots/tables
altruism_labels = c('Meagre','Unanswered','Generous')

### Part 1: fixed effects model

# create external text file with bugs model
model.file = 'bugs_donation_altruism_ordinal_fixed.txt'
type = 'fixed'
n_indep = length(altruism_labels)
n_categories = ifelse(categories=='four', 4, 5)
source('99_make_bugs.R')

# MCMC parameters
source('1_MCMC_parameters.R')

# prepare the random data
N = nrow(for.model) # 
X = as.matrix(select(for.model, starts_with('altruism_'))) # extract variables
X[,1] = 1 # reference group is meagre
P = ncol(X)
# ordinal donation as a matrix
formula = 'countrynum ~ -1 + factor(donate_ordinal)'
donate = model.matrix(as.formula(formula), data=for.model) # one answer per row
ncat = ncol(donate)
bdata = list(N = N, 
             X = X, 
             P = P,
             donate = donate, 
             ncat = ncat)
c_start = NA 
for(k in 1:(n_categories-2)){
  c_start = c(c_start, k*0.1)
}
inits = list(beta=rep(0,P), C=c_start)
inits = rep(list(inits), n_chains)

# run BUGS
parms = c('beta','C')
bugs.results.fixed = bugs(data=bdata, inits=inits, parameters=parms, model.file=model.file, bugs.seed=189,
                     n.chains=n_chains, n.iter=n_samples*2*n_thin, n.thin=n_thin, DIC=dic, debug=debug,
                     bugs.directory="c:/Program Files/WinBUGS14")
bugs.results.fixed$summary

### part 2: model with between country heterogeneity

# create external text file with bugs model
model.file = 'bugs_donation_altruism_ordinal_variable.txt'
type = 'random'
source('99_make_bugs.R')

# prepare the additional data needed data for the random model
R = toeplitz(c(1,0,0))
P = ncol(R)
bdata = list(N = N, 
             M = n.country,
             countrynum = for.model$countrynum,
             X = X, 
             R = R,
             P = P,
             donate = donate, 
             ncat = ncat)
C = matrix(rep(c_start, n.country), byrow=TRUE, ncol=ncat-1) # needs ordered start
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
save(countries, n.country, altruism_labels, bugs.results, bugs.results.fixed, file=outfile)

