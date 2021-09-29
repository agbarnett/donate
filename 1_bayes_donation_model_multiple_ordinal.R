# 1_bayes_donation_model_multiple_ordinal.R
# multiple variable model for ordinal version
# Sep 2021
library(dplyr)
library(R2WinBUGS)

# select response categories, five with missing data, four without
categories = 'four' 
# results are saved
outfile = 'results/bugs_model_multiple_ordinal.RData'
if(categories == 'four'){
  outfile = 'results/bugs_model_multiple_ordinal_four.RData'
}

# get the data
load('data/donate_analysis_ready.RData') # from 0_read_data_donation.R
# remove do not know / prefer not to say
if(categories == 'four'){
  data = filter(data, donate_ordinal != 2) 
}
## manage data
data = mutate(data, 
              countrynum = as.numeric(as.factor(country))) %>%
  group_by(country, countrynum) %>%
  mutate(lower = quantile(ideology, 0.33, na.rm=TRUE), # make tertiles of ideology per country
         upper = quantile(ideology, 0.67, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(politics = 1 + as.numeric(ideology>=lower) + as.numeric(ideology>upper), 
         politics = ifelse(is.na(politics)==TRUE, 4, politics), # separate group for missing
         age_group = 1 + as.numeric(age>=40) + as.numeric(age>=60), #
         education_level = ifelse(education_level=='Less than primary completed', 'Primary completed', education_level), # combine lower two categories
         male = gender %in% c('Homme','Male','Maschio','Masculino'), 
         risk_group = case_when(
           is.na(willing_risk) ~ 2,
           !is.na(willing_risk) & willing_risk <= 5 ~ 1, # unwilling = reference category
           !is.na(willing_risk) & willing_risk > 5 ~ 3),
          risk_group = factor(risk_group, levels=1:3, labels=c('unwilling','missing','willing')),
         altruism = relevel(altruism, 'Meagre'), # use Meagre as reference
         countrynum = as.numeric(as.factor(country)))

## set up the data as a matrix, do not use missing as reference category
small = select(data, donate_ordinal, country, countrynum) # other variables to not add to contrast matrix
formula = 'donate_ordinal ~ 1 + factor(male) + factor(education_level) + factor(age_group) + factor(politics) + risk_group + altruism + incomecat' # altruism and income are already factors
for.model = data.frame(model.matrix(as.formula(formula), data=data))
X_vars = c('intercept','gender_male','e_secondary','e_uni','e_unknown','age_mid','age_old','poll_central','poll_right','poll_missing','risk_missing','risk_willing','altruism_noanswer','altruism_generous','income_high','income_missing') # without reference groups
names(for.model) = X_vars
for.model = bind_cols(for.model, small)
# cbind(names(for.model), X_vars) # check

# create external text file with bugs model
model.file = 'bugs_donation_multiple_ordinal.txt'
type = 'random'
n_indep = 16
n_categories = ifelse(categories=='four', 4, 5)
source('99_make_bugs.R')

# MCMC parameters
source('1_MCMC_parameters.R')

## prepare the random data
# ordinal donation as a matrix
formula = 'countrynum ~ -1 + factor(donate_ordinal)'
donate = model.matrix(as.formula(formula), data=for.model) # one answer per row
ncat = ncol(donate)
#
N = nrow(for.model) # 
X = as.matrix(select(for.model, 'intercept', starts_with('gender_'), starts_with('e_'), starts_with('age_'), starts_with('poll_'), starts_with('risk_'), starts_with('altruism_'), starts_with('income_'))) # extract variables
P = ncol(X)
R = toeplitz(c(1,rep(0, P-1)))
n.country = n_distinct(for.model$country) # number of countries
# in bugs format
bdata = list(N = N, 
             M = n.country, 
             P = P, 
             X = X, 
             R = R, 
             ncat = ncat,
             donate = donate, 
             countrynum = for.model$countrynum)
# initial values
c_start = NA 
for(k in 1:(n_categories-2)){
  c_start = c(c_start, k*0.1)
}
C = matrix(rep(c_start, n.country), byrow=TRUE, ncol=ncat-1) # needs ordered start
inits = list(beta = matrix(rep(0, P*n.country), ncol=P), 
             C = C,
             tau.beta = rep(1,P), 
             mu.beta = rep(0,P))
inits = rep(list(inits), n_chains)

# run BUGS
parms = c('beta','C','mu.beta','tau.beta')
bugs.results =  bugs(data=bdata, inits=inits, parameters=parms, model.file=model.file, bugs.seed=89,
                     n.chains=n_chains, n.iter=n_samples*2*n_thin, n.thin=n_thin, DIC=FALSE, debug=TRUE,
                     bugs.directory="c:/Program Files/WinBUGS14")
bugs.results$summary

# save results
countries = unique(for.model$country)
save(bugs.results, X_vars, countries, file=outfile)

## export to jags - stuck with bugs instead
#save(inits, bdata, parms, model.file, file='//hpc-fs/barnetta/vaccine/jags_ready_ordinal.RData')
