# 1_bayes_donation_model_multiple_ordinal.R
# multiple variable model for ordinal version
# Sep 2021
library(dplyr)
library(R2WinBUGS)

# get the data
load('data/donate_analysis_ready.RData') # from 0_read_data_donation.R

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
bugs = file(model.file, 'w')
cat('model{
## Likelihood
  for (i in 1:N){
    donate[i, 1:ncat] ~ dmulti(pi[i, 1:ncat], 1)
    # Cumulative probability of > category k given cutpoint
    for (k in 1:(ncat-1)){
      logit(Q[i,k]) <- beta[countrynum[i],1]*X[i,1] +
  beta[countrynum[i],2]*X[i,2] +
  beta[countrynum[i],3]*X[i,3] +
  beta[countrynum[i],4]*X[i,4] +
  beta[countrynum[i],5]*X[i,5] +
  beta[countrynum[i],6]*X[i,6] +
  beta[countrynum[i],7]*X[i,7] +
  beta[countrynum[i],8]*X[i,8] +
  beta[countrynum[i],9]*X[i,9] +
  beta[countrynum[i],10]*X[i,10]+
  beta[countrynum[i],11]*X[i,11]+
  beta[countrynum[i],12]*X[i,12]+
  beta[countrynum[i],13]*X[i,13]+
  beta[countrynum[i],14]*X[i,14]+ 
  beta[countrynum[i],15]*X[i,15]+ 
  beta[countrynum[i],16]*X[i,16] - C[countrynum[i],k];
    }
    # Calculate probabilities
    pi[i,1] <- 1 - Q[i,1]; # Pr(cat=1) = 1 - Pr(cat>1);
    for (k in 2:(ncat-1)) {
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
C = matrix(rep(c(NA,0.1,0.2,0.3), n.country), byrow=TRUE, ncol=ncat-1) # needs ordered start
inits = list(beta = matrix(rep(0, P*n.country), ncol=P), 
             C = C,
             tau.beta = rep(1,P), 
             mu.beta = rep(0,P))
#inits = rep(list(inits), n_chains)

# run BUGS
parms = c('beta','C','mu.beta','tau.beta')
bugs.results =  bugs(data=bdata, inits=inits, parameters=parms, model.file=model.file, bugs.seed=89,
                     n.chains=n_chains, n.iter=n_samples*2*n_thin, n.thin=n_thin, DIC=FALSE, debug=TRUE,
                     bugs.directory="c:/Program Files/WinBUGS14")
bugs.results$summary

# save results
countries = unique(for.model$country)
save(bugs.results, X_vars, countries, file='results/bugs_model_multiple_ordinal.RData')

## export to jags - stuck with bugs
#save(inits, bdata, parms, model.file, file='//hpc-fs/barnetta/vaccine/jags_ready_ordinal.RData')
