# 1_bayes_donation_model_multiple.R
# multiple variable model
# June 2021
# question is does Income, political leaning or education predict the binary outcome any willingness to donate vaccine versus none
library(ggplot2)
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
         age_group = 1 + as.numeric(age>=40) + as.numeric(age>=60),
         education_level = ifelse(education_level=='Less than primary completed', 'Primary completed', education_level), 
         male = gender %in% c('Homme','Male','Maschio','Masculino'), 
         risk_group = case_when(
           is.na(willing_risk) ~ 2,
           !is.na(willing_risk) & willing_risk <= 5 ~ 1, # meagre = reference category
           !is.na(willing_risk) & willing_risk > 5 ~ 3),
          risk_group = factor(risk_group, levels=1:3, labels=c('unwilling','missing','willing')),
         altruism = relevel(altruism, 'Meagre'),
         countrynum = as.numeric(as.factor(country)))

## set up the data as a matrix, do not use missing as reference category
small = select(data, donateany, country, countrynum) # other variables to not add to contrast matrix
formula = 'donateany ~ 1 + factor(male) + factor(education_level) + factor(age_group) + factor(politics) + risk_group + altruism + incomecat' # altruism and income are already factors
for.model = data.frame(model.matrix(as.formula(formula), data=data))
X_vars = c('intercept','gender_male','e_secondary','e_uni','e_unknown','age_mid','age_old','poll_central','poll_right','poll_missing','risk_missing','risk_willing','altruism_noanswer','altruism_generous','income_high','income_missing') # without reference groups
names(for.model) = X_vars
for.model = bind_cols(for.model, small)

# create external text file with bugs model
model.file = 'bugs_donation_multiple.txt'
bugs = file(model.file, 'w')
cat('model{
for (j in 1:N) {
  donateany[j] ~ dbern(mu[j])
  logit(mu[j]) <- alpha[countrynum[j],1]*X[j,1] +
                  alpha[countrynum[j],2]*X[j,2] +
                  alpha[countrynum[j],3]*X[j,3] +
                  alpha[countrynum[j],4]*X[j,4] +
                  alpha[countrynum[j],5]*X[j,5] +
                  alpha[countrynum[j],6]*X[j,6] +
                  alpha[countrynum[j],7]*X[j,7] +
                  alpha[countrynum[j],8]*X[j,8] +
                  alpha[countrynum[j],9]*X[j,9] +
                  alpha[countrynum[j],10]*X[j,10]+
                  alpha[countrynum[j],11]*X[j,11]+
                  alpha[countrynum[j],12]*X[j,12]+
                  alpha[countrynum[j],13]*X[j,13]+
                  alpha[countrynum[j],14]*X[j,14]+ 
                  alpha[countrynum[j],15]*X[j,15]+ 
                  alpha[countrynum[j],16]*X[j,16] 
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
N = nrow(for.model) # 
X = as.matrix(select(for.model, 'intercept', starts_with('gender_'), starts_with('e_'), starts_with('age_'), starts_with('poll_'), starts_with('risk_'), starts_with('altruism_'), starts_with('income_'))) # extract variables
C = ncol(X)
R = toeplitz(c(1,rep(0, C-1)))
n.country = n_distinct(for.model$country) # number of countries
# in bugs format
bdata = list(N = N, 
             M = n.country, 
             C = C, 
             X = X, 
             R = R, 
             donateany = for.model$donateany, 
             countrynum = for.model$countrynum)
inits = list(alpha=matrix(rep(0, C*n.country), ncol=C), tau.alpha=rep(1,C), mu.alpha=rep(0,C))
inits = rep(list(inits), n_chains)

# run BUGS
parms = c('alpha','mu.alpha','tau.alpha')
bugs.results =  bugs(data=bdata, inits=inits, parameters=parms, model.file=model.file, bugs.seed=89,
                     n.chains=n_chains, n.iter=n_samples*2*n_thin, n.thin=n_thin, DIC=FALSE, debug=TRUE,
                     bugs.directory="c:/Program Files/WinBUGS14")
bugs.results$summary

# save results
countries = unique(for.model$country)
save(bugs.results, X_vars, countries, file='results/bugs_model_multiple.RData')


## check for smaller set of variables using elastic net
library(glmnet)
Xdash = X[,-1] # remove intercept
y = for.model$donateany
model = glmnet(y=y, x=Xdash, family='binomial', alpha=0.95)
plot(model)
cvfit = cv.glmnet(x=X, y=y, family='binomial', alpha=0.95)
plot(cvfit)


## check VIF (data has to be data frame)
data = select(for.model, 'donateany', starts_with('gender_'), starts_with('e_'), starts_with('age_'), starts_with('poll_'), starts_with('risk_'), starts_with('altruism_'), starts_with('income_')) # extract variables
names = names(data)
formula = paste(names[1], '~', paste(names[-1], collapse='+'))
full_model = glm(formula, data=data, family='binomial')
car::vif(full_model)
