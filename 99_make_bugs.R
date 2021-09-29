## 99_make_bugs.R
# make winbugs file for ordinal model depending on number of independent variables and categories
# September 2021

## a) fixed
if(type == 'fixed'){
  bugs = file(model.file, 'w')
cat('model{
  ## Likelihood
  for (i in 1:N){
    donate[i, 1:ncat] ~ dmulti(pi[i, 1:ncat], 1)
    # Cumulative probability of > category k given cutpoint
    for (k in 1:(ncat-1)){
      logit(Q[i,k]) <- ', sep='', file=bugs)
for(k in 1:(n_indep-1)){
    cat('beta[',k,']*X[i,',k,'] + ', sep='', file=bugs)
}
cat(' beta[',n_indep,']*X[i,',n_indep,'] - C[k];', sep='', file=bugs)
    cat('}
    # Calculate probabilities
    pi[i,1] <- 1 - Q[i,1]; # Pr(cat=1) = 1 - Pr(cat>1);
    for (k in 2:ncat-1) {
      pi[i,k] <- Q[i,k-1] - Q[i,k]; # Pr(cat>k-1) - Pr(cat>k);
    }
    pi[i,ncat] <- Q[i,ncat-1]; # Pr(cat=k) = Pr(cat>k-1);
  }
  ## Priors
  for (k in 1:P){beta[k] ~ dnorm(0, 0.00001)}
  # ordered cut-offs
  C[1] <- 0; # for identifiability
        ', sep='', file=bugs)
    for(k in 1:(n_categories-3)){  
  cat('C[',k+1,'] ~ dunif(C[',k,'], C[',k+2,'])
      ', sep='', file=bugs)
    }
  cat('C[',n_categories-1,'] ~ dunif(C[',n_categories-2,'], 10)
      ', sep='', file=bugs)
cat('}\n', file=bugs)
close(bugs)
} # end of type if

## b) random
if(type == 'random'){
bugs = file(model.file, 'w')
cat('model{
  ## Likelihood
  for (i in 1:N){
    donate[i, 1:ncat] ~ dmulti(pi[i, 1:ncat], 1)
    # Cumulative probability of > category k given cutpoint
    for (k in 1:(ncat-1)){
          logit(Q[i,k]) <- ', sep='', file=bugs)
for(k in 1:(n_indep-1)){
  cat('beta[countrynum[i],',k,']*X[i,',k,'] + ', sep='', file=bugs)
}
cat(' beta[countrynum[i],',n_indep,']*X[i,',n_indep,'] - C[countrynum[i],k];', sep='', file=bugs)
cat('}
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
        ', sep='', file=bugs)
for(k in 1:(n_categories-3)){  
  cat('C[i,',k+1,'] ~ dunif(C[i,',k,'], C[i,',k+2,'])
      ', sep='', file=bugs)
}
cat('C[i,',n_categories-1,'] ~ dunif(C[i,',n_categories-2,'], 10)
      ', sep='', file=bugs)
cat('}
    }\n', file=bugs)
close(bugs)

}
