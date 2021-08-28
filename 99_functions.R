# 99_functions.R
# functions for vaccine survey
# June 2021

# makes nice labels for table / plot
nice.label = function(intext){
  out = case_when(
    str_detect(intext, 'willing_risk') == TRUE ~ str_replace(intext, pattern='willing_risk', 'Willing to risk = '),
    str_detect(intext, '^ideology') == TRUE ~ str_replace(intext, pattern='^ideology', 'Ideology = '),
    str_detect(intext, '^income') == TRUE ~ str_replace(intext, pattern='income', 'Income = '),
    str_detect(intext, '^employment') == TRUE ~ str_remove(intext, pattern='employment'),
    str_detect(intext, '^country') == TRUE ~ str_remove(intext, pattern='country'),
    str_detect(intext, '^occupation') == TRUE ~ str_remove(intext, pattern='occupation'),
    str_detect(intext, '^education_level') == TRUE ~ str_remove(intext, pattern='education_level'),
    #str_detect(intext, 'occupation') == TRUE ~ str_remove(intext, pattern='occupation'),
    intext == 'education_levelUnknown/missing' ~ 'Unknown education',
    intext == 'employmentUnknown' ~ 'Unknown employment',
    intext == 'maleYes' ~ 'Male',
    intext == 'age' ~ 'Age +10 years',
    intext == 'int_lost_jobPrefer not to say' ~ 'Prefer not to say about lost income from COVID',
    intext == 'int_lost_jobYes' ~ 'Lost income from COVID',
    is.character(intext) ~ intext
  )
  return(out)
}

# function to process results from winbugs (version for multiple variable model with intercept)
process_results_multiple = function(in.results, max.category.num=3){

## construct probability of supporting vaccine donation from chains
names = names(in.results$sims.array[1,1,]) # names of all the stored results
all_results = NULL
for (countrynum in 1:n.country){
  index1 = grep(paste('^alpha\\[', countrynum, sep=''), names) # locations of the alpha's for this country
  indexi = grep(paste(',', 1 ,'\\]', sep=''), names) # intercept locations
  i = intersect(index1, indexi) # location of intercept for this country
  c1 = in.results$sims.array[,1,i] # chain 1  
  c2 = in.results$sims.array[,2,i] # chain 2
  frame = data.frame(countrynum = countrynum, num = 1, res=c(c1, c2), index=1:(length(c1)*2))
  all_results = bind_rows(all_results, frame)
  # other category levels
  for (num in 2:max.category.num){
    index2 = grep(paste(',', num ,'\\]', sep=''), names) # locations of these category levels
    u = intersect(index1, index2) # location of category (alpha + num)
    c1 = in.results$sims.array[,1,i] + in.results$sims.array[,1,u] # chain 1  
    c2 = in.results$sims.array[,2,i] + in.results$sims.array[,2,u] # chain 2
    frame = data.frame(countrynum = countrynum, num = num, res=c(c1, c2), index=1:(length(c1)*2))
    all_results = bind_rows(all_results, frame)
  }
}
# now add the overall means
indexi = grep('^mu.alpha\\[1\\]', names) # mean alpha for intercept
c1 = in.results$sims.array[,1,indexi] # chain 1  
c2 = in.results$sims.array[,2,indexi] # chain 2
frame = data.frame(countrynum = n.country+1, num = 1, res=c(c1, c2), index=1:(length(c1)*2)) # country number is + 1
all_results = bind_rows(all_results, frame)
for (num in 2:max.category.num){
  index2 = grep(paste('^mu.alpha\\[', num, '\\]', sep=''), names) # mean alpha for other questions
  c1 = in.results$sims.array[,1,indexi] + in.results$sims.array[,1,index2] # chain 1  
  c2 = in.results$sims.array[,2,indexi] + in.results$sims.array[,2,index2] # chain 2
  frame = data.frame(countrynum = n.country+1, num = num, res=c(c1, c2), index=1:(length(c1)*2))
  all_results = bind_rows(all_results, frame)
}
# 
all_results = mutate(all_results,
                     p = inv.logit(res), # inverse logit to get probability
                     overall = 1 + as.numeric(countrynum > n.country)) # binary for overall vs countries

return(all_results)

} # end of function


# function to process results from winbugs (version for random model)
process_results = function(in.results.random,
                           in.results.fixed,
                           max.category.num = 3){
  
  ### random estimates ###
  ## construct probability of supporting vaccine donation from chains
  names = names(in.results.random$sims.array[1,1,]) # names of all the stored results
  index = grep('^alpha\\[', names) # locations of the alpha's for countries
  c1 = in.results.random$sims.array[,1,index] # chain 1  
  c2 = in.results.random$sims.array[,2,index] # chain 2
  long1 = reshape2::melt(c1)
  long2 = reshape2::melt(c2)
  chains = bind_rows(long1, long2, .id='chain') %>%
    mutate(node = str_remove_all(Var2, "[^a-z]"),
           Var2 = str_remove_all(Var2, "[^0-9|,]")) %>%
    separate(Var2, c('countrynum','num'), sep=',') %>%
    mutate(index = as.numeric(Var1),
           countrynum = as.numeric(countrynum),
           num = as.numeric(num))
  # now add the overall means
  indexi = grep('^mu.alpha\\[', names) # mean alpha for intercept
  c1 = in.results.random$sims.array[,1,indexi] # chain 1  
  c2 = in.results.random$sims.array[,2,indexi] # chain 2
  long1 = reshape2::melt(c1)
  long2 = reshape2::melt(c2)
  o_chains = bind_rows(long1, long2, .id='chain') %>%
    mutate(node = str_remove_all(Var2, "[^a-z]"),
           index = 1:(length(c1)*2),
      countrynum = n.country+1,
      Var2 = as.numeric(str_remove_all(Var2, "[^0-9|,]"))) %>%
    rename('num' = 'Var2') %>%
    mutate(index = as.numeric(Var1),
           num = as.numeric(num))
  # 
  random_results = bind_rows(chains, o_chains) %>%
    rename('res' = 'value') %>% # to match other function
    mutate('type' = 'random',
           p = inv.logit(res), # inverse logit to get probability
           overall = 1 + as.numeric(countrynum > n.country)) # binary for overall vs countries
  
  ### fixed estimates ###
  ## construct probability of supporting vaccine donation from chains
  names = names(in.results.fixed$sims.array[1,1,]) # names of all the stored results
  index = grep('^alpha\\[', names) # locations of the alpha's for this country
  c1 = in.results.fixed$sims.array[,1,index] # chain 1  
  c2 = in.results.fixed$sims.array[,2,index] # chain 2
  long1 = reshape2::melt(c1)
  long2 = reshape2::melt(c2)
  chains = bind_rows(long1, long2, .id='chain') %>%
    mutate(num = 1, # intercept
           node = str_remove_all(Var2, "[^a-z]"),
           Var2 = str_remove_all(Var2, "[^0-9|,]"), 
           index = as.numeric(Var1),
           countrynum = as.numeric(Var2)) %>%
    select(-starts_with('Var'))
  # beta
  index = grep('^beta\\[', names) # locations of the beta's
  c1 = in.results.fixed$sims.array[,1,index] # chain 1  
  c2 = in.results.fixed$sims.array[,2,index] # chain 2
  long1 = reshape2::melt(c1)
  long2 = reshape2::melt(c2)
  chains_beta = bind_rows(long1, long2, .id='chain') %>%
    mutate(node = str_remove_all(Var2, "[^a-z]"),
           Var2 = str_remove_all(Var2, "[^0-9|,]"), 
           index = as.numeric(Var1),
           num = as.numeric(Var2)) %>%
    select(-starts_with('Var'))
  # add zero value for intercept
  zero = select(chains_beta, chain, index, node) %>%
    unique() %>%
    mutate(num = 0, value = 0)
  chains_beta = bind_rows(chains_beta, zero)
  # combine
  both = full_join(chains, chains_beta, by=c('chain','index')) %>% # repeats all values of beta across alpha
    mutate(num = num.x + num.y,
           value = value.x + value.y) %>%
    select(chain, index, num, countrynum, value)
  
  # now make the overall means
  o_chains = group_by(both, chain, index, num) %>%
    mutate(value = mean(value), # average across countries
           countrynum = n.country + 1) %>%
    ungroup()
  # 
  fixed_results = bind_rows(both, o_chains) %>%
    rename('res' = 'value') %>% # to match other function
    mutate('type' = 'fixed',
           p = inv.logit(res), # inverse logit to get probability
           overall = 1 + as.numeric(countrynum > n.country)) # binary for overall vs countries
  
  # combine fixed and random
  all_results = bind_rows(fixed_results, random_results) 
  
  return(all_results)
  
} # end of function

# lasso function to run multiple outcomes
lasso_this = function(indata, # data sets
                      expand = 0.1, # for creating space for labels on plot (just give right hand side); multiplier
                      family = 'binomial', # glm family
                      gtitle = '', # title for graph
                      outcome, # outcome variable 
                      o_levels, # outcome levels class as 'yes'
                      predictors # list of predictor variables
){
  # simple rename
  index = names(indata) == outcome
  names(indata)[index] = 'outcome'
  #
  indata = mutate(indata,
                  age = (age - 50)/10) %>% # scale age
    filter(!is.na(outcome)) # remove any missing the outcome
  if(family=='binomial'){indata = mutate(indata, outcome = as.numeric(outcome %in% o_levels)) } # binary outcome
  # make design matrix
  formula = paste('outcome ~ ', paste(predictors, collapse=' + '), sep='')
  design.matrix = model.matrix(as.formula(formula), data=indata)
  X = design.matrix[, -1] # no intercept
  y = indata$outcome
  # run lasso
  lasso = glmnet(x=X, y=y, family=family)
  # run cross-validation
  set.seed(23789)
  cv.lasso = cv.glmnet(x=X, y=y, family=family)
  # coefficients
  beta = coef(cv.lasso, s='lambda.1se')
  
  # refit best model using standard model to give CIs
  update_vars_to_include = row.names(beta)[which(beta!=0)]
  X.dash = X[, colnames(X) %in% update_vars_to_include]
  standard_model = glm(y~X.dash, family=family)
  # tidy estimates
  standard_model_ests= broom::tidy(standard_model, conf.int=TRUE) %>%
    mutate(term = str_remove(string=term, pattern='X.dash'))
  if(family == 'binomial'){ # extra step for binomial
    standard_model_ests = mutate(standard_model_ests,
           estimate = exp(estimate), # make odds ratios
           lower = exp(conf.low),
           upper = exp(conf.high))
    # stuff for plot
    yref = 1
    ylab = 'Odds ratio'
    xlabs = c('Lower odds','Higher odds')
  }
  if(family == 'gaussian'){
    standard_model_ests = rename(standard_model_ests,
                                 'lower' = 'conf.low',
                                 'upper' = 'conf.high')
    # stuff for plot
    yref = 0
    ylab = 'Estimated difference'
    xlabs = c('Disagree more','Agree more')
    
  }
  standard_model_ests = filter(standard_model_ests,
      !is.na(estimate), # small number missing
      term != '(Intercept)') %>% # do not need intercept
  rename('Variable' = 'term')
    
  ## refit best model with country-specific ideology
  # add country-by-ideology interaction to X.dash
  which_ideology = which(str_detect(colnames(X.dash), pattern='^ideology')) # which ideology variables were in the model
  which_country = which(str_detect(colnames(X.dash), pattern='^country')) # which countries variables were in the model
  X.dash.plus = X.dash
  for (i in which_ideology){
    for (j in which_country){
      if(str_detect(colnames(X.dash)[j], pattern='China') == FALSE){ # do not add China as ideology not collected
        X = as.matrix(X.dash[,i] * X.dash[,j]) # make interaction
        colnames(X) = paste(colnames(X.dash)[i], colnames(X.dash)[j], sep='')
        X.dash.plus = cbind(X.dash.plus, X)
      }
    }
  }
  standard_model_plus = glm(y~X.dash.plus, family=family)
  # extract means and confidence intervals for country effects
  ests = standard_model_plus$coefficients
  vcov = vcov(standard_model_plus)
  names = names(ests)
  names = str_remove_all(names, pattern='^X.dash.plus')
  which_ideology = names[str_detect(names, pattern='ideology')] # ideology names
  which_country = names[str_detect(names, pattern='country')] # country names
  countries = str_remove_all(names[str_detect(names,'^country')], '^country')
  countries = countries[countries !='China'] # ideology not measured in China
  ideologies = str_remove_all(names[str_detect(names,'^ideology')], '^ideology')
  ideologies = ideologies[!str_detect(ideologies, 'country')]
  interaction_estimates = NULL
  for (this_country in countries){
    for (this_ideology in ideologies){
      search1 = paste(c('^country', this_country, "$"), collapse='') 
      search2 = paste('ideology', this_ideology , '$', sep='') # needed because of 1 or 10
      search3 = paste(c('ideology', this_ideology, 'country', this_country), collapse='')
      to_detect = paste(c(search1, search2, search3), collapse='|')
      index = which(str_detect(names, to_detect))
      coef = sum(ests[index])
      var = sum(vcov[index,index])
      frame = data.frame(country =this_country, ideology=this_ideology, est = coef, var = var)
      interaction_estimates = bind_rows(interaction_estimates, frame)
    }
  }
  # make confidence intervals and prepare for plot
  interaction_estimates_plot = filter(interaction_estimates, ideology !='Missing') %>% # exclude missing
    mutate(ideology = as.numeric(ideology),
                                 z = qnorm(0.975),
                                 lower = est - (z*sqrt(var)),
                                 upper = est + (z*sqrt(var)),
                                 x = as.numeric(factor(country))) %>%
    arrange(x, ideology) %>%
    group_by(country) %>%
    mutate(jitter = (1:n() - 0.5 - (n()/2))/(2*n()), # jitter x-axis position to avoid overlap
           xj = x + jitter) %>%
    select(-z,  -jitter) %>%
    filter(var <= 100) %>% # exclude crazy estimates
    ungroup()
  if(family == 'binomial'){ # extra step for binomial
    interaction_estimates_plot = mutate(interaction_estimates_plot,
                                 est = exp(est), # make odds ratios
                                 lower = exp(lower),
                                 upper = exp(upper))
  }
  
  # colours:
  palette_int = brewer.pal(10, 'RdBu') # palette for ideology (red to blue)
  ideology_num = as.numeric(ideologies)
  ideology_num = ideology_num[order(ideology_num)]
  palette_int =palette_int[ideology_num]
  
  ## make plot for main estimates
  # prepare data
  for_plot = arrange(standard_model_ests, estimate) %>%
             mutate(colour = case_when(
              str_detect(Variable, pattern='country')==TRUE ~ 1,
              str_detect(Variable, pattern='education')==TRUE ~ 2,
              str_detect(Variable, pattern='occupation')==TRUE ~ 3,
              str_detect(Variable, pattern='employ')==TRUE ~ 4,
              str_detect(Variable, pattern='ideology')==TRUE ~ 5,
              str_detect(Variable, pattern='int_lost_job')==TRUE ~ 6,
              str_detect(Variable, pattern='willing_risk')==TRUE ~ 7,
              is.character(Variable)==TRUE ~ 8
            ),
            gtitle = gtitle,
            outcome = outcome,
            colour = factor(colour, levels=1:8),
            x=1:n()) %>%
    select(-std.error, -statistic, -p.value) # slim down
  N = nrow(for_plot)
  ## labels and colours
  # start with full colours and labels then slim down (gives consistent order and colour across plots)
  palette = brewer.pal(8, 'Dark2')
  legend_labels = c('Country','Education','Occupation','Employment','Ideology','Lost income from COVID','Take risks with your health','Other')
  to_use = as.numeric(unique(for_plot$colour)) # find just the colours that were used
  to_use = to_use[order(to_use)]
  palette = palette[to_use]
  legend_labels = legend_labels[to_use]
  # labels for each odds ratio
  est_labels = select(for_plot, x, estimate, lower, upper, colour, Variable) %>%
    mutate(label = nice.label(Variable)) # make nice labels
  # text for axis labels
  text1 = data.frame(x=0.7, y=yref, lower=yref, upper=yref, reference=1, label=xlabs[2])
  text2 = data.frame(x=0.7, y=yref, lower=yref, upper=yref, reference=1, label=xlabs[1])
  # plot
  plot = ggplot(data=for_plot, aes(x=x, y=estimate, ymin=lower, ymax=upper, col=colour))+
    scale_x_continuous(breaks=NULL, minor_breaks=1:N, labels=NULL, expand=expansion(mult = c(0.035, 0.015)))+ # narrow white space
    scale_y_continuous(expand = expansion(mult = c(0.005, expand)))+ # add space on right for labels
    scale_color_manual(NULL, values=palette, labels=legend_labels)+
    geom_hline(lty=2, yintercept=yref, col='gray33')+
    geom_point(size=2, shape=19)+
    geom_errorbar(width=0, size=1.02)+
    geom_text(data=text1, aes(x=x, y=y, label =label), adj=-0.1, vjust=1, col=grey(0.5))+
    geom_text(data=text2, aes(x=x, y=y, label =label), adj=1.1, vjust=1, col=grey(0.5))+
    geom_text(data=est_labels, aes(x=x, y=upper, label =label, col=colour), adj=0, show.legend = FALSE)+
    coord_flip()+
    ggtitle(gtitle)+
    theme_bw()+
    theme(
        legend.position = 'bottom',
        legend.justification = 'left',
        legend.box.spacing = unit(0, 'mm'), # reduce space between plot and legend
        legend.box.margin	= margin(t=0, r=0, b=0, l=0), # reduce space around legend
        legend.margin = margin(t=0, r=0, b=0, l=0, unit='mm'))+
    ylab(ylab)+
    xlab('')
  plot
  
  # just for binomial
  if(family=='binomial'){
    # rename for table
    standard_model_ests = rename(standard_model_ests, 'OR' = 'estimate')
  
    # truncate at upper limit and add arrows for upper limit
    upper_or_limit = 10
    interaction_estimates_plot = mutate(interaction_estimates_plot, 
                                        est = ifelse(est > upper_or_limit, upper_or_limit, est),
                                        upper = ifelse(upper > upper_or_limit, upper_or_limit, upper))
    add_arrows = filter(interaction_estimates_plot, upper == upper_or_limit)
    
  }
  
  ## make plot for interaction estimates
  plot_interaction = ggplot(data=interaction_estimates_plot, aes(x=xj, y=est, ymin=lower, ymax=upper, col=factor(ideology)))+
    scale_x_continuous(breaks=1:length(countries), labels=countries, expand=expansion(mult = c(0.035, 0.015)))+ # narrow white space
    scale_y_continuous(expand = expansion(mult = c(0.005, expand)))+ # add space on right for labels
    scale_color_manual('Ideology', values=palette_int, labels=ideology_num)+
    geom_hline(lty=2, yintercept=yref, col='gray33')+
    geom_point(size=2, shape=19)+
    geom_errorbar(width=0, size=1.02)+
    geom_text(data=text1, aes(x=x, y=y, label =label), adj=-0.1, vjust=1, col=grey(0.5))+
    geom_text(data=text2, aes(x=x, y=y, label =label), adj=1.1, vjust=1, col=grey(0.5))+
    coord_flip( )+ # truncate at OR of 10
    ggtitle(gtitle)+
    theme_bw()+
    theme(
      panel.grid.minor = element_blank(),
      legend.position = 'bottom',
      legend.justification = 'left',
      legend.box.spacing = unit(0, 'mm'), # reduce space between plot and legend
      legend.box.margin	= margin(t=0, r=0, b=0, l=0), # reduce space around legend
      legend.margin = margin(t=0, r=0, b=0, l=0, unit='mm'))+
    ylab(ylab)+
    xlab('')
  if(family == 'binomial'){ # only for binomial
    if(nrow(add_arrows)> 0){ # only if there are some over the limit
      plot_interaction = plot_interaction +
        geom_segment(data=add_arrows, arrow = arrow(length = unit(0.03, "npc")), aes(x=xj, xend=xj, y=upper_or_limit-1, yend=upper_or_limit, col=factor(ideology))) # arrow for upper limit
    }
  }
    
  ## return
  to.return = list()
  to.return$N = nrow(X)
  to.return$lasso = lasso
  to.return$cv.lasso = cv.lasso
  to.return$standard_model_ests = standard_model_ests
  to.return$for_plot = for_plot
  to.return$plot = plot
  to.return$plot_interaction = plot_interaction
  return(to.return)
  
}

# function to calculate median income; also used to add income to data
# also splits into high, low, missing
median_income = function(indata, 
                         what_country, 
                         lower_range_old,
                         lower_range_new,
                         upper_range_old,
                         upper_range_new,
                         what_income, 
                         extra = NULL, # extra stuff to clear at front
                         income_split = ' to ',
                         nil = NULL,
                         nil_income = NULL,
                         shy,
                         report = FALSE){
  # simple rename 
  index = names(indata) == what_income
  names(indata)[index] = 'Income'
  #
  all_respondents <- filter(indata, 
                        country == what_country) %>%
                 mutate(Income = ifelse(Income %in% shy, NA, Income), # turn shy into missing
   mid_range = str_squish(str_remove_all(string=Income, pattern='[^-|A-Z|a-z|0-9| ]'))) # remove anything bar letters, numbers and spaces
  # split into missing or not
  missing = filter(all_respondents, is.na(Income)) %>%
    mutate(middle = 3) # income is missing `3`
  one_country = filter(all_respondents, !is.na(Income))
  
  # extra fix for some countries
  if(what_country == 'China'){
    one_country <- mutate(one_country,
                          mid_range = str_replace_all(string=mid_range, pattern='4250084999', replacement = '42500 84999')) # one income that did not have space
  }
  # fix extra (if any)
  if(!is.null(extra)){
    one_country = mutate(one_country, mid_range = str_remove_all(string=mid_range, pattern=extra)) # 
  }

  # fix null (if any)
  if(!is.null(nil)){
    one_country = mutate(one_country, mid_range = ifelse(mid_range %in% nil, nil_income, mid_range)) # 
  }
  
  # fix lower and upper ranges (if any)
  if(!is.null(lower_range_old)){
    one_country = mutate(one_country, mid_range = ifelse(mid_range==lower_range_old, lower_range_new, mid_range)) # fix for range without 'to'
  }
  if(!is.null(upper_range_old)){
    one_country = mutate(one_country, mid_range = ifelse(mid_range==upper_range_old, upper_range_new, mid_range)) # fix for range without 'to'
  }
  print(table(one_country$mid_range)) # check
  one_country = mutate(one_country,  
                       tmp_chunks = str_split(mid_range, pattern=income_split)) %>%
    mutate(lower = map_chr(tmp_chunks, 1),
           upper = map_chr(tmp_chunks, 2),
           middle = (as.numeric(lower) + as.numeric(upper))/2) 
  # median income
  median = median(one_country$middle)
  if(report == TRUE){
    cat('25th income percentile = ', format(quantile(one_country$middle, 0.25), big.mark=','), '\n', sep='')
    cat('Median income = ', format(median(one_country$middle), big.mark=','), '\n', sep='')
    cat('75th income percentile = ', format(quantile(one_country$middle, 0.75), big.mark=','), '\n', sep='')
    # Gini
    gini = ineq(expanded, type="Gini")
    cat('Gini = ', format(gini, digits=2), '\n', sep='')
  }
  if(report==FALSE){
    # Transform to high low
    one_country = mutate(one_country, 
              middle = ifelse(middle < median, 1, 2)) 
    # Add back missing
    to_export = bind_rows(one_country, missing) %>%
      mutate(middle = factor(middle, levels=1:3, labels = c('Low','High','Missing'))) %>%
      select(id, country, middle) # only need merging variables as this is added to original data
    names(to_export)[names(to_export) == 'middle'] = what_income
    #
    return(to_export)
  }
}
