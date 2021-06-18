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
  c1 = bugs.results$sims.array[,1,i] # chain 1  
  c2 = bugs.results$sims.array[,2,i] # chain 2
  frame = data.frame(countrynum = countrynum, num = 1, res=c(c1, c2), index=1:(length(c1)*2))
  all_results = bind_rows(all_results, frame)
  # other category levels
  for (num in 2:max.category.num){
    index2 = grep(paste(',', num ,'\\]', sep=''), names) # locations of these category levels
    u = intersect(index1, index2) # location of category (alpha + num)
    c1 = bugs.results$sims.array[,1,i] + bugs.results$sims.array[,1,u] # chain 1  
    c2 = bugs.results$sims.array[,2,i] + bugs.results$sims.array[,2,u] # chain 2
    frame = data.frame(countrynum = countrynum, num = num, res=c(c1, c2), index=1:(length(c1)*2))
    all_results = bind_rows(all_results, frame)
  }
}
# now add the overall means
indexi = grep('^mu.alpha\\[1\\]', names) # mean alpha for intercept
c1 = bugs.results$sims.array[,1,indexi] # chain 1  
c2 = bugs.results$sims.array[,2,indexi] # chain 2
frame = data.frame(countrynum = n.country+1, num = 1, res=c(c1, c2), index=1:(length(c1)*2)) # country number is + 1
all_results = bind_rows(all_results, frame)
for (num in 2:max.category.num){
  index2 = grep(paste('^mu.alpha\\[', num, '\\]', sep=''), names) # mean alpha for other questions
  c1 = bugs.results$sims.array[,1,indexi] + bugs.results$sims.array[,1,index2] # chain 1  
  c2 = bugs.results$sims.array[,2,indexi] + bugs.results$sims.array[,2,index2] # chain 2
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
process_results_random = function(in.results, max.category.num=3){
  
  ## construct probability of supporting vaccine donation from chains
  names = names(in.results$sims.array[1,1,]) # names of all the stored results
  index = grep('^alpha\\[', names) # locations of the alpha's for this country
  c1 = bugs.results$sims.array[,1,index] # chain 1  
  c2 = bugs.results$sims.array[,2,index] # chain 2
  long1 = reshape2::melt(c1)
  long2 = reshape2::melt(c2)
  chains = bind_rows(long1, long2, .id='chain') %>%
    mutate(Var2 = str_remove_all(Var2, "[^0-9|,]")) %>%
    separate(Var2, c('countrynum','num'), sep=',') %>%
    mutate(index = as.numeric(Var1),
           countrynum = as.numeric(countrynum),
           num = as.numeric(num))
  # now add the overall means
  indexi = grep('^mu.alpha\\[', names) # mean alpha for intercept
  c1 = bugs.results$sims.array[,1,indexi] # chain 1  
  c2 = bugs.results$sims.array[,2,indexi] # chain 2
  long1 = reshape2::melt(c1)
  long2 = reshape2::melt(c2)
  o_chains = bind_rows(long1, long2, .id='chain') %>%
    mutate(
      index = 1:(length(c1)*2),
      countrynum = n.country+1,
      Var2 = as.numeric(str_remove_all(Var2, "[^0-9|,]"))) %>%
    rename('num' = 'Var2') %>%
    mutate(index = as.numeric(Var1),
           num = as.numeric(num))
  # 
  all_results = bind_rows(chains, o_chains) %>%
    rename('res' = 'value') %>% # to match other function
                mutate(p = inv.logit(res), # inverse logit to get probability
                       overall = 1 + as.numeric(countrynum > n.country)) # binary for overall vs countries
  
  return(all_results)
  
} # end of function

# lasso function to run lots of outcomes
lasso_this = function(indata, # data sets
                      expand = 0.1, # for creating space for labels on plot (just give right hand side); multiplier
                      family = 'binomial', # glm family
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
  
  # refit best model using standard model for CIs
  update_vars_to_include = row.names(beta)[which(beta!=0)]
  X.dash = X[, colnames(X) %in% update_vars_to_include]
  standard_model = glm(y~X.dash)
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
    
  ## make plot
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
            colour = factor(colour, levels=1:8),
            x=1:n())
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
    theme_bw()+
    theme(
        legend.position = 'top',
        legend.justification = 'left',
        legend.box.spacing = unit(0, 'mm'), # reduce space between plot and legend
        legend.box.margin	= margin(t=0, r=0, b=0, l=0), # reduce space around legend
        legend.margin = margin(t=0, r=0, b=0, l=0, unit='mm'))+
    ylab(ylab)+
    xlab('')
  plot
  
  # rename for table
  if(family=='binomial'){
    standard_model_ests = rename(standard_model_ests, 'OR' = 'estimate')
  }
  
  # return
  to.return = list()
  to.return$N = nrow(X)
  to.return$lasso = lasso
  to.return$cv.lasso = cv.lasso
  to.return$standard_model_ests = standard_model_ests
  to.return$plot = plot
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
