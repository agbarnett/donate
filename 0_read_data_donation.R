# 0_read_data_donation.R
# read the data for the charity donation analysis
# September 2021
library(dplyr)
library(janitor)
library(stringr)
library(purrr) # for map_chr
source('99_functions.R') # 

# list of countries for this analysis:
countries = c('Australia','Canada','France','Italy','Spain','UK','US')

# income data from larry
income = haven::read_dta('data/Income_from_Larry.dta') %>%
  select(id_country, country_string, age, hhequivinc) %>% # just keep household equivalent income
  rename('id' = 'id_country',
         'country' = 'country_string',
         'age_check' = 'age') # for checking merge
# make per country median for income
medians = group_by(income, country) %>%
  summarise(medinc = median(hhequivinc, na.rm=TRUE)) %>%
  ungroup()
income = left_join(income, medians, by='country') %>%
  mutate(incomecat = case_when(
    is.na(hhequivinc) == TRUE  ~ 3,
    hhequivinc < medinc ~ 1,
    hhequivinc >= medinc ~ 2
  ),
  incomecat = factor(incomecat, levels=1:3, labels=c('low','high','missing')))

## Survey data  
# from github https://github.com/CANDOUR-COVID/survey_data
raw = read.csv('data/CANDOUR.csv', stringsAsFactors = FALSE)
# used data management code in do file in stata folder
# make primary outcome `donateany`, any versus none donation; this groups together 'should not donate' with "Prefer not to say" and "Do not know" 
data = raw %>%
  clean_names() %>%
  filter(country %in% countries) %>% # just the selected countries
  mutate(
    # fix age
    age = ifelse(age > 1000, 2020 - age, age), # fixes for those who wrote a year
    #
    weights = ifelse(country %in% c('Canada','Spain'), 1, weights)) %>% # no weights needed
  rename('donate1'='geq_provision_1', # rename from Stata code
         'donate2'='geq_provision_2',
         'donate3'='geq_provision_3') %>%
  select(id, country, region_0, weights, age, gender, education_level, geq_donation, ideology, occupation, lottery_vignette, contains('income'), willing_risk, donation, donation_amount) %>%
  mutate(donate10 = as.numeric(geq_donation == 'Should donate 10%'),
         donateless = as.numeric(geq_donation == 'Should donate less than 10%'),
         donateany = as.numeric(geq_donation %in% c('Should donate more than 10%','Should donate 10%','Should donate less than 10%'))) %>%
  filter(age <= 105, # filter on age
         age >= 16 )

## add income data
data = left_join(data, income, by=c('country','id')) %>%
  select(-age_check) # not needed as merge is fine

## make categorical altruism variable
# first get median donation per country
donation_per_country = filter(data, !is.na(donation_amount)) %>% 
  group_by(country) %>%
  summarise(median  = median(donation_amount))
data = left_join(data, donation_per_country, by='country') %>%
  mutate(altruism = case_when(
    is.na(donation_amount) ~ 1, # missing, prefer not to say, and do not know
    !is.na(donation_amount) & donation_amount < median ~ 2, # less generous
    !is.na(donation_amount) & donation_amount >= median ~ 3 # more generous
  )) %>%
  mutate(altruism = factor(altruism, levels=1:3, labels=c('Not answered','Meagre','Generous')),
         altruism = relevel(altruism, ref='Meagre'))

# convert income by country
source('0_income_information_list.R') # contains all necessary income data in list
all_income = NULL
for (this_country in countries){
  # a) individual income
  with_income = median_income(
    indata = data,
    what_country= this_country, 
    what_income = 'income',
    lower_range_old = income_list[[this_country]]$lower_range_old,
    lower_range_new = income_list[[this_country]]$lower_range_new,
    upper_range_old = income_list[[this_country]]$upper_range_old,
    upper_range_new = income_list[[this_country]]$upper_range_new,
    extra = income_list[[this_country]]$extra,
    nil = income_list[[this_country]]$nil,
    nil_income = income_list[[this_country]]$nil_income,
    income_split = income_list[[this_country]]$income_split,
    shy = income_list[[this_country]]$shy,
    report = FALSE)
  all_income = bind_rows(all_income, with_income)

}
data = select(data, -income)#, -hh_income) # remove non-standardised variables from original data ...
data = left_join(data, all_income, by=c('country','id')) # ... add standardised variables (low/high/missing)

## make political ideology into relative tertiles per country
data = mutate(data, 
              countrynum = as.numeric(as.factor(country))) %>%
  group_by(country, countrynum) %>%
  mutate(lower = quantile(ideology, 0.33, na.rm=TRUE),
         upper = quantile(ideology, 0.67, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(ideologycat = 1 + as.numeric(ideology >= lower) + as.numeric(ideology > upper),
         ideologycat = ifelse(is.na(ideologycat)==TRUE, 4, ideologycat), # separate group for missing
         ideologycat = factor(ideologycat, levels=1:4, labels=c('Left','Central','Right','Missing'))) %>%
  select(-lower, -upper)
#with(data, table(group, country)) # check distribution

## adding variables
# fix-up gender, make risk and age groups, ordinal donation
data = mutate(data,
              donate_ordinal = case_when(
                geq_donation == 'Should not donate' ~ 1,
                geq_donation == 'Prefer not to say' ~ 2,
                geq_donation == 'Do not know' ~ 2,
                geq_donation == 'Should donate less than 10%' ~ 3,
                geq_donation == 'Should donate 10%' ~ 4,
                geq_donation == 'Should donate more than 10%' ~ 5
              ),
              risk_group = case_when(
                is.na(willing_risk) ~ 1,
                !is.na(willing_risk) & willing_risk <= 5 ~ 2,
                !is.na(willing_risk) & willing_risk > 5 ~ 3),
              risk_group = factor(risk_group, levels=1:3, labels=c('No answer','Unwilling','Willing')),
              age_group = 1 + as.numeric(age>=40) + as.numeric(age>=60),
              age_group = factor(age_group, levels=1:3, labels=c('Young','Mid','Old')),
              sex = case_when(
                gender == 'Male' ~ 1,
                gender == 'Maschio' ~ 1,
                gender == 'Masculino' ~ 1,
                gender == 'Homme' ~ 1,
                is.character(gender) ~ 0
              ),
              sex = factor(sex, levels=0:1, labels=c('Female','Male')))

## save the data ##
save(data, file='data/donate_analysis_ready.RData')

# weighted estimates, verifying I get the same as the Excel spreadsheet
summary = mutate(data,
                 wtdonate1 = donate10*weights,
                 wtdonate2 = donateless*weights) %>%
  group_by(country) %>%
  summarise(mean_10 = mean(wtdonate1, na.rm=TRUE),
            mean_less = mean(wtdonate2, na.rm=TRUE))
summary

# alternative weighted estimates 

library(survey)
s_design <- svydesign(weights=~weights, data=data, strata=~country, ids=~0)
svyciprop(formula = ~ donate10, s_design, method="logit")
svyciprop(formula = ~ donateless, s_design, method="logit")

confint(svymean(~ donate10+donateless, dclus1))
