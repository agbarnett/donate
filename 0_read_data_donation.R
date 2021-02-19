# 0_read_data_donation.R
# read the data for the charity donation analysis
# Feb 2021
library(dplyr)
library(janitor)

# does not work
#library(foreign)
#data = read.dta('stata/Donation.dta')

# list of countries for this analysis:
countries = c('Australia','Canada','France','Italy','Spain','UK','US')
# from github https://github.com/CANDOUR-COVID/survey_data
raw = read.csv('data/CANDOUR.csv', stringsAsFactors = FALSE)
# used data management code in do file in stata folder
data = raw %>%
  clean_names() %>%
  filter(country %in% countries) %>%
  mutate(weights = ifelse(country %in% c('Canada','Spain'), 1, weights)) %>% # no weights needed
  rename('donate1'='geq_provision_1', # rename from Stata code
         'donate2'='geq_provision_2',
         'donate3'='geq_provision_3') %>%
  select(id, country, weights, age, gender, education_level, geq_donation, ideology, occupation) %>%
  mutate(donate10 = as.numeric(geq_donation == 'Should donate 10%'),
         donateless = as.numeric(geq_donation == 'Should donate less than 10%'),
         donateany = as.numeric(geq_donation %in% c('Should donate more than 10%','Should donate 10%','Should donate less than 10%'))) 
# primary outcome, any versus none donation; this groups together 'should not donate' with "Prefer not to say" and "Do not know" 
 
# weighted estimates, verifying I get the same as the Excel spreadsheet
summary = mutate(data,
                 wtdonate1 = donate10*weights,
                 wtdonate2 = donateless*weights) %>%
  group_by(country) %>%
  summarise(mean_10 = mean(wtdonate1, na.rm=TRUE),
            mean_less = mean(wtdonate2, na.rm=TRUE))
summary

# save the data
save(data, file='data/donate_analysis_ready.RData')
