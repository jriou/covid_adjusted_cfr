# Setup ----
source("setup.R")

# Controls ----

# day_start is the day before first data point of interest (before exponential growth)
# day_data is the first day with data
# day_max is the last day with data
# day_quarantine is the day they introduced quarantine

day_start = as.Date("2020-02-17") #
day_data = as.Date("2020-02-18")
day_max = as.Date("2020-03-11") # FOR GITHUB DATA THIS IS 2020-03-11
day_quarantine = as.Date("2020-02-28") #alert level orange to red.


## Age distribution in South Korea for 9 age classes ----

age_dist = read_excel("data/age_structure.xlsx") %>%
  filter(country=="Korea, Republic of") %>% #name?
  gather("age","n",3:23) %>%
  mutate(n=as.numeric(n),age2=c(0,0,10,10,20,20,30,30,40,40,50,50,60,60,70,70,80,80,80,80,80)) %>%
  group_by(age2) %>%
  summarise(n=sum(n)) %>%
  pull(n)
age_dist = age_dist/sum(age_dist)

## Population in South Korea? ----
# (https://data.worldbank.org/country/korea-rep, 2018)
pop_t = 51.635e6

## Case incidence in Korea up to 2020-03-02 ----
#source:	https://doi.org/10.3346/jkms.2020.35.e112
#data:CUMULATIVE INCIDENCE	figure 2
#data:cases	  calculated CUM INC
#data:deaths	table3
#age:cases    figure 2C
#age:deaths	  figure 4

# korea_data_03march = read.csv("data/korea/korea_pub_cases.csv") %>%
#   tbl_df() %>%
#   mutate(date=ymd(paste("2020",month,day,sep="-"))) %>%
#   filter(date>=day_data)
# incidence_cases = pull(korea_data_03march,cases)

#incidence cases from github: https://raw.githubusercontent.com/jihoo-kim/Coronavirus-Dataset/master/patient.csv


korea_data_11march = read.csv("data/south-korea/korea_github_cases.csv") %>%
  tbl_df() %>%
  mutate(date=ymd(paste("2020",month,day,sep="-"))) %>%
  filter(date>=day_data) %>%
  mutate(deaths=ifelse(is.na(deaths),0,deaths))
incidence_cases = pull(korea_data_11march,cases)


## Deaths incidence in Northern Korea up to 2020-03-01 ----
# https://doi.org/10.3346/jkms.2020.35.e112
# incidence_deaths = pull(korea_data_03march,deaths)

# https://raw.githubusercontent.com/jihoo-kim/Coronavirus-Dataset/master/patient.csv
incidence_deaths = pull(korea_data_11march,deaths)



## Age distribution of cases in Korea up to 2020-03-02 ----
# https://doi.org/10.3346/jkms.2020.35.e112

age_distributions_cases_deaths_3march = read.csv("data/south-korea/korea_pub_ages.csv") %>%
  tbl_df()

cases_tmax = pull(age_distributions_cases_deaths_3march,cases)
prop_cases_tmax = cases_tmax / sum(cases_tmax)

mort_tmax = pull(age_distributions_cases_deaths_3march,deaths)
prop_mort_tmax = mort_tmax / sum(mort_tmax)



