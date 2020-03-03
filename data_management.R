# Set-up
library(tidyverse)
library(lubridate)
library(rstan)
library(cowplot)
theme_set(theme_bw())
source("setup.R")

### Controls ################################################

day_start = as.Date("2019-12-30")
day_data = as.Date("2019-12-08")
if(day_data<=day_start) day_data = day_start + 2
day_tmax = as.Date("2020-02-11")
day_tswitch = as.Date("2020-01-23")
time_lag = 100

t0 = 0
tmax = as.numeric(day_tmax-day_start)
S = tmax+time_lag
ts = 1:S
t_data = as.numeric(day_data-day_start)
D = tmax-t_data+1
tswitch = as.numeric(day_tswitch-day_start)

### Data ####################################################

# Age distribution in China for 9 age classes ---------
# (https://www.worldometers.info/demographics/china-demographics/) 

age_dist<- c(
  83932437 + 86735183,
  84262751 + 82341859,
  87158167+ 97989003,
  128738970 + 100091455, 
  96274146 + 119837617,
  123445382 + 98740491,
  77514139 + 74149766,
  44949689 + 26544616,
  16181417 + 7581777 + 2305024 + 475193 + 74692)
age_dist<-age_dist/sum(age_dist)

# Population in Hubei ---------
# (National Bureau of Statistics of China) 

pop_t = 59020000

# Case incidence in Hubei up to 2020-02-11 ---------
# (Chinese CDC Weekly, The epidemiological characteristics of an outbreak...)

confirmed_cases = read.csv("data/confirmed_cases.csv") %>%
  tbl_df() %>%
  mutate(date=ymd(paste(year,month,day,sep="-"))) %>%
  filter(date>=day_data)
incidence_cases = pull(confirmed_cases,confirmed_cases_hubei)

# Deaths incidence in Hubei up to 2020-02-11 (GuangchuangYu) -------

# Death incidence in China
# remotes::install_github("GuangchuangYu/nCov2019")
# library(nCov2019)
# x = get_nCov2019(lang='en')
# china_incidence = summary(x) %>%
#   tbl_df() %>%
#   mutate(date=as.Date(date,"%m.%d")) %>%
#   filter(date<=day_tmax,date>=day_data) %>%
#   mutate(incidence=confirm-lag(confirm,default=0))%>%
#   mutate(mort=dead-lag(dead,default=0)) 
# sum(china_incidence$incidence)
# 
# # ggplot(china_incidence,aes(date, mort)) +
# #   geom_col(fill='firebrick') + theme_minimal(base_size = 14) +
# #   xlab(NULL) + ylab(NULL)
# 
# death = china_incidence$mort
# incidence_deaths_tot = c(rep(0,tmax-t_data+1-length(death)),death)
# 
# # Correct for Hubei only
# incidence_deaths = round(979/sum(incidence_deaths_tot)*incidence_deaths_tot)

incidence_deaths = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 2, 3, 7, 7, 14, 13, 21, 23, 23, 33, 38, 40, 40, 50, 56, 58, 64,
                     64, 76, 78, 85, 95, 85)

# Age distribution of cases in mainland China as of 2020-02-11 --------- 
# (Chinese CDC Weekly, The epidemiological characteristics of an outbreak...)

cases_tmax = c(416,549,3619,7600,8571,10008,8583,3918,1408)
prop_cases_tmax = cases_tmax / sum(cases_tmax)

# Age distribution of deaths in mainland China as of 2020-02-11 --------- 
# (Chinese CDC Weekly, The epidemiological characteristics of an outbreak...)

mort_tmax = c(0,1,7,18,38,130,309,312,208)
prop_mort_tmax = mort_tmax / sum(mort_tmax)









