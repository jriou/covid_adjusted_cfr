#Switzerland: lockdown: 17.03.2020
#source: OFSP

# Setup ----
source("setup.R")
library("lubridate")
library("dplyr")
age_class=function(x,min_age,age_range,max_age){
  age_lim=seq(min_age,max_age,age_range)
  return(sapply(as.list(x),function(x) sum(age_lim<=x)))
}

# Controls ----

day_start = as.Date("2020-03-01")
day_data = as.Date("2020-03-02")
day_max = as.Date("2020-04-23") #10 days before last date of data
day_quarantine = as.Date("2020-03-16")

## Population in Switzerland ----
# (Source OFS 2018) https://www.bfs.admin.ch/asset/fr/349-1800
pop_t = 8.545e6

#******************************************************************************************************************************
#age dist from: http://data.un.org/Data.aspx?d=POP&f=tableCode%3A22
excel_age_dist=read.csv("data/contact_matrix/age_distribution.csv",header=T)

age_dist = as.data.frame(excel_age_dist) %>%
  filter(Country.or.Area=="Switzerland") %>%
  filter(Year == 2018) %>%
  filter(Age %in% as.character(0:120)) %>% 
  mutate(age_group=age_class(as.numeric(levels(Age))[Age],10,10,80)) %>% 
  group_by(Sex,age_group) %>%
  summarise(n=sum(Value)) %>%
  pull(n)
age_dist = age_dist/sum(age_dist)

## Case incidence by date of symptoms and deaths in Switzerland, situation 24.04.2020
cases_deaths=read.csv("data/switzerland/agg_data.csv",header=T) %>%
  transmute(date=ymd(as.character(date)),new_cases=count_pos,new_cases_onset=onset_dt,new_deaths=death_dt) %>%
  filter(date>=ymd(day_data),date<=ymd(day_max))

incidence_cases = pull(cases_deaths,new_cases_onset)
incidence_deaths = pull(cases_deaths,new_deaths)
p_underreport_cases=sum(incidence_cases)/sum(cases_deaths$new_cases)
p_underreport_deaths=1

plot(cases_deaths$date,incidence_cases,type="l",col="firebrick")
plot(cases_deaths$date,incidence_deaths,type="l",col="seagreen")

#Age distribution of cases and deaths, situation 24.04.2020
agedistr_cases_deaths=read.csv("data/switzerland/age_cases_mort.csv",header=T) %>%
  select(age_class,cases,deaths)
agedistr_cases = pull(agedistr_cases_deaths,cases)
agedistr_deaths = pull(agedistr_cases_deaths,deaths)

plot(agedistr_cases/sum(agedistr_cases),type="b",col="firebrick",ylim=c(0,1))
points(agedistr_deaths/sum(agedistr_deaths),type="b",col="seagreen")