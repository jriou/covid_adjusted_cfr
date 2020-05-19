#Germany: lockdown: 20.03.2020, school closure: 14.03.2020, social distancing: 12.03.2020
# https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Situationsberichte/2020-04-26-en.pdf?__blob=publicationFile
# https://experience.arcgis.com/experience/478220a4c454480e823b17327b2bf1d4

# Setup ----
source("setup.R")
library("lubridate")
library("dplyr")
age_class=function(x,age_range,max_age){
  age_lim=seq(0,max_age,age_range)
  return(sapply(as.list(x),function(x) sum(age_lim<=x)))
}

# Controls ----

day_start = as.Date("2020-03-02")
day_data = as.Date("2020-03-03")
day_max = as.Date("2020-04-16")
day_quarantine = as.Date("2020-03-12")

## Population in badenw ----
# (Source Eurostat 2019)

pop_t = 11.07e6

#******************************************************************************************************************************
#Data 1: date of symptoms
#age dist from: http://data.un.org/Data.aspx?d=POP&f=tableCode%3A22
excel_age_dist=read.csv("data/germany/age_distribution_germany.csv",header=T)

age_dist100 = as.data.frame(excel_age_dist) %>%
  filter(Age %in% as.character(0:120)) %>% 
  mutate(age_group=age_class(as.numeric(levels(Age))[Age],10,100)) %>% 
  group_by(Sex,age_group) %>%
  summarise(n=sum(Value))
age_dist = as.data.frame(excel_age_dist) %>%
  filter(Age %in% as.character(0:120)) %>% 
  mutate(age_group=age_class(as.numeric(levels(Age))[Age],10,80)) %>% 
  group_by(age_group) %>%
  summarise(n=sum(Value)) %>%
  pull(n)
age_dist = age_dist/sum(age_dist)

## Case incidence by date of symptoms in Southern Germany ----
#https://experience.arcgis.com/experience/478220a4c454480e823b17327b2bf1d4
data_germany=read.csv("data/germany/cases_mort_germany.csv",header=T)


## Deaths incidence in Southern Germany,
#https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Situationsberichte/2020-04-26-en.pdf?__blob=publicationFile
#do it for each date, to have mortality over time in Bavaria, from 18-03-2020 (no deaths before) to 26-04-2020
data_badenw=transmute(data_germany,date=dmy(as.character(date)),new_cases=cases_badenw,new_deaths=c(NA,diff(deaths_badenw))) %>%
  filter(date>=ymd(day_data),date<=ymd(day_max))

incidence_cases=data_badenw$new_cases
incidence_deaths=data_badenw$new_deaths
p_underreport_cases=sum(incidence_cases)/31196
p_underreport_deaths=1

plot(data_badenw$date,incidence_cases,type="l",col="firebrick")
plot(data_badenw$date,incidence_deaths,type="l",col="seagreen")

## Age distribution of cases in Germany, situation 26 April 2020
#https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Situationsberichte/2020-04-26-en.pdf?__blob=publicationFile
age_cases = read.csv("data/germany/age_cases.csv",header = F) %>%
  tbl_df() %>%
  transmute(inc=V2) %>%
  mutate(age=rep(seq(0,100,10),2),age2=rep(c(seq(0,80,10),80,80),2),sex=rep(c("Female","Male"),each=11),pop=age_dist100$n,
         cases=inc*pop/100000) %>%
  group_by(age2) %>%
  summarise(cases=round(sum(cases)))
agedistr_cases = pull(age_cases,cases)

## Age distribution of cases in Germany, situation 26 April 2020
#https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Situationsberichte/2020-04-26-en.pdf?__blob=publicationFile
#age group 0-4, 5-14, 15-34, 35-59, 60-79, >=80
age_cases_march = read.csv("data/germany/age_cases17032020.csv",header = F) %>%
  tbl_df() %>%
  transmute(cases=V2) %>%
  mutate(age=rep(c(0,5,15,35,60,80),2),sex=rep(c("Male","Female"),each=6)) %>%
  group_by(age) %>%
  summarise(cases=round(sum(cases))) %>%
  pull(cases)

#Age distribution of deaths in Germany, situation 26 April 2020
#https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Situationsberichte/2020-04-26-en.pdf?__blob=publicationFile
age_mort = read.csv("data/germany/age_mortality.csv",header = F) %>%
  tbl_df() %>%
  transmute(mort=V2) %>%
  mutate(age=rep(seq(30,100,10),2),age2=rep(c(seq(30,80,10),80,80),2),sex=rep(c("Female","Male"),each=8)) %>%
  group_by(age2) %>%
  summarise(mort=round(sum(mort))) %>%
  pull(mort)
agedistr_deaths = c(0,0,0,age_mort)

plot(agedistr_cases/sum(agedistr_cases),type="b",col="firebrick",ylim=c(0,1))
points(agedistr_deaths/sum(agedistr_deaths),type="b",col="seagreen")