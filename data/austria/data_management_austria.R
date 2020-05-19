#Austria: lockdown: 16.03.2020, school closure: 14.03.2020, social distancing: 16.03.2020
#source: https://info.gesundheitsministerium.at/

# Setup ----
source("setup.R")
library("lubridate")
library("dplyr")
age_class=function(x,min_age,age_range,max_age){
  age_lim=seq(min_age,max_age,age_range)
  return(sapply(as.list(x),function(x) sum(age_lim<=x)))
}

# Controls ----

day_start = as.Date("2020-03-10")
day_data = as.Date("2020-03-11")
day_max = as.Date("2020-04-14") #10 days before last date of data
day_quarantine = as.Date("2020-03-16")

## Population in Austria ----
# (Source Eurostat 2019)
pop_t = 8.859e6

#******************************************************************************************************************************
#age dist from: http://data.un.org/Data.aspx?d=POP&f=tableCode%3A22
excel_age_dist=read.csv("data/contact_matrix/age_distribution.csv",header=T)

age_dist = as.data.frame(excel_age_dist) %>%
  filter(Country.or.Area=="Austria") %>%
  filter(Source.Year == 2018) %>%
  filter(Age %in% as.character(0:120)) %>% 
  mutate(age_group=age_class(as.numeric(levels(Age))[Age],5,10,75)) %>% 
  group_by(Sex,age_group) %>%
  summarise(n=sum(Value)) %>%
  pull(n)
age_dist = age_dist/sum(age_dist)

## Case incidence by date of symptoms and deaths in Austria, situation 25.04.2020
cases=read.table("data/austria/Epikurve.csv",sep=";",header=T) %>%
  transmute(date=dmy(time),new_cases=tag_Erkrankungen) %>%
filter(date>=ymd(day_data),date<=ymd(day_max))

deaths=read.table("data/austria/TodesfaelleTimeline.csv",sep=";",header=T) %>%
  transmute(date=dmy(time),new_deaths=c(NA,diff(Todesfalle))) %>%
filter(date>=ymd(day_data),date<=ymd(day_max))

incidence_cases = pull(cases,new_cases)
incidence_deaths = pull(deaths,new_deaths)
p_underreport_deaths=1
p_underreport_cases=1

plot(cases$date,incidence_cases,type="l",col="firebrick")
plot(deaths$date,incidence_deaths,type="l",col="seagreen")

#Age distribution of cases and deaths, situation 25.04.2020
agedistr_cases=read.table("data/austria/Altersverteilung.csv",sep=";",header=T) %>%
  transmute(age=age, cases=number, age2=c(1:9,9)) %>%
  group_by(age2) %>%
  summarise(n=sum(cases)) %>%
  pull(n)
  
agedistr_deaths=read.table("data/austria/AltersverteilungTodesfaelle.csv",sep=";",header=T) %>%
  transmute(age=age, new_deaths=number, age2=c(1:9,9)) %>%
  group_by(age2) %>%
  summarise(n=sum(new_deaths)) %>%
  pull(n)

plot(agedistr_cases/sum(agedistr_cases),type="b",col="firebrick",ylim=c(0,1))
points(agedistr_deaths/sum(agedistr_deaths),type="b",col="seagreen")