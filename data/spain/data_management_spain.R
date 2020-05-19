#Spain: lockdown: 14.03.2020, school closure: 13.03.2020, social distancing: 09.03.2020
#https://info.gesundheitsministerium.at/

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
day_max = as.Date("2020-04-16") #10 days before last day of data about number of new symptomatic cases
day_quarantine = as.Date("2020-03-09")

## Population in Spain ----
# (Source Eurostat 2019)

pop_t = 46.94e6

#******************************************************************************************************************************
#Data 1: date of symptoms
#age dist from: http://data.un.org/Data.aspx?d=POP&f=tableCode%3A22
excel_age_dist=read.csv("data/contact_matrix/age_distribution.csv",header=T)

age_dist = as.data.frame(excel_age_dist) %>%
  filter(Country.or.Area=="Spain") %>%
  filter(Year == 2018) %>%
  filter(Age %in% as.character(0:120)) %>% 
  mutate(age_group=age_class(as.numeric(levels(Age))[Age],10,10,80)) %>% 
  group_by(Sex,age_group) %>%
  summarise(n=sum(Value)) %>%
  pull(n)
age_dist = age_dist/sum(age_dist)

#Cases by day of symptom onsets, Spain, situation 28-04-2020
#https://www.mscbs.gob.es/profesionales/saludPublica/ccayes/alertasActual/nCov-China/documentos/Actualizacion_89_COVID-19.pdf
cases=read.csv("data/spain/spain_cases.csv",header = F) %>%
  transmute(date=seq(dmy("20-02-2020"),dmy("24-04-2020"),by="day"),new_cases=round(V2)) %>%
  filter(date>=ymd(day_data),date<=ymd(day_max))

## Death by day in Spain, dataset loaded from https://covid19.isciii.es/ (link at the bottom left), situation 28-04-2020
mort=read.csv("data/spain/spain_mortality.csv",header=T) %>%
  transmute(date=dmy(FECHA),deaths=round(Fallecidos)) %>%
  group_by(date) %>%
  summarise(deaths=sum(deaths,na.rm=T)) %>%
  mutate(new_deaths=c(NA,diff(deaths))) %>%
filter(date>=ymd(day_data),date<=ymd(day_max))

incidence_cases = cases$new_cases
incidence_deaths = mort$new_deaths
p_underreport_deaths=1
p_underreport_cases=142343/180689

plot(cases$date,incidence_cases,type="l",col="firebrick")
plot(mort$date,incidence_deaths,type="l",col="seagreen")

#number of cases reported on 16.04.2020 (reporte 17.04.2020): 188'068
#https://www.mscbs.gob.es/profesionales/saludPublica/ccayes/alertasActual/nCov-China/documentos/Actualizacion_78_COVID-19.pdf

#Age distribution of cases and deaths, situation 28-04-2020
#https://www.mscbs.gob.es/profesionales/saludPublica/ccayes/alertasActual/nCov-China/documentos/Actualizacion_89_COVID-19.pdf
agedistr=read.csv("data/spain/agedistr_spain.csv",header=T) %>%
  mutate(age2=c(1:9,9))  %>%
  group_by(age2) %>%
  summarise(cases=sum(cases,na.rm=T),deaths=sum(deaths,na.rm=T))
agedistr_cases=pull(agedistr,cases)
agedistr_deaths=pull(agedistr,deaths)

plot(agedistr_cases/sum(agedistr_cases),type="b",col="firebrick",ylim=c(0,1))
points(agedistr_deaths/sum(agedistr_deaths),type="b",col="seagreen")