library(nCov2019)
library(dplyr)
library(ggplot2)
library(lubridate)

source("data_new/china/data_management_china.R")
  
#China, Hubei
x <- load_nCov2019(lang='en')
d<-x[1]
d<-filter(d,province=="Hubei") %>%
  group_by(time) %>%
  summarise(deaths=sum(cum_dead),cases=sum(cum_confirm)) %>%
  mutate(new_deaths=c(0,diff(deaths)),new_cases=c(0,diff(cases)))
ggplot(data=d,aes(x=time,y=new_deaths))+geom_bar(fill="grey",stat="identity")
ggplot(data=d,aes(x=time,y=new_cases))+geom_bar(fill="grey",stat="identity")

incidence_cases_report<-filter(d,time<ymd("2020-02-12")) %>%
  pull(cases)

d_old<-select(confirmed_cases,time=date,new_cases=confirmed_cases_hubei)
d_old$new_deaths=incidence_deaths

ggplot()+geom_line(data=d,aes(x=time,y=new_cases))+
 geom_line(data=d_old,aes(x=time,y=new_cases),colour="blue")
sum_old<-filter(d_old,time<=ymd("2020-02-12")) %>%
  summarise(sum=sum(new_cases))
sum_update<-filter(d,time<=ymd("2020-03-10")) %>%
  summarise(sum=sum(new_cases))
sum_old/sum_update

ggplot()+geom_line(data=d,aes(x=time,y=new_deaths))+
  geom_line(data=d_old,aes(x=time,y=new_deaths),colour="blue")
sum_update<-summarise(d,sum=sum(new_deaths))
sum_old<-filter(d,time<ymd("2020-04-17")) %>%
  summarise(sum=sum(new_deaths))
sum_old/sum_update

#ECDC
# https://www.ecdc.europa.eu/en/publications-data/download-todays-data-geographic-distribution-covid-19-cases-worldwide

#China, Hubei
# https://github.com/GuangchuangYu/nCov2019
# http://www.bcloud.org/e/
list_china=list(K=9,
                day_max,
                day_start,
                day_data,
                day_quarantine,
                pop_t,
                age_dist,
                contact,
                incidence_cases,
                incidence_deaths,
                agedistr_cases,
                agedistr_deaths,
                p_underreport_deaths=1,
                p_underreport_cases=1)

#Italy:  lockdown: 11.03.2020, school closure: 05.03.2020, social distancing: 09.03.2020
# https://github.com/pcm-dpc
# https://www.epicentro.iss.it/coronavirus/sars-cov-2-sorveglianza-dati
# https://www.epicentro.iss.it/coronavirus/bollettino/Infografica_15marzo%20eng.pdf

#Germany: lockdown: 20.03.2020, school closure: 14.03.2020, social distancing: 12.03.2020
# https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Situationsberichte/2020-04-17-en.pdf?__blob=publicationFile
# https://experience.arcgis.com/experience/478220a4c454480e823b17327b2bf1d4

#Sweden: no lockdown, school closure: 18.03.2020, social distancing: 16.03.2020
# https://en.wikipedia.org/wiki/2020_coronavirus_pandemic_in_Sweden#cite_note-FHM_Aktuellt_läge-14
  
  
#Switerland:  lockdown: 17.03.2020, school closure: 17.03.2020, social distancing: 17.03.2020


d<-read.table(file = 'C:/Users/ahauser/Downloads/coronavirus_histrical_2020-03-14.tsv', sep = '\t', header = TRUE) 
d<-filter(d,province=="Hubei")%>%
  group_by(time) %>%
  summarise(death=sum(cum_dead))
plot(d$time[-1],diff(d$death))