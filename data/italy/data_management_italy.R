# Setup ----
source("setup.R")

# Controls ----

day_start = as.Date("2020-02-07")
day_data = as.Date("2020-02-08")
day_max = as.Date("2020-03-03")
day_quarantine = as.Date("2020-03-08")

## Age distribution in Italy for 9 age classes ----

age_dist = read_excel("data/age_structure.xlsx") %>%
  filter(country=="Italy") %>%
  gather("age","n",3:23) %>%
  mutate(n=as.numeric(n),age2=c(0,0,10,10,20,20,30,30,40,40,50,50,60,60,70,70,80,80,80,80,80)) %>%
  group_by(age2) %>%
  summarise(n=sum(n)) %>%
  pull(n)
age_dist = age_dist/sum(age_dist)

## Population in Lombardia + Emilia Romagna + Veneto + Marche + Piemonte ----
# (Source Eurostat 2018)

pop_t = sum(c(10.04e6,4.453e6,4.905e6,1.522e6,4.356e6))

## Case incidence in Northern Italy up to 2020-03-12 ----
# (https://www.epicentro.iss.it/coronavirus/bollettino/Bollettino-sorveglianza-integrata-COVID-19_09-marzo-2020.pdf)
# (https://github.com/pcm-dpc/COVID-19/blob/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale.csv)

italy_data_15march = read.csv("data/italy/italy_data_15march.csv") %>%
  tbl_df() %>%
  mutate(date=ymd(paste("2020",month,day,sep="-"))) %>%
  filter(date>=day_data,date<=day_max)
incidence_cases = pull(italy_data_15march,cases)

ggplot(italy_data_15march) +
  geom_col(aes(x=date,y=cases))

## Deaths incidence in Northern Italy up to 2020-03-12 ----
# (https://github.com/pcm-dpc/COVID-19/blob/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale.csv)

incidence_deaths = pull(italy_data_15march,deaths)

ggplot(italy_data_15march) +
  geom_col(aes(x=date,y=deaths))


## Age distribution of cases in Northern Italy up to 2020-03-08 ----
# (https://www.epicentro.iss.it/coronavirus/bollettino/Bollettino-sorveglianza-integrata-COVID-19_09-marzo-2020.pdf)

age_distributions_cases_deaths_15march = read.csv("data/italy/age_distributions_cases_deaths_15march.csv") %>%
  tbl_df()

cases_tmax = pull(age_distributions_cases_deaths_15march,cases)
prop_cases_tmax = cases_tmax / sum(cases_tmax)

mort_tmax = pull(age_distributions_cases_deaths_15march,deaths)
prop_mort_tmax = mort_tmax / sum(mort_tmax)



