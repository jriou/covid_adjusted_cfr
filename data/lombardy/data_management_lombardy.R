# Setup ----
source("setup.R")

# Controls ----

day_start = as.Date("2020-02-10")
day_data = as.Date("2020-02-11")
day_max = as.Date("2020-04-25")
day_quarantine = as.Date("2020-02-24")

## Age distribution in lombardy for 9 age classes ----

age_dist = read_excel("data/age_structure.xlsx") %>%
  filter(country=="Italy") %>%
  gather("age","n",3:23) %>%
  mutate(n=as.numeric(n),age2=c(0,0,10,10,20,20,30,30,40,40,50,50,60,60,70,70,80,80,80,80,80)) %>%
  group_by(age2) %>%
  summarise(n=sum(n)) %>%
  pull(n)
age_dist = age_dist/sum(age_dist)

## Population in Lombardia ----
# (Source Eurostat 2018)

pop_t = 10.04e6

#******************************************************************************************************************************
#Data 1: date of symptoms, up to 2020-03-08

## Case incidence by date of symptoms in lombardy  ----
# https://www.epicentro.iss.it/coronavirus/bollettino/Bolletino-sorveglianza-integrata-COVID-19_28-aprile-2020_appendix.pdf

lombardy_data_4may = read.csv("data/lombardy/lombardy_data_4may.csv") %>%
  tbl_df() %>%
  mutate(date=ymd(paste("2020",month,day,sep="-"))) %>% 
  filter(date>=day_data,date<=day_max)
incidence_cases = pull(lombardy_data_4may,cases)

ggplot(lombardy_data_4may) +
  geom_col(aes(x=date,y=cases))

## Deaths incidence in lombardy ----
# (https://github.com/pcm-dpc/COVID-19/blob/master/dati-regioni/dpc-covid19-ita-regioni.csv)

incidence_deaths = pull(lombardy_data_4may,deaths)

ggplot(lombardy_data_4may) +
  geom_col(aes(x=date,y=cases),fill="seagreen") +
  geom_col(aes(x=date,y=deaths),fill="firebrick") +
  geom_vline(xintercept = day_quarantine,linetype=2) +
  geom_vline(xintercept = day_max,linetype=2)


## Age distribution of cases ----
# (https://www.epicentro.iss.it/coronavirus/bollettino/Bolletino-sorveglianza-integrata-COVID-19_28-aprile-2020_appendix.pdf)

age_distributions_cases_deaths_4may = read.csv("data/lombardy/age_distributions_cases_deaths_4may.csv") %>%
  tbl_df()

cases_tmax = pull(age_distributions_cases_deaths_4may,cases)
prop_cases_tmax = cases_tmax / sum(cases_tmax)

mort_tmax = pull(age_distributions_cases_deaths_4may,deaths)
prop_mort_tmax = mort_tmax / sum(mort_tmax)

agedistr_cases = cases_tmax
agedistr_deaths = mort_tmax

plot(prop_mort_tmax,type="b",col="firebrick")
points(prop_cases_tmax,type="b",col="seagreen")

# Underreporting in all of Italy

p_underreport_cases = sum(incidence_cases)/74346
