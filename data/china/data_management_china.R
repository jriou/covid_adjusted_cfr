# Setup ----
source("setup.R")

# Controls ----

day_start = as.Date("2019-12-31")
day_data = as.Date("2020-01-01")
day_max = as.Date("2020-02-11")
day_quarantine = as.Date("2020-01-20")


## Age distribution in China for 9 age classes ----

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

## Population in Hubei ----

# (National Bureau of Statistics of China) 
pop_t = 59020000

## Case incidence in Hubei up to 2020-02-11 ----

# (Chinese CDC Weekly, The epidemiological characteristics of an outbreak...)
confirmed_cases = read.csv("data/china/confirmed_cases.csv") %>%
  tbl_df() %>%
  mutate(date=ymd(paste(year,month,day,sep="-"))) %>%
  filter(date>=day_data)
incidence_cases = pull(confirmed_cases,confirmed_cases_hubei)


## Deaths incidence in Hubei up to 2020-02-11 (GuangchuangYu) -------

# Death incidence in China
# remotes::install_github("GuangchuangYu/nCov2019")
# library(nCov2019)
# x = get_nCov2019(lang='en')
# china_incidence = summary(x) %>%
#   tbl_df() %>%
#   mutate(date=as.Date(date,"%m.%d")) %>%
#   filter(date<=day_max,date>=day_data) %>%
#   mutate(incidence=confirm-lag(confirm,default=0))%>%
#   mutate(mort=dead-lag(dead,default=0))
# save(china_incidence,file="data/china/china_incidence.Rdata")
load("data/china/china_incidence.Rdata")
sum(china_incidence$incidence)

# ggplot(china_incidence,aes(date, incidence)) +
#   geom_col(fill='chartreuse4') + theme_minimal(base_size = 14) +
#   xlab(NULL) + ylab(NULL)
# ggplot(china_incidence,aes(date, mort)) +
#   geom_col(fill='firebrick') + theme_minimal(base_size = 14) +
#   xlab(NULL) + ylab(NULL)

death = china_incidence$mort
incidence_deaths_tot = c(rep(0,day_max-day_data+1-length(death)),death)

# Correct for confirmed cases in Hubei only (979 total)
incidence_deaths = round(979/sum(incidence_deaths_tot)*incidence_deaths_tot)



## Age distribution of cases in mainland China as of 2020-02-11 ----

# (Chinese CDC Weekly, The epidemiological characteristics of an outbreak...)
cases_tmax = c(416,549,3619,7600,8571,10008,8583,3918,1408)
prop_cases_tmax = cases_tmax / sum(cases_tmax)

## Age distribution of deaths in mainland China as of 2020-02-11 ----

# (Chinese CDC Weekly, The epidemiological characteristics of an outbreak...)
mort_tmax = c(0,1,7,18,38,130,309,312,208)
prop_mort_tmax = mort_tmax / sum(mort_tmax)


## Comorbidities among cases and deaths ----

comorbidities_total = c(44672-23690,1023-617)
hypertension_tmax = c(2683,161)
diabetes_tmax = c(1102,80)
cvd_tmax = c(873,92)
crd_tmax = c(511,32)
cancer_tmax = c(107,6)

# IHME 2017 http://ghdx.healthdata.org/gbd-results-tool
dist_diabetes_pop =
  read.csv("data/china/IHME_data_china_2017.csv",stringsAsFactors = F) %>% 
  tbl_df() %>%
  filter(measure_name=="Prevalence",metric_name %in% c("Percent","Number"),year==2017,cause_name=="Diabetes mellitus") %>%
  mutate(age_group=recode_factor(age_name,
                                 "1 to 4"="1","5 to 9"="1",
                                 "10 to 14"="10","15 to 19"="10",
                                 "20 to 24"="20","25 to 29"="20",
                                 "30 to 34"="30","35 to 39"="30",
                                 "40 to 44"="40","45 to 49"="40",
                                 "50 to 54"="50","55 to 59"="50",
                                 "60 to 64"="60","65 to 69"="60",
                                 "70 to 74"="70","75 to 79"="70",
                                 "80 plus"="80")) %>%
  select(sex_id,age_name,age_group,metric_name,val) %>%
  spread(metric_name,val) %>%
  mutate(sample_size=Number/Percent) %>%
  group_by(age_group) %>%
  summarise(Number=sum(Number),sample_size=sum(sample_size)) %>%
  mutate(p=Number/sample_size) %>%
  pull(p)

dist_crd_pop =
  read.csv("data/china/IHME_data_china_2017.csv",stringsAsFactors = F) %>% 
  tbl_df() %>% 
  filter(measure_name=="Prevalence",metric_name %in% c("Percent","Number"),year==2017,cause_name=="Chronic respiratory diseases") %>%
  mutate(age_group=recode_factor(age_name,
                                 "1 to 4"="1","5 to 9"="1",
                                 "10 to 14"="10","15 to 19"="10",
                                 "20 to 24"="20","25 to 29"="20",
                                 "30 to 34"="30","35 to 39"="30",
                                 "40 to 44"="40","45 to 49"="40",
                                 "50 to 54"="50","55 to 59"="50",
                                 "60 to 64"="60","65 to 69"="60",
                                 "70 to 74"="70","75 to 79"="70",
                                 "80 plus"="80")) %>%
  select(sex_id,age_name,age_group,metric_name,val) %>%
  spread(metric_name,val) %>%
  mutate(sample_size=Number/Percent) %>%
  group_by(age_group) %>%
  summarise(Number=sum(Number),sample_size=sum(sample_size)) %>%
  mutate(p=Number/sample_size) %>%
  pull(p)

dist_cvd_pop =
  read.csv("data/china/IHME_data_china_2017.csv",stringsAsFactors = F) %>% 
  tbl_df() %>% 
  filter(measure_name=="Prevalence",metric_name %in% c("Percent","Number"),year==2017,cause_name=="Cardiovascular diseases") %>%
  mutate(age_group=recode_factor(age_name,
                                 "1 to 4"="1","5 to 9"="1",
                                 "10 to 14"="10","15 to 19"="10",
                                 "20 to 24"="20","25 to 29"="20",
                                 "30 to 34"="30","35 to 39"="30",
                                 "40 to 44"="40","45 to 49"="40",
                                 "50 to 54"="50","55 to 59"="50",
                                 "60 to 64"="60","65 to 69"="60",
                                 "70 to 74"="70","75 to 79"="70",
                                 "80 plus"="80")) %>%
  select(sex_id,age_name,age_group,metric_name,val) %>%
  spread(metric_name,val) %>%
  mutate(sample_size=Number/Percent) %>%
  group_by(age_group) %>%
  summarise(Number=sum(Number),sample_size=sum(sample_size)) %>%
  mutate(p=Number/sample_size) %>%
  pull(p)

dist_hypertension_pop = c(.01,.01,0.087,0.123,0.199,0.338,0.480,.611,.611)

# Later deaths
laterdeaths = read.csv("data/china/time_series_19-covid-Deaths.csv") %>%
  tbl_df() %>%
  filter(Country.Region=="China",Province.State=="Hubei") %>%
  gather("date","mort",5:57) %>%
  mutate(date2=gsub("X","",date),
         date2=mdy(date2),
         deaths=mort-lag(mort,1,default = 0)) 
  
  
  
  