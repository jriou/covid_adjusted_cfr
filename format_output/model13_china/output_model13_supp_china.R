# Setup ----
library(rstan)
library(tidyverse)
library(lubridate)
library(cowplot)
library(ggpubr)
theme_set(theme_bw())
library(xtable)

# Settings ----
source("data/china/data_management_china.R")
source("format_output/functions_model13.R")
figpath = "/home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/format_output/figures/"
col.cases1 = "darkcyan"
col.cases2 = "chartreuse3"
col.cases3 = "deepskyblue2"
col.deaths1 = "firebrick"
col.deaths2 =  "gold"
col.deaths3 = "darkorange"

# Load posterior samples (first extract chains and save them in run_model13_china.R) ----
l = load("posterior_samples/sensitivity_2020-03-18.Rdata")


# RHO 50 % -----------------------------------------------------------



## total symptomatic cases
print(S_model13B,"predicted_total_overall_symptomatic_cases")
print(S_model13B_rho50,"predicted_total_overall_symptomatic_cases")
summary(S_model13B_rho50,"predicted_total_overall_symptomatic_cases")[[1]] / 
  summary(S_model13B,"predicted_total_overall_symptomatic_cases")[[1]]
## reporting rate by age
print(S_model13B,"rho")
print(S_model13B_rho50,"rho")
## total cases
print(S_model13B,"predicted_total_overall_all_cases")
print(S_model13B_rho50,"predicted_total_overall_all_cases")
summary(S_model13B_rho50,"predicted_total_overall_all_cases")[[1]]/
  summary(S_model13B,"predicted_total_overall_all_cases")[[1]]
## CFR among symptomatics
print(S_model13B,"cfr_D_symptomatic",digits_summary = 4)
print(S_model13B_rho50,"cfr_D_symptomatic",digits_summary = 4)
summary(S_model13B_rho50,"cfr_D_symptomatic",digits_summary = 4)[[1]]/
  summary(S_model13B,"cfr_D_symptomatic",digits_summary = 4)[[1]]
## CFR among all
print(S_model13B,"cfr_D_all",digits_summary = 4)
print(S_model13B_rho50,"cfr_D_all",digits_summary = 4)


# Figure 2: model fit ----

i1 = plot_incidence_cases(S_model13B_rho50,D_S_model13B,start_date=as.Date("2019-12-31"),end_date=as.Date("2020-03-20"),
                          col1=col.cases1,col2=col.cases2,col3=col.cases3)
i2 = plot_total_cases(S_model13B_rho50,D_S_model13B,
                          col1=col.cases1,col2=col.cases2,col3=col.cases3)
i3 = plot_agedist_cases(S_model13B_rho50,D_S_model13B,
                      col1=col.cases1,col2=col.cases2,col3=col.cases3)

j1 = plot_incidence_deaths2(S_model13B_rho50,D_S_model13B,start_date=as.Date("2019-12-31"),end_date=as.Date("2020-03-20"),
                           col1=col.deaths1,col2=col.deaths2)
j2 = plot_total_deaths(S_model13B_rho50,D_S_model13B,start_date=as.Date("2019-12-31"),end_date=as.Date("2020-03-20"),
                           col1=col.deaths1,col2=col.deaths3)
j3 = plot_agedist_deaths(S_model13B_rho50,D_S_model13B,
                        col1=col.deaths1,col2=col.deaths3)

ggarrange(ggarrange(i1,i2,i3,j1,j2,j3,nrow=2,ncol=3,widths = c(3,1,2),labels=LETTERS),
          legend,nrow=2,heights=c(6,1))
ggsave(file=paste0(figpath,"modelfit_china_rho50.pdf"),width=10,height=6)

  

# Figure 3: age specific CFR

k1 = plot_agedist_cfr(S_model13B_rho50,D_S_model13B_rho50,col1="tomato",col2=col.deaths2,col3=col.cases3,
                      insert=c(.7,7,.1,.22))

k2 = plot_cfr(S_model13B_rho50,D_S_model13B_rho50,col1="tomato",col2=col.deaths2,col3=col.cases3)

ggarrange(k1,k2,widths = c(3.5,1),labels=LETTERS)

ggsave(file=paste0(figpath,"cfr_china_rho50.pdf"),width=8,height=4)

# supp table
rest = summary(S_model13B,pars=c("beta","psi","pi","eta","nu","xi"))[[1]][,c(6,4,8)] %>%
  as.data.frame() %>%
  tbl_df() %>%
  mutate(par=c("beta","psi","pi","eta","nu","xi"),
    prev=paste0(round(`50%`,3)," [",round(`2.5%`,3),"-",round(`97.5%`,3),"]"),
    age_group="General") %>%
  select(par,prev,age_group) %>%
  spread(par,prev)
rho = summary(S_model13B,pars="rho")[[1]][,c(6,4,8)] %>%
  as.data.frame() %>%
  tbl_df() %>%
  mutate(age_group=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),
         rho=paste0(round(`50%`,2)," [",round(`2.5%`,2),"-",round(`97.5%`,2),"]")) %>%
  select(age_group,rho)
epsilon = summary(S_model13B,pars="epsilon")[[1]][,c(6,4,8)] %>%
  as.data.frame() %>%
  tbl_df() %>%
  mutate(age_group=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),
         epsilon=paste0(round(`50%`,3)," [",round(`2.5%`,3),"-",round(`97.5%`,3),"]")) %>%
  select(epsilon)

bind_rows(data.frame(age_group="General",rho="-",epsilon="-"),
          bind_cols(rho,epsilon)) %>%
  left_join(rest) %>%
  xtable(.) %>%
  print(.,include.rownames=FALSE)






# _gamma15 -----------------------------------------------------------



## total symptomatic cases
print(S_model13B,"predicted_total_overall_symptomatic_cases")
print(S_model13B_gamma15,"predicted_total_overall_symptomatic_cases")
summary(S_model13B_gamma15,"predicted_total_overall_symptomatic_cases")[[1]] / 
  summary(S_model13B,"predicted_total_overall_symptomatic_cases")[[1]]
## reporting rate by age
print(S_model13B,"rho")
print(S_model13B_gamma15,"rho")
## total cases
print(S_model13B,"predicted_total_overall_all_cases")
print(S_model13B_gamma15,"predicted_total_overall_all_cases")
summary(S_model13B_gamma15,"predicted_total_overall_all_cases")[[1]]/
  summary(S_model13B,"predicted_total_overall_all_cases")[[1]]
## CFR among symptomatics
print(S_model13B,"cfr_D_symptomatic",digits_summary = 4)
print(S_model13B_gamma15,"cfr_D_symptomatic",digits_summary = 4)
summary(S_model13B_gamma15,"cfr_D_symptomatic",digits_summary = 4)[[1]]/
  summary(S_model13B,"cfr_D_symptomatic",digits_summary = 4)[[1]]
## CFR among all
print(S_model13B,"cfr_D_all",digits_summary = 4)
print(S_model13B_gamma15,"cfr_D_all",digits_summary = 4)


# Figure 2: model fit ----

i1 = plot_incidence_cases(S_model13B_gamma15,D_S_model13B,start_date=as.Date("2019-12-31"),end_date=as.Date("2020-03-20"),
                          col1=col.cases1,col2=col.cases2,col3=col.cases3)
i2 = plot_total_cases(S_model13B_gamma15,D_S_model13B,
                      col1=col.cases1,col2=col.cases2,col3=col.cases3)
i3 = plot_agedist_cases(S_model13B_gamma15,D_S_model13B,
                        col1=col.cases1,col2=col.cases2,col3=col.cases3)

j1 = plot_incidence_deaths2(S_model13B_gamma15,D_S_model13B,start_date=as.Date("2019-12-31"),end_date=as.Date("2020-03-20"),
                            col1=col.deaths1,col2=col.deaths2)
j2 = plot_total_deaths(S_model13B_gamma15,D_S_model13B,start_date=as.Date("2019-12-31"),end_date=as.Date("2020-03-20"),
                       col1=col.deaths1,col2=col.deaths3)
j3 = plot_agedist_deaths(S_model13B_gamma15,D_S_model13B,
                         col1=col.deaths1,col2=col.deaths3)

hh = filter(laterdeaths,date2>day_max) %>%
  mutate(deaths=ifelse(deaths==0,3,deaths))
j1bis = j1 + 
  geom_point(data=hh,aes(x=date2,y=deaths),shape=24,fill="grey")
ggarrange(ggarrange(i1,i2,i3,j1bis,j2,j3,nrow=2,ncol=3,widths = c(3,1,2),labels=LETTERS),
          legend,nrow=2,heights=c(6,1))
ggsave(file=paste0(figpath,"modelfit_china_gamma15.pdf"),width=10,height=6)




# _gamma15 -----------------------------------------------------------



## total symptomatic cases
print(S_model13B,"predicted_total_overall_symptomatic_cases")
print(S_model13_shenzhen,"predicted_total_overall_symptomatic_cases")
summary(S_model13_shenzhen,"predicted_total_overall_symptomatic_cases")[[1]] / 
  summary(S_model13B,"predicted_total_overall_symptomatic_cases")[[1]]
## reporting rate by age
print(S_model13B,"rho")
print(S_model13_shenzhen,"rho")
## total cases
print(S_model13B,"predicted_total_overall_all_cases")
print(S_model13_shenzhen,"predicted_total_overall_all_cases")
summary(S_model13_shenzhen,"predicted_total_overall_all_cases")[[1]]/
  summary(S_model13B,"predicted_total_overall_all_cases")[[1]]
## CFR among symptomatics
print(S_model13B,"cfr_D_symptomatic",digits_summary = 4)
print(S_model13_shenzhen,"cfr_D_symptomatic",digits_summary = 4)
summary(S_model13_shenzhen,"cfr_D_symptomatic",digits_summary = 4)[[1]]/
  summary(S_model13B,"cfr_D_symptomatic",digits_summary = 4)[[1]]
## CFR among all
print(S_model13B,"cfr_D_all",digits_summary = 4)
print(S_model13_shenzhen,"cfr_D_all",digits_summary = 4)


# Figure 2: model fit ----

i1 = plot_incidence_cases(S_model13_shenzhen,D_S_model13B,start_date=as.Date("2019-12-31"),end_date=as.Date("2020-03-20"),
                          col1=col.cases1,col2=col.cases2,col3=col.cases3)
i2 = plot_total_cases(S_model13_shenzhen,D_S_model13B,
                      col1=col.cases1,col2=col.cases2,col3=col.cases3)
i3 = plot_agedist_cases(S_model13_shenzhen,D_S_model13B,
                        col1=col.cases1,col2=col.cases2,col3=col.cases3)

j1 = plot_incidence_deaths2(S_model13_shenzhen,D_S_model13B,start_date=as.Date("2019-12-31"),end_date=as.Date("2020-03-20"),
                            col1=col.deaths1,col2=col.deaths2)
j2 = plot_total_deaths(S_model13_shenzhen,D_S_model13B,start_date=as.Date("2019-12-31"),end_date=as.Date("2020-03-20"),
                       col1=col.deaths1,col2=col.deaths3)
j3 = plot_agedist_deaths(S_model13_shenzhen,D_S_model13B,
                         col1=col.deaths1,col2=col.deaths3)

hh = filter(laterdeaths,date2>day_max) %>%
  mutate(deaths=ifelse(deaths==0,3,deaths))
j1bis = j1 + 
  geom_point(data=hh,aes(x=date2,y=deaths),shape=24,fill="grey")
ggarrange(ggarrange(i1,i2,i3,j1bis,j2,j3,nrow=2,ncol=3,widths = c(3,1,2),labels=LETTERS),
          legend,nrow=2,heights=c(6,1))
ggsave(file=paste0(figpath,"modelfit_china_gamma15.pdf"),width=10,height=6)



