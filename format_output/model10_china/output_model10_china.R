# Setup
library(rstan)
library(tidyverse)
library(lubridate)
library(cowplot)
theme_set(theme_bw())

# Settings
col.cases1 = "darkcyan"
col.cases2 = "chartreuse3"
col.cases3 = "deepskyblue2"
col.deaths1 = "firebrick"
col.deaths2 =  "darkgoldenrod1"

# Load posterior samples
setwd("posterior_samples/")
SIM_model13 = read_stan_csv(dir(".",pattern = 'SIM_model13_2020-03-09-19-10-33_[[:digit:]]+.csv'))
D_SIM_model13 = read_rdump("data_SIM_model13.R")
S_model13 = read_stan_csv(dir(".",pattern = 'S_model13_2020-03-09-19-10-33_[[:digit:]]+.csv'))
D_S_model13 = read_rdump("data_S_model13.R")

# Format data and model outputs
data_incidence_cases = data.frame(time=D_SIM_model13$t_data:D_SIM_model13$S,incidence_data=D_SIM_model13$incidence_cases)
data_incidence_deaths = data.frame(time=D_SIM_model13$t_data:D_SIM_model13$S,incidence=D_SIM_model13$incidence_deaths)
pp1 = c("predicted_reported_incidence_symptomatic_cases",
       "predicted_overall_incidence_symptomatic_cases",
       "predicted_overall_incidence_all_cases")
SIM_pred_cases = summary(SIM_model13,pars=pp1)[[1]] %>%
  tbl_df() %>%
  mutate(time=rep(1:D_SIM_model13$S,3),
         date=time+as.Date("2019-12-30"),
         type=rep(pp1,each=D_SIM_model13$S)) %>%
  left_join(data_incidence_cases)
S_pred_cases = summary(S_model13,pars=pp1)[[1]] %>%
  tbl_df() %>%
  mutate(time=rep(1:D_S_model13$S,3),
         date=time+as.Date("2019-12-30"),
         type=rep(pp1,each=D_S_model13$S)) %>%
  left_join(data_incidence_cases)

pp2 = c("predicted_overall_incidence_deaths")
SIM_pred_deaths = summary(SIM_model13,pars=pp2)[[1]] %>%
  tbl_df() %>%
  mutate(time=1:D_SIM_model13$S,
         date=time+as.Date("2019-12-30"),
         sim="prior",
         type=rep(pp2,each=D_SIM_model13$S)) %>%
  left_join(data_incidence_deaths)
S_pred_deaths = summary(S_model13,pars=pp2)[[1]] %>%
  tbl_df() %>%
  mutate(time=1:D_S_model13$S,
         date=time+as.Date("2019-12-30"),
         sim="prior",
         type=rep(pp2,each=D_S_model13$S)) %>%
  left_join(data_incidence_deaths)

# Prior predictive checks
ggplot(SIM_pred_cases) +
  geom_ribbon(aes(x=date,ymin=`2.5%`,ymax=`97.5%`,group=type,fill=type),alpha=.5) +
  geom_line(aes(x=date,y=`50%`,group=type,colour=type))

# Posterior predictive checks
ggplot(S_pred_cases) +
  geom_ribbon(aes(x=date,ymin=`25%`,ymax=`75%`,group=type,fill=type),alpha=.5) +
  geom_line(aes(x=date,y=`50%`,group=type,colour=type))

