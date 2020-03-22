# Setup ----
library(rstan)
library(tidyverse)
library(lubridate)
library(cowplot)
library(ggpubr)
theme_set(theme_bw())

# Settings ----
source("format_output/functions_model13.R")
figpath = "/home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/format_output/figures/"
col.cases1 = "darkcyan"
col.cases2 = "chartreuse3"
col.cases3 = "deepskyblue2"
col.deaths1 = "firebrick"
col.deaths2 =  "darkgoldenrod1"

# Load posterior samples (first extract chains and save them in run_model13_china.R) ----
l = load("posterior_samples/model13_2020-03-13.Rdata")
l = load("posterior_samples/model13IT_2020-03-16.Rdata")

# Figure 1: data ----
# China
source("data/china/data_management_china.R")
confirmed_cases$deaths= D_S_model13B$incidence_deaths
h1 = ggplot(confirmed_cases) +
  geom_col(aes(x=date,y=confirmed_cases_hubei),fill=col.cases1,width=1,col="black") +
  geom_col(aes(x=date,y=deaths),fill=col.deaths1,width=1,col="black") +
  geom_vline(xintercept=day_quarantine+.5,linetype=2) +
  geom_vline(xintercept=day_max+.5,linetype=2) +
  scale_y_continuous(expand=expand_scale(mult=c(0,.05)),
                     labels=function(x) paste0(x/1000,"K")) +
  scale_x_date(limits=c(as.Date("2019-12-06"),day_max+1)) +
  labs(x="Date",y="N") +
  annotate("segment",x=as.Date("2019-12-15"),xend=as.Date("2019-12-08"),y=900,yend=5,arrow=arrow(length=unit(.2,"cm"))) +
  annotate("label",x=as.Date("2019-12-15"),y=900,label="First case") +
  annotate("segment",x=as.Date("2020-01-01"),xend=as.Date("2020-01-01"),y=900,yend=5,arrow=arrow(length=unit(.2,"cm"))) +
  annotate("label",x=as.Date("2020-01-01"),y=900,label="Market closure") +
  annotate("segment",x=as.Date("2020-01-03"),xend=as.Date("2020-01-20"),y=2150,yend=2150,arrow=arrow(length=unit(.2,"cm"))) +
  annotate("label",x=as.Date("2020-01-03"),y=2150,label="Control measures\nimplemented") +
  annotate("segment",x=as.Date("2020-01-03"),xend=as.Date("2020-02-11"),y=2750,yend=2750,arrow=arrow(length=unit(.2,"cm"))) +
  annotate("label",x=as.Date("2020-01-03"),y=2750,label="End of data collection") 
h1

dd2 = data.frame(age_class=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),
                 pop = age_dist,
                 cases_tmax=cases_tmax,
                 mort_tmax=mort_tmax,
                 prop_cases_tmax=prop_cases_tmax) %>%
  mutate(prop_mort_tmax=prop_mort_tmax,
         rawcfr=mort_tmax/cases_tmax,
         rawcfr2=rawcfr*prop_cases_tmax)
llab = c(pop="Chinese population",
         prop_cases_tmax="Reported cases",
         prop_mort_tmax="Reported deaths")
g1 = select(dd2,age_class,pop,prop_cases_tmax,prop_mort_tmax) %>%
  gather("var","value",2:4) %>%
  ggplot() +
  geom_col(aes(x=age_class,y=value,fill=var),colour="black",width=1) +
  facet_wrap(~var,labeller = labeller(var=llab)) +
  geom_hline(yintercept=0) +
  scale_fill_manual(values=c("seagreen",col.cases1,col.deaths1),guide=FALSE) +
  coord_flip() +
  scale_y_continuous(labels=scales::percent,breaks=c(0,.1,.2,.3),expand = expand_scale(mult=c(0,.11))) +
  labs(x="Age group",y="Proportion") +
  theme(axis.text.x=element_text(angle=45,hjust=1))
g1

# Italy
source("data/italy/data_management_italy.R")
itcases = data.frame(cases=D_S_model13ITB$incidence_cases,
                     time=1:D_S_model13ITB$D,
                     deaths=D_S_model13ITB$incidence_deaths) %>%
  mutate(date=day_data+time)

h2 = ggplot(itcases) +
  geom_col(aes(x=date,y=cases),fill=col.cases1,width=1,col="black") +
  geom_col(aes(x=date,y=deaths),fill=col.deaths1,width=1,col="black") +
  geom_vline(xintercept=day_max+.5,linetype=2) +
  scale_y_continuous(expand=expand_scale(mult=c(0,.10))) +
  scale_x_date(limits=c(as.Date("2020-02-06"),day_max+.5)) +
  labs(x="Date",y="N") +
  annotate("segment",x=as.Date("2020-02-20"),xend=as.Date("2020-02-20"),y=500,yend=100,arrow=arrow(length=unit(.2,"cm"))) +
  annotate("label",x=as.Date("2020-02-20"),y=500,label="Local transmission\nidentified") +
  annotate("segment",x=as.Date("2020-02-20"),xend=as.Date("2020-03-03")+.5,y=900,yend=900,arrow=arrow(length=unit(.2,"cm"))) +
  annotate("label",x=as.Date("2020-02-20"),y=900,label="End of data collection") 
h2

dd2 = data.frame(age_class=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),
                 pop = D_S_model13ITB$age_dist,
                 cases_tmax=D_S_model13ITB$agedistr_cases/sum(D_S_model13ITB$agedistr_cases),
                 mort_tmax=D_S_model13ITB$agedistr_deaths/sum(D_S_model13ITB$agedistr_deaths)) 
g2 = select(dd2,age_class,pop,cases_tmax,mort_tmax) %>%
  gather("var","value",2:4) %>%
  mutate(var2=factor(var,levels=c("pop","cases_tmax","mort_tmax"),labels=c("Italian population","Reported cases","Reported deaths"))) %>%
  ggplot() +
  geom_col(aes(x=age_class,y=value,fill=var2),colour="black",width=1) +
  facet_wrap(~var2) +
  geom_hline(yintercept=0) +
  scale_fill_manual(values=c("seagreen",col.cases1,col.deaths1),guide=FALSE) +
  coord_flip() +
  scale_y_continuous(labels=scales::percent,breaks=c(0,.2,.4,.6),expand = expand_scale(mult=c(0,.11))) +
  labs(x="Age group",y="Proportion") +
  theme(axis.text.x=element_text(angle=45,hjust=1))
g2

ggarrange(h1,g1,h2,g2,nrow=2,ncol=2,widths=c(1.3,1),labels=LETTERS)

ggsave(file=paste0(figpath,"data_fig.pdf"),width=10,height=5.5)


# Figure comorbidities ----
# China
source("data/china/data_management_china.R")

com = data.frame(age_group=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),
                 diabetes=dist_diabetes_pop,
                 crd=dist_crd_pop,
                 cvd=dist_cvd_pop,
                 hypertension=dist_hypertension_pop) %>%
  tbl_df() %>%
  gather("comorbidity","prev",2:5) %>%
  mutate(com2=factor(comorbidity,labels=c("Chronic respiratory disease","Cardio-vascular disease","Diabetes","Hypertension")))
ggplot(com) +
  geom_col(aes(x=age_group,y=prev,fill=com2),colour="black",width=1) +
  facet_wrap(~com2) +
  scale_fill_discrete(guide=F) +
  scale_y_continuous(expand=expand_scale(mult=c(0,0.05)),labels=scales::percent) +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  labs(x=NULL,y="General population prevalence in China")
ggsave(file=paste0(figpath,"comorbidity_china.pdf"),width=8,height=4.5)


# Trace plots ----

# China
stan_trace(S_model13B,pars=c("beta","eta","xi","nu","psi","pi","rho","epsilon"))
ggsave(file=paste0(figpath,"traceplot_china.pdf"),width=10,height=6)

# Italy
stan_trace(S_model13ITB,pars=c("beta","psi","pi","rho","epsilon"))
ggsave(file=paste0(figpath,"traceplot_italy.pdf"),width=10,height=6)
