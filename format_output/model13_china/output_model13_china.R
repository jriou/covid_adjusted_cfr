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
l = load("posterior_samples/model13_2020-03-13.Rdata")

print(S_model13B,"rho")

# Text

#crude CFR
sum(D_S_model13B$incidence_deaths)/sum(D_S_model13B$incidence_cases) *100
## reduction in transmissibility after 20 Jan
1 - summary(S_model13B,"eta",digits_summary = 5)[[1]]
## total symptomatic cases
print(S_model13B,"predicted_total_overall_symptomatic_cases")
summary(S_model13B,"predicted_total_overall_symptomatic_cases")[[1]] / sum(D_S_model13B$incidence_cases) 
## reporting rate by age
print(S_model13B,"beta")
print(S_model13B,"rho")
## proportion symptomatic (data)
print(S_model13B,"psi")
## total cases
print(S_model13B,"predicted_total_overall_all_cases")
## CFR among symptomatics
print(S_model13B,"cfr_D_symptomatic",digits_summary = 4)
## CFR among all
print(S_model13B,"cfr_D_all",digits_summary = 4)
## CFR among all by age
print(S_model13B,"cfr_D_all_by_age",digits_summary = 4)
## CFR among symptomatics by age
print(S_model13B,"cfr_D_symptomatic_by_age",digits_summary = 4)
summary(S_model13B,"cfr_D_symptomatic_by_age",digits_summary = 4)[[1]] * 100

# Figure 1: data ----

h1 = ggplot(confirmed_cases) +
  geom_col(aes(x=date,y=confirmed_cases_hubei),fill=col.cases1,width=1,col="black") +
  geom_vline(xintercept=day_quarantine+.5,linetype=2) +
  geom_vline(xintercept=day_max+.5,linetype=2) +
  scale_y_continuous(expand=expand_scale(mult=c(0,.05)),
                     labels=function(x) paste0(x/1000,"K")) +
  scale_x_date(limits=c(as.Date("2020-01-01"),day_max+1)) +
  labs(x="Date",y="N") #+
  # annotate("segment",x=as.Date("2019-12-15"),xend=as.Date("2019-12-08"),y=900,yend=5,arrow=arrow(length=unit(.2,"cm"))) +
  # annotate("label",x=as.Date("2019-12-15"),y=900,label="First case") +
  # annotate("segment",x=as.Date("2020-01-01"),xend=as.Date("2020-01-01"),y=900,yend=5,arrow=arrow(length=unit(.2,"cm"))) +
  # annotate("label",x=as.Date("2020-01-01"),y=900,label="Market closure") +
  # annotate("segment",x=as.Date("2020-01-03"),xend=as.Date("2020-01-20"),y=2150,yend=2150,arrow=arrow(length=unit(.2,"cm"))) +
  # annotate("label",x=as.Date("2020-01-03"),y=2150,label="Control measures\nimplemented") +
  # annotate("segment",x=as.Date("2020-01-03"),xend=as.Date("2020-02-11"),y=2750,yend=2750,arrow=arrow(length=unit(.2,"cm"))) +
  # annotate("label",x=as.Date("2020-01-03"),y=2750,label="End of data collection") 
h1
confirmed_cases$deaths= incidence_deaths
h2 = ggplot(confirmed_cases) +
  geom_col(aes(x=date,y=deaths),fill=col.deaths1,col="black",width=1) +
  geom_vline(xintercept=day_quarantine+.5,linetype=2) +
  geom_vline(xintercept=day_max+.5,linetype=2) +
  scale_y_continuous(expand=expand_scale(mult=c(0,.05))) +
  scale_x_date(limits=c(as.Date("2020-01-01"),day_max+1)) +
  labs(x="Date",y="N") 
h2

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

contacts = matrix(D_S_model13A$contact,nrow=9,byrow=T) %>%
  data.frame() %>%
  tbl_df()
names(contacts) = c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+")
contacts$i = c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+")
g2 = gather(contacts,"j","c",1:9) %>%
  ggplot() +
  geom_tile(aes(x=i,y=j,fill=c)) +
  scale_fill_gradient(low="darkblue",high="darkgoldenrod1") +
  labs(x="Age group",y="Age group",fill="Contacts")+
  theme(axis.text.x=element_text(angle=45,hjust=1))

ggarrange(h1,g1,h2,g2,nrow=2,ncol=2,widths=c(1.3,1),labels=LETTERS)

ggsave(file=paste0(figpath,"data_china2.pdf"),width=10,height=5.5)


# Figure 2: model fit ----

i1 = plot_incidence_cases(S_model13B,D_S_model13B,start_date=as.Date("2019-12-31"),end_date=as.Date("2020-03-20"),
                          col1=col.cases1,col2=col.cases2,col3=col.cases3)
i2 = plot_total_cases(S_model13B,D_S_model13B,
                          col1=col.cases1,col2=col.cases2,col3=col.cases3)
i3 = plot_agedist_cases(S_model13B,D_S_model13B,
                      col1=col.cases1,col2=col.cases2,col3=col.cases3)

j1 = plot_incidence_deaths2(S_model13B,D_S_model13B,start_date=as.Date("2019-12-31"),end_date=as.Date("2020-03-20"),
                           col1=col.deaths1,col2=col.deaths2)
j2 = plot_total_deaths(S_model13B,D_S_model13B,start_date=as.Date("2019-12-31"),end_date=as.Date("2020-03-20"),
                           col1=col.deaths1,col2=col.deaths3)
j3 = plot_agedist_deaths(S_model13B,D_S_model13B,
                        col1=col.deaths1,col2=col.deaths3)

# Build legend
data_legend = data.frame(type=c("Reported cases","Symptomatic cases","Symptomatic and\nasymptomatic cases","Reported deaths","Projected deaths","Overall deaths"),
                         col=c(col.cases1,col.deaths1,col.cases2,col.deaths2,col.cases3,col.deaths3),
                         x=1:6,y=1:6,stringsAsFactors = FALSE)
leg = ggplot(data_legend) +
  geom_point(aes(x=x,y=y,fill=factor(type,levels=c("Reported cases","Reported deaths","Symptomatic cases","Projected deaths","Symptomatic and\nasymptomatic cases","Overall deaths"))),
             colour="black",shape=22,size=6) +
  scale_fill_manual(values=data_legend$col,guide=guide_legend()) +
  labs(fill=NULL) +  
  theme(legend.position="bottom",
        legend.text=element_text(size=12),
        legend.key.height=unit(12,"pt"),
        legend.margin=margin(4,5,4,4,"pt"),
        legend.box.background = element_rect(fill="white",colour = "white",size=.4)
  )
legend <- cowplot::get_legend(leg)
gnull = ggplot() +  geom_blank() + theme_minimal()


ggarrange(ggarrange(i1,i2,i3,j1,j2,j3,nrow=2,ncol=3,widths = c(3,1,2),labels=LETTERS),
          legend,nrow=2,heights=c(6,1))
ggsave(file=paste0(figpath,"modelfit_china.pdf"),width=10,height=5)

  
hh = filter(laterdeaths,date2>day_max) %>%
  mutate(deaths=ifelse(deaths==0,3,deaths))
j1bis = j1 + 
  geom_point(data=hh,aes(x=date2,y=deaths),shape=24,fill="grey")
ggarrange(ggarrange(i1,i2,i3,j1bis,j2,j3,nrow=2,ncol=3,widths = c(3,1,2),labels=LETTERS),
          legend,nrow=2,heights=c(6,1))
ggsave(file=paste0(figpath,"modelfit_china_bis.pdf"),width=10,height=5)



# Figure 3: age specific CFR

k1 = plot_agedist_cfr(S_model13B,D_S_model13B,col1="tomato",col2=col.deaths2,col3=col.cases3,
                      insert=c(.7,7,.18,.45))
k1

k2 = plot_cfr(S_model13B,D_S_model13B,col1="tomato",col2=col.deaths2,col3=col.cases3)


ggarrange(k1,k2,widths = c(3.5,1),labels=LETTERS)

ggsave(file=paste0(figpath,"cfr_china.pdf"),width=8,height=4)

# Figure 3bis: with comorbidities

diab = data.frame(type="Diabetes",
                  expected_cases=sum(dist_diabetes_pop*D_S_model13B$agedistr_cases/sum(D_S_model13B$agedistr_cases)),
                  expected_deaths=sum(dist_diabetes_pop*D_S_model13B$agedistr_deaths/sum(D_S_model13B$agedistr_deaths)),
                  cases_k=diabetes_tmax[[1]],
                  cases_n= comorbidities_total[[1]],
                  deaths_k=diabetes_tmax[[2]],
                  deaths_n=comorbidities_total[[2]]
                  )
crd = data.frame(type="CRD",
                 expected_cases=sum(dist_crd_pop*D_S_model13B$agedistr_cases/sum(D_S_model13B$agedistr_cases)),
                 expected_deaths=sum(dist_crd_pop*D_S_model13B$agedistr_deaths/sum(D_S_model13B$agedistr_deaths)),
                  cases_k=crd_tmax[[1]],
                  cases_n= comorbidities_total[[1]],
                  deaths_k=crd_tmax[[2]],
                  deaths_n=comorbidities_total[[2]]
)
cvd = data.frame(type="CVD",
                 expected_cases=sum(dist_cvd_pop*D_S_model13B$agedistr_cases/sum(D_S_model13B$agedistr_cases)),
                 expected_deaths=sum(dist_cvd_pop*D_S_model13B$agedistr_deaths/sum(D_S_model13B$agedistr_deaths)),
                 cases_k=cvd_tmax[[1]],
                 cases_n= comorbidities_total[[1]],
                 deaths_k=cvd_tmax[[2]],
                 deaths_n=comorbidities_total[[2]]
)
hpt = data.frame(type="Hypertension",
                 expected_cases=sum(dist_hypertension_pop*D_S_model13B$agedistr_cases/sum(D_S_model13B$agedistr_cases)),
                 expected_deaths=sum(dist_hypertension_pop*D_S_model13B$agedistr_deaths/sum(D_S_model13B$agedistr_deaths)),
                 cases_k=hypertension_tmax[[1]],
                 cases_n= comorbidities_total[[1]],
                 deaths_k=hypertension_tmax[[2]],
                 deaths_n=comorbidities_total[[2]]
)
cc = bind_rows(diab,crd,cvd,hpt) %>%
  tbl_df() %>%
  mutate(p_cases=qbeta(0.5,cases_k+1,cases_n-cases_k+1),
         pinf_cases=qbeta(0.025,cases_k+1,cases_n-cases_k+1),
         psup_cases=qbeta(0.975,cases_k+1,cases_n-cases_k+1),
         p_deaths=qbeta(0.5,deaths_k+1,deaths_n-deaths_k+1),
         pinf_deaths=qbeta(0.025,deaths_k+1,deaths_n-deaths_k+1),
         psup_deaths=qbeta(0.975,deaths_k+1,deaths_n-deaths_k+1)
  )

k3 = ggplot(cc) +
  geom_col(aes(x=type,y=expected_deaths),fill="white",colour="black") +
  geom_pointrange(aes(x=type,y=p_deaths,ymin=pinf_deaths,ymax=psup_deaths),colour=col.deaths1) +
  scale_y_continuous(expand=expand_scale(mult=c(0,0.05)),labels=scales::percent) +
  scale_x_discrete(labels=c("Chronic respiratory dis.","Cardio-vascular dis.","Diabetes","Hypertension")) +
  labs(x=NULL,y="Proportion") +
  theme(axis.text.x=element_text(angle=45,hjust=1))
k3

ggarrange(k1,k2,k3,nrow=1,ncol=3,widths = c(3.5,1.4,1.8),labels=LETTERS)

ggsave(file=paste0(figpath,"cfr_china_bis.pdf"),width=9,height=3.8)



# supp table
summary(S_model13B,pars=c("beta","psi","pi","eta","nu","xi"))[[1]][,c(6,4,8)] %>%
  as.data.frame() %>%
  tbl_df() %>%
  mutate(par=c("beta","psi","pi","eta","nu","xi"),
    prev=paste0(round(`50%`,3)," [",round(`2.5%`,3),"-",round(`97.5%`,3),"]")) %>%
  select(par,prev)  %>%
  xtable(.) %>%
  print(.,include.rownames=FALSE)


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
  xtable(.) %>%
  print(.,include.rownames=FALSE)
