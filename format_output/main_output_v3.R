# Setup ----
library(rstan)
library(tidyverse)
library(lubridate)
library(cowplot)
library(ggpubr)
theme_set(theme_bw())

# Settings ----
source("format_output/functions_model16.R")
figpath = "/home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/format_output/figures_v3/"
col.cases1 = "darkcyan"
col.cases2 = "chartreuse3"
col.cases3 = "deepskyblue2"
col.deaths1 = "firebrick"
col.deaths2 =  "gold"
col.deaths3 = "darkorange"

# Load posterior samples (first extract chains and save them in run_model13_china.R) ----
l = load("posterior_samples/posterior_samples_HU_2020-05-04.Rdata")
load("posterior_samples/posterior_samples_CH_2020-05-04.Rdata")
load("posterior_samples/posterior_samples_LO_2020-05-05.Rdata")
load("posterior_samples/posterior_samples_AU_2020-05-04.Rdata")
load("posterior_samples/posterior_samples_BW_2020-05-04.Rdata")
load("posterior_samples/posterior_samples_BA_2020-05-04.Rdata")
load("posterior_samples/posterior_samples_SP_2020-05-04.Rdata")
load("posterior_samples/posterior_samples_supp_HU_2020-05-04.Rdata")

# Text ----

#crude CFR
sum(D_S_model16A$incidence_deaths)/sum(D_S_model16A$incidence_cases) *100
## reduction in transmissibility after 20 Jan
1 - summary(S_model16A,"eta",digits_summary = 5,probs = c(.5,0.025,0.975))[[1]]
## total cases
print(S_model16A,"predicted_total_overall_all_cases",probs = c(.5,0.025,0.975))
print(S_model16A,"predicted_total_overall_symptomatic_cases",probs = c(.5,0.025,0.975))
summary(S_model16A,"predicted_total_overall_symptomatic_cases",probs = c(.5,0.025,0.975))[[1]] / sum(D_S_model16A$incidence_cases) 
print(S_model16B,"predicted_total_overall_all_cases",probs = c(.5,0.025,0.975))
print(S_model16B,"predicted_total_overall_symptomatic_cases",probs = c(.5,0.025,0.975)) 
print(S_model16D,"predicted_total_overall_all_cases",probs = c(.5,0.025,0.975))
print(S_model16B,"predicted_total_overall_symptomatic_cases",probs = c(.5,0.025,0.975))
## reporting rate by age
print(S_model16A,"beta")
print(S_model16A,"rho",probs = c(.5,0.025,0.975))
print(S_model16B,"rho",probs = c(.5,0.025,0.975))

## proportion symptomatic (data)
print(S_model16A,"psi")
## total cases
print(S_model16A,"predicted_total_overall_all_cases",probs = c(.5,0.025,0.975))
## total deaths
print(S_model16A,"predicted_total_overall_deaths_delay",probs = c(.5,0.025,0.975))
print(S_model16B,"predicted_total_overall_deaths_delay",probs = c(.5,0.025,0.975))
## CFR among symptomatics
print(S_model16A,"cfr_D_symptomatic",digits_summary = 4,probs = c(.5,0.025,0.975))
print(S_model16B,"cfr_D_symptomatic",digits_summary = 4,probs = c(.5,0.025,0.975))
## CFR among all
print(S_model16A,"cfr_D_all",digits_summary = 4,probs = c(.5,0.025,0.975))
print(S_model16B,"cfr_D_all",digits_summary = 4,probs = c(.5,0.025,0.975))
## CFR among all by age
summary(S_model16A,"cfr_D_symptomatic_by_age",digits_summary = 4,probs = c(.5,0.025,0.975))[[1]]*100
## CFR among symptomatics by age
summary(S_model16A,"cfr_D_symptomatic_by_age",digits_summary = 4,probs = c(.5,0.025,0.975))[[1]] * 100
summary(S_model16B,"cfr_D_symptomatic_by_age",digits_summary = 4,probs = c(.5,0.025,0.975))[[1]] * 100
summary(S_model16A,"epsilon",digits_summary = 4,probs = c(.5,0.025,0.975))[[1]] * 100
# later correction of cases and deaths
1/D_S_model16B$p_underreport_cases
1/D_S_model16B$p_underreport_deaths
corr = 1/D_S_model16B$p_underreport_deaths/(1/D_S_model16B$p_underreport_cases)
summary(S_model16A,"cfr_D_all",digits_summary = 4,probs = c(.5,0.025,0.975))[[1]]*corr
print(S_model16B,"cfr_D_all",digits_summary = 4,probs = c(.5,0.025,0.975))

# Figure 1: data ----
# China
source("data/china/data_management_china.R")
confirmed_cases$deaths= D_S_model16A$incidence_deaths
h1 = ggplot(confirmed_cases) +
  geom_col(aes(x=date,y=confirmed_cases_hubei),fill=col.cases1,width=1,col="black") +
  geom_vline(xintercept=day_quarantine+.5,linetype=2) +
  geom_vline(xintercept=day_max+.5,linetype=2) +
  scale_y_continuous(expand=expand_scale(mult=c(0,.05)),
                     labels=function(x) paste0(x/1000,"K")) +
  scale_x_date(limits=c(as.Date("2019-12-06"),day_max+3)) +
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
h2 = ggplot(confirmed_cases) +
  geom_col(aes(x=date,y=deaths),fill=col.deaths1,width=1,col="black") +
  geom_vline(xintercept=day_quarantine+.5,linetype=2) +
  geom_vline(xintercept=day_max+.5,linetype=2) +
  scale_y_continuous(expand=expand_scale(mult=c(0,.05))) +
  scale_x_date(limits=c(as.Date("2019-12-06"),day_max+3)) +
  labs(x="Date",y="N") +
  annotate("segment",x=as.Date("2020-01-01"),xend=as.Date("2020-01-13"),y=50,yend=1,arrow=arrow(length=unit(.2,"cm"))) +
  annotate("label",x=as.Date("2020-01-01"),y=50,label="First death") 
h2

dd2 = data.frame(age_class=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),
                 pop = age_dist,
                 cases_tmax=cases_tmax,
                 mort_tmax=mort_tmax,
                 prop_cases_tmax=prop_cases_tmax) %>%
  mutate(prop_mort_tmax=prop_mort_tmax,
         rawcfr=mort_tmax/cases_tmax,
         rawcfr2=rawcfr*prop_cases_tmax)
llab = c(pop="General population",
         prop_cases_tmax="Reported cases",
         prop_mort_tmax="Reported deaths")
g1 = select(dd2,age_class,pop,prop_cases_tmax,prop_mort_tmax) %>%
  gather("var","value",2:4) %>%
  ggplot() +
  geom_col(aes(x=age_class,y=value,fill=var),colour="black",width=1) +
  facet_wrap(~var,nrow=1,labeller = labeller(var=llab)) +
  geom_hline(yintercept=0) +
  scale_fill_manual(values=c("seagreen",col.cases1,col.deaths1),guide=FALSE) +
  coord_flip() +
  scale_y_continuous(labels=scales::percent,breaks=c(0,.1,.2,.3),expand = expand_scale(mult=c(0,.11))) +
  labs(x="Age group",y="Proportion") +
  theme(axis.text.x=element_text(angle=45,hjust=1))
g1

g2 = data.frame(contacts=D_S_model16A$contact) %>%
  tbl_df() %>%
  mutate(age1=rep(c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),9),age2=rep(c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),each=9)) %>%
  ggplot() +
  geom_tile(aes(x=age2,y=age1,fill=contacts)) +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.key.height = unit(30,"pt"),
        legend.margin = margin(0,0,0,0,"pt")) +
  scale_fill_gradient(low="dodgerblue4",high="yellow",breaks=c(0,1,2,3)) +
  labs(x="Age group",y="Age group",fill=NULL)


ggarrange(h1,g1,h2,g2,nrow=2,ncol=2,labels=LETTERS,widths=c(1.4,1))

ggsave(file=paste0(figpath,"fig1.pdf"),width=10,height=5.5)


# Figure 2: fit ----

i1 = plot_incidence_cases(S_model16A,D_S_model16A,start_date = day_start,end_date = day_max)
i2 = plot_total_cases(S_model16A,D_S_model16A)
i3 = plot_agedist_cases(S_model16A,D_S_model16A)
i4 = plot_incidence_deaths(S_model16A,D_S_model16A,start_date = day_start,end_date = day_max+50)
i5 = plot_total_deaths(S_model16A,D_S_model16A)
i6 = plot_agedist_deaths(S_model16A,D_S_model16A)


data_legend = data.frame(type=c("Reported cases","Symptomatic cases","Infected cases","Reported deaths","Projected deaths","Total deaths"),
                         col=c(col.cases1,col.deaths1,col.cases2,col.deaths2,col.cases3,col.deaths3),
                         x=1:6,y=1:6,stringsAsFactors = FALSE)
leg = ggplot(data_legend) +
  geom_point(aes(x=x,y=y,fill=factor(type,levels=c("Reported cases","Reported deaths","Symptomatic cases","Projected deaths","Infected cases","Total deaths"))),
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


ggarrange(ggarrange(i1,i2,i3,i4,i5,i6,nrow=2,ncol=3,widths = c(3,1,2),labels=LETTERS),
          legend,nrow=2,heights=c(6,1))

ggsave(file=paste0(figpath,"fig2.pdf"),width=8,height=5)

hh = filter(laterdeaths,date2>day_max) %>%
  mutate(deaths=ifelse(deaths==0,3,deaths))
i4bis = i4 + 
  geom_point(data=hh,aes(x=date2,y=deaths),shape=24,fill="grey")
ggarrange(ggarrange(i1,i2,i3,i4bis,i5,i6,nrow=2,ncol=3,widths = c(3,1,2),labels=LETTERS),
          legend,nrow=2,heights=c(6,1))

ggsave(file=paste0(figpath,"external_validation_hubei.pdf"),width=8,height=5)



# Figure 3: additional results ----

j1_tt = summary(S_model16A,"rho")[[1]] %>%
  as.data.frame() %>%
  tbl_df() %>%
  mutate(age_class=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),
         age_class2=1:9-.2,
         type="Baseline",
         shape=c(rep("Estimated",8),"Fixed"))
j1_tt = summary(S_model16B,"rho")[[1]] %>%
  as.data.frame() %>%
  tbl_df() %>%
  mutate(age_class=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),
         age_class2=1:9,
         type="After correction",
         shape=c(rep("Estimated",8),"Fixed")) %>%
  bind_rows(j1_tt)
j1_tt = summary(S_model16D,"rho")[[1]] %>%
  as.data.frame() %>%
  tbl_df() %>%
  mutate(age_class=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),
         age_class2=1:9+.2,
         type="Lower susceptibility\nof children",
         shape=c(rep("Estimated",8),"Fixed"))%>%
  bind_rows(j1_tt)
j1_tt = mutate(j1_tt,type=factor(type,levels=c("Baseline","After correction","Lower susceptibility\nof children")))
j1 = ggplot(j1_tt) +
  geom_line(aes(x=age_class2,y=`50%`,colour=type)) +
  geom_pointrange(aes(x=age_class2,y=`50%`,ymin=`2.5%`,ymax=`97.5%`,colour=type,fill=type,shape=shape)) +
  scale_x_continuous(breaks=1:9,labels=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+")) +
  scale_y_continuous(labels=scales::percent,limits=c(0,1)) +
  scale_fill_manual(values=c("mediumorchid1","mediumorchid4","pink")) +
  scale_color_manual(values=c("mediumorchid1","mediumorchid4","pink")) +
  scale_shape_manual(values=c(21,22)) +
  labs(x="Age group",y="Ascertainment",fill=NULL,colour=NULL,shape=NULL) +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position = c(.75,.35),
        legend.margin = margin(0.1,0.1,0.1,0.1,unit="pt"),
        legend.spacing = unit(0,"pt"))
j1

j2_tt = summary(S_model16A,"epsilon")[[1]] %>%
  as.data.frame() %>%
  tbl_df() %>%
  mutate(age_class=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),
         age_class2=1:9-.2,
         type="Hubei")
j2_tt = summary(S_model16B,"epsilon")[[1]] %>%
  as.data.frame() %>%
  tbl_df() %>%
  mutate(age_class=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),
         age_class2=1:9,
         type="Hubei*") %>%
  bind_rows(j2_tt)
j2_tt = summary(S_model16D,"epsilon")[[1]] %>%
  as.data.frame() %>%
  tbl_df() %>%
  mutate(age_class=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),
         age_class2=1:9+.2,
         type="Hubei**") %>%
  bind_rows(j2_tt)
j2 = ggplot() +
  geom_line(data=j2_tt,aes(x=age_class2,y=`50%`,color=type)) +
  geom_pointrange(data=j2_tt,aes(x=age_class2,y=`50%`,ymin=`2.5%`,ymax=`97.5%`,color=type)) +
  scale_x_continuous(breaks=1:9,labels=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+")) +
  scale_y_continuous(labels=scales::percent,limits=c(0,0.5)) +
  scale_color_manual(values=c(col.cases2,"darkgreen","lightgreen"),labels=c("Baseline","After correction","Lower susceptibility\nof children","Fixed")) +
  labs(x="Age group",y="sCFR",color=NULL) +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position = c(.4,.7),
        legend.margin = margin(0.1,0.1,0.1,0.1,unit="pt") )
j2


mort_by_rho_80p = NULL
for(i in c(10,20,30,40,50,60,70,80,90)) {
  tcfr = sum(get(paste0("D_S_model16E",i))$incidence_deaths)/sum(get(paste0("D_S_model16E",i))$incidence_cases)
  dcfr = data.frame(rowname='crude_cfr',`2.5%`=tcfr,`50%`=tcfr,`97.5%`=tcfr,rho_80p=i)
  names(dcfr) = c("rowname","2.5%","50%","97.5%","rho_80p")
  mort_by_rho_80p = bind_rows(mort_by_rho_80p,
                              summary(get(paste0("S_model16E",i)),pars=c("cfr_A_symptomatic","cfr_B_symptomatic","cfr_C_symptomatic","cfr_D_symptomatic","cfr_C_all","cfr_D_all"))[[1]] %>%
                                as.data.frame() %>%
                                rownames_to_column() %>%
                                tbl_df() %>%
                                mutate(rho_80p=i),
                              dcfr)
}
tcfr = sum(D_S_model16A$incidence_deaths)/sum(D_S_model16A$incidence_cases)
dcfr = data.frame(rowname='crude_cfr',`2.5%`=tcfr,`50%`=tcfr,`97.5%`=tcfr,rho_80p=100)
names(dcfr) = c("rowname","2.5%","50%","97.5%","rho_80p")
mort_by_rho_80p = bind_rows(mort_by_rho_80p,
                            summary(S_model16A,pars=c("cfr_A_symptomatic","cfr_B_symptomatic","cfr_C_symptomatic","cfr_D_symptomatic","cfr_C_all","cfr_D_all"))[[1]] %>%
                              as.data.frame() %>%
                              rownames_to_column() %>%
                              tbl_df() %>%
                              mutate(rho_80p=100),
                            dcfr)
j3_tt = filter(mort_by_rho_80p,rowname %in% c("crude_cfr","cfr_D_symptomatic","cfr_D_all")) %>%
  mutate(type=factor(rowname,levels=c("crude_cfr","cfr_D_symptomatic","cfr_D_all"),
                     labels=c("CFR","sCFR","IFR"))) %>%
  mutate(rho2=ifelse(type=="CFR",rho_80p-1,rho_80p),
         rho2=ifelse(type=="IFR",rho_80p+1,rho2))
j3 = ggplot(j3_tt) +
  geom_line(aes(x=rho2,y=`50%`,color=type)) +
  geom_pointrange(aes(x=rho2,y=`50%`,ymin=`2.5%`,ymax=`97.5%`,color=type)) +
  scale_x_continuous(breaks=c(10,20,30,40,50,60,70,80,90,100),labels=function(x) paste0(x,"%")) +
  scale_y_continuous(labels=scales::percent,limits=c(0,.05)) +
  scale_color_manual(values=c(col.deaths1,col.cases2,col.cases3)) +
  labs(x="Ascertainment of individuals aged 80+",y="Proportion",fill=NULL,colour=NULL)+
  theme(legend.position = c(.2,.8),
        legend.margin = margin(0.1,0.1,0.1,0.1,unit="pt") )
j3

most_by_date = NULL
for(i in  42-c(5,10,15,20,25,30)) {
  tcfr = sum(get(paste0("D_S_model16F",i))$incidence_deaths)/sum(get(paste0("D_S_model16F",i))$incidence_cases)
  dcfr = data.frame(rowname='crude_cfr',`2.5%`=tcfr,`50%`=tcfr,`97.5%`=tcfr,date=i)
  names(dcfr) = c("rowname","2.5%","50%","97.5%","date")
  most_by_date = bind_rows(most_by_date,
                           summary(get(paste0("S_model16F",i)),pars=c("cfr_A_symptomatic","cfr_B_symptomatic","cfr_C_symptomatic","cfr_D_symptomatic","cfr_C_all","cfr_D_all"))[[1]] %>%
                             as.data.frame() %>%
                             rownames_to_column() %>%
                             tbl_df() %>%
                             mutate(date=i),
                           dcfr)
}
tcfr = sum(D_S_model16A$incidence_deaths)/sum(D_S_model16A$incidence_cases)
dcfr = data.frame(rowname='crude_cfr',`2.5%`=tcfr,`50%`=tcfr,`97.5%`=tcfr,date=42)
names(dcfr) = c("rowname","2.5%","50%","97.5%","date")
most_by_date = bind_rows(most_by_date,
                         summary(S_model16A,pars=c("cfr_A_symptomatic","cfr_B_symptomatic","cfr_C_symptomatic","cfr_D_symptomatic","cfr_C_all","cfr_D_all"))[[1]] %>%
                           as.data.frame() %>%
                           rownames_to_column() %>%
                           tbl_df() %>%
                           mutate(date=42),
                         dcfr)
j4_tt = filter(most_by_date,rowname %in% c("crude_cfr","cfr_D_symptomatic","cfr_D_all")) %>%
  mutate(type=factor(rowname,levels=c("crude_cfr","cfr_D_symptomatic","cfr_D_all"),
                     labels=c("CFR","sCFR","IFR"))) %>%
  mutate(date2=day_start+date)  %>%
  mutate(date3=ymd(ifelse(type=="CFR",as.character(date2-1),as.character(date2))),
         date3=ymd(ifelse(type=="IFR",as.character(date2+1),as.character(date3))))

j4 = ggplot(j4_tt) +
  # geom_ribbon(aes(x=date2,ymin=`2.5%`,ymax=`97.5%`,fill=type),alpha=.5) +
  # geom_line(aes(x=date2,y=`50%`,colour=type)) +
  geom_line(aes(x=date3,y=`50%`,color=type)) +
  geom_pointrange(aes(x=date3,y=`50%`,ymin=`2.5%`,ymax=`97.5%`,color=type)) +
  scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(values=c(col.deaths1,col.cases2,col.cases3)) +
  scale_x_date(breaks=day_max-c(0,5,10,15,20,25,30),
               date_labels = "%b %d") +
  labs(x="Date of analysis",y="Proportion",fill=NULL,colour=NULL)+
  scale_color_manual(values=c(col.deaths1,col.cases2,col.cases3)) +
  theme(legend.position = c(.8,.8),
        legend.margin = margin(0.1,0.1,0.1,0.1,unit="pt") )
j4

ggarrange(j1,j2,j3,j4,labels=LETTERS)

ggsave(file=paste0(figpath,"fig3.pdf"),width=8,height=6)

# Figure 4: other countries ----

S_model16AHU = S_model16A
D_S_model16AHU = D_S_model16A

coun_tt = c("AU","BA","BW","CH","HU","LO","SP")
count_txt = c("Austria","Bavaria (G)","Baden-\nWürttemberg (G)","Switzerland","Hubei (C)","Lombardy (I)","Spain")

k1_tt = NULL
for(i in 1:length(coun_tt)) {
  tcfr = sum(get(paste0("D_S_model16A",coun_tt[i]))$incidence_deaths)/sum(get(paste0("D_S_model16A",coun_tt[i]))$incidence_cases/get(paste0("D_S_model16A",coun_tt[i]))$p_underreport_cases)
  dcfr = data.frame(rowname='crude_cfr',`2.5%`=tcfr,`50%`=tcfr,`97.5%`=tcfr,countryID=coun_tt[i],countryname=count_txt[i],type="CFR")
  names(dcfr) = c("rowname","2.5%","50%","97.5%","countryID","countryname","type")
  k1_tt = bind_rows(k1_tt,summary(get(paste0("S_model16A",coun_tt[i])),pars=c("cfr_D_symptomatic","cfr_D_all"))[[1]] %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    tbl_df() %>%
    mutate(countryID=coun_tt[i],
           countryname=count_txt[i],
           type=c("sCFR","IFR")),
    dcfr)
}
k1_tt = mutate(k1_tt,type=factor(type,levels=c("CFR","sCFR","IFR")))

k1 = ggplot(k1_tt) +
  geom_pointrange(aes(x=countryname,y=`50%`,ymin=`2.5%`,ymax=`97.5%`,colour=type),position=position_dodge2(.4),size=.4) +
  scale_color_manual(values=c(col.deaths1,col.cases2,col.cases3)) +
  scale_y_continuous(labels=scales::percent) +
  labs(x="Location",y="Fatality rate",colour=NULL)+
  theme(legend.position = c(.15,.7),
        legend.margin = margin(0.1,0.1,0.1,0.1,unit="pt") ,
        axis.text.x=element_text(angle=45,hjust=1))
k1

k2_tt = NULL
for(i in 1:length(coun_tt)) {
  k2_tt = rbind(k2_tt,summary(get(paste0("S_model16A",coun_tt[i])),pars=c("cfr_D_all_by_age"))[[1]] %>%
                  as.data.frame() %>%
                  rownames_to_column() %>%
                  tbl_df() %>%
                  mutate(countryID=coun_tt[i],
                         countryname=count_txt[i],
                         age_class2=1:9))
}
# adapt to intervals for Austria
k2_tt = mutate(k2_tt,age_class2=ifelse(countryID=="AU" & age_class2!=1,age_class2-.5,age_class2))
filter(k2_tt,age_class2==9)
k2 = ggplot(k2_tt) +
  geom_line(aes(x=age_class2,y=`50%`,colour=countryname)) +
  geom_point(aes(x=age_class2,y=`50%`,colour=countryname)) +
  # geom_pointrange(aes(x=age_class2,y=`50%`,ymin=`2.5%`,ymax=`97.5%`,colour=countryname)) +
  scale_x_continuous(breaks=1:9,labels=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+")) +
  scale_y_continuous(labels=scales::percent) +
  labs(x="Age group",y="Infection fatality rate",colour=NULL)+
  theme(legend.position = c(.25,.55),
        legend.margin = margin(0.1,0.1,0.1,0.1,unit="pt") ,
        axis.text.x=element_text(angle=45,hjust=1))
k2


k3_tt = NULL
for(i in 1:length(coun_tt)) {
  k3_tt = rbind(k3_tt,summary(get(paste0("S_model16A",coun_tt[i])),pars=c("rho"))[[1]] %>%
                  as.data.frame() %>%
                  rownames_to_column() %>%
                  tbl_df() %>%
                  mutate(countryID=coun_tt[i],
                         countryname=count_txt[i],
                         age_class2=1:9,
                         type=c(rep("Estimated",8),"Fixed")))
}
# adapt to intervals for Austria
k3_tt = mutate(k3_tt,age_class2=ifelse(countryID=="AU" & age_class2!=1,age_class2-.5,age_class2))

k3 = ggplot(k3_tt) +
  geom_line(aes(x=age_class2,y=`50%`,colour=countryname)) +
  geom_point(aes(x=age_class2,y=`50%`,colour=countryname,fill=countryname,shape=type)) +
  scale_x_continuous(breaks=1:9,labels=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+")) +
  scale_y_continuous(labels=scales::percent) +
  scale_shape_manual(values=c(21,22)) +
  scale_colour_discrete(guide=FALSE) +
  scale_fill_discrete(guide=FALSE) +
  labs(x="Age group",y="Ascertainment",colour=NULL,shape=NULL)+
  theme(legend.position = c(.2,.8),
        legend.margin = margin(0.1,0.1,0.1,0.1,unit="pt") ,
        axis.text.x=element_text(angle=45,hjust=1))
k3

k4_tt = NULL
for(i in 1:length(coun_tt)) {
  k4_tt = rbind(k4_tt,
                data.frame(cases=get(paste0("D_S_model16A",coun_tt[i]))$agedistr_cases,
                           deaths=get(paste0("D_S_model16A",coun_tt[i]))$agedistr_deaths) %>%
                  mutate(countryID=coun_tt[i],
                         countryname=count_txt[i],
                         age_class2=1:9))
}
k4_tt = mutate(k4_tt,age_class2=ifelse(countryID=="AU" & age_class2!=1,age_class2-.5,age_class2))
k4_tt = group_by(k4_tt,countryID) %>%
  mutate(p_cases=cases/sum(cases),
         p_deaths=deaths/sum(deaths))
k4_tt = mutate(k4_tt,p_cases=ifelse(countryID=="BW",p_cases+.002,p_cases))

k4 = ggplot(k4_tt) +
  geom_line(aes(x=age_class2,y=p_cases,group=countryname,colour=countryname)) +
  geom_point(aes(x=age_class2,y=p_cases,group=countryname,colour=countryname)) +
  scale_fill_discrete(guide=F) +
  scale_x_continuous(breaks=1:9,labels=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+")) +
  scale_y_continuous(labels=scales::percent) +
  scale_colour_discrete(guide=FALSE) +
  labs(x="Age group",y="Proportion of cases",colour=NULL)+
  theme(legend.position = c(.2,.8),
        legend.margin = margin(0.1,0.1,0.1,0.1,unit="pt") ,
        axis.text.x=element_text(angle=45,hjust=1))
k4


ggarrange(k1,k2,k4,k3,nrow=2,ncol=2,labels=LETTERS)

ggsave(file=paste0(figpath,"fig4.pdf"),width=8,height=6)

  

# Table ----


S_model16AHU2 = S_model16B
D_S_model16AHU2 = D_S_model16B
S_model16AHU3 = S_model16D
D_S_model16AHU3 = D_S_model16D
coun_tt = c("HU","HU2","HU3","AU","BA","BW","CH","LO","SP")
count_txt = c("Hubei (C)","Hubei (C) after correction","Hubei children","Austria","Bavaria (G)","Baden-\nWürttemberg (G)","Switzerland","Lombardy (I)","Spain")


k1_tt = NULL
for(i in 1:length(coun_tt)) {
  tcfr = sum(get(paste0("D_S_model16A",coun_tt[i]))$incidence_deaths)/sum(get(paste0("D_S_model16A",coun_tt[i]))$incidence_cases/get(paste0("D_S_model16A",coun_tt[i]))$p_underreport_cases)
  dcfr = data.frame(rowname='crude_cfr',output=paste0(sprintf("%.1f",tcfr*100)),countryID=coun_tt[i],countryname=count_txt[i],type="CFR")
  tar = summary(get(paste0("S_model16A",coun_tt[i])),pars="predicted_total_overall_all_cases")[[1]] %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    tbl_df() %>%
    mutate(rowname="Attack rate",
           `2.5%`=`2.5%`/get(paste0("D_S_model16A",coun_tt[i]))$pop_t,
           `50%`=`50%`/get(paste0("D_S_model16A",coun_tt[i]))$pop_t,
           `97.5%`=`97.5%`/get(paste0("D_S_model16A",coun_tt[i]))$pop_t,
           countryID=coun_tt[i],
           countryname=count_txt[i],
           type="Attack rate",
           output= paste0(sprintf("%.1f",`50%`*100),"% (",sprintf("%.1f",`2.5%`*100),"-",sprintf("%.1f",`97.5%`*100),")"))
  k1_tt = bind_rows(k1_tt,summary(get(paste0("S_model16A",coun_tt[i])),pars=c("predicted_total_overall_all_cases","predicted_total_overall_deaths_delay","cfr_D_symptomatic","cfr_D_all"))[[1]] %>%
                  as.data.frame() %>%
                  rownames_to_column() %>%
                  tbl_df() %>%
                  mutate(countryID=coun_tt[i],
                         countryname=count_txt[i],
                         type=c("Total cases","Total deaths","SFR","IFR"),
                         output=ifelse(type %in% c("Total cases","Total deaths"),
                                paste0(signif(`50%`,3)," (",signif(`2.5%`,3),"-",signif(`97.5%`,3),")"),
                                paste0(sprintf("%.1f",`50%`*100),"% (",sprintf("%.1f",`2.5%`*100),"-",sprintf("%.1f",`97.5%`*100),")"))),
                tar,dcfr)
}
k1_tt = mutate(k1_tt,type=factor(type,levels=c("Total cases","Attack rate","Total deaths","CFR","SFR","IFR"))) 

k1_table = k1_tt %>%
  select(countryname,type,output) %>%
  spread(type,output) %>%
  xtable()
print(k1_table, include.rownames=F)

# Figure supplementary Fits -------------------------------

source('format_output/functions_model16.R')
data_legend = data.frame(type=c("Reported cases","Symptomatic cases","Infected cases","Reported deaths","Projected deaths","Total deaths"),
                         col=c(col.cases1,col.deaths1,col.cases2,col.deaths2,col.cases3,col.deaths3),
                         x=1:6,y=1:6,stringsAsFactors = FALSE)
leg = ggplot(data_legend) +
  geom_point(aes(x=x,y=y,fill=factor(type,levels=c("Reported cases","Reported deaths","Symptomatic cases","Projected deaths","Infected cases","Total deaths"))),
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

# Model fit in Austria ----
source("data/austria/data_management_austria.R")
ggarrange(ggarrange(plot_incidence_cases(S_model16AAU,D_S_model16AAU,start_date = day_start,end_date = day_max) + labs(title="Austria"),
                    plot_total_cases(S_model16AAU,D_S_model16AAU),
                    plot_agedist_cases(S_model16AAU,D_S_model16AAU),
                    plot_incidence_deaths(S_model16AAU,D_S_model16AAU,start_date = day_start,end_date = day_max+50),
                    plot_total_deaths(S_model16AAU,D_S_model16AAU),
                    plot_agedist_deaths(S_model16AAU,D_S_model16AAU),
                    nrow=2,ncol=3,widths = c(3,1,2),labels=LETTERS),
          legend,
          stan_trace(S_model16AAU,pars=c("beta","eta","xi","nu","psi","pi","rho","epsilon","phi")) + scale_x_continuous(breaks=c(0,250,500)),
          nrow=3,heights=c(6,1,8),labels=c("","","G"))
ggsave(file="format_output/figures_v3/supp_fit_austria.pdf",width=8,height=10)

# Model fit in BW ----
source("data/germany/data_management_badenw.R")
ggarrange(ggarrange(plot_incidence_cases(S_model16ABW,D_S_model16ABW,start_date = day_start,end_date = day_max) + labs(title="Baden-Württemberg"),
                    plot_total_cases(S_model16ABW,D_S_model16ABW),
                    plot_agedist_cases(S_model16ABW,D_S_model16ABW),
                    plot_incidence_deaths(S_model16ABW,D_S_model16ABW,start_date = day_start,end_date = day_max+50),
                    plot_total_deaths(S_model16ABW,D_S_model16ABW),
                    plot_agedist_deaths(S_model16ABW,D_S_model16ABW),
                    nrow=2,ncol=3,widths = c(3,1,2),labels=LETTERS),
          legend,
          stan_trace(S_model16ABW,pars=c("beta","eta","xi","nu","psi","pi","rho","epsilon","phi")) + scale_x_continuous(breaks=c(0,250,500)),
          nrow=3,heights=c(6,1,8),labels=c("","","G"))
ggsave(file="format_output/figures_v3/supp_fit_badenw.pdf",width=8,height=10)


# Model fit in BA ----
source("data/germany/data_management_bavaria.R")
ggarrange(ggarrange(plot_incidence_cases(S_model16ABA,D_S_model16ABA,start_date = day_start,end_date = day_max) + labs(title="Bavaria") ,
                    plot_total_cases(S_model16ABA,D_S_model16ABA),
                    plot_agedist_cases(S_model16ABA,D_S_model16ABA),
                    plot_incidence_deaths(S_model16ABA,D_S_model16ABA,start_date = day_start,end_date = day_max+50),
                    plot_total_deaths(S_model16ABA,D_S_model16ABA),
                    plot_agedist_deaths(S_model16ABA,D_S_model16ABA),
                    nrow=2,ncol=3,widths = c(3,1,2),labels=LETTERS),
          legend,
          stan_trace(S_model16ABA,pars=c("beta","eta","xi","nu","psi","pi","rho","epsilon","phi")) + scale_x_continuous(breaks=c(0,250,500)),
          nrow=3,heights=c(6,1,8),labels=c("","","G"))
ggsave(file="format_output/figures_v3/supp_fit_bavaria.pdf",width=8,height=10)

# Model fit in LO ----
source("data/lombardy/data_management_lombardy.R")
ggarrange(ggarrange(plot_incidence_cases(S_model16ALO,D_S_model16ALO,start_date = day_start,end_date = day_max) + labs(title="Lombardy"),
                    plot_total_cases(S_model16ALO,D_S_model16ALO),
                    plot_agedist_cases(S_model16ALO,D_S_model16ALO),
                    plot_incidence_deaths(S_model16ALO,D_S_model16ALO,start_date = day_start,end_date = day_max+50),
                    plot_total_deaths(S_model16ALO,D_S_model16ALO),
                    plot_agedist_deaths(S_model16ALO,D_S_model16ALO),
                    nrow=2,ncol=3,widths = c(3,1,2),labels=LETTERS),
          legend,
          stan_trace(S_model16ALO,pars=c("beta","eta","xi","nu","psi","pi","rho","epsilon","phi")) + scale_x_continuous(breaks=c(0,250,500)),
          nrow=3,heights=c(6,1,8),labels=c("","","G"))
ggsave(file="format_output/figures_v3/supp_fit_lombardy.pdf",width=8,height=10)

# Model fit in SP ----
source("data/spain/data_management_spain.R")
ggarrange(ggarrange(plot_incidence_cases(S_model16ASP,D_S_model16ASP,start_date = day_start,end_date = day_max) + labs(title="Spain") ,
                    plot_total_cases(S_model16ASP,D_S_model16ASP),
                    plot_agedist_cases(S_model16ASP,D_S_model16ASP),
                    plot_incidence_deaths(S_model16ASP,D_S_model16ASP,start_date = day_start,end_date = day_max+50),
                    plot_total_deaths(S_model16ASP,D_S_model16ASP),
                    plot_agedist_deaths(S_model16ASP,D_S_model16ASP),
                    nrow=2,ncol=3,widths = c(3,1,2),labels=LETTERS),
          legend,
          stan_trace(S_model16ASP,pars=c("beta","eta","xi","nu","psi","pi","rho","epsilon","phi")) + scale_x_continuous(breaks=c(0,250,500)),
          nrow=3,heights=c(6,1,8),labels=c("","","G"))
ggsave(file="format_output/figures_v3/supp_fit_spain.pdf",width=8,height=10)

# Model fit in CH ----
source("data/switzerland/data_management_switzerland.R")
ggarrange(ggarrange(plot_incidence_cases(S_model16ACH,D_S_model16ACH,start_date = day_start,end_date = day_max) + labs(title="Switzerland") ,
                    plot_total_cases(S_model16ACH,D_S_model16ACH),
                    plot_agedist_cases(S_model16ACH,D_S_model16ACH),
                    plot_incidence_deaths(S_model16ACH,D_S_model16ACH,start_date = day_start,end_date = day_max+50),
                    plot_total_deaths(S_model16ACH,D_S_model16ACH),
                    plot_agedist_deaths(S_model16ACH,D_S_model16ACH),
                    nrow=2,ncol=3,widths = c(3,1,2),labels=LETTERS),
          legend,
          stan_trace(S_model16ACH,pars=c("beta","eta","xi","nu","psi","pi","rho","epsilon","phi")) + scale_x_continuous(breaks=c(0,250,500)),
          nrow=3,heights=c(6,1,8),labels=c("","","G"))
ggsave(file="format_output/figures_v3/supp_fit_switzerland.pdf",width=8,height=10)



# Sensitivity analyses ----------------------------------------------------------------

l = load("posterior_samples/posterior_samples_supp_HU_2020-05-04.Rdata")



# Model fit in 16B ----
# corrected data

source("data/china/data_management_china.R")
plot_incidence_deaths_bis = function(samples,data_list,col1="firebrick",col2="gold",col3="darkorange",start_date,end_date) {
    t0 = data_list$t0
    tmax2 = data_list$S
    D = data_list$D
    G = data_list$G
    tswitch = data_list$tswitch
    S = data_list$S
    y = rstan::extract(samples,"y")[[1]]
    data_incidence_deaths = data.frame(time=1:D,incidence=data_list$incidence_deaths)
    predicted_overall_incidence_deaths = rstan::summary(samples,"predicted_overall_incidence_deaths")[[1]] %>%
      tbl_df() %>%
      mutate(time=1:(D+G),
             date=time+start_date) %>%
      left_join(data_incidence_deaths)
    predicted_overall_incidence_deaths_tmax = filter(predicted_overall_incidence_deaths,time<=D)
    ggplot() +
      geom_ribbon(data=predicted_overall_incidence_deaths,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill=col2,alpha=1) +
      geom_line(data=predicted_overall_incidence_deaths,aes(x=date,y=`50%`),linetype=3) +
      geom_ribbon(data=predicted_overall_incidence_deaths_tmax,aes(x=date,ymin=`2.5%`*data_list$p_underreport_death,ymax=`97.5%`*data_list$p_underreport_death),fill=col1,alpha=1) +
      geom_line(data=predicted_overall_incidence_deaths_tmax,aes(x=date,y=`50%`*data_list$p_underreport_death)) +
      geom_point(data=predicted_overall_incidence_deaths,aes(x=date,y=incidence),shape=21,fill="white") +
      coord_cartesian(xlim=c(start_date,end_date)) +
      labs(x="Time",y="Deaths per day") +
      scale_y_continuous(expand=expand_scale(mult=c(0,.05))) +
      geom_vline(xintercept=tswitch+start_date,linetype=2) +
      geom_vline(xintercept=tmax2+start_date,linetype=2) 
  }
  
ggarrange(ggarrange(plot_incidence_cases(S_model16B,D_S_model16B,start_date = day_start,end_date = day_max),
                    plot_total_cases(S_model16B,D_S_model16B),
                    plot_agedist_cases(S_model16B,D_S_model16B),
                    plot_incidence_deaths_bis(S_model16B,D_S_model16B,start_date = day_start,end_date = day_max+50),
                    plot_total_deaths(S_model16B,D_S_model16B),
                    plot_agedist_deaths(S_model16B,D_S_model16B),
                    nrow=2,ncol=3,widths = c(3,1,2),labels=LETTERS),
          legend,nrow=2,heights=c(6,1))
ggsave(file="format_output/figures_v3/supp_fit_16B.pdf",width=8,height=5)


# Model fit in 16D ----
# 50% susceptibility in children
ggarrange(ggarrange(plot_incidence_cases(S_model16D,D_S_model16D,start_date = day_start,end_date = day_max),
                    plot_total_cases(S_model16D,D_S_model16D),
                    plot_agedist_cases(S_model16D,D_S_model16D),
                    plot_incidence_deaths(S_model16D,D_S_model16D,start_date = day_start,end_date = day_max+50),
                    plot_total_deaths(S_model16D,D_S_model16D),
                    plot_agedist_deaths(S_model16D,D_S_model16D),
                    nrow=2,ncol=3,widths = c(3,1,2),labels=LETTERS),
          legend,nrow=2,heights=c(6,1))
ggsave(file="format_output/figures_v3/supp_fit_16D.pdf",width=8,height=5)

# Model fit in 16E ----
# lower ascertainment of 80+
ggarrange(ggarrange(plot_incidence_cases(S_model16E50,D_S_model16E50,start_date = day_start,end_date = day_max),
                    plot_total_cases(S_model16E50,D_S_model16E50),
                    plot_agedist_cases(S_model16E50,D_S_model16E50),
                    plot_incidence_deaths(S_model16E50,D_S_model16E50,start_date = day_start,end_date = day_max+50),
                    plot_total_deaths(S_model16E50,D_S_model16E50),
                    plot_agedist_deaths(S_model16E50,D_S_model16E50),
                    nrow=2,ncol=3,widths = c(3,1,2),labels=LETTERS),
          legend,nrow=2,heights=c(6,1))
ggsave(file="format_output/figures_v3/supp_fit_16E50.pdf",width=8,height=5)

ggarrange(ggarrange(plot_incidence_cases(S_model16E10,D_S_model16E10,start_date = day_start,end_date = day_max),
                    plot_total_cases(S_model16E10,D_S_model16E10),
                    plot_agedist_cases(S_model16E10,D_S_model16E10),
                    plot_incidence_deaths(S_model16E10,D_S_model16E10,start_date = day_start,end_date = day_max+50),
                    plot_total_deaths(S_model16E10,D_S_model16E10),
                    plot_agedist_deaths(S_model16E10,D_S_model16E10),
                    nrow=2,ncol=3,widths = c(3,1,2),labels=LETTERS),
          legend,nrow=2,heights=c(6,1))
ggsave(file="format_output/figures_v3/supp_fit_16E10.pdf",width=8,height=5)


coun_tt = c(10,20,30,40,50,60,70,80,90)
count_txt = paste0(coun_tt,"%")
k1_tt = NULL
for(i in 1:length(coun_tt)) {
  k1_tt = rbind(k1_tt,summary(get(paste0("S_model16E",coun_tt[i])),pars=c("predicted_total_overall_all_cases","predicted_total_overall_deaths_delay","cfr_A_symptomatic","cfr_D_symptomatic","cfr_D_all"))[[1]] %>%
                  as.data.frame() %>%
                  rownames_to_column() %>%
                  tbl_df() %>%
                  mutate(countryID=coun_tt[i],
                         countryname=count_txt[i],
                         p_underreport=c(1,1,get(paste0("D_S_model16E",coun_tt[i]))$p_underreport_cases,1,1),
                         median2=`50%`*p_underreport,
                         lower2=`2.5%`*p_underreport,
                         upper2=`97.5%`*p_underreport,
                         type=c("Total cases","Total deaths","CFR","SFR","IFR")))
}
k1_tt = mutate(k1_tt,type=factor(type,levels=c("Total cases","Total deaths","CFR","SFR","IFR")))

k1_table = k1_tt %>%
  mutate(output=
           ifelse(type %in% c("Total cases","Total deaths"),
                  paste0(signif(median2,3)," (",signif(lower2,3),"-",signif(upper2,3),")"),
                  paste0(sprintf("%.1f",median2*100),"% (",sprintf("%.1f",lower2*100),"-",sprintf("%.1f",upper2*100),")"))) %>%
  select(countryname,type,output) %>%
  spread(type,output) %>%
  xtable()
print(k1_table, include.rownames=F)





# Model fit in 16F ----
# shorter period




coun_tt = c(12,17,22,27,32,37)
count_txt = paste0(coun_tt,"%")
k1_tt = NULL
for(i in 1:length(coun_tt)) {
  k1_tt = rbind(k1_tt,summary(get(paste0("S_model16F",coun_tt[i])),pars=c("predicted_total_overall_all_cases","predicted_total_overall_deaths_delay","cfr_A_symptomatic","cfr_D_symptomatic","cfr_D_all"))[[1]] %>%
                  as.data.frame() %>%
                  rownames_to_column() %>%
                  tbl_df() %>%
                  mutate(countryID=coun_tt[i],
                         countryname=count_txt[i],
                         p_underreport=c(1,1,get(paste0("D_S_model16F",coun_tt[i]))$p_underreport_cases,1,1),
                         median2=`50%`*p_underreport,
                         lower2=`2.5%`*p_underreport,
                         upper2=`97.5%`*p_underreport,
                         type=c("Total cases","Total deaths","CFR","SFR","IFR")))
}
k1_tt = mutate(k1_tt,type=factor(type,levels=c("Total cases","Total deaths","CFR","SFR","IFR")))

k1_table = k1_tt %>%
  mutate(output=
           ifelse(type %in% c("Total cases","Total deaths"),
                  paste0(signif(median2,3)," (",signif(lower2,3),"-",signif(upper2,3),")"),
                  paste0(sprintf("%.1f",median2*100),"% (",sprintf("%.1f",lower2*100),"-",sprintf("%.1f",upper2*100),")"))) %>%
  select(countryname,type,output) %>%
  spread(type,output) %>%
  xtable()
print(k1_table, include.rownames=F)

ggarrange(ggarrange(plot_incidence_cases(S_model16F12,D_S_model16F12,start_date = day_start,end_date = day_max),
                    plot_total_cases(S_model16F12,D_S_model16F12),
                    plot_agedist_cases(S_model16F12,D_S_model16F12),
                    plot_incidence_deaths(S_model16F12,D_S_model16F12,start_date = day_start,end_date = day_max+50),
                    plot_total_deaths(S_model16F12,D_S_model16F12),
                    plot_agedist_deaths(S_model16F12,D_S_model16F12),
                    nrow=2,ncol=3,widths = c(3,1,2),labels=LETTERS),
          legend,nrow=2,heights=c(6,1))
ggsave(file="format_output/figures_v3/supp_fit_16F12.pdf",width=8,height=5)

ggarrange(ggarrange(plot_incidence_cases(S_model16F17,D_S_model16F17,start_date = day_start,end_date = day_max),
                    plot_total_cases(S_model16F17,D_S_model16F17),
                    plot_agedist_cases(S_model16F17,D_S_model16F17),
                    plot_incidence_deaths(S_model16F17,D_S_model16F17,start_date = day_start,end_date = day_max+50),
                    plot_total_deaths(S_model16F17,D_S_model16F17),
                    plot_agedist_deaths(S_model16F17,D_S_model16F17),
                    nrow=2,ncol=3,widths = c(3,1,2),labels=LETTERS),
          legend,nrow=2,heights=c(6,1))
ggsave(file="format_output/figures_v3/supp_fit_16F17.pdf",width=8,height=5)

ggarrange(ggarrange(plot_incidence_cases(S_model16F22,D_S_model16F22,start_date = day_start,end_date = day_max),
                    plot_total_cases(S_model16F22,D_S_model16F22),
                    plot_agedist_cases(S_model16F22,D_S_model16F22),
                    plot_incidence_deaths(S_model16F22,D_S_model16F22,start_date = day_start,end_date = day_max+50),
                    plot_total_deaths(S_model16F22,D_S_model16F22),
                    plot_agedist_deaths(S_model16F22,D_S_model16F22),
                    nrow=2,ncol=3,widths = c(3,1,2),labels=LETTERS),
          legend,nrow=2,heights=c(6,1))
ggsave(file="format_output/figures_v3/supp_fit_16F22.pdf",width=8,height=5)
ggarrange(ggarrange(plot_incidence_cases(S_model16F27,D_S_model16F27,start_date = day_start,end_date = day_max),
                    plot_total_cases(S_model16F27,D_S_model16F27),
                    plot_agedist_cases(S_model16F27,D_S_model16F27),
                    plot_incidence_deaths(S_model16F27,D_S_model16F27,start_date = day_start,end_date = day_max+50),
                    plot_total_deaths(S_model16F27,D_S_model16F27),
                    plot_agedist_deaths(S_model16F27,D_S_model16F27),
                    nrow=2,ncol=3,widths = c(3,1,2),labels=LETTERS),
          legend,nrow=2,heights=c(6,1))
ggsave(file="format_output/figures_v3/supp_fit_16F27.pdf",width=8,height=5)

ggarrange(ggarrange(plot_incidence_cases(S_model16F32,D_S_model16F32,start_date = day_start,end_date = day_max),
                    plot_total_cases(S_model16F32,D_S_model16F32),
                    plot_agedist_cases(S_model16F32,D_S_model16F32),
                    plot_incidence_deaths(S_model16F32,D_S_model16F32,start_date = day_start,end_date = day_max+50),
                    plot_total_deaths(S_model16F32,D_S_model16F32),
                    plot_agedist_deaths(S_model16F32,D_S_model16F32),
                    nrow=2,ncol=3,widths = c(3,1,2),labels=LETTERS),
          legend,nrow=2,heights=c(6,1))
ggsave(file="format_output/figures_v3/supp_fit_16F32.pdf",width=8,height=5)

ggarrange(ggarrange(plot_incidence_cases(S_model16F37,D_S_model16F37,start_date = day_start,end_date = day_max),
                    plot_total_cases(S_model16F37,D_S_model16F37),
                    plot_agedist_cases(S_model16F37,D_S_model16F37),
                    plot_incidence_deaths(S_model16F37,D_S_model16F37,start_date = day_start,end_date = day_max+50),
                    plot_total_deaths(S_model16F37,D_S_model16F37),
                    plot_agedist_deaths(S_model16F37,D_S_model16F37),
                    nrow=2,ncol=3,widths = c(3,1,2),labels=LETTERS),
          legend,nrow=2,heights=c(6,1))
ggsave(file="format_output/figures_v3/supp_fit_16F37.pdf",width=8,height=5)




# Model fit in 16G ----
# lower contribution presymptomatic to 30%
ggarrange(ggarrange(plot_incidence_cases(S_model16G,D_S_model16G,start_date = day_start,end_date = day_max),
                    plot_total_cases(S_model16G,D_S_model16G),
                    plot_agedist_cases(S_model16G,D_S_model16G),
                    plot_incidence_deaths(S_model16G,D_S_model16G,start_date = day_start,end_date = day_max+50),
                    plot_total_deaths(S_model16G,D_S_model16G),
                    plot_agedist_deaths(S_model16G,D_S_model16G),
                    nrow=2,ncol=3,widths = c(3,1,2),labels=LETTERS),
          legend,nrow=2,heights=c(6,1))
ggsave(file="format_output/figures_v3/supp_fit_16G.pdf",width=8,height=5)

# lower contribution presymptomatic to 0%
ggarrange(ggarrange(plot_incidence_cases(S_model16H,D_S_model16H,start_date = day_start,end_date = day_max),
                    plot_total_cases(S_model16H,D_S_model16H),
                    plot_agedist_cases(S_model16H,D_S_model16H),
                    plot_incidence_deaths(S_model16H,D_S_model16H,start_date = day_start,end_date = day_max+50),
                    plot_total_deaths(S_model16H,D_S_model16H),
                    plot_agedist_deaths(S_model16H,D_S_model16H),
                    nrow=2,ncol=3,widths = c(3,1,2),labels=LETTERS),
          legend,nrow=2,heights=c(6,1))
ggsave(file="format_output/figures_v3/supp_fit_16H.pdf",width=8,height=5)


# Model fit in 16I ---- 
# 15 days (+/- 10 from onset to death)
ggarrange(ggarrange(plot_incidence_cases(S_model16I,D_S_model16I,start_date = day_start,end_date = day_max),
                    plot_total_cases(S_model16I,D_S_model16I),
                    plot_agedist_cases(S_model16I,D_S_model16I),
                    plot_incidence_deaths(S_model16I,D_S_model16I,start_date = day_start,end_date = day_max+50),
                    plot_total_deaths(S_model16I,D_S_model16I),
                    plot_agedist_deaths(S_model16I,D_S_model16I),
                    nrow=2,ncol=3,widths = c(3,1,2),labels=LETTERS),
          legend,nrow=2,heights=c(6,1))
ggsave(file="format_output/figures_v3/supp_fit_16I.pdf",width=8,height=5)


# Additional results -----

additional_res = function(samples,place,age_group=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+")) {
  r1 = summary(samples, pars=c("beta","psi","pi","eta","nu","xi","phi"))[[1]]  %>%
    as.data.frame() %>%
    tbl_df() %>%
    mutate(par=c("$\\beta$","$\\psi$","$\\pi$","$\\eta$","$\\nu$","$\\xi$","$\\phi_1$","$\\phi_2$"),
           int=c("Probability of transmission upon contact",
                 "Proportion of symptomatic infections",
                 "Initial proportion of infected in the population",
                 "Reduction of transmission due to control measures",
                 "Delay of implementation of control measures (days)",
                 "Slope of implementation of control measures",
                 "Overdispersion parameter for cases",
                 "Overdispersion parameter for deaths"),
           nn = c(2,2,6,3,1,2,1,1),
           prev=paste0(sprintf(paste0("%.",nn,"f"),`50%`)," [",sprintf(paste0("%.",nn,"f"),`2.5%`),"-",sprintf(paste0("%.",nn,"f"),`97.5%`),"]")) %>%
    select(par,int,prev)
  names(r1) = c("Parameter","Interpretation","Posterior median (95\\%CrI)")
  print(xtable(r1,caption = paste0("Posterior distributions of the general parameters in ",place)),include.rownames =FALSE,sanitize.text.function=identity)
  r2 = summary(samples,pars="rho")[[1]]%>%
    as.data.frame() %>%
    tbl_df() %>%
    mutate(age_group=age_group,
           nn = 2,
           rho=paste0(sprintf(paste0("%.",nn,"f"),`50%`*100),"\\% [",sprintf(paste0("%.",nn,"f"),`2.5%`*100),"-",sprintf(paste0("%.",nn,"f"),`97.5%`*100),"]")) %>%
    select(age_group,rho)
  r3 = summary(samples,pars="epsilon")[[1]]%>%
    as.data.frame() %>%
    tbl_df() %>%
    mutate(age_group=age_group,
           nn = c(3,3,3,3,3,2,2,2,2),
           epsilon=paste0(sprintf(paste0("%.",nn,"f"),`50%`*100),"\\% [",sprintf(paste0("%.",nn,"f"),`2.5%`*100),"-",sprintf(paste0("%.",nn,"f"),`97.5%`*100),"]")) %>%
    select(epsilon)
  r4 = summary(samples,pars="cfr_D_all_by_age")[[1]]%>%
    as.data.frame() %>%
    tbl_df() %>%
    mutate(age_group=age_group,
           nn = c(3,3,3,3,3,2,2,2,2),
           ifr=paste0(sprintf(paste0("%.",nn,"f"),`50%`*100),"\\% [",sprintf(paste0("%.",nn,"f"),`2.5%`*100),"-",sprintf(paste0("%.",nn,"f"),`97.5%`*100),"]")) %>%
    select(ifr)
  rr = bind_cols(r2,r3,r4)
  names(rr) = c("Age group","Ascertainment ($\\rho_k$)","SFR ($\\epsilon_k$)","IFR ($\\epsilon_k\\cdot\\psi$)")
  print(xtable(rr,caption = paste0("Posterior distributions of the age-specific parameters in ",place)),include.rownames =FALSE,sanitize.text.function=identity)
}

additional_res(S_model16A,place="Hubei, China")
additional_res(S_model16AAU,place="Austria")
additional_res(S_model16ABW,place="Baden-Württemberg, Germany")
additional_res(S_model16ABA,place="Bavaria, Germany")
additional_res(S_model16ALO,place="Lombardy, Italy")
additional_res(S_model16ASP,place="Spain")
additional_res(S_model16ACH,place="Switzerland")


