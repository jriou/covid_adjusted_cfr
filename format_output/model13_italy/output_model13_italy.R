# Setup
library(rstan)
library(tidyverse)
library(lubridate)
library(cowplot)
theme_set(theme_bw())
library(ggpubr)

# Settings
source("format_output/functions_model13.R")
figpath = "/home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/format_output/figures/"
col.cases1 = "darkcyan"
col.cases2 = "chartreuse3"
col.cases3 = "deepskyblue2"
col.deaths1 = "firebrick"
col.deaths2 =  "gold"
col.deaths3 = "darkorange"

# Load posterior samples (first extract chains and save them in run_model13IT_ITBly.R)
l = load("posterior_samples/model13IT_2020-03-18.Rdata")


sum(D_S_model13ITB$incidence_deaths)/sum(D_S_model13ITB$incidence_cases) *100
## total symptomatic cases
sum(D_S_model13ITB$incidence_cases)
print(S_model13ITB,"predicted_total_overall_symptomatic_cases")
summary(S_model13ITB,"predicted_total_overall_symptomatic_cases")[[1]] / sum(D_S_model13ITB$incidence_cases) 
## reporting rate by age
print(S_model13ITB,"rho")
## proportion symptomatic (data)
print(S_model13ITB,"psi")
## total cases
print(S_model13ITB,"predicted_total_overall_all_cases")
## total cases
print(S_model13ITB,"predicted_total_overall_deaths_delay")
## reported deaths
sum(D_S_model13ITB$incidence_deaths) 
## CFR among symptomatics
print(S_model13ITB,"cfr_D_symptomatic",digits_summary = 4)
## CFR among all
print(S_model13ITB,"cfr_D_all",digits_summary = 4)
## CFR among all by age
summary(S_model13ITB,"epsilon",digits_summary = 4)[[1]][,c(6,4,8)]*100
## comorbidities
exp(summary(S_model15B1,c("kappa","psi"))[[1]][,c("50%","2.5%","97.5%")])

# Figure 1: data ----
source("data/italy/data_management_italy.R")

itcases = data.frame(cases=D_S_model13ITB$incidence_cases,
                     time=1:D_S_model13ITB$D) %>%
  mutate(date=day_data+time)
itdeaths = data.frame(deaths=D_S_model13ITB$incidence_deaths,
                     time=1:D_S_model13ITB$D) %>%
  mutate(date=day_data+time)

h1 = ggplot(itcases) +
  geom_col(aes(x=date,y=cases),fill=col.cases1,width=1,col="black") +
  geom_vline(xintercept=day_max+.5,linetype=2) +
  scale_y_continuous(expand=expand_scale(mult=c(0,.05))) +
  scale_x_date(limits=c(as.Date("2020-02-06"),day_max+.5)) +
  labs(x="Date",y="N")
h1
h2 = ggplot(itdeaths) +
  geom_col(aes(x=date,y=deaths),fill=col.deaths1,col="black",width=1) +
  geom_vline(xintercept=day_max+.5,linetype=2) +
  scale_y_continuous(expand=expand_scale(mult=c(0,.05))) +
  scale_x_date(limits=c(as.Date("2020-02-06"),day_max+.5)) +
  labs(x="Date",y="N") 
h2

dd2 = data.frame(age_class=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),
                 pop = D_S_model13ITB$age_dist,
                 cases_tmax=D_S_model13ITB$agedistr_cases/sum(D_S_model13ITB$agedistr_cases),
                 mort_tmax=D_S_model13ITB$agedistr_deaths/sum(D_S_model13ITB$agedistr_deaths)) 
g1 = select(dd2,age_class,pop,cases_tmax,mort_tmax) %>%
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
g1

contacts = matrix(D_S_model13ITB$contact,nrow=9,byrow=T) %>%
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

ggsave(file=paste0(figpath,"data_italy.pdf"),width=10,height=5.5)


# Figure 2: model fit ----
D_S_model13ITB$tswitch = D_S_model13ITB$D

i1 = plot_incidence_cases(S_model13ITB,D_S_model13ITB,start_date=day_start,end_date=day_max+60,
                          col1=col.cases1,col2=col.cases2,col3=col.cases3)
i1
i2 = plot_total_cases(S_model13ITB,D_S_model13ITB,
                      col1=col.cases1,col2=col.cases2,col3=col.cases3)
i2
i3 = plot_agedist_cases(S_model13ITB,D_S_model13ITB,
                        col1=col.cases1,col2=col.cases2,col3=col.cases3)
i3

j1 = plot_incidence_deaths2(S_model13ITB,D_S_model13ITB,start_date=day_start,end_date=day_max+60,
                           col1=col.deaths1,col2=col.deaths2)
j1
j2 = plot_total_deaths(S_model13ITB,D_S_model13ITB,
                       col1=col.deaths1,col2=col.deaths3)
j2
j3 = plot_agedist_deaths(S_model13ITB,D_S_model13ITB,
                         col1=col.deaths1,col2=col.deaths3)
j3

# Build legend
data_legend = data.frame(type=c("Reported cases","Symptomatic cases","All infections","Reported deaths","Projected deaths","Overall deaths"),
                         col=c(col.cases1,col.deaths1,col.cases2,col.deaths2,col.cases3,col.deaths3),
                         x=1:6,y=1:6,stringsAsFactors = FALSE)
leg = ggplot(data_legend) +
  geom_point(aes(x=x,y=y,fill=factor(type,levels=c("Reported cases","Reported deaths","Symptomatic cases","Projected deaths","All infections","Overall deaths"))),
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

ggsave(file=paste0(figpath,"modelfit_italy.pdf"),width=10,height=5)


# Figure 3: age specific CFR

k1 = plot_agedist_cfr(S_model13ITB,D_S_model13ITB,col1="tomato",col2=col.deaths2,col3=col.cases3)

k2 = plot_cfr(S_model13ITB,D_S_model13ITB,col1="tomato",col2=col.deaths2,col3=col.cases3)
k3 = ggplot() + geom_blank() + theme_minimal()

ggarrange(k1,k2,k3,nrow=1,ncol=3,widths = c(3.5,1.4,1.8),labels=c("A","B"))
ggarrange(k1,k2,widths = c(3.5,1.2),labels=c("A","B"))

ggsave(file=paste0(figpath,"cfr_italy.pdf"),width=9,height=3.8)





# supp table
summary(S_model13ITB,pars=c("beta","psi","pi"))[[1]][,c(6,4,8)] %>%
  as.data.frame() %>%
  tbl_df() %>%
  mutate(par=c("beta","psi","pi"),
         prev=paste0(round(`50%`,3)," [",round(`2.5%`,3),"-",round(`97.5%`,3),"]")) %>%
  select(par,prev)  %>%
  xtable(.) %>%
  print(.,include.rownames=FALSE)


rho = summary(S_model13ITB,pars="rho")[[1]][,c(6,4,8)] %>%
  as.data.frame() %>%
  tbl_df() %>%
  mutate(age_group=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),
         rho=paste0(round(`50%`,2)," [",round(`2.5%`,2),"-",round(`97.5%`,2),"]")) %>%
  select(age_group,rho)
epsilon = summary(S_model13ITB,pars="epsilon")[[1]][,c(6,4,8)] %>%
  as.data.frame() %>%
  tbl_df() %>%
  mutate(age_group=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),
         epsilon=paste0(round(`50%`,3)," [",round(`2.5%`,3),"-",round(`97.5%`,3),"]")) %>%
  select(epsilon)

bind_rows(data.frame(age_group="General",rho="-",epsilon="-"),
          bind_cols(rho,epsilon)) %>%
  xtable(.) %>%
  print(.,include.rownames=FALSE)
