source("data_management.R")
source("setup.R")

## Load chosen model ----
load("model/model10_2020-02-28.Rdata")

col.cases = "darkcyan"
col.cases.adj = "chartreuse3"
col.deaths = "firebrick"
col.deaths.adj =  "darkgoldenrod1"

## Description ----

dd = data.frame(death=incidence_deaths,
                date=day_data+0:(length(incidence_deaths)-1)) %>%
  tbl_df()

h1 = read.csv("data/confirmed_cases.csv") %>%
  tbl_df() %>%
  mutate(date=ymd(paste(year,month,day,sep="-"))) %>%
  left_join(dd) %>% 
  ggplot() +
  geom_col(aes(x=date,y=confirmed_cases_hubei),stat="identity",fill=col.cases,colour="black") +
  geom_col(aes(x=date,y=death),stat="identity",fill=col.deaths,colour="black") +
  geom_vline(xintercept=c(day_tswitch-3,day_tmax),linetype=2) +
  # annotate("label",x=c(day_tswitch,day_tmax),y=c(2900,2900),label=c("a","b")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(0,3050)) +
  # annotate("segment",x=as.Date("2019-11-20"),xend=as.Date("2019-12-04"),y=50,yend=50) +
  annotate("segment",x=as.Date("2019-12-15"),xend=as.Date("2019-12-08"),y=900,yend=5,arrow=arrow(length=unit(.2,"cm"))) +
  annotate("label",x=as.Date("2019-12-15"),y=900,label="First case") +
  annotate("segment",x=as.Date("2020-01-01"),xend=as.Date("2020-01-01"),y=900,yend=5,arrow=arrow(length=unit(.2,"cm"))) +
  annotate("label",x=as.Date("2020-01-01"),y=900,label="Market closure") +
  annotate("segment",x=as.Date("2020-01-03"),xend=as.Date("2020-01-20"),y=2150,yend=2150,arrow=arrow(length=unit(.2,"cm"))) +
  annotate("label",x=as.Date("2020-01-03"),y=2150,label="Control measures\nimplemented") +
  annotate("segment",x=as.Date("2020-01-03"),xend=as.Date("2020-02-11"),y=2750,yend=2750,arrow=arrow(length=unit(.2,"cm"))) +
  annotate("label",x=as.Date("2020-01-03"),y=2750,label="End of data collection") +
  labs(x="Time",y="N") 
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
         prop_cases_tmax="Confirmed cases",
         prop_mort_tmax="Reported deaths")
h2 = select(dd2,age_class,pop,prop_cases_tmax,prop_mort_tmax) %>%
  gather("var","value",2:4) %>%
  ggplot() +
  geom_col(aes(x=age_class,y=value,fill=var),colour="black",width=1) +
  facet_wrap(~var,labeller = labeller(var=llab)) +
  geom_hline(yintercept=0) +
  scale_fill_manual(values=c("seagreen",col.cases,col.deaths),guide=FALSE) +
  coord_flip() +
  scale_y_continuous(labels=scales::percent,breaks=c(0,.15,.3),expand = c(0,.02)) +
  labs(x=NULL,y="Proportion")
h2

diamond_asympto = data.frame(symp=c(0,2,25,27,19,28,76,95,29),
           asymp=c(1,3,3,7,8,31,101,139,25),
           age_class=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+")) %>%
  mutate(prop=symp/(symp+asymp),
         prop_inf=qbeta(0.025,symp+1,asymp+1),
         prop_sup=qbeta(0.975,symp+1,asymp+1)) %>%
  mutate(name="Symptomatics")

h3 = ggplot(diamond_asympto) +
  geom_col(aes(x=age_class,y=prop),fill="deepskyblue2",colour="black",width=1) +
  geom_errorbar(aes(x=age_class,ymin=prop_inf,ymax=prop_sup),width=0.2) +
  geom_hline(yintercept=.486,linetype=2) +
  coord_flip() +
  scale_y_continuous(labels=scales::percent,breaks=c(0,.5,1),expand = c(0,.06),limits=c(0,1)) +
  facet_wrap(~name) +
  labs(x=NULL,y="Proportion")
h3

plot_grid(h1,h2,h3,labels=c("A","B","C"),rel_widths = c(1.2,1,.4),nrow=1)
ggsave(file="figures/fig_desc.pdf",width=10.5,height=3.2)
ggsave(file="figures/fig_desc.png",width=10.5,height=3.2)

## Fit ----
check_hmc_diagnostics(S_model10)
print(S_model10,pars=c("beta","eta","epsilon","xi","nu","rho_K","pi","phi"),digits_summary=4)

round((1-summary(S_model10,pars=c("eta"),digits_summary=4)[[1]])*100,1)

g1 = plot_incidence_cases(S_model10,D_S_model10,col1=col.cases,col2=col.cases.adj,show.asympto=TRUE) +
  # annotate("label",x=c(day_tswitch,day_tmax),y=c(13000,13000),label=c("a","b")) +
  coord_cartesian(xlim=c(as.Date("2019-12-30"),as.Date("2020-03-31")),ylim=c(0,13000*1.03)) +
  scale_y_continuous(breaks=c(0,5000,10000),labels=c("0","5K","10K"))
g1
g3 = plot_incidence_deaths(S_model10,D_S_model10,col1=col.deaths,col2=col.deaths.adj,show.after.tmax = TRUE)+
  # annotate("label",x=c(day_tswitch,day_tmax),y=c(125,125),label=c("a","b")) +
  coord_cartesian(xlim=c(as.Date("2019-12-30"),as.Date("2020-03-31")),ylim=c(0,125*1.03))
g3
g4 = plot_agedist_cases(S_model10,D_S_model10,col1=col.cases,col2=col.cases.adj)
g4
g5 = plot_agedist_deaths(S_model10,D_S_model10,col1=col.deaths,col2=col.deaths.adj)

data_legend = data.frame(type=c("Reported cases","Symptomatic cases","Overall cases","Reported deaths","Overall deaths"),
                         col=c(col.cases,col.cases.adj,"deepskyblue2",col.deaths,col.deaths.adj),
                         x=1:5,y=1:5,stringsAsFactors = FALSE)
leg = ggplot(data_legend) +
  geom_point(aes(x=x,y=y,fill=factor(type,levels=c("Reported cases","Symptomatic cases","Overall cases","Reported deaths","Overall deaths"))),
             colour="black",shape=22,size=6) +
  scale_fill_manual(values=data_legend$col,guide=guide_legend()) +
  labs(fill=NULL) +  
  theme(legend.position="bottom",
        legend.text=element_text(size=12),
        legend.key.height=unit(12,"pt"),
        legend.margin=margin(4,5,4,4,"pt"),
        legend.box.background = element_rect(fill="white",colour = "grey60",size=.4)
  )
legend <- cowplot::get_legend(leg)
gnull = ggplot() +  geom_blank() + theme_minimal()

plot_grid(
  plot_grid(g1,g4,g3,g5,nrow=2,rel_widths = c(2,1),labels=LETTERS[1:4]),
  plot_grid(gnull,legend,rel_widths = c(.1,1)),
  ncol=1,rel_heights = c(10,1))
ggsave(file="figures/fig_fit.pdf",width=10.5,height=6.4)
ggsave(file="figures/fig_fit.png",width=10.5,height=6.4)


sum(D_S_model10$incidence_cases)
print(S_model10,pars="predicted_total_reported_cases")
print(S_model10,pars="predicted_total_overall_cases")

summary(S_model10,pars="predicted_total_overall_cases")[[1]]/.486

sum(D_S_model10$incidence_deaths)
print(S_model10,pars="predicted_total_overall_deaths_tmax")
print(S_model10,pars="predicted_total_overall_deaths_delay")

# CFR ----------

summary(S_model10,pars=c("cfr_A","cfr_B","cfr_C","cfr_D"))[[1]]
  
summary(S_model10,pars=c("cfr_A_by_age"))[[1]]


sum_cfr = summary(S_model10,pars=c("cfr_A_by_age","cfr_B_by_age","cfr_C_by_age","cfr_D_by_age","cfr_A","cfr_B","cfr_C","cfr_D"))[[1]] %>%
  tbl_df() %>%
  mutate(age_class=c(rep(c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),4),rep("Overall",4)),
         type=c(rep(c("Crude","Adjusted for delayed mortality","Adjusted for unidentified young cases","Adjusted for both"),each=9),
                c("Crude","Adjusted for delayed mortality","Adjusted for unidentified young cases","Adjusted for both"))) %>%
  mutate(type=factor(type,levels=c("Crude","Adjusted for delayed mortality","Adjusted for unidentified young cases","Adjusted for both")),
         age_class=factor(age_class,levels=c("Overall","0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"))) %>%
  mutate(overall=ifelse(age_class=="Overall",1,0)) %>%
  mutate(output_cfr_sympto=paste0(signif(`50%`*100,2),"% (",signif(`2.5%`*100,2),"-",signif(`97.5%`*100,2),")")) %>%
  mutate(output_cfr_asympto=paste0(signif(`50%`*100*.486,2),"% (",signif(`2.5%`*100*.486,2),"-",signif(`97.5%`*100*.486,2),")")) %>%
  mutate(output_cfodds_sympto=paste0(sprintf("%.0f",1/(`50%`))," (",sprintf("%.0f",1/`97.5%`),"-",sprintf("%.0f",1/`2.5%`),")")) %>%
  mutate(output_cfodds_asympto=paste0(sprintf("%.0f",1/(`50%`*.486))," (",sprintf("%.0f",1/(`97.5%`*.486)),"-",sprintf("%.0f",1/(`2.5%`*.486)),")"))




ggplot() +
  # geom_hline(data=filter(sum_cfr,overall==1),aes(yintercept=`50%`)) +
  geom_pointrange(data=filter(sum_cfr,overall==0),aes(x=age_class,y=`50%`,ymin=`2.5%`,ymax=`97.5%`,colour=type)) +
  scale_colour_manual(values=c("grey30",col.deaths.adj,col.cases.adj,"blue"),guide=FALSE) +
  facet_wrap(~type,nrow=2) +
  scale_y_continuous(labels=scales::percent) +
  labs(x="Age class",y="Case fatality ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


ggplot() +
  geom_pointrange(data=filter(sum_cfr,overall==0),
                  aes(x=age_class,y=`50%`,ymin=`2.5%`,ymax=`97.5%`,colour=type),
                  position=position_dodge(.5)) +
  scale_colour_manual(values=c("grey30",col.deaths.adj,col.cases.adj,"blue")) +
  scale_y_continuous(labels=scales::percent) +
  labs(x="Age class",y="Case fatality ratio") +
  theme_bw()

# Table CFR among symptomatics
select(sum_cfr,age_class,type,output_cfr_sympto) %>%
  spread(age_class,output_cfr_sympto) %>%
  write.csv(file="tables/output_cfr_sympto.csv")

# Table CFR among asymptomatics
select(sum_cfr,age_class,type,output_cfr_asympto) %>%
  spread(age_class,output_cfr_asympto) %>%
  write.csv(file="tables/output_cfr_asympto.csv")

# Table CFodds among symptomatics
select(sum_cfr,age_class,type,output_cfodds_sympto) %>%
  spread(age_class,output_cfodds_sympto) %>%
  write.csv(file="tables/output_cfodds_sympto.csv")

# Table CFodds among asymptomatics
select(sum_cfr,age_class,type,output_cfodds_asympto) %>%
  spread(age_class,output_cfodds_asympto) %>%
  write.csv(file="tables/output_cfodds_asympto.csv")

select(sum_cfr,age_class,type,output_cfodds_asympto) %>%
  spread(age_class,output_cfodds_asympto) %>%
  write.csv(file="tables/output_cfodds_asympto.csv")
