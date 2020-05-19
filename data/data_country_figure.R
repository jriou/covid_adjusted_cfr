# Setup ----
library(rstan)
library(tidyverse)
library(lubridate)
library(cowplot)
library(ggpubr)
theme_set(theme_bw())
library(xtable)
library(socialmixr)


# Settings ----
country_name="Lombardy"
region_name="Lombardy"

list_surveys()
survey_china <- get_survey("https://doi.org/10.5281/zenodo.3366396")
survey_polymod <- get_survey("https://doi.org/10.5281/zenodo.1043437")

plot_country =function(country_name,region_name){
  country=casefold(country_name)
  region=casefold(region_name)
  source(paste("data/",country,"/data_management_",region,".R",sep=""))
  col.cases1 = "darkcyan"
  col.cases2 = "chartreuse3"
  col.cases3 = "deepskyblue2"
  col.deaths1 = "firebrick"
  col.deaths2 =  "gold"
  col.deaths3 = "darkorange"
  
  d_country=data.frame(date=seq(day_data,day_max,by="day"),confirmed_cases=incidence_cases,confirmed_deaths=incidence_deaths)
  # Figure 1: data ----
  
  h1 = ggplot(d_country) +
    geom_col(aes(x=date,y=confirmed_cases),fill=col.cases1,width=1,col="black") +
    geom_vline(xintercept=day_data-.5,linetype=2) +
    geom_vline(xintercept=day_max+.5,linetype=2) +
    scale_y_continuous(expand=expand_scale(mult=c(0,.05)),
                       labels=function(x) paste0(x/1000,"K")) +
    #scale_x_date(limits=c(as.Date("2020-01-01"),day_max+1)) +
    labs(x="Date",y="N") #+
  # annotate("segment",x=as.Date("2019-12-15"),xend=as.Date("2019-12-08"),y=900,yend=5,arrow=arrow(length=unit(.2,"cm"))) +
  # annotate("label",x=as.Date("2019-12-15"),y=900,label="First case") +
  # annotate("segment",x=as.Date("2020-01-01"),xend=as.Date("2020-01-01"),y=900,yend=5,arrow=arrow(length=unit(.2,"cm"))) +
  # annotate("label",x=as.Date("2020-01-01"),y=900,label="Market closure") +
  # annotate("segment",x=as.Date("2020-01-03"),xend=as.Date("2020-01-20"),y=2150,yend=2150,arrow=arrow(length=unit(.2,"cm"))) +
  # annotate("label",x=as.Date("2020-01-03"),y=2150,label="Control measures\nimplemented") +
  # annotate("segment",x=as.Date("2020-01-03"),xend=as.Date("2020-02-11"),y=2750,yend=2750,arrow=arrow(length=unit(.2,"cm"))) +
  # annotate("label",x=as.Date("2020-01-03"),y=2750,label="End of data collection") 
  
  h2 = ggplot(d_country) +
    geom_col(aes(x=date,y=confirmed_deaths),fill=col.deaths1,col="black",width=1) +
    geom_vline(xintercept=day_data-.5,linetype=2) +
    geom_vline(xintercept=day_max+.5,linetype=2) +
    scale_y_continuous(expand=expand_scale(mult=c(0,.05))) +
    #scale_x_date(limits=c(as.Date("2020-01-01"),day_max+1)) +
    labs(x="Date",y="N") 
  
  if(country_name=="Austria"){
    age_class_country=c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75+")
  }else{
    age_class_country=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+")
  }

  dd2 = data.frame(age_class=age_class_country,
                   pop = age_dist,
                   cases_tmax=agedistr_cases,
                   mort_tmax=agedistr_deaths,
                   prop_cases_tmax=agedistr_cases/sum(agedistr_cases),
                   prop_mort_tmax=agedistr_deaths/sum(agedistr_deaths))
  llab = c(pop="General population",
           prop_cases_tmax="Reported cases",
           prop_mort_tmax="Reported deaths")
  g1 = select(dd2,age_class,pop,prop_cases_tmax,prop_mort_tmax) %>%
    gather("var","value",2:4) %>%
    ggplot() +
    geom_col(aes(x=factor(age_class,levels=age_class_country),y=value,fill=var),colour="black",width=1) +
    facet_wrap(~var,labeller = labeller(var=llab)) +
    geom_hline(yintercept=0) +
    scale_fill_manual(values=c("seagreen",col.cases1,col.deaths1),guide=FALSE) +
    coord_flip() +
    scale_y_continuous(labels=scales::percent,expand = expand_scale(mult=c(0,.11))) +
    labs(x="Age group",y="Proportion") +
    theme(axis.text.x=element_text(angle=45,hjust=1))
  
  if(country_name=="China"){
    survey <- survey_china
  } else{
    survey <- survey_polymod
  }
  if(country_name=="Austria"){
    contact <- contact_matrix(survey = survey, age.limits=c(0,5,15,25,35,45,55,65,75),symmetric = TRUE)
  }else{
    contact <- contact_matrix(survey = survey, age.limits=c(0,10,20,30,40,50,60,70,80),symmetric = TRUE)
  }
 
  contacts = matrix(contact$matrix,nrow=9,byrow=T) %>%
    data.frame() %>%
    tbl_df()
  names(contacts) = age_class_country
  contacts$i = age_class_country
  g2 = gather(contacts,"j","c",1:9) %>%
    ggplot() +
    geom_tile(aes(x=factor(i,levels=age_class_country),y=factor(j,levels=age_class_country),fill=c)) +
    scale_fill_gradient(low="darkblue",high="darkgoldenrod1") +
    labs(x="Age group",y="Age group",fill="Contacts")+
    theme(axis.text.x=element_text(angle=45,hjust=1))
  g_final <- ggarrange(h1,g1,h2,g2,nrow=2,ncol=2,widths=c(1.3,1),labels=LETTERS)
  ggsave(file=paste0("format_output/figures/data_",region,".pdf"),width=10,height=5.5)
  return(g_final)
}

plot_country("Lombardy","Lombardy")
plot_country("Austria","Austria")
plot_country("Germany","Bavaria")
plot_country("Germany","Badenw")
plot_country("Spain","Spain")
plot_country("Switzerland","Switzerland")
plot_country("China","China")
