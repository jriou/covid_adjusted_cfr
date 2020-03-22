# general functions
qsum = function(x) c(`50%`=median(x),quantile(x,c(0.025,0.975)))
inv_logit = function(x) exp(x)/(1+exp(x))
logit = function(x) log(x/(1-x))

# plot functions
plot_incidence_cases = function(samples,data_list,col1="darkcyan",col2="chartreuse3",col3="deepskyblue2",start_date,end_date) {
  t0 = data_list$t0
  tmax2 = data_list$S
  D = data_list$D
  if(!is.null(data_list$D1)) D =  data_list$D1
  G = data_list$G
  tswitch = data_list$tswitch
  S = data_list$S
  y = rstan::extract(samples,"y")[[1]]
  data_incidence_cases = data.frame(time=1:D,incidence=data_list$incidence_cases) %>%
    mutate(date=time+start_date)
  predicted_reported_incidence_symptomatic_cases = rstan::summary(samples,"predicted_reported_incidence_symptomatic_cases")[[1]] %>%
    tbl_df() %>%
    mutate(time=1:S,
           date=time+start_date) %>%
    left_join(data_incidence_cases) %>%
    filter(time<=tmax2)
  predicted_overall_incidence_symptomatic_cases = rstan::summary(samples,"predicted_overall_incidence_symptomatic_cases")[[1]] %>%
    tbl_df() %>%
    mutate(time=1:S,
           date=time+start_date) %>%
    left_join(data_incidence_cases) %>%
    filter(time<=tmax2)
  predicted_overall_incidence_all_cases = rstan::summary(samples,"predicted_overall_incidence_all_cases")[[1]] %>%
    tbl_df() %>%
    mutate(time=1:S,
           date=time+start_date) %>%
    left_join(data_incidence_cases) %>%
    filter(time<=tmax2)
  ggplot() +
    geom_ribbon(data=predicted_overall_incidence_all_cases,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill=col3,alpha=1) +
    geom_line(data=predicted_overall_incidence_all_cases,aes(x=date,y=`50%`),linetype=3) +
    geom_ribbon(data=predicted_overall_incidence_symptomatic_cases,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill=col2,alpha=1) +
    geom_line(data=predicted_overall_incidence_symptomatic_cases,aes(x=date,y=`50%`),linetype=3) +
    geom_ribbon(data=predicted_reported_incidence_symptomatic_cases,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill=col1,alpha=1) +
    geom_line(data=predicted_reported_incidence_symptomatic_cases,aes(x=date,y=`50%`)) +
    geom_point(data=data_incidence_cases,aes(x=date,y=incidence),shape=21,fill="white") +
    coord_cartesian(xlim=c(start_date,end_date)) +
    labs(x="Time",y="Cases per day") +
    scale_y_continuous(expand=expand_scale(mult=c(0,.05)),
                       labels=function(x) paste0(x/1000,"K")) +
    geom_vline(xintercept=tswitch+start_date,linetype=2) +
    geom_vline(xintercept=tmax2+start_date,linetype=2) 
}

plot_total_cases = function(samples,data_list,col1="darkcyan",col2="chartreuse3",col3="deepskyblue2") {
  t0 = data_list$t0
  tmax2 = data_list$S
  D = data_list$D
  G = data_list$G
  tswitch = data_list$tswitch
  S = data_list$S
  y = rstan::extract(samples,"y")[[1]]
  data_incidence_cases = data.frame(time=1:D,incidence=data_list$incidence_cases)
  totals = rstan::summary(samples,c("predicted_total_reported_symptomatic_cases","predicted_total_overall_symptomatic_cases","predicted_total_overall_all_cases"))[[1]] %>%
    tbl_df() %>%
    mutate(type=c("predicted_total_reported_symptomatic_cases","predicted_total_overall_symptomatic_cases","predicted_total_overall_all_cases"),
           type2=c("R","S","A"),
           type3="") %>%
    mutate(type2=factor(type2,levels=c("R","S","A")))
  ggplot(totals) +
    annotate("col",x=1,y=sum(data_incidence_cases$incidence),colour="black",fill="white",width=.4) +
    geom_pointrange(aes(x=type3,ymin=`2.5%`,y=`50%`,ymax=`97.5%`,colour=type),stat="identity",
                 position=position_dodge(-0.3)) +
    scale_colour_manual(values=c(col3,col2,col1),guide=FALSE) +
    scale_y_continuous(expand=expand_scale(mult=c(0,.05)),limits=c(0,max(totals$`97.5%`)),
                       labels=function(x) paste0(x/1000,"K")) +
    scale_x_discrete(expand=expand_scale(add=.4)) +
    labs(x="",y="Total cases") 
}
plot_agedist_cases = function(samples,data_list,col1="darkcyan",col2="chartreuse3",col3="deepskyblue2") {
  t0 = data_list$t0
  tmax2 = data_list$S
  D = data_list$D
  G = data_list$G
  tswitch = data_list$tswitch
  S = data_list$S
  y = rstan::extract(samples,"y")[[1]]
  cases_agedist = data_list$agedistr_cases / sum(data_list$agedistr_cases) * sum(data_list$incidence_cases)
  deaths_agedist = data_list$agedistr_deaths / sum(data_list$agedistr_deaths) * sum(data_list$incidence_deaths)
  agedist_data = data.frame(
    age_group=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),
    genpop=data_list$age_dist,
    cases=cases_agedist,
    deaths=deaths_agedist,
    totpop=data_list$pop_t,
    age_group2=1:9)
  pp = c("predicted_total_reported_symptomatic_cases_by_age",
         "predicted_total_overall_symptomatic_cases_by_age",
         "predicted_total_overall_all_cases_by_age")
  pred = rstan::summary(samples,pp)[[1]] %>%
    tbl_df() %>%
    mutate(age_group=rep(c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),3),
           type=rep(pp,each=9),
           type2=rep(c("R","S","A"),each=9)) %>%
    mutate(type2=factor(type2,levels=c("R","S","A"))) %>%
    left_join(agedist_data) %>%
    group_by(type) %>%
    mutate(med_scale=`50%`/sum(`50%`),low_scale=`2.5%`/sum(`2.5%`),high_scale=`97.5%`/sum(`97.5%`))
  ggplot(pred,aes(x=age_group)) +
    geom_col(data=agedist_data,aes(y=cases),fill="white",colour="black") +
    geom_pointrange(aes(ymin=`2.5%`,y=`50%`,ymax=`97.5%`,colour=type),stat="identity",position=position_dodge(-0.3)) +
    scale_colour_manual(values=c(col3,col2,col1),guide=F) +
    scale_y_continuous(expand=expand_scale(mult=c(0,.05)),
                     labels=function(x) paste0(x/1000,"K")) +
    labs(x="Age group",y="Total cases") +
    theme(axis.text.x=element_text(angle=45,hjust=1))
}

plot_incidence_deaths = function(samples,data_list,col1="darkcyan",col2="chartreuse3",col3="deepskyblue2",start_date,end_date) {
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
    geom_ribbon(data=predicted_overall_incidence_deaths_tmax,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill=col1,alpha=1) +
    geom_line(data=predicted_overall_incidence_deaths_tmax,aes(x=date,y=`50%`)) +
    geom_point(data=predicted_overall_incidence_deaths,aes(x=date,y=incidence),shape=21,fill="white") +
    coord_cartesian(xlim=c(start_date,end_date)) +
    labs(x="Time",y="Deaths per day") +
    scale_y_continuous(expand=expand_scale(mult=c(0,.05))) +
    geom_vline(xintercept=tswitch+start_date,linetype=2) +
    geom_vline(xintercept=tmax2+start_date,linetype=2) 
}
plot_total_deaths = function(samples,data_list,col1="darkcyan",col2="chartreuse3",col3="deepskyblue2",start_date,end_date) {
  t0 = data_list$t0
  tmax2 = data_list$S
  D = data_list$D
  G = data_list$G
  tswitch = data_list$tswitch
  S = data_list$S
  y = rstan::extract(samples,"y")[[1]]
  data_incidence_deaths = data.frame(time=1:D,incidence=data_list$incidence_deaths)
  totals = rstan::summary(samples,c("predicted_total_overall_deaths_tmax","predicted_total_overall_deaths_delay"))[[1]] %>%
    tbl_df() %>%
    mutate(type=c("predicted_total_reported_symptomatic_cases","predicted_total_overall_symptomatic_cases"),
           type2=c("R","A"),
           type3="") %>%
    mutate(type2=factor(type2,levels=c("R","A")))
  ggplot(totals) +
    annotate("col",x=1,y=sum(data_incidence_deaths$incidence),colour="black",fill="white",width=.4) +
    geom_pointrange(aes(x=type3,ymin=`2.5%`,y=`50%`,ymax=`97.5%`,colour=type),stat="identity",
                 position=position_dodge(-0.3)) +
    scale_colour_manual(values=c(col2,col1),guide=FALSE) +
    scale_y_continuous(expand=expand_scale(mult=c(0,.05))) +
    scale_x_discrete(expand=expand_scale(add=.5)) +
    labs(x="",y="Total deaths") 
}


plot_agedist_deaths = function(samples,data_list,col1="darkcyan",col2="chartreuse3") {
  t0 = data_list$t0
  tmax2 = data_list$S
  D = data_list$D
  G = data_list$G
  tswitch = data_list$tswitch
  S = data_list$S
  y = rstan::extract(samples,"y")[[1]]
  cases_agedist = data_list$agedistr_cases / sum(data_list$agedistr_cases) * sum(data_list$incidence_cases)
  deaths_agedist = data_list$agedistr_deaths / sum(data_list$agedistr_deaths) * sum(data_list$incidence_deaths)
  agedist_data = data.frame(
    age_group=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),
    genpop=data_list$age_dist,
    cases=cases_agedist,
    deaths=deaths_agedist,
    totpop=data_list$pop_t,
    age_group2=1:9)
  pp = c("predicted_total_overall_deaths_tmax_by_age",
         "predicted_total_overall_deaths_delay_by_age")
  pred = rstan::summary(samples,pp)[[1]] %>%
    tbl_df() %>%
    mutate(age_group=rep(c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),2),
           type=rep(pp,each=9),
           type2=rep(c("R","A"),each=9)) %>%
    mutate(type2=factor(type2,levels=c("R","A"))) %>%
    left_join(agedist_data) %>%
    group_by(type) %>%
    mutate(med_scale=`50%`/sum(`50%`),low_scale=`2.5%`/sum(`2.5%`),high_scale=`97.5%`/sum(`97.5%`))
  ggplot(pred,aes(x=age_group)) +
    geom_col(data=agedist_data,aes(y=deaths),fill="white",colour="black") +
    geom_pointrange(aes(ymin=`2.5%`,y=`50%`,ymax=`97.5%`,colour=type),stat="identity",position=position_dodge(-0.3)) +
    scale_colour_manual(values=c(col2,col1),guide=F) +
    scale_fill_manual(values=c(col2,col1),guide=F) +
    scale_y_continuous(expand=expand_scale(mult=c(0,.05))) +
    labs(x="Age group",y="Total deaths") +
    theme(axis.text.x=element_text(angle=45,hjust=1))
}


plot_cfr = function(samples,data_list,col1,col2,col3) {
  t0 = data_list$t0
  tmax2 = data_list$S
  D = data_list$D
  G = data_list$G
  tswitch = data_list$tswitch
  S = data_list$S
  y = rstan::extract(samples,"y")[[1]]
  pp = c("cfr_A_symptomatic",
                "cfr_D_symptomatic",
                "cfr_D_all"
  )
  pp2 = c("Crude","Among symptomatics","Among all infected")
  pred = rstan::summary(samples,pp)[[1]] %>%
    tbl_df() %>%
    mutate(type=pp,
           type2=factor(pp2,levels=pp2)) 
  ggplot(pred) +
    geom_pointrange(aes(x=type2,ymin=`2.5%`,y=`50%`,ymax=`97.5%`,colour=type2),
                    stat="identity") +
    scale_colour_manual(values=c(col1,col2,col3),guide=FALSE) +
    scale_y_continuous(expand=expand_scale(mult=c(0,.05)),labels = scales::percent) +
    coord_cartesian(ylim=c(0,max(pred$`97.5%`))) +
    labs(x=NULL,y="Mortality",colour=NULL) +
    theme(axis.text.x=element_text(angle=45,hjust=1))
}


plot_agedist_cfr= function(samples,data_list,col1,col2,col3,insert=c(.7,7,.35,.95)) {
  t0 = data_list$t0
  tmax2 = data_list$S
  D = data_list$D
  G = data_list$G
  tswitch = data_list$tswitch
  S = data_list$S
  y = rstan::extract(samples,"y")[[1]]
  xshift = .15
  A = rstan::summary(samples,"cfr_A_symptomatic_by_age")[[1]] %>%
    tbl_df() %>%
    mutate(age_group=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),
           age_group2=1:9,
           type2="Crude",
           x2=age_group2-xshift) 
  B = rstan::summary(samples,"epsilon")[[1]] %>%
    tbl_df() %>%
    mutate(age_group=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),
           age_group2=1:9,
           type2="Among symptomatics",
           x2=age_group2) 
  epsilon2 = rstan::extract(samples,"epsilon")[[1]]
  psi = rstan::extract(samples,"psi")[[1]]
  C = t(apply(epsilon2,2,function(x) qsum(x*psi))) %>%
    as.data.frame(.) %>%
    tbl_df() %>%
    mutate(age_group=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),
           age_group2=1:9,
           type2="Among all infected",
           x2=age_group2+xshift) 
  pred_by_age = bind_rows(A,B,C) %>%
    mutate(type2=factor(type2,levels=c("Crude","Among symptomatics","Among all infected")))
  wind = ggplot(filter(pred_by_age,age_group2<=6)) +
    geom_pointrange(aes(x=x2,ymin=`2.5%`,y=`50%`,ymax=`97.5%`,colour=type2),
                    stat="identity") +
    # geom_line(aes(x=x2,y=`50%`,colour=type2)) +
    scale_colour_manual(values=c(col1,col2,col3),guide=FALSE) +
    scale_x_continuous(breaks=1:9,labels=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+")) +
    scale_y_continuous(expand=expand_scale(mult=c(0,.05)),labels = scales::percent) +
    labs(x=NULL,y=NULL) +
    theme(axis.text.x=element_text(angle=45,hjust=1),
          panel.border=element_rect(linetype=2,colour="grey30"),
          panel.grid.minor=element_blank())
  g1 = ggplot(pred_by_age) +
    geom_pointrange(aes(x=x2,ymin=`2.5%`,y=`50%`,ymax=`97.5%`,colour=type2),
                    stat="identity") +
    scale_colour_manual(values=c(col1,col2,col3)) +
    scale_x_continuous(breaks=1:9,labels=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+")) +
    scale_y_continuous(expand=expand_scale(mult=c(0,.05)),labels = scales::percent) +
    labs(x=NULL,y="Mortality",colour=NULL) +
    theme(axis.text.x=element_text(angle=45,hjust=1),
          legend.position="bottom")
  pred_by_age %>%
    mutate(res=paste0(round(`50%`,3)," [",round(`2.5%`,3),"-",round(`97.5%`,3),"]")) %>%
    select(type2,age_group,res) %>%
    spread(type2,res)  %>%
    xtable(.) %>%
    print(.,include.rownames=FALSE)
  g1 + 
    annotation_custom(
      ggplotGrob(wind), 
      xmin = insert[1], xmax = insert[2], ymin = insert[3], ymax = insert[4]
    )  +
    annotate("rect",xmin=.5,xmax=6.5,ymin=0,ymax=.037,fill=NA,colour="grey30",linetype=2) + 
    annotate("segment",x=3.5,xend=mean(insert[1:2]),y=.037,yend=insert[3]*1.2,colour="grey30",linetype=2) 
}


plot_incidence_deaths2 = function(samples,data_list,col1="darkcyan",col2="chartreuse3",col3="deepskyblue2",start_date,end_date) {
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
    geom_ribbon(data=predicted_overall_incidence_deaths_tmax,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill=col1,alpha=1) +
    geom_line(data=predicted_overall_incidence_deaths_tmax,aes(x=date,y=`50%`)) +
    geom_point(data=predicted_overall_incidence_deaths,aes(x=date,y=incidence),shape=21,fill="white") +
    coord_cartesian(xlim=c(start_date,end_date)) +
    labs(x="Time",y="Deaths per day") +
    scale_y_continuous(expand=expand_scale(mult=c(0,.05))) +
    geom_vline(xintercept=tswitch+start_date,linetype=2) +
    geom_vline(xintercept=tmax2+start_date,linetype=2) 
}