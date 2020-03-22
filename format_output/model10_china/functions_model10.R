

bashfile_rdump = function(model_name,data_file,warmup=1000,iter=1000,adapt_delta=0.8,max_depth=12,init=0.5,timelimit=24,chains=4,priorpredcheck=TRUE) {
  data_file$inference=0
  with(data_file,stan_rdump(ls(data_file),file=paste0("model/data_SIM_",model_name,".R")))
  data_file$inference=1
  with(data_file,stan_rdump(ls(data_file),file=paste0("model/data_S_",model_name,".R")))
  if(priorpredcheck) {
    tt_launch = c("#!/bin/bash",
                  paste0("#SBATCH --job-name='",model_name,"'"),
                  "#SBATCH --partition=all",
                  paste0("#SBATCH --array=1-",chains),
                  paste0("#SBATCH --cpus-per-task=4"),
                  paste0("#SBATCH --mem-per-cpu=32G"),
                  paste0("#SBATCH --time=",timelimit,":00:00"),
                  "",
                  "module load Boost/1.66.0-foss-2018a",
                  paste0("./",model_name," sample num_warmup=",warmup," num_samples=",iter," \\"),
                  paste0("      adapt delta=",adapt_delta," algorithm=hmc engine=nuts max_depth=",max_depth," init=",init," \\"),
                  paste0("      data file=data_SIM_",model_name,".R output file=SIM_",model_name,"_",gsub(" |:","-",Sys.time()),"_${SLURM_ARRAY_TASK_ID}.csv refresh=10"),
                  paste0("./",model_name," sample num_warmup=",warmup," num_samples=",iter," \\"),
                  paste0("      adapt delta=",adapt_delta," algorithm=hmc engine=nuts max_depth=",max_depth," init=",init," \\"),
                  paste0("      data file=data_S_",model_name,".R output file=S_",model_name,"_",gsub(" |:","-",Sys.time()),"_${SLURM_ARRAY_TASK_ID}.csv refresh=10")
    )
  } else {
    tt_launch = c("#!/bin/bash",
                  paste0("#SBATCH --job-name='",model_name,"'"),
                  "#SBATCH --partition=all",
                  paste0("#SBATCH --array=1-",chains),
                  paste0("#SBATCH --cpus-per-task=4"),
                  paste0("#SBATCH --mem-per-cpu=8G"),
                  paste0("#SBATCH --time=",timelimit,":00:00"),
                  "",
                  "module load Boost/1.66.0-foss-2018a",
                  paste0("./",model_name," sample num_warmup=",warmup," num_samples=",iter," \\"),
                  paste0("      adapt delta=",adapt_delta," algorithm=hmc engine=nuts max_depth=",max_depth," init=",init," \\"),
                  paste0("      data file=data_S_",model_name,".R output file=S_",model_name,"_",gsub(" |:","-",Sys.time()),"_${SLURM_ARRAY_TASK_ID}.csv refresh=10")
    )
  }
  writeLines(tt_launch,paste0("model/sb_",model_name,".sh"))
}

qsum = function(x) c(`50%`=median(x),quantile(x,c(0.025,0.975)))

# plot functions
plot_incidence_cases = function(samples,data_list,col1="grey",col2="red",show.asympto=FALSE) {
  t0 = data_list$t0
  tmax2 = data_list$tmax2
  D = data_list$D
  tswitch = data_list$tswitch
  S = data_list$S
  y = rstan::extract(samples,"y")[[1]]
  
  data_incidence_cases = data.frame(time=t_data:tmax,incidence=data_list$incidence_cases)
  output_incidence_cases = rstan::summary(samples,"predicted_reported_incidence_cases")[[1]] %>%
    tbl_df() %>%
    mutate(time=1:S,
           date=time+as.Date("2019-12-30")) %>%
    left_join(data_incidence_cases) %>%
    filter(time<=tmax2)
  predicted_overall_incidence_cases = rstan::summary(samples,"predicted_overall_incidence_cases")[[1]] %>%
    tbl_df() %>%
    mutate(time=1:S,
           date=time+as.Date("2019-12-30")) %>%
    left_join(data_incidence_cases) %>%
    filter(time<=tmax2)
  
  ggplot() +
    geom_ribbon(data=predicted_overall_incidence_cases,aes(x=date,ymin=`2.5%`/.48,ymax=`97.5%`/.48),fill="deepskyblue2",alpha=1) +
    geom_line(data=predicted_overall_incidence_cases,aes(x=date,y=`50%`/.48),linetype=3) +
    geom_ribbon(data=predicted_overall_incidence_cases,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill=col2,alpha=1) +
    geom_line(data=predicted_overall_incidence_cases,aes(x=date,y=`50%`),linetype=3) +
    geom_ribbon(data=output_incidence_cases,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill=col1,alpha=1) +
    geom_line(data=output_incidence_cases,aes(x=date,y=`50%`)) +
    geom_point(data=output_incidence_cases,aes(x=date,y=incidence)) +
    coord_cartesian(xlim=c(as.Date("2019-12-30"),tmax2-1+as.Date("2019-12-30"))) +
    labs(x="Time",y="Cases per day") +
    geom_vline(xintercept=tswitch-3+as.Date("2019-12-30"),linetype=2) +
    geom_vline(xintercept=tmax2+as.Date("2019-12-30"),linetype=2) 
}

plot_incidence_deaths = function(samples,data_list,col1="grey",col2="red",show.after.tmax=FALSE) {
  t0 = data_list$t0
  tmax2 = data_list$tmax2
  D = data_list$D
  tswitch = data_list$tswitch
  S = data_list$S
  y = rstan::extract(samples,"y")[[1]]
  data_incidence_cases = data.frame(time=t_data:tmax,incidence=data_list$incidence_deaths)
  output_incidence_deaths = rstan::summary(samples,"predicted_overall_incidence_deaths")[[1]] %>%
    tbl_df() %>%
    mutate(time=1:S,
           date=time+as.Date("2019-12-30")) %>%
    left_join(data_incidence_cases)
  output_incidence_deaths2 = filter(output_incidence_deaths,time>=tmax2) 
  
  g = ggplot() +
    geom_ribbon(data=filter(output_incidence_deaths,time<=tmax2),aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill=col1,alpha=1) +
    geom_point(data=filter(output_incidence_deaths,time<=tmax2),aes(x=date,y=incidence)) +
    geom_ribbon(data=output_incidence_deaths2,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill=col2,alpha=1) +
    geom_line(data=filter(output_incidence_deaths,time<=tmax2),aes(x=date,y=`50%`)) +
    geom_line(data=output_incidence_deaths2,aes(x=date,y=`50%`),linetype=3) +
    labs(x="Time",y="Deaths per day") +
    geom_vline(xintercept=tswitch-3+as.Date("2019-12-30"),linetype=2) +
    geom_vline(xintercept=tmax2+as.Date("2019-12-30"),linetype=2) +
    coord_cartesian(xlim=c(as.Date("2019-12-30"),as.Date("2020-04-05")))
  if(!show.after.tmax) g = g + coord_cartesian(xlim=c(as.Date("2019-12-30"),tmax2-1+as.Date("2019-12-30")))
  return(g)
}

plot_agedist_cases = function(samples,data_list,col1="grey",col2="red") {
  t0 = data_list$t0
  tmax2 = data_list$tmax2
  D = data_list$D
  tswitch = data_list$tswitch
  S = data_list$S
  y = rstan::extract(samples,"y")[[1]]
  
  predicted_total_reported_cases_by_age = rstan::summary(samples,"predicted_total_reported_cases_by_age")[[1]] %>%
    tbl_df() %>%
    mutate(age_class=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+")) %>%
    bind_cols(data=data_list$agedistr_cases,pop=data_list$age_dist*data_list$pop_t/100000) %>%
    mutate(p_mean=mean/pop,
           p_low=`2.5%`/pop,
           p_high=`97.5%`/pop,
           p_data=data/pop,
           p_mean2=p_mean/129,
           p_low2=p_low/129,
           p_high2=p_high/129,
           p_data2=p_data/129)
  predicted_total_overall_cases_by_age = rstan::summary(samples,"predicted_total_overall_cases_by_age")[[1]] %>%
    tbl_df() %>%
    mutate(age_class=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+")) %>%
    bind_cols(data=data_list$agedistr_cases,pop=data_list$age_dist*data_list$pop_t/100000) %>%
    mutate(p_mean=mean/pop,
           p_low=`2.5%`/pop,
           p_high=`97.5%`/pop,
           p_data=data/pop,
           p_mean2=p_mean/129,
           p_low2=p_low/129,
           p_high2=p_high/129,
           p_data2=p_data/129)
  # Raw numbers
  # g = ggplot() +
  #   geom_bar(data=predicted_total_reported_cases_by_age,aes(x=age_class,y=data),stat="identity",fill=NA,colour="black") +
  #   geom_pointrange(data=predicted_total_reported_cases_by_age,aes(x=age_class,y=mean,ymin=`2.5%`,ymax=`97.5%`),colour=col1) +
  #   geom_pointrange(data=predicted_total_overall_cases_by_age,aes(x=age_class,y=mean,ymin=`2.5%`,ymax=`97.5%`),colour=col2) +
  #   labs(x="Age class",y="N")
  # Scaled by population
  # g = ggplot() +
  #   geom_bar(data=predicted_total_reported_cases_by_age,aes(x=age_class,y=p_data),stat="identity",fill=NA,colour="black") +
  #   geom_pointrange(data=predicted_total_reported_cases_by_age,aes(x=age_class,y=p_mean,ymin=p_low,ymax=p_high),colour=col1) +
  #   geom_pointrange(data=predicted_total_overall_cases_by_age,aes(x=age_class,y=p_mean,ymin=p_low,ymax=p_high),colour=col2) +
  #   labs(x="Age class",y="Total cases per 100,000") +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #   scale_y_continuous(expand = c(0,0))
  # Scaled by population and max
  g = ggplot() +
    geom_bar(data=predicted_total_reported_cases_by_age,aes(x=age_class,y=p_data2),stat="identity",fill=NA,colour="black") +
    geom_pointrange(data=predicted_total_reported_cases_by_age,aes(x=age_class,y=p_mean2,ymin=p_low2,ymax=p_high2),colour=col1) +
    geom_pointrange(data=predicted_total_overall_cases_by_age,aes(x=age_class,y=p_mean2,ymin=p_low2,ymax=p_high2),colour=col2) +
    geom_pointrange(data=predicted_total_overall_cases_by_age,aes(x=age_class,y=p_mean2/.48,ymin=p_low2/.48,ymax=p_high2/.48),colour="deepskyblue2") +
    labs(x="Age class",y="Relative total cases") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_hline(yintercept=1,linetype=3) +
    scale_y_continuous(expand = c(0,0))
  return(g)
}

plot_agedist_deaths = function(samples,data_list,col1="grey",col2="red") {
  t0 = data_list$t0
  tmax2 = data_list$tmax2
  D = data_list$D
  tswitch = data_list$tswitch
  S = data_list$S
  y = rstan::extract(samples,"y")[[1]]
  
  predicted_total_overall_deaths_tmax_by_age = rstan::summary(samples,"predicted_total_overall_deaths_tmax_by_age")[[1]] %>%
    tbl_df() %>%
    mutate(age_class=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+")) %>%
    bind_cols(data=data_list$agedistr_deaths,pop=data_list$age_dist*data_list$pop_t/100000) %>%
    mutate(p_mean=mean/pop,
           p_low=`2.5%`/pop,
           p_high=`97.5%`/pop,
           p_data=data/pop,
           p_mean2=p_mean/19.1,
           p_low2=p_low/19.1,
           p_high2=p_high/19.1,
           p_data2=p_data/19.1)
  predicted_total_overall_deaths_delay_by_age = rstan::summary(samples,"predicted_total_overall_deaths_delay_by_age")[[1]] %>%
    tbl_df() %>%
    mutate(age_class=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+")) %>%
    bind_cols(data=data_list$agedistr_deaths,pop=data_list$age_dist*data_list$pop_t/100000) %>%
    mutate(p_mean=mean/pop,
           p_low=`2.5%`/pop,
           p_high=`97.5%`/pop,
           p_data=data/pop,
           p_mean2=p_mean/19.1,
           p_low2=p_low/19.1,
           p_high2=p_high/19.1,
           p_data2=p_data/19.1)
  # Raw numbers
  # g = ggplot() +
  #   geom_bar(data=predicted_total_overall_deaths_tmax_by_age,aes(x=age_class,y=data),stat="identity",fill=NA,colour="black") +
  #   geom_pointrange(data=predicted_total_overall_deaths_tmax_by_age,aes(x=age_class,y=mean,ymin=`2.5%`,ymax=`97.5%`),colour=col1) +
  #   geom_pointrange(data=predicted_total_overall_deaths_delay_by_age,aes(x=age_class,y=mean,ymin=`2.5%`,ymax=`97.5%`),colour=col2) +
  #   labs(x="Age class",y="Total deaths") +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #   scale_y_continuous(expand = c(0,0))
  # Scaled by population
  # g = ggplot() +
  #   geom_bar(data=predicted_total_overall_deaths_tmax_by_age,aes(x=age_class,y=p_data),stat="identity",fill=NA,colour="black") +
  #   geom_pointrange(data=predicted_total_overall_deaths_tmax_by_age,aes(x=age_class,y=p_mean,ymin=p_low,ymax=p_high),colour=col1) +
  #   geom_pointrange(data=predicted_total_overall_deaths_delay_by_age,aes(x=age_class,y=p_mean,ymin=p_low,ymax=p_high),colour=col2) +
  #   labs(x="Age class",y="Proportion") +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #   scale_y_continuous(expand = c(0,0))
  # Scaled by population and max
  g = ggplot() +
    geom_bar(data=predicted_total_overall_deaths_tmax_by_age,aes(x=age_class,y=p_data2),stat="identity",fill=NA,colour="black") +
    geom_pointrange(data=predicted_total_overall_deaths_tmax_by_age,aes(x=age_class,y=p_mean2,ymin=p_low2,ymax=p_high2),colour=col1) +
    geom_pointrange(data=predicted_total_overall_deaths_delay_by_age,aes(x=age_class,y=p_mean2,ymin=p_low2,ymax=p_high2),colour=col2) +
    labs(x="Age class",y="Relative total deaths") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_hline(yintercept=1,linetype=3) +
    scale_y_continuous(expand = c(0,0)) 
  g
  return(g)
}