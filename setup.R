# Set-up
library(tidyverse)
library(lubridate)
library(rstan)
library(cowplot)
theme_set(theme_bw())
library(readxl)
library(xtable)

bashfile_rdump = function(model_name,id="",data_file,warmup=1000,iter=1000,adapt_delta=0.8,max_depth=12,init=0.5,timelimit=24,chains=4) {
  tt = gsub(" |:","-",Sys.time())
  data_file$inference=0
  with(data_file,stan_rdump(ls(data_file),file=paste0("run_models/data_SIM_",model_name,id,"_",tt,".R")))
  data_file$inference=1
  with(data_file,stan_rdump(ls(data_file),file=paste0("run_models/data_S_",model_name,id,"_",tt,".R")))
  tt_launch = c("#!/bin/bash",
                  paste0("#SBATCH --job-name='",model_name,id,"'"),
                  "#SBATCH --partition=all",
                  paste0("#SBATCH --array=1-",chains),
                  paste0("#SBATCH --cpus-per-task=4"),
                  paste0("#SBATCH --mem-per-cpu=8G"),
                  paste0("#SBATCH --time=",timelimit,":00:00"),
                  paste0("#SBATCH --qos=job_highprio_covid19"),
                  "",
                  "module load Boost/1.66.0-foss-2018a",
                  paste0("./",model_name," sample num_warmup=",warmup," num_samples=",iter," \\"),
                  paste0("      adapt delta=",adapt_delta," algorithm=hmc engine=nuts max_depth=",max_depth," init=",init," \\"),
                  paste0("      data file=data_S_",model_name,id,"_",tt,".R output file=S_",model_name,id,"_",tt,"_${SLURM_ARRAY_TASK_ID}.csv refresh=10")
    )
  writeLines(tt_launch,paste0("run_models/sb_",model_name,id,".sh"))
}

qsum = function(x) c(`50%`=median(x),quantile(x,c(0.025,0.975)))
