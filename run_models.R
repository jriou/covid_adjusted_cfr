source("data_management.R")
source("setup.R")

# Model 1 ----------------------------------------------------------

# Format data and set priors
data_list_model1 = list(
  day_start=day_start,
  day_data=day_data,
  day_tmax=day_tmax,
  dat_tswitch=day_tswitch,
  time_lag=time_lag,
  K=9,
  age_dist=age_dist,
  pop_t = pop_t,
  t0=0,
  t_data = t_data,
  tmax =tmax,
  tmax2 = tmax,
  tswitch = tswitch,
  S=S,
  ts = 1:S,
  D = D,
  incidence_cases = incidence_cases,
  incidence_deaths = incidence_deaths,
  agedistr_cases = cases_tmax,
  agedistr_deaths = mort_tmax,
  p_beta=10,
  p_eta=c(1,1),
  p_pi=c(1,999),
  p_incubation=c(5.3,3.2), # Lindon 2020
  p_onset_to_death =c(15.0,6.9), # Lindon 2020
  p_epsilon=c(1,1),
  p_phi = 1/10,
  p_rho=c(1,1),
  inference=1,
  doprint=0
)

# Test
M_model1 = stan_model("model/model1.stan")
T_model1 = sampling(M_model1,
                    data = data_list_model1,
                    iter = 100,
                    chains = 1,
                    init=0,
                    control=list(max_treedepth=10,adapt_delta=0.8))
print(T_model1,pars=c("beta","eta","tau","mu","epsilon","rho_K","pi"))
summary(T_model1,"output_incidence_cases")[[1]] %>%
  tbl_df() %>%
  mutate(t=data_list_model1$t_data:data_list_model1$tmax) %>%
  bind_cols(inc=data_list_model1$incidence_cases) %>%
  ggplot() +
  geom_point(aes(x=t,y=inc)) +
  geom_line(aes(x=t,y=`50%`))

# Create data and bash files
bashfile_rdump("model1",data_list_model1,warmup=500,iter=500,adapt_delta=0.8,max_depth=10,init=0,timelimit=24,chains=4,priorpredcheck=FALSE)

# Copy on cluster
system("scp /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/model1.stan /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/sb_model1.sh /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/data_SIM_model1.R /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/data_S_model1.R UBELIX:projects/COVID_age/model/.")

# Run on cluster
# >> R
# >> library(rstan)
SIM_model1 = read_stan_csv(dir(".",pattern = 'SIM_model1_2020-02-21-11-02-02_[[:digit:]]+.csv'))
D_SIM_model1 = read_rdump("data_SIM_model1.R")
check_hmc_diagnostics(SIM_model1)
print(SIM_model1,pars=c("beta","eta","tau","mu","epsilon","rho_K","pi"),digits_summary=4)
S_model1 = read_stan_csv(dir(".",pattern = 'S_model1_2020-02-21-11-02-02_[[:digit:]]+.csv'))
D_S_model1 = read_rdump("data_S_model1.R")
check_hmc_diagnostics(S_model1)
print(S_model1,pars=c("beta","eta","tau","mu","epsilon","rho_K","pi"),digits_summary=4)
save(SIM_model1,D_SIM_model1,
     S_model1,D_S_model1,
     file="model1_2020-02-21.Rdata")

# Retrieve from cluster
system("scp  UBELIX:projects/COVID_age/model/model1_2020-02-21.Rdata /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/.")

# Load
load("model/model1_2020-02-21.Rdata")

# Prior predictive check
check_hmc_diagnostics(SIM_model1)
print(SIM_model1,pars=c("beta","eta","tau","mu","epsilon","rho_K","pi","phi"),digits_summary=4)
plot_incidence_cases(SIM_model1,D_SIM_model1)
plot_incidence_deaths(SIM_model1,D_SIM_model1)
pairs(SIM_model1,pars=c("beta","eta","tau","mu","output_incidence_cases[50]"))


for(i in 101:200) {
g = extract(SIM_model1,"y")[[1]][i,,37:45] %>%
  apply(.,1,function(x) sum(x) *D_SIM_model1$pop_t ) %>%
  data.frame() %>%
  tbl_df() %>%
  mutate(time=1:S,
         incidence=.-lag(.,default = 0)) %>%
  filter(time<D_SIM_model1$tmax) %>%
  mutate(data_incidence=D_SIM_model1$incidence_cases) %>%
  ggplot() +
  geom_point(aes(x=time,y=data_incidence)) +
  geom_line(aes(x=time,y=incidence)) +
  geom_vline(xintercept=D_SIM_model1$tswitch,linetype=2)
  print(g)
}



# Posterior predictive check
check_hmc_diagnostics(S_model1)
print(S_model1,pars=c("beta","eta","tau","mu","epsilon","rho_K","pi","phi"),digits_summary=4)
print(S_model1,pars=c("incubation","onset_to_death"),digits_summary=4)
print(SIM_model1,pars=c("incubation","onset_to_death"),digits_summary=4)
plot_incidence_cases(S_model1,D_S_model1) + geom_vline(aes(xintercept=as.Date("2020-01-23")))
plot_incidence_deaths(S_model1,D_S_model1) + geom_vline(aes(xintercept=as.Date("2020-01-23")))

print(S_model1,pars=c("cfr_A","cfr_B","cfr_C","cfr_D","cfr_overall"),digits_summary = 4)


summary(S_model1,pars=c("predicted_reported_incidence_cases","predicted_overall_incidence_cases","predicted_overall_incidence_deaths"))[[1]] %>%
  data.frame() %>%
  rownames_to_column() %>%
  tbl_df() %>%
  mutate(time=rep(2:D_S_model1$S,each=3*D_S_model1$K),
         date=time+as.Date("2019-12-01"),
         age_class=rep(c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),3*(D_S_model1$S-1)),
         type=rep(c("predicted_reported_incidence_cases","predicted_overall_incidence_cases","predicted_overall_incidence_deaths"),each=(D_S_model1$S-1)*D_S_model1$K)) %>% 
  ggplot(aes(x=time)) +
  geom_line(aes(y=mean)) +
  facet_grid(age_class ~ type)



extract(S_model1,pars=c("cfr_D"))[[1]]








# Model 2 ----------------------------------------------------------

# Format data and set priors
data_list_model2 = list(
  time_lag=100,
  K=9,
  age_dist=age_dist,
  pop_t = pop_t,
  t0=0,
  t_data = t_data,
  tmax =tmax,
  tmax2 = tmax,
  tswitch = tswitch,
  S=S,
  ts = 1:S,
  D = D,
  incidence_cases = incidence_cases,
  incidence_deaths = incidence_deaths,
  agedistr_cases = cases_tmax,
  agedistr_deaths = mort_tmax,
  p_beta=1,
  p_eta=c(1,1),
  p_pi=c(1,999),
  p_incubation=5.2,
  p_infectious=2.3,
  p_recovery_or_death=20.2,
  p_asymptomatic_age80plus=1,
  p_epsilon=c(1,1),
  p_phi = 1/100,
  p_rho=c(1,1),
  inference=0,
  doprint=0
)

# Test
M_model2 = stan_model("model/model2.stan")
T_model2 = sampling(M_model2,
                    data = data_list_model2,
                    iter = 10,
                    chains = 1,
                    init=0,
                    control=list(max_treedepth=10,adapt_delta=0.8))
print(T_model2,pars=c("beta","eta","epsilon","rho_K","pi"))

for(i in 1:5) {
  g = extract(T_model2,"comp_diffC")[[1]][i,,] %>%
    apply(.,1,sum) %>%
    data.frame() %>%
    tbl_df() %>%
    mutate(time=1:S) %>%
    filter(time<data_list_model2$tmax) %>%
    mutate(data_incidence=data_list_model2$incidence_cases) %>%
    ggplot() +
    geom_point(aes(x=time,y=data_incidence)) +
    geom_line(aes(x=time,y=.)) +
    geom_vline(xintercept=data_list_model2$tswitch,linetype=2)
  print(g)
}

# Create data and bash files
bashfile_rdump("model2",data_list_model2,warmup=500,iter=500,adapt_delta=0.8,max_depth=10,init=0.5,timelimit=96,chains=4,priorpredcheck=FALSE)

# Copy on cluster
system("scp /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/model2.stan /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/sb_model2.sh /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/data_SIM_model2.R /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/data_S_model2.R UBELIX:projects/COVID_age/model/.")

# Run on cluster
# >> R
# >> library(rstan)
SIM_model2 = read_stan_csv(dir(".",pattern = 'SIM_model2_2020-02-24-16-54-12_[[:digit:]]+.csv'))
D_SIM_model2 = read_rdump("data_SIM_model2.R")
check_hmc_diagnostics(SIM_model2)
print(SIM_model2,pars=c("beta","eta","epsilon","rho_K","pi"),digits_summary=4)
S_model2 = read_stan_csv(dir(".",pattern = 'S_model2_2020-02-24-16-54-12_[[:digit:]]+.csv'))
D_S_model2 = read_rdump("data_S_model2.R")
check_hmc_diagnostics(S_model2)
print(S_model2,pars=c("beta","eta","epsilon","rho_K","pi"),digits_summary=4)
save(SIM_model2,D_SIM_model2,
     S_model2,D_S_model2,
     file="model2_2020-02-24.Rdata")

# Retrieve from cluster
system("scp  UBELIX:projects/COVID_age/model/model2_2020-02-24.Rdata /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/.")

# Load
load("model/model2_2020-02-24.Rdata")

# Prior predictive check
check_hmc_diagnostics(SIM_model2)
print(SIM_model2,pars=c("beta","eta","epsilon","rho_K","pi","phi"),digits_summary=4)
plot_incidence_cases(SIM_model2,D_SIM_model2)
plot_incidence_deaths(SIM_model2,D_SIM_model2)
pairs(SIM_model2,pars=c("beta","eta","output_incidence_cases[50]"))

for(i in 1:10) {
  g = extract(SIM_model2,"y")[[1]][i,,37:45] %>%
    apply(.,1,function(x) sum(x) *D_SIM_model2$pop_t ) %>%
    data.frame() %>%
    tbl_df() %>%
    mutate(time=1:S,
           incidence=.-lag(.,default = 0)) %>%
    filter(time<D_SIM_model2$tmax) %>%
    mutate(data_incidence=D_SIM_model2$incidence_cases) %>%
    ggplot() +
    geom_point(aes(x=time,y=data_incidence)) +
    geom_line(aes(x=time,y=incidence)) +
    geom_vline(xintercept=D_SIM_model2$tswitch,linetype=2)
  print(g)
}



# Posterior predictive check
check_hmc_diagnostics(S_model2)
print(S_model2,pars=c("beta","eta","epsilon","rho_K","pi","phi"),digits_summary=4)
plot_incidence_cases(S_model2,D_S_model2)
plot_incidence_deaths(S_model2,D_S_model2)

print(S_model2,pars=c("cfr_A_by_age","cfr_B_by_age","cfr_C_by_age","cfr_D_by_age"),digits_summary = 4)
print(S_model2,pars=c("cfr_A","cfr_B","cfr_C","cfr_D"),digits_summary = 4)

print(S_model2,pars=c("output_agedistr_cases"),digits_summary = 4)

summary(S_model2,pars=c("predicted_reported_incidence_cases","predicted_overall_incidence_cases","predicted_overall_incidence_deaths"))[[1]] %>%
  data.frame() %>%
  rownames_to_column() %>%
  tbl_df() %>%
  mutate(time=rep(2:D_S_model2$S,each=3*D_S_model2$K),
         date=time+as.Date("2019-12-01"),
         age_class=rep(c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),3*(D_S_model2$S-1)),
         type=rep(c("predicted_reported_incidence_cases","predicted_overall_incidence_cases","predicted_overall_incidence_deaths"),each=(D_S_model2$S-1)*D_S_model2$K)) %>% 
  ggplot(aes(x=time)) +
  geom_line(aes(y=mean)) +
  facet_grid(age_class ~ type)




# Model 3 ----------------------------------------------------------

# Format data and set priors
data_list_model3 = list(
  time_lag=100,
  K=9,
  age_dist=age_dist,
  pop_t = pop_t,
  t0=0,
  t_data = t_data,
  tmax =tmax,
  tmax2 = tmax,
  tswitch = tswitch,
  S=S,
  ts = 1:S,
  D = D,
  incidence_cases = incidence_cases,
  incidence_deaths = incidence_deaths,
  agedistr_cases = cases_tmax,
  agedistr_deaths = mort_tmax,
  p_beta=1,
  p_eta=c(1,1),
  p_pi=c(1,999),
  p_incubation=5.2,
  p_infectious=2.3,
  p_recovery_or_death=12.7,
  p_asymptomatic_age80plus=0.5,
  p_epsilon=c(1,1),
  p_phi = 1/100,
  p_rho=c(1,1),
  inference=0,
  doprint=0
)

# Test
M_model3 = stan_model("model/model3.stan")
T_model3 = sampling(M_model3,
                    data = data_list_model3,
                    iter = 100,
                    chains = 1,
                    init=0,
                    control=list(max_treedepth=10,adapt_delta=0.8))
print(T_model3,pars=c("beta","eta","tau","mu","epsilon","rho_K","pi"))
summary(T_model3,"output_incidence_cases")[[1]] %>%
  tbl_df() %>%
  mutate(t=data_list_model3$t_data:data_list_model3$tmax) %>%
  bind_cols(inc=data_list_model3$incidence_cases) %>%
  ggplot() +
  geom_point(aes(x=t,y=inc)) +
  geom_line(aes(x=t,y=`50%`))

# Create data and bash files
bashfile_rdump("model3",data_list_model3,warmup=500,iter=500,adapt_delta=0.8,max_depth=10,init=0.5,timelimit=96,chains=4,priorpredcheck=FALSE)

# Copy on cluster
system("scp /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/model3.stan /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/sb_model3.sh /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/data_SIM_model3.R /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/data_S_model3.R UBELIX:projects/COVID_age/model/.")

# Run on cluster
# >> R
# >> library(rstan)
SIM_model3 = read_stan_csv(dir(".",pattern = 'SIM_model3_2020-02-24-16-57-58_[[:digit:]]+.csv'))
D_SIM_model3 = read_rdump("data_SIM_model3.R")
check_hmc_diagnostics(SIM_model3)
print(SIM_model3,pars=c("beta","eta","epsilon","rho_K","pi"),digits_summary=4)
S_model3 = read_stan_csv(dir(".",pattern = 'S_model3_2020-02-24-16-57-58_[[:digit:]]+.csv'))
D_S_model3 = read_rdump("data_S_model3.R")
check_hmc_diagnostics(S_model3)
print(S_model3,pars=c("beta","eta","epsilon","rho_K","pi"),digits_summary=4)
save(SIM_model3,D_SIM_model3,
     S_model3,D_S_model3,
     file="model3_2020-02-24.Rdata")

# Retrieve from cluster
system("scp  UBELIX:projects/COVID_age/model/model3_2020-02-24.Rdata /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/.")

# Load
load("model/model3_2020-02-24.Rdata")

# Prior predictive check
check_hmc_diagnostics(SIM_model3)
print(SIM_model3,pars=c("beta","eta","epsilon","rho_K","pi","phi"),digits_summary=4)
plot_incidence_cases(SIM_model3,D_SIM_model3)
plot_incidence_deaths(SIM_model3,D_SIM_model3)
pairs(SIM_model3,pars=c("beta","eta","output_incidence_cases[50]"))


for(i in 101:200) {
  g = extract(SIM_model3,"y")[[1]][i,,37:45] %>%
    apply(.,1,function(x) sum(x) *D_SIM_model3$pop_t ) %>%
    data.frame() %>%
    tbl_df() %>%
    mutate(time=1:S,
           incidence=.-lag(.,default = 0)) %>%
    filter(time<D_SIM_model3$tmax) %>%
    mutate(data_incidence=D_SIM_model3$incidence_cases) %>%
    ggplot() +
    geom_point(aes(x=time,y=data_incidence)) +
    geom_line(aes(x=time,y=incidence)) +
    geom_vline(xintercept=D_SIM_model3$tswitch,linetype=2)
  print(g)
}



# Posterior predictive check
check_hmc_diagnostics(S_model3)
print(S_model3,pars=c("beta","eta","epsilon","rho_K","pi","phi"),digits_summary=4)
plot_incidence_cases(S_model3,D_S_model3)
plot_incidence_deaths(S_model3,D_S_model3)

print(S_model3,pars=c("cfr_A","cfr_B","cfr_C","cfr_D"),digits_summary = 4)
print(S_model3,pars=c("cfr_A_by_age","cfr_B_by_age","cfr_C_by_age","cfr_D_by_age"),digits_summary = 4)



summary(S_model3,pars=c("predicted_reported_incidence_cases","predicted_overall_incidence_cases","predicted_overall_incidence_deaths"))[[1]] %>%
  data.frame() %>%
  rownames_to_column() %>%
  tbl_df() %>%
  mutate(time=rep(2:D_S_model3$S,each=3*D_S_model3$K),
         date=time+as.Date("2019-12-01"),
         age_class=rep(c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),3*(D_S_model3$S-1)),
         type=rep(c("predicted_reported_incidence_cases","predicted_overall_incidence_cases","predicted_overall_incidence_deaths"),each=(D_S_model3$S-1)*D_S_model3$K)) %>% 
  ggplot(aes(x=time)) +
  geom_line(aes(y=mean)) +
  facet_grid(age_class ~ type)







# Model 4 ----------------------------------------------------------

# Format data and set priors
data_list_model4 = list(
  time_lag=100,
  K=9,
  age_dist=age_dist,
  pop_t = pop_t,
  t0=0,
  t_data = t_data,
  tmax =tmax,
  tmax2 = tmax,
  tswitch = tswitch,
  S=S,
  ts = 1:S,
  D = D,
  incidence_cases = incidence_cases,
  incidence_deaths = incidence_deaths,
  agedistr_cases = cases_tmax,
  agedistr_deaths = mort_tmax,
  p_beta=1,
  p_eta=c(1,1),
  p_delta=20.2,
  p_pi=c(1,999),
  p_incubation=5.2,
  p_infectious=2.3,
  p_asymptomatic_age80plus=1,
  p_epsilon=c(1,1),
  p_phi = 1/100,
  p_rho=c(1,1),
  inference=0,
  doprint=0
)

# Test
M_model4 = stan_model("model/model4.stan")
T_model4 = sampling(M_model4,
                    data = data_list_model4,
                    iter = 10,
                    chains = 1,
                    init=0,
                    control=list(max_treedepth=10,adapt_delta=0.8))
print(T_model4,pars=c("beta","eta","delta","epsilon","rho_K","pi"))

# Create data and bash files
bashfile_rdump("model4",data_list_model4,warmup=500,iter=500,adapt_delta=0.8,max_depth=10,init=0.5,timelimit=96,chains=4,priorpredcheck=FALSE)

# Copy on cluster
system("scp /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/model4.stan /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/sb_model4.sh /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/data_SIM_model4.R /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/data_S_model4.R UBELIX:projects/COVID_age/model/.")

# Run on cluster
# >> R
# >> library(rstan)
# SIM_model4 = read_stan_csv(dir(".",pattern = 'SIM_model4_2020-02-24-16-54-12_[[:digit:]]+.csv'))
# D_SIM_model4 = read_rdump("data_SIM_model4.R")
# check_hmc_diagnostics(SIM_model4)
# print(SIM_model4,pars=c("beta","eta","epsilon","rho_K","pi"),digits_summary=4)
S_model4 = read_stan_csv(dir(".",pattern = 'S_model4_2020-02-26-12-18-45_[[:digit:]]+.csv'))
D_S_model4 = read_rdump("data_S_model4.R")
check_hmc_diagnostics(S_model4)
print(S_model4,pars=c("beta","eta","epsilon","rho_K","pi"),digits_summary=4)
save(
  # SIM_model4,D_SIM_model4,
     S_model4,D_S_model4,
     file="model4_2020-02-26.Rdata")

# Retrieve from cluster
system("scp  UBELIX:projects/COVID_age/model/model4_2020-02-26.Rdata /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/.")

# Load
load("model/model4_2020-02-26.Rdata")

# Prior predictive check
check_hmc_diagnostics(SIM_model4)
print(SIM_model4,pars=c("beta","eta","epsilon","rho_K","pi","phi"),digits_summary=4)
plot_incidence_cases(SIM_model4,D_SIM_model4)
plot_incidence_deaths(SIM_model4,D_SIM_model4)
pairs(SIM_model4,pars=c("beta","eta","output_incidence_cases[50]"))



# Posterior predictive check
check_hmc_diagnostics(S_model4)
print(S_model4,pars=c("beta","eta","epsilon","delta","rho_K","pi","phi"),digits_summary=4)
plot_incidence_cases(S_model4,D_S_model4)
plot_incidence_deaths(S_model4,D_S_model4)
plot_incidence_deaths(S_model4,D_S_model4,show.after.tmax = TRUE)

print(S_model4,pars=c("cfr_A_by_age","cfr_B_by_age","cfr_C_by_age","cfr_D_by_age"),digits_summary = 4)
print(S_model4,pars=c("cfr_A","cfr_B","cfr_C","cfr_D"),digits_summary = 4)

print(S_model4,pars=c("output_agedistr_cases"),digits_summary = 4)

summary(S_model4,pars=c("predicted_reported_incidence_cases","predicted_overall_incidence_cases","predicted_overall_incidence_deaths"))[[1]] %>%
  data.frame() %>%
  rownames_to_column() %>%
  tbl_df() %>%
  mutate(time=rep(2:D_S_model4$S,each=3*D_S_model4$K),
         date=time+as.Date("2019-12-01"),
         age_class=rep(c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),3*(D_S_model4$S-1)),
         type=rep(c("predicted_reported_incidence_cases","predicted_overall_incidence_cases","predicted_overall_incidence_deaths"),each=(D_S_model4$S-1)*D_S_model4$K)) %>% 
  ggplot(aes(x=time)) +
  geom_line(aes(y=mean)) +
  facet_grid(age_class ~ type)




# Model 5 ----------------------------------------------------------

# Format data and set priors
data_list_model5 = list(
  time_lag=100,
  K=9,
  age_dist=age_dist,
  pop_t = pop_t,
  t0=0,
  t_data = t_data,
  tmax =tmax,
  tmax2 = tmax,
  tswitch = tswitch,
  S=S,
  ts = 1:S,
  D = D,
  incidence_cases = incidence_cases,
  incidence_deaths = incidence_deaths,
  agedistr_cases = cases_tmax,
  agedistr_deaths = mort_tmax,
  p_beta=1,
  p_eta=c(1,1),
  p_delta=12.7,
  p_pi=c(1,999),
  p_incubation=5.2,
  p_infectious=2.3,
  p_asymptomatic_age80plus=1,
  p_epsilon=c(1,1),
  p_phi = 1/100,
  p_rho=c(1,1),
  inference=0,
  doprint=0
)

# Test
M_model5 = stan_model("model/model5.stan")
T_model5 = sampling(M_model5,
                    data = data_list_model5,
                    iter = 10,
                    chains = 1,
                    init=0,
                    control=list(max_treedepth=10,adapt_delta=0.8))
print(T_model5,pars=c("beta","eta","delta","epsilon","rho_K","pi"))


# Create data and bash files
bashfile_rdump("model5",data_list_model5,warmup=500,iter=500,adapt_delta=0.8,max_depth=10,init=0.5,timelimit=96,chains=4,priorpredcheck=FALSE)

# Copy on cluster
system("scp /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/model5.stan /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/sb_model5.sh /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/data_SIM_model5.R /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/data_S_model5.R UBELIX:projects/COVID_age/model/.")

# Run on cluster
# >> R
# >> library(rstan)
# SIM_model5 = read_stan_csv(dir(".",pattern = 'SIM_model5_2020-02-24-16-54-12_[[:digit:]]+.csv'))
# D_SIM_model5 = read_rdump("data_SIM_model5.R")
# check_hmc_diagnostics(SIM_model5)
# print(SIM_model5,pars=c("beta","eta","epsilon","rho_K","pi"),digits_summary=4)
S_model5 = read_stan_csv(dir(".",pattern = 'S_model5_2020-02-26-12-51-49_[[:digit:]]+.csv'))
D_S_model5 = read_rdump("data_S_model5.R")
check_hmc_diagnostics(S_model5)
print(S_model5,pars=c("beta","eta","epsilon","rho_K","pi"),digits_summary=4)
save(
  # SIM_model5,D_SIM_model5,
     S_model5,D_S_model5,
     file="model5_2020-02-26.Rdata")

# Retrieve from cluster
system("scp  UBELIX:projects/COVID_age/model/model5_2020-02-26.Rdata /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/.")

# Load
load("model/model5_2020-02-26.Rdata")

# Prior predictive check
check_hmc_diagnostics(SIM_model5)
print(SIM_model5,pars=c("beta","eta","epsilon","rho_K","pi","phi"),digits_summary=4)
plot_incidence_cases(SIM_model5,D_SIM_model5)
plot_incidence_deaths(SIM_model5,D_SIM_model5)
pairs(SIM_model5,pars=c("beta","eta","output_incidence_cases[50]"))



# Posterior predictive check
check_hmc_diagnostics(S_model5)
print(S_model5,pars=c("beta","eta","epsilon","rho_K","pi","phi"),digits_summary=4)
plot_incidence_cases(S_model5,D_S_model5)
plot_incidence_deaths(S_model5,D_S_model5)

print(S_model5,pars=c("cfr_A_by_age","cfr_B_by_age","cfr_C_by_age","cfr_D_by_age"),digits_summary = 4)
print(S_model5,pars=c("cfr_A","cfr_B","cfr_C","cfr_D"),digits_summary = 4)

print(S_model5,pars=c("output_agedistr_cases"),digits_summary = 4)

summary(S_model5,pars=c("predicted_reported_incidence_cases","predicted_overall_incidence_cases","predicted_overall_incidence_deaths"))[[1]] %>%
  data.frame() %>%
  rownames_to_column() %>%
  tbl_df() %>%
  mutate(time=rep(2:D_S_model5$S,each=3*D_S_model5$K),
         date=time+as.Date("2019-12-01"),
         age_class=rep(c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),3*(D_S_model5$S-1)),
         type=rep(c("predicted_reported_incidence_cases","predicted_overall_incidence_cases","predicted_overall_incidence_deaths"),each=(D_S_model5$S-1)*D_S_model5$K)) %>% 
  ggplot(aes(x=time)) +
  geom_line(aes(y=mean)) +
  facet_grid(age_class ~ type)






# Model 4 ----------------------------------------------------------

# Format data and set priors
data_list_model6 = list(
  time_lag=100,
  K=9,
  age_dist=age_dist,
  pop_t = pop_t,
  t0=0,
  t_data = t_data,
  tmax =tmax,
  tmax2 = tmax,
  tswitch = tswitch,
  S=S,
  ts = 1:S,
  D = D,
  incidence_cases = incidence_cases,
  incidence_deaths = incidence_deaths,
  agedistr_cases = cases_tmax,
  agedistr_deaths = mort_tmax,
  p_beta=1,
  p_eta=c(1,1),
  p_delta=12.7,
  p_pi=c(1,999),
  p_incubation=5.2,
  p_infectious=2.3,
  p_asymptomatic_age80plus=1,
  p_epsilon=c(1,1),
  p_phi = 1/100,
  p_rho=c(1,1),
  p_xi=1,
  p_nu=1/5,
  inference=0,
  doprint=0
)

# Test
M_model6 = stan_model("model/model6.stan")
T_model6 = sampling(M_model6,
                    data = data_list_model6,
                    iter = 10,
                    chains = 1,
                    init=0,
                    control=list(max_treedepth=10,adapt_delta=0.8))
print(T_model6,pars=c("beta","eta","delta","epsilon","rho_K","pi"))

for(i in 1:5) {
  g = extract(T_model6,"comp_diffC")[[1]][i,,] %>%
    apply(.,1,sum) %>%
    data.frame() %>%
    tbl_df() %>%
    mutate(time=1:S) %>%
    filter(time<data_list_model6$tmax) %>%
    mutate(data_incidence=data_list_model6$incidence_cases) %>%
    ggplot() +
    geom_point(aes(x=time,y=data_incidence)) +
    geom_line(aes(x=time,y=.)) +
    geom_vline(xintercept=data_list_model6$tswitch,linetype=2)
  print(g)
}

# Create data and bash files
bashfile_rdump("model6",data_list_model6,warmup=500,iter=500,adapt_delta=0.8,max_depth=10,init=0.5,timelimit=96,chains=4,priorpredcheck=FALSE)

# Copy on cluster
system("scp /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/model6.stan /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/sb_model6.sh /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/data_SIM_model6.R /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/data_S_model6.R UBELIX:projects/COVID_age/model/.")

# Run on cluster
# >> R
# >> library(rstan)
# SIM_model6 = read_stan_csv(dir(".",pattern = 'SIM_model6_2020-02-24-16-54-12_[[:digit:]]+.csv'))
# D_SIM_model6 = read_rdump("data_SIM_model6.R")
# check_hmc_diagnostics(SIM_model6)
# print(SIM_model6,pars=c("beta","eta","epsilon","rho_K","pi"),digits_summary=4)
S_model6 = read_stan_csv(c(dir(".",pattern = 'S_model6_2020-02-26-13-38-32_[[:digit:]]+.csv')[c(1,2,4)],
                           dir(".",pattern = 'S_model7_2020-02-26-13-38-44_[[:digit:]]+.csv')[2]))
D_S_model6 = read_rdump("data_S_model6.R")
check_hmc_diagnostics(S_model6)
print(S_model6,pars=c("beta","eta","epsilon","rho_K","pi"),digits_summary=4)
save(
  # SIM_model6,D_SIM_model6,
     S_model6,D_S_model6,
     file="model6_2020-02-26.Rdata")

# Retrieve from cluster
system("scp  UBELIX:projects/COVID_age/model/model6_2020-02-26.Rdata /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/.")

# Load
load("model/model6_2020-02-26.Rdata")

# Prior predictive check
check_hmc_diagnostics(SIM_model6)
print(SIM_model6,pars=c("beta","eta","epsilon","delta","rho_K","pi","phi"),digits_summary=4)
plot_incidence_cases(SIM_model6,D_SIM_model6)
plot_incidence_deaths(SIM_model6,D_SIM_model6)
pairs(SIM_model6,pars=c("beta","eta","output_incidence_cases[50]"))

# Posterior predictive check
check_hmc_diagnostics(S_model6)
print(S_model6,pars=c("beta","eta","epsilon","delta","rho_K","pi","phi"),digits_summary=4)

g1 = plot_incidence_cases(S_model6,D_S_model6,col=col.cases)
g2 = plot_incidence_deaths(S_model6,D_S_model6,col1=col.deaths,col2=col.deaths.adj)
g3 = plot_incidence_deaths(S_model6,D_S_model6,col1=col.deaths,col2=col.deaths.adj,show.after.tmax = TRUE)
g4 = plot_agedist_cases(S_model6,D_S_model6,col1=col.cases,col2=col.cases.adj)
g5 = plot_agedist_deaths(S_model6,D_S_model6,col1=col.deaths,col2=col.deaths.adj)

data_legend = data.frame(type=c("Reported cases","Overall cases","Reported deaths","Overall deaths"),
                         col=c(col.cases,col.cases.adj,col.deaths,col.deaths.adj),
                         x=1:4,y=1:4,stringsAsFactors = FALSE)
leg = ggplot(data_legend) +
  geom_point(aes(x=x,y=y,fill=factor(type,levels=c("Reported cases","Overall cases","Reported deaths","Overall deaths"))),
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
  plot_grid(gnull,legend,rel_widths = c(.3,1)),
  ncol=1,rel_heights = c(10,1))








print(S_model6,pars=c("cfr_A_by_age","cfr_B_by_age","cfr_C_by_age","cfr_D_by_age"),digits_summary = 4)
print(S_model6,pars=c("cfr_A","cfr_B","cfr_C","cfr_D"),digits_summary = 4)

print(S_model6,pars=c("output_agedistr_cases"),digits_summary = 4)

summary(S_model6,pars=c("predicted_reported_incidence_cases","predicted_overall_incidence_cases","predicted_overall_incidence_deaths"))[[1]] %>%
  data.frame() %>%
  rownames_to_column() %>%
  tbl_df() %>%
  mutate(time=rep(2:D_S_model6$S,each=3*D_S_model6$K),
         date=time+as.Date("2019-12-01"),
         age_class=rep(c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),3*(D_S_model6$S-1)),
         type=rep(c("predicted_reported_incidence_cases","predicted_overall_incidence_cases","predicted_overall_incidence_deaths"),each=(D_S_model6$S-1)*D_S_model6$K)) %>% 
  ggplot(aes(x=time)) +
  geom_line(aes(y=mean)) +
  facet_grid(age_class ~ type)




# Model 7 ----------------------------------------------------------

# Format data and set priors
data_list_model7 = list(
  time_lag=100,
  K=9,
  age_dist=age_dist,
  pop_t = pop_t,
  t0=0,
  t_data = t_data,
  tmax =tmax,
  tmax2 = tmax,
  tswitch = tswitch,
  S=S,
  ts = 1:S,
  D = D,
  incidence_cases = incidence_cases,
  incidence_deaths = incidence_deaths,
  agedistr_cases = cases_tmax,
  agedistr_deaths = mort_tmax,
  p_beta=1,
  p_eta=c(1,1),
  p_delta=12.7,
  p_pi=c(1,999),
  p_incubation=5.2,
  p_infectious=2.3,
  p_asymptomatic_age80plus=1,
  p_epsilon=c(1,1),
  p_phi = 1/100,
  p_rho=c(1,1),
  p_xi=1,
  p_nu=1/5,
  inference=0,
  doprint=0
)

# Test
M_model7 = stan_model("model/model7.stan")
T_model7 = sampling(M_model7,
                    data = data_list_model7,
                    iter = 10,
                    chains = 1,
                    init=0,
                    control=list(max_treedepth=10,adapt_delta=0.8))
print(T_model7,pars=c("beta","eta","delta","epsilon","rho_K","pi"))

for(i in 1:5) {
  g = extract(T_model7,"comp_diffC")[[1]][i,,] %>%
    apply(.,1,sum) %>%
    data.frame() %>%
    tbl_df() %>%
    mutate(time=1:S) %>%
    filter(time<data_list_model7$tmax) %>%
    mutate(data_incidence=data_list_model7$incidence_cases) %>%
    ggplot() +
    geom_point(aes(x=time,y=data_incidence)) +
    geom_line(aes(x=time,y=.)) +
    geom_vline(xintercept=data_list_model7$tswitch,linetype=2)
  print(g)
}

# Create data and bash files
bashfile_rdump("model7",data_list_model7,warmup=500,iter=500,adapt_delta=0.8,max_depth=10,init=0.5,timelimit=96,chains=4,priorpredcheck=FALSE)

# Copy on cluster
system("scp /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/model7.stan /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/sb_model7.sh /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/data_SIM_model7.R /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/data_S_model7.R UBELIX:projects/COVID_age/model/.")

# Run on cluster
# >> R
# >> library(rstan)
SIM_model7 = read_stan_csv(dir(".",pattern = 'SIM_model7_2020-02-24-16-54-12_[[:digit:]]+.csv'))
D_SIM_model7 = read_rdump("data_SIM_model7.R")
check_hmc_diagnostics(SIM_model7)
print(SIM_model7,pars=c("beta","eta","epsilon","rho_K","pi"),digits_summary=4)
S_model7 = read_stan_csv(dir(".",pattern = 'S_model7_2020-02-24-16-54-12_[[:digit:]]+.csv'))
D_S_model7 = read_rdump("data_S_model7.R")
check_hmc_diagnostics(S_model7)
print(S_model7,pars=c("beta","eta","epsilon","rho_K","pi"),digits_summary=4)
save(SIM_model7,D_SIM_model7,
     S_model7,D_S_model7,
     file="model7_2020-02-24.Rdata")

# Retrieve from cluster
system("scp  UBELIX:projects/COVID_age/model/model7_2020-02-24.Rdata /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/.")

# Load
load("model/model7_2020-02-24.Rdata")

# Prior predictive check
check_hmc_diagnostics(SIM_model7)
print(SIM_model7,pars=c("beta","eta","epsilon","rho_K","pi","phi"),digits_summary=4)
plot_incidence_cases(SIM_model7,D_SIM_model7)
plot_incidence_deaths(SIM_model7,D_SIM_model7)
pairs(SIM_model7,pars=c("beta","eta","output_incidence_cases[50]"))

for(i in 1:10) {
  g = extract(SIM_model7,"y")[[1]][i,,37:45] %>%
    apply(.,1,function(x) sum(x) *D_SIM_model7$pop_t ) %>%
    data.frame() %>%
    tbl_df() %>%
    mutate(time=1:S,
           incidence=.-lag(.,default = 0)) %>%
    filter(time<D_SIM_model7$tmax) %>%
    mutate(data_incidence=D_SIM_model7$incidence_cases) %>%
    ggplot() +
    geom_point(aes(x=time,y=data_incidence)) +
    geom_line(aes(x=time,y=incidence)) +
    geom_vline(xintercept=D_SIM_model7$tswitch,linetype=2)
  print(g)
}



# Posterior predictive check
check_hmc_diagnostics(S_model7)
print(S_model7,pars=c("beta","eta","epsilon","rho_K","pi","phi"),digits_summary=4)
plot_incidence_cases(S_model7,D_S_model7)
plot_incidence_deaths(S_model7,D_S_model7)

print(S_model7,pars=c("cfr_A_by_age","cfr_B_by_age","cfr_C_by_age","cfr_D_by_age"),digits_summary = 4)
print(S_model7,pars=c("cfr_A","cfr_B","cfr_C","cfr_D"),digits_summary = 4)

print(S_model7,pars=c("output_agedistr_cases"),digits_summary = 4)

summary(S_model7,pars=c("predicted_reported_incidence_cases","predicted_overall_incidence_cases","predicted_overall_incidence_deaths"))[[1]] %>%
  data.frame() %>%
  rownames_to_column() %>%
  tbl_df() %>%
  mutate(time=rep(2:D_S_model7$S,each=3*D_S_model7$K),
         date=time+as.Date("2019-12-01"),
         age_class=rep(c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),3*(D_S_model7$S-1)),
         type=rep(c("predicted_reported_incidence_cases","predicted_overall_incidence_cases","predicted_overall_incidence_deaths"),each=(D_S_model7$S-1)*D_S_model7$K)) %>% 
  ggplot(aes(x=time)) +
  geom_line(aes(y=mean)) +
  facet_grid(age_class ~ type)








# Model 8 ----------------------------------------------------------

# Format data and set priors
data_list_model8 = list(
  time_lag=100,
  K=9,
  age_dist=age_dist,
  pop_t = pop_t,
  t0=0,
  t_data = t_data,
  tmax =tmax,
  tmax2 = tmax,
  tswitch = tswitch,
  S=S,
  ts = 1:S,
  D = D,
  incidence_cases = incidence_cases,
  incidence_deaths = incidence_deaths,
  agedistr_cases = cases_tmax,
  agedistr_deaths = mort_tmax,
  p_beta=1,
  p_eta=c(1,1),
  p_pi=c(1,999),
  p_incubation=5.6, # Incubation including WR (Linton)
  p_infectious=2.9, # Onset to hospital admission (Liu)
  p_isolation_to_death=c(34.6,2), # Onset to death - 2.9 (Linton)
  p_symptomatic=0.48, # Diamond princess
  p_epsilon=c(1,1),
  p_phi = 1/100,
  p_rho=c(1,1),
  p_xi=1,
  p_nu=1/5,
  inference=0,
  doprint=0
)

# Test
M_model8 = stan_model("model/model8.stan")
T_model8 = sampling(M_model8,
                    data = data_list_model8,
                    iter = 10,
                    chains = 1,
                    init=0,
                    control=list(max_treedepth=10,adapt_delta=0.8))
print(T_model8,pars=c("beta","eta","delta","epsilon","rho_K","pi"))

# Create data and bash files
bashfile_rdump("model8",data_list_model8,warmup=500,iter=500,adapt_delta=0.8,max_depth=10,init=0.5,timelimit=96,chains=4,priorpredcheck=FALSE)

# Copy on cluster
system("scp /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/model8.stan /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/sb_model8.sh /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/data_SIM_model8.R /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/data_S_model8.R UBELIX:projects/COVID_age/model/.")

# Run on cluster
# >> R
# >> library(rstan)
SIM_model8 = read_stan_csv(dir(".",pattern = 'SIM_model8_2020-02-24-16-54-12_[[:digit:]]+.csv'))
D_SIM_model8 = read_rdump("data_SIM_model8.R")
check_hmc_diagnostics(SIM_model8)
print(SIM_model8,pars=c("beta","eta","epsilon","rho_K","pi"),digits_summary=4)
S_model8 = read_stan_csv(dir(".",pattern = 'S_model8_2020-02-24-16-54-12_[[:digit:]]+.csv'))
D_S_model8 = read_rdump("data_S_model8.R")
check_hmc_diagnostics(S_model8)
print(S_model8,pars=c("beta","eta","epsilon","rho_K","pi"),digits_summary=4)
save(SIM_model8,D_SIM_model8,
     S_model8,D_S_model8,
     file="model8_2020-02-24.Rdata")

# Retrieve from cluster
system("scp  UBELIX:projects/COVID_age/model/model8_2020-02-24.Rdata /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/.")

# Load
load("model/model8_2020-02-24.Rdata")

# Prior predictive check
check_hmc_diagnostics(SIM_model8)
print(SIM_model8,pars=c("beta","eta","epsilon","rho_K","pi","phi"),digits_summary=4)
plot_incidence_cases(SIM_model8,D_SIM_model8)
plot_incidence_deaths(SIM_model8,D_SIM_model8)
pairs(SIM_model8,pars=c("beta","eta","output_incidence_cases[50]"))

for(i in 1:10) {
  g = extract(SIM_model8,"y")[[1]][i,,37:45] %>%
    apply(.,1,function(x) sum(x) *D_SIM_model8$pop_t ) %>%
    data.frame() %>%
    tbl_df() %>%
    mutate(time=1:S,
           incidence=.-lag(.,default = 0)) %>%
    filter(time<D_SIM_model8$tmax) %>%
    mutate(data_incidence=D_SIM_model8$incidence_cases) %>%
    ggplot() +
    geom_point(aes(x=time,y=data_incidence)) +
    geom_line(aes(x=time,y=incidence)) +
    geom_vline(xintercept=D_SIM_model8$tswitch,linetype=2)
  print(g)
}



# Posterior predictive check
check_hmc_diagnostics(S_model8)
print(S_model8,pars=c("beta","eta","epsilon","rho_K","pi","phi"),digits_summary=4)
plot_incidence_cases(S_model8,D_S_model8)
plot_incidence_deaths(S_model8,D_S_model8)

print(S_model8,pars=c("cfr_A_by_age","cfr_B_by_age","cfr_C_by_age","cfr_D_by_age"),digits_summary = 4)
print(S_model8,pars=c("cfr_A","cfr_B","cfr_C","cfr_D"),digits_summary = 4)

print(S_model8,pars=c("output_agedistr_cases"),digits_summary = 4)

summary(S_model8,pars=c("predicted_reported_incidence_cases","predicted_overall_incidence_cases","predicted_overall_incidence_deaths"))[[1]] %>%
  data.frame() %>%
  rownames_to_column() %>%
  tbl_df() %>%
  mutate(time=rep(2:D_S_model8$S,each=3*D_S_model8$K),
         date=time+as.Date("2019-12-01"),
         age_class=rep(c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),3*(D_S_model8$S-1)),
         type=rep(c("predicted_reported_incidence_cases","predicted_overall_incidence_cases","predicted_overall_incidence_deaths"),each=(D_S_model8$S-1)*D_S_model8$K)) %>% 
  ggplot(aes(x=time)) +
  geom_line(aes(y=mean)) +
  facet_grid(age_class ~ type)








# Model 9 ----------------------------------------------------------

# Format data and set priors
data_list_model9 = list(
  time_lag=100,
  K=9,
  age_dist=age_dist,
  pop_t = pop_t,
  t0=0,
  t_data = t_data,
  tmax =tmax,
  tmax2 = tmax,
  tswitch = tswitch,
  S=S,
  ts = 1:S,
  D = D,
  incidence_cases = incidence_cases,
  incidence_deaths = incidence_deaths,
  agedistr_cases = cases_tmax,
  agedistr_deaths = mort_tmax,
  p_beta=1,
  p_eta=c(1,1),
  p_pi=c(1,999),
  p_incubation=5.6, # Incubation including WR (Linton)
  p_infectious=2.9, # Onset to hospital admission (Liu)
  p_isolation_to_death=17.3, # Onset to death - 2.9 (Linton)
  p_symptomatic=0.48, # Diamond princess
  p_epsilon=c(1,1),
  p_phi = 1/100,
  p_rho=c(1,1),
  p_xi=1,
  p_nu=1/5,
  inference=0,
  doprint=0
)

# Test
M_model9 = stan_model("model/model9.stan")
T_model9 = sampling(M_model9,
                    data = data_list_model9,
                    iter = 10,
                    chains = 1,
                    init=0,
                    control=list(max_treedepth=10,adapt_delta=0.8))
print(T_model9,pars=c("beta","eta","delta","epsilon","rho_K","pi"))

# Create data and bash files
bashfile_rdump("model9",data_list_model9,warmup=500,iter=500,adapt_delta=0.8,max_depth=10,init=0.5,timelimit=96,chains=4,priorpredcheck=FALSE)

# Copy on cluster
system("scp /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/model9.stan /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/sb_model9.sh /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/data_SIM_model9.R /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/data_S_model9.R UBELIX:projects/COVID_age/model/.")

# Run on cluster
# >> R
# >> library(rstan)
SIM_model9 = read_stan_csv(dir(".",pattern = 'SIM_model9_2020-02-24-16-54-12_[[:digit:]]+.csv'))
D_SIM_model9 = read_rdump("data_SIM_model9.R")
check_hmc_diagnostics(SIM_model9)
print(SIM_model9,pars=c("beta","eta","epsilon","rho_K","pi"),digits_summary=4)
S_model9 = read_stan_csv(dir(".",pattern = 'S_model9_2020-02-24-16-54-12_[[:digit:]]+.csv'))
D_S_model9 = read_rdump("data_S_model9.R")
check_hmc_diagnostics(S_model9)
print(S_model9,pars=c("beta","eta","epsilon","rho_K","pi"),digits_summary=4)
save(SIM_model9,D_SIM_model9,
     S_model9,D_S_model9,
     file="model9_2020-02-24.Rdata")

# Retrieve from cluster
system("scp  UBELIX:projects/COVID_age/model/model9_2020-02-24.Rdata /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/.")

# Load
load("model/model9_2020-02-24.Rdata")

# Prior predictive check
check_hmc_diagnostics(SIM_model9)
print(SIM_model9,pars=c("beta","eta","epsilon","rho_K","pi","phi"),digits_summary=4)
plot_incidence_cases(SIM_model9,D_SIM_model9)
plot_incidence_deaths(SIM_model9,D_SIM_model9)
pairs(SIM_model9,pars=c("beta","eta","output_incidence_cases[50]"))

for(i in 1:10) {
  g = extract(SIM_model9,"y")[[1]][i,,37:45] %>%
    apply(.,1,function(x) sum(x) *D_SIM_model9$pop_t ) %>%
    data.frame() %>%
    tbl_df() %>%
    mutate(time=1:S,
           incidence=.-lag(.,default = 0)) %>%
    filter(time<D_SIM_model9$tmax) %>%
    mutate(data_incidence=D_SIM_model9$incidence_cases) %>%
    ggplot() +
    geom_point(aes(x=time,y=data_incidence)) +
    geom_line(aes(x=time,y=incidence)) +
    geom_vline(xintercept=D_SIM_model9$tswitch,linetype=2)
  print(g)
}



# Posterior predictive check
check_hmc_diagnostics(S_model9)
print(S_model9,pars=c("beta","eta","epsilon","rho_K","pi","phi"),digits_summary=4)
plot_incidence_cases(S_model9,D_S_model9)
plot_incidence_deaths(S_model9,D_S_model9)

print(S_model9,pars=c("cfr_A_by_age","cfr_B_by_age","cfr_C_by_age","cfr_D_by_age"),digits_summary = 4)
print(S_model9,pars=c("cfr_A","cfr_B","cfr_C","cfr_D"),digits_summary = 4)

print(S_model9,pars=c("output_agedistr_cases"),digits_summary = 4)

summary(S_model9,pars=c("predicted_reported_incidence_cases","predicted_overall_incidence_cases","predicted_overall_incidence_deaths"))[[1]] %>%
  data.frame() %>%
  rownames_to_column() %>%
  tbl_df() %>%
  mutate(time=rep(2:D_S_model9$S,each=3*D_S_model9$K),
         date=time+as.Date("2019-12-01"),
         age_class=rep(c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"),3*(D_S_model9$S-1)),
         type=rep(c("predicted_reported_incidence_cases","predicted_overall_incidence_cases","predicted_overall_incidence_deaths"),each=(D_S_model9$S-1)*D_S_model9$K)) %>% 
  ggplot(aes(x=time)) +
  geom_line(aes(y=mean)) +
  facet_grid(age_class ~ type)









# Model 10 ----------------------------------------------------------

# get gamma from discretizing onset_to_death in Linton
linton_mean = 20.2
linton_sd = 11.6

get_par_lnorm = function(m,s) {
  mu = log(m) - 1/2*log((s/m)^2+1)
  sigma = sqrt(log((s/m)^2+1))
  return(list(mu=mu,sigma=sigma))
}
linton_pars = get_par_lnorm(linton_mean,linton_sd)

gamma = numeric(60)
for(i in 1:60) {
  gamma[i] = plnorm(i+.5,linton_pars$mu,linton_pars$sigma)-plnorm(i-.5,linton_pars$mu,linton_pars$sigma)
}
gamma = gamma/sum(gamma)
plot(density(rlnorm(1000,linton_pars$mu,linton_pars$sigma)))
lines(gamma,col="red")

# Format data and set priors
data_list_model10 = list(
  time_lag=100,
  K=9,
  age_dist=age_dist,
  pop_t = pop_t,
  t0=0,
  t_data = t_data,
  tmax =tmax,
  tmax2 = tmax,
  tswitch = tswitch,
  S=S,
  ts = 1:S,
  D = D,
  incidence_cases = incidence_cases,
  incidence_deaths = incidence_deaths,
  agedistr_cases = cases_tmax,
  agedistr_deaths = mort_tmax,
  p_beta=1,
  p_eta=c(1,1),
  p_pi=c(1,999),
  p_incubation=5.6, # Incubation including WR (Linton)
  p_infectious=2.9, # Infectious period (Liu)
  p_psi=0.49, # Diamond princess
  p_epsilon=c(1,1),
  p_phi = 1/100,
  p_rho=c(1,1),
  p_xi=1,
  p_nu=1/5,
  G=60,
  p_gamma=gamma,
  inference=0,
  doprint=0
)

# Test
M_model10 = stan_model("model/model10.stan")
T_model10 = sampling(M_model10,
                    data = data_list_model10,
                    iter = 5,
                    chains = 1,
                    init=0,
                    control=list(max_treedepth=10,adapt_delta=0.8))
print(T_model10,pars=c("beta","eta","delta","epsilon","rho_K","pi"))

# Create data and bash files
bashfile_rdump("model10",data_list_model10,warmup=500,iter=500,adapt_delta=0.8,max_depth=10,init=0.5,timelimit=96,chains=4,priorpredcheck=FALSE)

# Copy on cluster
system("scp /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/model10.stan /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/sb_model10.sh /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/data_SIM_model10.R /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/data_S_model10.R UBELIX:projects/COVID_age/model/.")

# Run on cluster
# >> R
# >> library(rstan)
SIM_model10 = read_stan_csv(dir(".",pattern = 'SIM_model10_2020-02-28-17-37-21_[[:digit:]]+.csv'))
D_SIM_model10 = read_rdump("data_SIM_model10.R")
check_hmc_diagnostics(SIM_model10)
print(SIM_model10,pars=c("beta","eta","epsilon","rho_K","pi"),digits_summary=4)
S_model10 = read_stan_csv(dir(".",pattern = 'S_model10_2020-02-28-17-37-21_[[:digit:]]+.csv'))
D_S_model10 = read_rdump("data_S_model10.R")
check_hmc_diagnostics(S_model10)
print(S_model10,pars=c("beta","eta","epsilon","rho_K","pi"),digits_summary=4)
save(
  # SIM_model10,D_SIM_model10,
     S_model10,D_S_model10,
     file="model10_2020-02-28.Rdata")

# Retrieve from cluster
system("scp  UBELIX:projects/COVID_age/model/model10_2020-02-28.Rdata /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/.")

# Load
load("model/model10_2020-02-28.Rdata")

# Prior predictive check
check_hmc_diagnostics(SIM_model10)
print(SIM_model10,pars=c("beta","eta","epsilon","rho_K","pi","phi"),digits_summary=4)
plot_incidence_cases(SIM_model10,D_SIM_model10)
plot_incidence_deaths(SIM_model10,D_SIM_model10)
pairs(SIM_model10,pars=c("beta","eta","output_incidence_cases[50]"))

# Posterior predictive check
check_hmc_diagnostics(S_model10)
print(S_model10,pars=c("beta","eta","epsilon","rho_K","pi","phi"),digits_summary=4)
plot_incidence_cases(S_model10,D_S_model10)
plot_incidence_deaths(S_model10,D_S_model10,show.after.tmax = TRUE)

print(S_model10,pars=c("cfr_A_by_age","cfr_B_by_age","cfr_C_by_age","cfr_D_by_age"),digits_summary = 4)
print(S_model10,pars=c("cfr_A","cfr_B","cfr_C","cfr_D"),digits_summary = 4)

print(S_model10,pars=c("output_agedistr_cases"),digits_summary = 4)


g1 = plot_incidence_cases(S_model10,D_S_model10,col=col.cases)
g2 = plot_incidence_deaths(S_model10,D_S_model10,col1=col.deaths,col2=col.deaths.adj)
g3 = plot_incidence_deaths(S_model10,D_S_model10,col1=col.deaths,col2=col.deaths.adj,show.after.tmax = TRUE)
g4 = plot_agedist_cases(S_model10,D_S_model10,col1=col.cases,col2=col.cases.adj)
g5 = plot_agedist_deaths(S_model10,D_S_model10,col1=col.deaths,col2=col.deaths.adj)


data_legend = data.frame(type=c("Reported cases","Overall cases","Reported deaths","Overall deaths"),
                         col=c(col.cases,col.cases.adj,col.deaths,col.deaths.adj),
                         x=1:4,y=1:4,stringsAsFactors = FALSE)
leg = ggplot(data_legend) +
  geom_point(aes(x=x,y=y,fill=factor(type,levels=c("Reported cases","Overall cases","Reported deaths","Overall deaths"))),
             colour="black",shape=22,size=6) +
  scale_fill_manual(values=data_legend$col,guide=guide_legend()) +
  labs(fill=NULL) +  
  theme(legend.position="bottom",
        legend.text=element_text(size=12),
        legend.key.height=unit(12,"pt"),
        legend.margin=margin(4,5,4,4,"pt"),
        legend.box.background = element_rect(fill="white",colour = "grey60",size=.4)
  )
leg
legend <- cowplot::get_legend(leg)
gnull = ggplot() +  geom_blank() + theme_minimal()

plot_grid(
  plot_grid(g1,g4,g3,g5,nrow=2,rel_widths = c(2,1),labels=LETTERS[1:4]),
  plot_grid(gnull,legend,rel_widths = c(.3,1)),
  ncol=1,rel_heights = c(10,1))







# Model 11 ----------------------------------------------------------

# get gamma from discretizing onset_to_death in Linton
linton_mean = 20.2
linton_sd = 11.6

get_par_lnorm = function(m,s) {
  mu = log(m) - 1/2*log((s/m)^2+1)
  sigma = sqrt(log((s/m)^2+1))
  return(list(mu=mu,sigma=sigma))
}
linton_pars = get_par_lnorm(linton_mean,linton_sd)

gamma = numeric(60)
for(i in 1:60) {
  gamma[i] = plnorm(i+.5,linton_pars$mu,linton_pars$sigma)-plnorm(i-.5,linton_pars$mu,linton_pars$sigma)
}
gamma = gamma/sum(gamma)
plot(density(rlnorm(1000,linton_pars$mu,linton_pars$sigma)))
lines(gamma,col="red")

# Format data and set priors
data_list_model11 = list(
  time_lag=100,
  K=9,
  age_dist=age_dist,
  pop_t = pop_t,
  t0=0,
  t_data = t_data,
  tmax =tmax,
  tmax2 = tmax,
  tswitch = tswitch,
  S=S,
  ts = 1:S,
  D = D,
  incidence_cases = incidence_cases,
  incidence_deaths = incidence_deaths,
  agedistr_cases = cases_tmax,
  agedistr_deaths = mort_tmax,
  p_beta=1,
  p_eta=c(1,1),
  p_pi=c(1,999),
  p_incubation=5.6, # Incubation including WR (Linton)
  p_infectious=2.9, # Infectious period (Liu)
  p_psi=0.48, # Diamond princess
  p_epsilon=c(1,1),
  p_phi = 1/100,
  p_rho=c(1,1),
  p_xi=1,
  p_nu=1/5,
  G=60,
  p_gamma=gamma,
  inference=0,
  doprint=0,
  contact=c(0.810810810810811, 0.0908559523547268, 0.372736406439194, 1.27360250772, 0.200569529052988, 0.375083342749019, 0.60252680195839, 0.0934189610338407, 0.0225225225225225, 0.0904095466183592, 2.4392523364486, 0.140093983348316, 0.706545801082683, 0.942937990573664, 0.27920963239528, 0.326366336169345, 0.196893495540358, 0.106045179398683, 0.289504965045504, 0.109348487445688, 1.76086956521739, 0.923069180041088, 0.93772012267962, 0.965186137047983, 0.274120168579709, 0.116564256844925, 0.0773400190233669, 0.91820215964926, 0.511898453884688, 0.85680985412458, 2.70542635658915, 1.41323192857305, 0.993399938008648, 0.719603621821669, 0.146103509716984, 0.07633130138862, 0.13119227828341, 0.619819944222649, 0.789700390093264, 1.28218991206025, 2.17699115044248, 1.1199461877854, 0.514253349451317, 0.496466649026704, 0.101504389707241, 0.259078294801222, 0.193808465356441, 0.858341528544101, 0.951750199084178, 1.18265232149625, 2.31730769230769, 0.977037933291252, 0.606164987575222, 0.4393566902894, 0.552747314447092, 0.300880970126328, 0.323770338070664, 0.915670466885606, 0.721247101248993, 1.29765260904839, 2.76146788990826, 0.959867553314515, 0.340125585278128, 0.121809161946423, 0.25799743320884, 0.1956843612527, 0.264241585561661, 0.989672909331423, 1.14428055461657, 1.36428769674242, 1.96363636363636, 1.0266513139522, 0.0447824174066075, 0.211894911958445, 0.197988778289041, 0.210517772531686, 0.308554588199316, 1.26474943927563, 0.737190168823191, 1.56555579008225, 2.0625)
)

# Test
M_model11 = stan_model("model/model11.stan")
T_model11 = sampling(M_model11,
                     data = data_list_model11,
                     iter = 5,
                     chains = 1,
                     init=0,
                     control=list(max_treedepth=10,adapt_delta=0.8))
print(T_model11,pars=c("beta","eta","delta","epsilon","rho_K","pi"))

# Create data and bash files
bashfile_rdump("model11",data_list_model11,warmup=500,iter=500,adapt_delta=0.8,max_depth=10,init=0.5,timelimit=96,chains=4,priorpredcheck=FALSE)

# Copy on cluster
system("scp /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/model11.stan /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/sb_model11.sh /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/data_SIM_model11.R /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/data_S_model11.R UBELIX:projects/COVID_age/model/.")

# Run on cluster
# >> R
# >> library(rstan)
SIM_model11 = read_stan_csv(dir(".",pattern = 'SIM_model11_2020-02-28-17-37-21_[[:digit:]]+.csv'))
D_SIM_model11 = read_rdump("data_SIM_model11.R")
check_hmc_diagnostics(SIM_model11)
print(SIM_model11,pars=c("beta","eta","epsilon","rho_K","pi"),digits_summary=4)
S_model11 = read_stan_csv(dir(".",pattern = 'S_model11_2020-02-28-17-37-21_[[:digit:]]+.csv'))
D_S_model11 = read_rdump("data_S_model11.R")
check_hmc_diagnostics(S_model11)
print(S_model11,pars=c("beta","eta","epsilon","rho_K","pi"),digits_summary=4)
save(
  # SIM_model11,D_SIM_model11,
  S_model11,D_S_model11,
  file="model11_2020-02-28.Rdata")

# Retrieve from cluster
system("scp  UBELIX:projects/COVID_age/model/model11_2020-02-28.Rdata /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/.")

# Load
load("model/model11_2020-02-28.Rdata")

# Prior predictive check
check_hmc_diagnostics(SIM_model11)
print(SIM_model11,pars=c("beta","eta","epsilon","rho_K","pi","phi"),digits_summary=4)
plot_incidence_cases(SIM_model11,D_SIM_model11)
plot_incidence_deaths(SIM_model11,D_SIM_model11)
pairs(SIM_model11,pars=c("beta","eta","output_incidence_cases[50]"))

# Posterior predictive check
check_hmc_diagnostics(S_model11)
print(S_model11,pars=c("beta","eta","epsilon","rho_K","pi","phi"),digits_summary=4)
plot_incidence_cases(S_model11,D_S_model11)
plot_incidence_deaths(S_model11,D_S_model11,show.after.tmax = TRUE)

print(S_model11,pars=c("cfr_A_by_age","cfr_B_by_age","cfr_C_by_age","cfr_D_by_age"),digits_summary = 4)
print(S_model11,pars=c("cfr_A","cfr_B","cfr_C","cfr_D"),digits_summary = 4)

print(S_model11,pars=c("output_agedistr_cases"),digits_summary = 4)


g1 = plot_incidence_cases(S_model11,D_S_model11,col=col.cases)
g2 = plot_incidence_deaths(S_model11,D_S_model11,col1=col.deaths,col2=col.deaths.adj)
g3 = plot_incidence_deaths(S_model11,D_S_model11,col1=col.deaths,col2=col.deaths.adj,show.after.tmax = TRUE)
g4 = plot_agedist_cases(S_model11,D_S_model11,col1=col.cases,col2=col.cases.adj)
g5 = plot_agedist_deaths(S_model11,D_S_model11,col1=col.deaths,col2=col.deaths.adj)


data_legend = data.frame(type=c("Reported cases","Overall cases","Reported deaths","Overall deaths"),
                         col=c(col.cases,col.cases.adj,col.deaths,col.deaths.adj),
                         x=1:4,y=1:4,stringsAsFactors = FALSE)
leg = ggplot(data_legend) +
  geom_point(aes(x=x,y=y,fill=factor(type,levels=c("Reported cases","Overall cases","Reported deaths","Overall deaths"))),
             colour="black",shape=22,size=6) +
  scale_fill_manual(values=data_legend$col,guide=guide_legend()) +
  labs(fill=NULL) +  
  theme(legend.position="bottom",
        legend.text=element_text(size=12),
        legend.key.height=unit(12,"pt"),
        legend.margin=margin(4,5,4,4,"pt"),
        legend.box.background = element_rect(fill="white",colour = "grey60",size=.4)
  )
leg
legend <- cowplot::get_legend(leg)
gnull = ggplot() +  geom_blank() + theme_minimal()

plot_grid(
  plot_grid(g1,g4,g3,g5,nrow=2,rel_widths = c(2,1),labels=LETTERS[1:4]),
  plot_grid(gnull,legend,rel_widths = c(.3,1)),
  ncol=1,rel_heights = c(10,1))
