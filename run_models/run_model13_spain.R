# Setup ----
source("setup.R")
source("data/spain/data_management_spain.R")

## Compute discrete distribution of time from onset to death ----

# Fom Linton et al
linton_mean = 20.2
linton_sd = 11.6
# Get lognormal parameters from mean and sd
get_par_lnorm = function(m,s) {
  mu = log(m) - 1/2*log((s/m)^2+1)
  sigma = sqrt(log((s/m)^2+1))
  return(list(mu=mu,sigma=sigma))
}
linton_pars = get_par_lnorm(linton_mean,linton_sd)
# Discretize
gamma = numeric(60)
for(i in 1:60) {
  gamma[i] = plnorm(i+.5,linton_pars$mu,linton_pars$sigma)-plnorm(i-.5,linton_pars$mu,linton_pars$sigma)
}
# Normalize
gamma = gamma/sum(gamma)
# plot(density(rlnorm(1000,linton_pars$mu,linton_pars$sigma)))
# lines(gamma,col="red")


## Format for stan

t0 = 0
S = as.numeric(day_max-day_start)
t_data = as.numeric(day_data-day_start)
ts = t_data:S
D = S-t_data+1
tswitch = as.numeric(day_quarantine-day_start)
data_list_model13 = list(
  K=9,
  age_dist=age_dist,
  pop_t = pop_t,
  t0=0,
  t_data = t_data,
  tswitch = tswitch,
  S=S,
  ts = ts,
  D = length(incidence_cases),
  incidence_cases = incidence_cases,
  incidence_deaths = incidence_deaths,
  agedistr_cases = cases_tmax,
  agedistr_deaths = mort_tmax,
  p_beta=1,
  p_eta=c(1,1),
  p_pi=c(1,999),
  p_incubation=5.95, # Bi et al
  p_infectious=2.4, # Bi et al
  p_epsilon=c(1,1),
  p_phi = 1/100,
  p_rho=c(1,1),
  p_xi=1,
  p_nu=1/5,
  p_chi=5,
  G=60,
  p_gamma=gamma,
  inference=0,
  doprint=0,
  contact=c(5.35148784013064, 0.792403698555513, 0.597888237674624, 1.81141261889314, 1.15880446573503, 0.426179823985014, 0.272566579444233, 0.0667965764245607, 0.0524292496334014, 0.962410237985733, 6.96132818205631, 1.05406336733019, 1.13101315162061, 1.62168368448698, 0.541175903065237, 0.131491288694527, 0.0540877655519872, 0.0424539866267057, 0.452019159060351, 1.34005151271159, 5.2217886689935, 3.01397324511203, 2.22431927441854, 1.43787568086573, 0.197174914122326, 0.049091677623549, 0.038532511077935, 1.86647113340635, 1.27397280273288, 2.63210519557764, 5.87106096844647, 3.36191516449852, 1.4530452744515, 0.523291006732789, 0.10038463556598, 0.0787928274047232, 1.10821166176017, 1.96582937687123, 2.05373971751821, 3.3943240797736, 4.28881981759823, 1.48129197804708, 0.305549086413716, 0.108848089031009, 0.0854358701807222, 0.894549948516997, 1.57636342878353, 1.99319025424205, 2.30845452446609, 2.57838327243766, 2.46844694468242, 0.521659898101079, 0.123368775343464, 0.0968332909509898, 0.616712123753189, 0.467315717289364, 0.722196642987352, 1.53068487644559, 0.963926067358351, 0.932445331806381, 1.58924678766151, 0.22742940762708, 0.178511442123438, 0.49296743101773, 0.761168134312711, 0.392685418560048, 0.922634583854662, 1.19274137966985, 0.916395267151485, 0.868894131016553, 0.5953456757922, 0.46729231833502, 0.0825597022962617, 0.0695377130823058, 0.063939746734617, 0.176136931824655, 0.229871638980727, 0.228790005278394, 0.317695213030174, 0.5953456757922, 0.5953456757922))

# test model13 for spain 
M_model13 = stan_model("models/model13.stan")
T_model13 = sampling(M_model13,
                     data = data_list_model13,
                     iter = 5,
                     chains = 1,
                     init=0,
                     control=list(max_treedepth=10,adapt_delta=0.8))
print(T_model13,pars=c("beta","eta","epsilon","rho","pi"))


# Create data and bash files ----

# symptomatic = 50% (B) (Diamond princess https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.10.2000180) 
data_list_model13$p_psi=c(522,114) 
bashfile_rdump("model13",id="ES-B",data_list_model13,warmup=500,iter=500,adapt_delta=0.8,max_depth=12,init=0.5,timelimit=2,chains=4)

data_list_model13$inference=1
options(mc.cores = parallel::detectCores())
S_model13ITB = sampling(M_model13,
                     data = data_list_model13,
                     iter = 1000,
                     chains = 4,
                     init=0.5,
                     control=list(max_treedepth=12,adapt_delta=0.8))
D_S_model13ITB = data_list_model13
# Copy on cluster
# system("scp /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/run_models/sb_model13ES* /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/run_models/data_S_model13ES* UBELIX:projects/COVID_age/model/.")

# Copy back posterior samples
system("scp  UBELIX:projects/COVID_age/model/S_model13ES-B_2020-03-26-10-36-33*  /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/posterior_samples/.")
system("scp  UBELIX:projects/COVID_age/model/data_S_model13ES-B_2020-03-26-10-36-33.R /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/posterior_samples/.")

# Load posterior samples 
D_S_model13ESB = read_rdump("posterior_samples/data_S_model13ES-B_2020-03-26-10-36-33.R")
S_model13ESB = read_stan_csv(paste0("posterior_samples/",dir("posterior_samples",pattern = 'S_model13ES-B_2020-03-26-10-36-33_[[:digit:]]+.csv')))

# Checks
check_hmc_diagnostics(S_model13ESB)
print(S_model13ESB,pars=c("beta","epsilon","rho","pi","psi"),digits_summary=4)
print(S_model13ESB,pars=c("cfr_A_symptomatic","cfr_B_symptomatic","cfr_C_symptomatic","cfr_D_symptomatic","cfr_C_all","cfr_D_all"),digits_summary=5)

# Save
save(
     S_model13ESB,D_S_model13ESB,
     file="posterior_samples/model13ES_2020-03-25.Rdata")




