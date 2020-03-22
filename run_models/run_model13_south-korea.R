# Setup ----
source("setup.R")
source("data/south-korea/data_management_south_korea.R")

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
  D = D,
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
  contact=c(3.521460976678, 0.807606635584979, 0.58017221932137, 1.45700400214389, 0.935832062414259, 0.404107489354786, 0.205352168108353, 0.0486140785410915, 0.0363597064909721, 1.13965890032514, 13.3263555794776, 1.37335196151099, 1.31217092938662, 2.01921605371974, 0.727844341372875, 0.157232043396022, 0.0541469767195827, 0.0404979018420231, 0.524125317641722, 1.77138293677653, 5.37342282955803, 2.6119530994001, 2.12498433097382, 1.49061176889747, 0.25082030298294, 0.0475140570215849, 0.0355369724027551, 1.73548563162998, 1.37863424962868, 2.33028799573755, 5.01693062997484, 3.00626105376023, 1.51976579426414, 0.479233855478861, 0.0777661367602361, 0.0581632727060529, 1.06279703666294, 2.49460796812014, 1.96780746460473, 3.08801652748333, 4.3284390359967, 1.70659785329017, 0.357728965065922, 0.0984010009252109, 0.0735966127391337, 1.01852479190416, 2.2415190271648, 2.27873878669667, 2.47625738362958, 2.98171697267116, 3.09444915679377, 0.606931880875396, 0.117623609620315, 0.0879736909666683, 0.6337638954766, 0.609156467872957, 0.857177966163359, 1.51627761423349, 1.18954954052422, 1.10638728069481, 1.23464881046495, 0.159781597736603, 0.119504723131815, 0.452002800912675, 0.708904909045733, 0.347302041274133, 0.710632862584685, 0.964118814556132, 0.721873750011454, 0.624025387196449, 0.426281556119421, 0.318826824001527, 0.0432872868178882, 0.0530447730927251, 0.0501000809302181, 0.0959421918690952, 0.158160905518649, 0.187119413095164, 0.200487698765865, 0.426281556119421, 0.426281556119421)
  )


# test
M_model13 = stan_model("models/model13.stan")
T_model13 = sampling(M_model13,
                     data = data_list_model13,
                     iter = 5,
                     chains = 1,
                     init=0,
                     control=list(max_treedepth=10,adapt_delta=0.8))
print(T_model13,pars=c("beta","eta","epsilon","rho_K","pi"))


# Create data and bash files ----
# symptomatic = 80% (A) (Bi et al) 
data_list_model13$p_psi=c(71,18) 
bashfile_rdump("model13",id="SK-A",data_list_model13,warmup=500,iter=500,adapt_delta=0.9,max_depth=11,init=0.5,timelimit=24,chains=4)

# symptomatic = 50% (B) (Diamond princess https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.10.2000180) 
data_list_model13$p_psi=c(522,114) 
bashfile_rdump("model13",id="SK-B",data_list_model13,warmup=500,iter=500,adapt_delta=0.9,max_depth=11,init=0.5,timelimit=24,chains=4)

# Copy on cluster
# system("scp /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/models/model13.stan /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/run_models/sb_model13* /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/run_models/data_SIM_model13* /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/run_models/data_S_model13* UBELIX:projects/COVID_age/model/.")

# Copy back posterior samples
system("scp  UBELIX:projects/COVID_age/model/S_model13SK-A_2020-03-14-22-20-44* /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/posterior_samples/.")
system("scp  UBELIX:projects/COVID_age/model/data_S_model13SK-B_2020-03-14-22-20-45* /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/posterior_samples/.")

# Load posterior samples 
D_S_model13SKA = read_rdump("posterior_samples/data_S_model13SK-A_2020-03-14-22-20-44.R")
S_model13SKA = read_stan_csv(paste0("posterior_samples/",dir("posterior_samples",pattern = 'S_model13SK-A_2020-03-14-22-20-44_[[:digit:]]+.csv')))
D_S_model13SKB = read_rdump("posterior_samples/data_S_model13SK-B_2020-03-14-22-20-45.R")
S_model13SKB = read_stan_csv(paste0("posterior_samples/",dir("posterior_samples",pattern = 'S_model13SK-B_2020-03-14-22-20-45_[[:digit:]]+.csv')))

# Checks
check_hmc_diagnostics(S_model13SKA)
print(S_model13SKA,pars=c("beta","eta","epsilon","rho","pi","xi","nu","psi"),digits_summary=4)
print(S_model13SKA,pars=c("cfr_A_symptomatic","cfr_B_symptomatic","cfr_C_symptomatic","cfr_D_symptomatic","cfr_C_all","cfr_D_all"),digits_summary=5)
check_hmc_diagnostics(S_model13SKB)
print(S_model13SKB,pars=c("beta","eta","epsilon","rho","pi","xi","nu","psi"),digits_summary=4)
print(S_model13SKB,pars=c("cfr_A_symptomatic","cfr_B_symptomatic","cfr_C_symptomatic","cfr_D_symptomatic","cfr_C_all","cfr_D_all"),digits_summary=5)


# Save
save(
  # S_model13SKA,D_S_model13SKA,
     S_model13SKB,D_S_model13SKB,
     file="posterior_samples/model13SK_2020-03-13.Rdata")




