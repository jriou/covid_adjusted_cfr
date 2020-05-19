# Setup ----
source("setup.R")
source("data/lombardy/data_management_lombardy.R")

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
  contact=c(5.13567073170732, 1.17274819632136, 0.982359525171638, 2.21715890088845, 1.29666356906914, 0.828866413937242, 0.528700773224482, 0.232116187961884, 0.0975205061876398, 1.01399087153423, 10.420788530466, 1.5084165224448, 1.46323525034693, 2.30050630727188, 1.0455742822567, 0.396916593664865, 0.276112578159939, 0.0867321859134207, 0.787940961549209, 1.39931415327149, 4.91448118586089, 2.39551550152373, 2.08291844616138, 1.67353143324194, 0.652483430981848, 0.263165822550241, 0.107498717856296, 1.53454251726848, 1.17129688889679, 2.06708280469829, 3.91165644171779, 2.74588910732349, 1.66499320847473, 1.02145416818956, 0.371633336270256, 0.112670158106901, 0.857264438638371, 1.7590640625625, 1.71686658407219, 2.62294018855816, 3.45916114790287, 1.87635185962704, 0.862205884832066, 0.523958801433231, 0.205791955532149, 0.646645383952458, 0.943424739130445, 1.62776721065554, 1.87677409215498, 2.21415705015835, 2.5920177383592, 1.10525460534109, 0.472961105423521, 0.282448363507455, 0.504954014454259, 0.438441714821823, 0.77694120330432, 1.40954408148402, 1.24556204828388, 1.35307720400585, 1.70385674931129, 0.812686154912104, 0.270111273681845, 0.305701280434649, 0.420580126969344, 0.432113761275257, 0.707170907986224, 1.04376196943771, 0.798427737704416, 1.12065725135372, 1.33035714285714, 0.322575366839763, 0.237578345845701, 0.24437789962337, 0.326505855457376, 0.396586297530862, 0.758318763302674, 0.881999483055259, 0.688988121391528, 0.596692087603768, 0.292682926829268)
)

# test model13 for lombardy because incidence_case is not known for the whole period
M_model13 = stan_model("models/model13.stan")
T_model13 = sampling(M_model13,
                     data = data_list_model13,
                     iter = 5,
                     chains = 1,
                     init=0,
                     control=list(max_treedepth=10,adapt_delta=0.8))
print(T_model13,pars=c("beta","eta","epsilon","rho","pi"))


# Create data and bash files ----
# symptomatic = 80% (A) (Bi et al) 
data_list_model13$p_psi=c(71,18) 
bashfile_rdump("model13",id="LO-A",data_list_model13,warmup=500,iter=500,adapt_delta=0.8,max_depth=10,init=0.5,timelimit=2,chains=4)

# symptomatic = 50% (B) (Diamond princess https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.10.2000180) 
data_list_model13$p_psi=c(522,114) 
bashfile_rdump("model13",id="LO-B",data_list_model13,warmup=500,iter=500,adapt_delta=0.99,max_depth=10,init=0.5,timelimit=2,chains=4)

data_list_model13$inference=1
options(mc.cores = parallel::detectCores())
S_model13ITB = sampling(M_model13,
                     data = data_list_model13,
                     iter = 1000,
                     chains = 4,
                     init=0.5,
                     control=list(max_treedepth=10,adapt_delta=0.8))
D_S_model13ITB = data_list_model13
# Copy on cluster
# system("scp /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/models/model13.stan /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/run_models/sb_model13LO* /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/run_models/data_S_model13LO* UBELIX:projects/COVID_age/model/.")

# Copy back posterior samples
system("scp  UBELIX:projects/COVID_age/model/S_model13LO-B_2020-03-25-08-45-07_*  /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/posterior_samples/.")
system("scp  UBELIX:projects/COVID_age/model/data_S_model13LO-B_2020-03-25-08-45-07.R /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/posterior_samples/.")

# Load posterior samples 
D_S_model13LOB = read_rdump("posterior_samples/data_S_model13LO-B_2020-03-25-08-45-07.R")
S_model13LOB = read_stan_csv(paste0("posterior_samples/",dir("posterior_samples",pattern = 'S_model13LO-B_2020-03-25-08-45-07_[[:digit:]]+.csv')))

# Checks
check_hmc_diagnostics(S_model13LOB)
print(S_model13LOB,pars=c("beta","epsilon","rho","pi","psi"),digits_summary=4)
print(S_model13LOB,pars=c("cfr_A_symptomatic","cfr_B_symptomatic","cfr_C_symptomatic","cfr_D_symptomatic","cfr_C_all","cfr_D_all"),digits_summary=5)

# Save
save(
     S_model13LOB,D_S_model13LOB,
     file="posterior_samples/model13LO_2020-03-25.Rdata")




