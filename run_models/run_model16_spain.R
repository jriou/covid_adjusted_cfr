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

## Compute proportion of symptomatics ----
# uncertainty on symptomatic rate from a systematic review (pooled mean: 19% (prediction interval 11 – 29%))
m = 1-0.19
low = 1-0.29
high = 1-0.11
v = (mean(c(m-low,high-m))/qnorm(0.975))^2

(p_psi_alpha = m*(m/v*(1-m)-1))
(p_psi_beta = (1-m)/m*p_psi_alpha)
round(1-qbeta(c(0.5,0.025,0.975),p_psi_alpha,p_psi_beta),2)

## Compute reduced transmissibility of pre- and a-symptomatics kappa ----
# He et al: 44% (25–69%) of secondary cases were infected during the index cases’ presymptomatic stage http://nature.com/articles/s41591-020-0869-5
# Liu estimated 46% (21 – 46%) (S. Funk group) https://wellcomeopenresearch.org/articles/5-58
# Ganyani (Wallinga en Hens) 48% (32-67%) https://medrxiv.org/content/10.1101/2020.03.05.20031815v1

# porameters
gt = 5.2
q_P = 0.46
incubation = 5
tau_2 = 1/2.3
tau_1 = 1/(incubation - 1/tau_2)
psi = rbeta(1000,p_psi_alpha,p_psi_beta)

# infectious duration
mu = (1-q_P)/(gt-1/tau_1-1/tau_2)
1/mu
kappa = (q_P*tau_2*psi)/((1-q_P)*mu-(1-psi)*q_P*tau_2)
qsum(kappa)

# check gp
1/tau_1 + q_P * 1/tau_2 +(1-q_P)*(1/tau_2+1/mu)

# check q_P
(kappa/tau_2)/(kappa/tau_2 + (1-psi)*kappa/mu + psi*1/mu)

## Contact matrix ----

contact_matrix_europe = c(5.13567073170732, 1.17274819632136, 0.982359525171638, 2.21715890088845, 1.29666356906914, 0.828866413937242, 0.528700773224482, 0.232116187961884, 0.0975205061876398, 1.01399087153423, 10.420788530466, 1.5084165224448, 1.46323525034693, 2.30050630727188, 1.0455742822567, 0.396916593664865, 0.276112578159939, 0.0867321859134207, 0.787940961549209, 1.39931415327149, 4.91448118586089, 2.39551550152373, 2.08291844616138, 1.67353143324194, 0.652483430981848, 0.263165822550241, 0.107498717856296, 1.53454251726848, 1.17129688889679, 2.06708280469829, 3.91165644171779, 2.74588910732349, 1.66499320847473, 1.02145416818956, 0.371633336270256, 0.112670158106901, 0.857264438638371, 1.7590640625625, 1.71686658407219, 2.62294018855816, 3.45916114790287, 1.87635185962704, 0.862205884832066, 0.523958801433231, 0.205791955532149, 0.646645383952458, 0.943424739130445, 1.62776721065554, 1.87677409215498, 2.21415705015835, 2.5920177383592, 1.10525460534109, 0.472961105423521, 0.282448363507455, 0.504954014454259, 0.438441714821823, 0.77694120330432, 1.40954408148402, 1.24556204828388, 1.35307720400585, 1.70385674931129, 0.812686154912104, 0.270111273681845, 0.305701280434649, 0.420580126969344, 0.432113761275257, 0.707170907986224, 1.04376196943771, 0.798427737704416, 1.12065725135372, 1.33035714285714, 0.322575366839763, 0.237578345845701, 0.24437789962337, 0.326505855457376, 0.396586297530862, 0.758318763302674, 0.881999483055259, 0.688988121391528, 0.596692087603768, 0.292682926829268)

data.frame(contacts=contact_matrix_europe) %>%
  tbl_df() %>%
  mutate(age1=rep(1:9,9),age2=rep(1:9,each=9)) %>%
  ggplot() +
  geom_tile(aes(x=age2,y=age1,fill=contacts))

## Format for stan

t0 = 0
S = as.numeric(day_max-day_start)
t_data = as.numeric(day_data-day_start)
ts = t_data:S
D = as.numeric(day_max-day_data+1)
tswitch = as.numeric(day_quarantine-day_start)
data_list_model16A = list(
  # structure
  K=9,
  age_dist=age_dist,
  pop_t = pop_t,
  t0=0,
  # Controls
  t_data = t_data,
  tswitch = tswitch,
  S=S,
  ts = ts,
  inference=0,
  doprint=0,
  # Data to fit
  D = D,
  incidence_cases = incidence_cases,
  incidence_deaths = incidence_deaths,
  agedistr_cases = agedistr_cases,
  agedistr_deaths = agedistr_deaths,
  # Priors
  p_beta=1,
  p_eta=c(1,1),
  p_pi=c(1,999),
  p_psi=c(p_psi_alpha,p_psi_beta),
  p_epsilon=c(1,1),
  p_phi = 1/100,
  p_rho=c(1,1),
  p_xi=1,
  p_nu=1/5,
  # Fixed parameters
  contact=contact_matrix_europe,
  p_q_P=q_P,
  p_incubation=incubation,
  p_preclinical=1/tau_2,
  p_generation_time=gt,
  p_children_trans=1,
  # Fixed corrections
  p_report_80plus=1,
  p_underreport_deaths=1,
  p_underreport_cases=p_underreport_cases,
  # Fixed delays
  G=60,
  p_gamma=gamma
)

# test
M_model16 = stan_model("models/model16.stan")
T_model16 = sampling(M_model16,data = data_list_model16A,iter = 5,chains = 1,init=0.5,control=list(max_treedepth=10,adapt_delta=0.8))
# print(T_model16,pars=c("beta","eta","epsilon","rho","pi"))

# Create data and bash files ----
bashfile_rdump("model16",id="SP-A",data_list_model16A,warmup=500,iter=500,adapt_delta=0.8,max_depth=10,init=0.5,timelimit=12,chains=4)

# system("scp /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/models/model16.stan /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/run_models/sb_model16SP* /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/run_models/data_S_model16SP* UBELIX:projects/COVID_age/model/.")

# Copy back posterior samples
system("scp  UBELIX:projects/COVID_age/model/S_model16SP-A_2020-05-04-16-26-01_*  /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/posterior_samples/.")
system("scp  UBELIX:projects/COVID_age/model/data_S_model16SP-A_2020-05-04-16-26-01.R /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/posterior_samples/.")

# SPad posterior samples 
D_S_model16ASP = read_rdump("posterior_samples/data_S_model16SP-A_2020-05-04-16-26-01.R")
S_model16ASP = read_stan_csv(paste0("posterior_samples/",dir("posterior_samples",pattern = 'S_model16SP-A_2020-05-04-16-26-01_[[:digit:]]+.csv')))

# Checks
check_hmc_diagnostics(S_model16ASP)
print(S_model16ASP,pars=c("beta","eta","epsilon","rho","pi","psi"),digits_summary=4)
print(S_model16ASP,pars=c("cfr_A_symptomatic","cfr_B_symptomatic","cfr_C_symptomatic","cfr_D_symptomatic","cfr_C_all","cfr_D_all"),digits_summary=5)


# PSPts
source('format_output/functions_model16.R')
plot_incidence_cases(S_model16ASP,D_S_model16ASP,start_date = day_start,end_date = day_max)
plot_total_cases(S_model16ASP,D_S_model16ASP)
plot_agedist_cases(S_model16ASP,D_S_model16ASP)
plot_incidence_deaths(S_model16ASP,D_S_model16ASP,start_date = day_start,end_date = day_max+50)
plot_total_deaths(S_model16ASP,D_S_model16ASP)
plot_agedist_deaths(S_model16ASP,D_S_model16ASP)

# Save
save(S_model16ASP,D_S_model16ASP,file="posterior_samples/posterior_samples_SP_2020-05-04.Rdata")
