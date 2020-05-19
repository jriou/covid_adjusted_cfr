# Setup ----
source("setup.R")
source("data/china/data_management_china.R")

#### Model 16 A with data up to 11 February -----------------------------------------------

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
round(qbeta(c(0.5,0.025,0.975),p_psi_alpha,p_psi_beta),2)

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

## Contact matrix in China ----

contact_matrix_china = c(0.810810810810811, 0.0908559523547268, 0.372736406439194, 1.27360250772, 0.200569529052988, 0.375083342749019, 0.60252680195839, 0.0934189610338407, 0.0225225225225225, 0.0904095466183592, 2.4392523364486, 0.140093983348316, 0.706545801082683, 0.942937990573664, 0.27920963239528, 0.326366336169345, 0.196893495540358, 0.106045179398683, 0.289504965045504, 0.109348487445688, 1.76086956521739, 0.923069180041088, 0.93772012267962, 0.965186137047983, 0.274120168579709, 0.116564256844925, 0.0773400190233669, 0.91820215964926, 0.511898453884688, 0.85680985412458, 2.70542635658915, 1.41323192857305, 0.993399938008648, 0.719603621821669, 0.146103509716984, 0.07633130138862, 0.13119227828341, 0.619819944222649, 0.789700390093264, 1.28218991206025, 2.17699115044248, 1.1199461877854, 0.514253349451317, 0.496466649026704, 0.101504389707241, 0.259078294801222, 0.193808465356441, 0.858341528544101, 0.951750199084178, 1.18265232149625, 2.31730769230769, 0.977037933291252, 0.606164987575222, 0.4393566902894, 0.552747314447092, 0.300880970126328, 0.323770338070664, 0.915670466885606, 0.721247101248993, 1.29765260904839, 2.76146788990826, 0.959867553314515, 0.340125585278128, 0.121809161946423, 0.25799743320884, 0.1956843612527, 0.264241585561661, 0.989672909331423, 1.14428055461657, 1.36428769674242, 1.96363636363636, 1.0266513139522, 0.0447824174066075, 0.211894911958445, 0.197988778289041, 0.210517772531686, 0.308554588199316, 1.26474943927563, 0.737190168823191, 1.56555579008225, 2.0625)

data.frame(contacts=contact_matrix_china) %>%
  tbl_df() %>%
  mutate(age1=rep(1:9,9),age2=rep(1:9,each=9)) %>%
  ggplot() +
  geom_tile(aes(x=age2,y=age1,fill=contacts))


## Format for stan ----
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
  agedistr_cases = cases_tmax,
  agedistr_deaths = mort_tmax,
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
  contact=contact_matrix_china,
  p_q_P=q_P,
  p_incubation=incubation,
  p_preclinical=1/tau_2,
  p_generation_time=gt,
  p_children_trans=1,
  # Fixed corrections
  p_report_80plus=1,
  p_underreport_deaths=1,
  p_underreport_cases=1,
  # Fixed delays
  G=60,
  p_gamma=gamma
)

# test
M_model16 = stan_model("models/model16.stan")
T_model16 = sampling(M_model16,data = data_list_model16A,iter = 5,chains = 1,init=0.5,control=list(max_treedepth=10,adapt_delta=0.8))
# print(T_model16,pars=c("beta","eta","epsilon","rho","pi"))


# Create data and bash files ----
bashfile_rdump("model16",id="A",data_list_model16A,warmup=500,iter=500,adapt_delta=0.8,max_depth=10,init=0.5,timelimit=12,chains=4)

# Copy on cluster
# system("scp /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/models/model16.stan /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/run_models/sb_model16* /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/run_models/data_SIM_model16A* /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/run_models/data_S_model16* UBELIX:projects/COVID_age/model/.")
# module load Boost/1.66.0-foss-2018a

# Copy back posterior samples
system("scp UBELIX:projects/COVID_age/model/S_model16A_2020-05-04-11-01-46* /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/posterior_samples/.")
system("scp UBELIX:projects/COVID_age/model/data_S_model16A_2020-05-04-11-01-46.R /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/posterior_samples/.")

# Load posterior samples 
D_S_model16A = read_rdump("posterior_samples/data_S_model16A_2020-05-04-11-01-46.R")
S_model16A = read_stan_csv(paste0("posterior_samples/",dir("posterior_samples",pattern = 'S_model16A_2020-05-04-11-01-46_[[:digit:]]+.csv')))

# Checks
check_hmc_diagnostics(S_model16A)
print(S_model16A,pars=c("beta","eta","epsilon","rho","pi","xi","nu","psi","phi","mu","kappa"),digits_summary=4)
print(S_model16A,pars=c("cfr_A_symptomatic","cfr_B_symptomatic","cfr_C_symptomatic","cfr_D_symptomatic","cfr_C_all","cfr_D_all"),digits_summary=5)

# Plots
source('format_output/functions_model16.R')
plot_incidence_cases(S_model16A,D_S_model16A,start_date = day_start,end_date = day_max)
plot_total_cases(S_model16A,D_S_model16A)
plot_agedist_cases(S_model16A,D_S_model16A)
plot_incidence_deaths(S_model16A,D_S_model16A,start_date = day_start,end_date = day_max+50)
plot_total_deaths(S_model16A,D_S_model16A)
plot_agedist_deaths(S_model16A,D_S_model16A)

# Save
save(S_model16A,D_S_model16A,file="posterior_samples/posterior_samples_HU_2020-05-04.Rdata")


#### Model 16 B:  model16A + later correction of deaths and cases -----------------------------------------------

data_list_model16B = data_list_model16A

## Later correction of cases and deaths ----
# While only 41092 cases were reported by date of disease onset until Feb 11 in China CDC Weekly, representing the first wave of the epidemic
# in Hubei, an additional 26667 cases were reported 

data_list_model16B$p_underreport_cases = 41092/67759
data_list_model16B$p_underreport_deaths = 3222/4512

# Create data and bash files ----
bashfile_rdump("model16",id="B",data_list_model16B,warmup=500,iter=500,adapt_delta=0.8,max_depth=10,init=0.5,timelimit=12,chains=4)

# Copy on cluster
# system("scp /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/run_models/sb_model16B.sh /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/run_models/data_S_model16B_* UBELIX:projects/COVID_age/model/.")
# module load Boost/1.66.0-foss-2018a

# Copy back posterior samples
system("scp  UBELIX:projects/COVID_age/model/S_model16B_2020-05-04-16-57-04* /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/posterior_samples/.")
system("scp  UBELIX:projects/COVID_age/model/data_S_model16B_2020-05-04-16-57-04* /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/posterior_samples/.")

# Load posterior samples 
D_S_model16B = read_rdump("posterior_samples/data_S_model16B_2020-05-04-16-57-04.R")
S_model16B = read_stan_csv(paste0("posterior_samples/",dir("posterior_samples",pattern = 'S_model16B_2020-05-04-16-57-04_[[:digit:]]+.csv')))

# Checks
check_hmc_diagnostics(S_model16B)
print(S_model16B,pars=c("beta","eta","epsilon","rho","pi","xi","nu","psi","mu","kappa"),digits_summary=4)
print(S_model16B,pars=c("cfr_A_symptomatic","cfr_B_symptomatic","cfr_C_symptomatic","cfr_D_symptomatic","cfr_C_all","cfr_D_all"),digits_summary=5)


print(S_model16B,pars=c("epsilon"),digits_summary=4)
print(S_model16A,pars=c("epsilon"),digits_summary=4)

# Plots
source('format_output/functions_model16.R')
plot_incidence_cases(S_model16B,D_S_model16B,start_date = day_start,end_date = day_max)
plot_total_cases(S_model16B,D_S_model16B)
plot_agedist_cases(S_model16B,D_S_model16B)
plot_incidence_deaths(S_model16B,D_S_model16B,start_date = day_start,end_date = day_max+50)
plot_total_deaths(S_model16B,D_S_model16B)
plot_agedist_deaths(S_model16B,D_S_model16B)





#### Model 16 C lower susceptibility of children -----------------------------------------------

data_list_model16C = data_list_model16A

data_list_model16C$p_children_trans = 0.1

# Create data and bash files ----
bashfile_rdump("model16",id="C",data_list_model16C,warmup=500,iter=500,adapt_delta=0.8,max_depth=10,init=0.5,timelimit=12,chains=4)

# Copy on cluster
# system("scp /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/run_models/sb_model16C.sh /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/run_models/data_S_model16C_* UBELIX:projects/COVID_age/model/.")
# module load Boost/1.66.0-foss-2018a

# Copy back posterior samples
system("scp UBELIX:projects/COVID_age/model/S_model16C_2020-05-07-10-19-49* /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/posterior_samples/.")
system("scp UBELIX:projects/COVID_age/model/data_S_model16C_2020-05-07-10-19-49* /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/posterior_samples/.")

# Load posterior samples 
D_S_model16C = read_rdump("posterior_samples/data_S_model16C_2020-05-07-10-19-49.R")
S_model16C = read_stan_csv(paste0("posterior_samples/",dir("posterior_samples",pattern = 'S_model16C_2020-05-07-10-19-49_[[:digit:]]+.csv')))

# Checks
check_hmc_diagnostics(S_model16C)
print(S_model16C,pars=c("beta","eta","epsilon","rho","pi","xi","nu","psi","phi"),digits_summary=4)
print(S_model16C,pars=c("cfr_A_symptomatic","cfr_B_symptomatic","cfr_C_symptomatic","cfr_D_symptomatic","cfr_C_all","cfr_D_all"),digits_summary=5)

# Plots
source('format_output/functions_model16.R')
plot_incidence_cases(S_model16C,D_S_model16C,start_date = day_start,end_date = day_max)
plot_total_cases(S_model16C,D_S_model16C)
plot_agedist_cases(S_model16C,D_S_model16C)
plot_incidence_deaths(S_model16C,D_S_model16C,start_date = day_start,end_date = day_max+50)
plot_total_deaths(S_model16C,D_S_model16C)
plot_agedist_deaths(S_model16C,D_S_model16C)







#### Model 16 D lower susceptibility of children -----------------------------------------------

data_list_model16D = data_list_model16A

data_list_model16D$p_children_trans = 0.5

# Create data and bash files ----
bashfile_rdump("model16",id="D",data_list_model16D,warmup=500,iter=500,adapt_delta=0.8,max_depth=10,init=0.5,timelimit=12,chains=4)

# Copy on cluster
# system("scp /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/run_models/sb_model16D.sh /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/run_models/data_S_model16D_* UBELIX:projects/COVID_age/model/.")
# module load Boost/1.66.0-foss-2018a

# Copy back posterior samples
system("scp UBELIX:projects/COVID_age/model/S_model16D_2020-05-04-16-59-33* /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/posterior_samples/.")
system("scp UBELIX:projects/COVID_age/model/data_S_model16D_2020-05-04-16-59-33* /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/posterior_samples/.")

# Load posterior samples 
D_S_model16D = read_rdump("posterior_samples/data_S_model16D_2020-05-04-16-59-33.R")
S_model16D = read_stan_csv(paste0("posterior_samples/",dir("posterior_samples",pattern = 'S_model16D_2020-05-04-16-59-33_[[:digit:]]+.csv')))

# Checks
check_hmc_diagnostics(S_model16D)
print(S_model16D,pars=c("beta","eta","epsilon","rho","pi","xi","nu","psi","phi"),digits_summary=4)
print(S_model16D,pars=c("cfr_A_symptomatic","cfr_B_symptomatic","cfr_C_symptomatic","cfr_D_symptomatic","cfr_C_all","cfr_D_all"),digits_summary=5)

# Plots
source('format_output/functions_model16.R')
plot_incidence_cases(S_model16D,D_S_model16D,start_date = day_start,end_date = day_max)
plot_total_cases(S_model16D,D_S_model16D)
plot_agedist_cases(S_model16D,D_S_model16D)
plot_incidence_deaths(S_model16D,D_S_model16D,start_date = day_start,end_date = day_max+50)
plot_total_deaths(S_model16D,D_S_model16D)
plot_agedist_deaths(S_model16D,D_S_model16D)




#### Model 16 E varying rho for 80+ -----------------------------------------------

data_list_model16E = data_list_model16A

# Create data and bash files ----
for(i in c(10,20,30,40,50,60,70,80,90)) {
  data_list_model16E$p_report_80plus = i/100
  bashfile_rdump("model16",id=paste0("E-",i),data_list_model16E,warmup=500,iter=500,adapt_delta=0.8,max_depth=10,init=0.5,timelimit=12,chains=4)
}
# Copy on cluster
# system("scp /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/run_models/sb_model16E-* /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/run_models/data_S_model16E-* UBELIX:projects/COVID_age/model/.")
# module load Boost/1.66.0-foss-2018a

# Copy back posterior samples
system("scp UBELIX:projects/COVID_age/model/S_model16E-*_2020-05-04-16-59-54* /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/posterior_samples/.")
system("scp UBELIX:projects/COVID_age/model/data_S_model16E*2020-05-04-16-59-54* /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/posterior_samples/.")

# Load posterior samples 
for(i in c(10,20,30,40,50,60,70,80,90)) {
  D = read_rdump(paste0("posterior_samples/data_S_model16E-",i,"_2020-05-04-16-59-54.R"))
  S = read_stan_csv(paste0("posterior_samples/",dir("posterior_samples",pattern = paste0('S_model16E-',i,'_2020-05-04-16-59-54_[[:digit:]]+.csv'))))
  assign(paste0("D_S_model16E",i),D)
  assign(paste0("S_model16E",i),S)
  cat(i)
}
# Checks
check_hmc_diagnostics(S_model16E10)
print(S_model16E10,pars=c("beta","eta","epsilon","rho","pi","xi","nu","psi","phi"),digits_summary=4)
print(S_model16E10,pars=c("cfr_A_symptomatic","cfr_B_symptomatic","cfr_C_symptomatic","cfr_D_symptomatic","cfr_C_all","cfr_D_all"),digits_summary=5)

# Plots
source('format_output/functions_model16.R')

mort_by_rho_80p = NULL
for(i in c(10,20,30,40,50,60,70,80,90)) {
  mort_by_rho_80p = bind_rows(mort_by_rho_80p,
  summary(get(paste0("S_model16E",i)),pars=c("cfr_A_symptomatic","cfr_B_symptomatic","cfr_C_symptomatic","cfr_D_symptomatic","cfr_C_all","cfr_D_all"))[[1]] %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    tbl_df() %>%
    mutate(rho_80p=i))
}
mort_by_rho_80p = bind_rows(mort_by_rho_80p,
                            summary(S_model16A,pars=c("cfr_A_symptomatic","cfr_B_symptomatic","cfr_C_symptomatic","cfr_D_symptomatic","cfr_C_all","cfr_D_all"))[[1]] %>%
                              as.data.frame() %>%
                              rownames_to_column() %>%
                              tbl_df() %>%
                              mutate(rho_80p=100))
filter(mort_by_rho_80p,rowname %in% c("cfr_A_symptomatic","cfr_D_symptomatic","cfr_D_all")) %>%
  mutate(type=factor(rowname,levels=c("cfr_A_symptomatic","cfr_D_symptomatic","cfr_D_all"),
                     labels=c("Crude","Adjusted among symptomatics","Adjusted among infected"))) %>%
ggplot() +
  geom_ribbon(aes(x=rho_80p,ymin=`2.5%`,ymax=`97.5%`,fill=type),alpha=.5) +
  geom_line(aes(x=rho_80p,y=`50%`,colour=type)) +
  scale_x_continuous(breaks=c(10,20,30,40,50,60,70,80,90,100),labels=function(x) paste0(x,"%")) +
  scale_y_continuous(labels=scales::percent) +
  labs(x="Ascertainment of 80+ year old",y="Proportion",fill="Mortality",colour="Mortality")
ggsave("format_output/figures/sst_mortality_by_rho_80p.pdf",width=7,height=5)




#### Model 16 F with increasing times -----------------------------------------------

for(i in c(5,10,15,20,25,30)) {
  data_list_model16F = data_list_model16A
  data_list_model16F$t0 = 0
  data_list_model16F$S = as.numeric(day_max-i-day_start)
  data_list_model16F$t_data = as.numeric(day_data-day_start)
  data_list_model16F$ts = t_data:data_list_model16F$S
  data_list_model16F$D = as.numeric(day_max-i-day_data+1)
  data_list_model16F$tswitch = as.numeric(day_quarantine-day_start)
  data_list_model16F$incidence_cases = incidence_cases[1:data_list_model16F$D]
  data_list_model16F$incidence_deaths = incidence_deaths[1:data_list_model16F$D]
  bashfile_rdump("model16",id=paste0("F-",42-i),data_list_model16F,warmup=500,iter=500,adapt_delta=0.8,max_depth=10,init=0.5,timelimit=12,chains=4)
}

# Copy on cluster
# system("scp /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/run_models/sb_model16F-* /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/run_models/data_S_model16F-* UBELIX:projects/COVID_age/model/.")
# # module load Boost/1.66.0-foss-2018a

# Copy back posterior samples

system("scp UBELIX:projects/COVID_age/model/S_model16F*_2020-05-05-16-16-40* /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/posterior_samples/.")
system("scp UBELIX:projects/COVID_age/model/data_S_model16F*_2020-05-05-16-16-40* /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/posterior_samples/.")

# Load posterior samples 
for(i in 42-c(5,10,15,20,25,30)) {
  D = read_rdump(paste0("posterior_samples/data_S_model16F-",i,"_2020-05-05-16-16-40.R"))
  S = read_stan_csv(paste0("posterior_samples/",dir("posterior_samples",pattern = paste0('S_model16F-',i,'_2020-05-05-16-16-40_[[:digit:]]+.csv'))))
  assign(paste0("D_S_model16F",i),D)
  assign(paste0("S_model16F",i),S)
  cat(i)
}

# Checks
check_hmc_diagnostics(S_model16F12)
print(S_model16F12,pars=c("beta","eta","epsilon","rho","pi","xi","nu","psi","phi"),digits_summary=4)
print(S_model16F12,pars=c("cfr_A_symptomatic","cfr_B_symptomatic","cfr_C_symptomatic","cfr_D_symptomatic","cfr_C_all","cfr_D_all"),digits_summary=5)

# Plots
source('format_output/functions_model16.R')
plot_incidence_cases(S_model16F12,D_S_model16F12,start_date = day_start,end_date = day_max)
plot_total_cases(S_model16F12,D_S_model16F12)
plot_agedist_cases(S_model16F12,D_S_model16F12)
plot_incidence_deaths(S_model16F12,D_S_model16F12,start_date = day_start,end_date = day_max+50)
plot_total_deaths(S_model16F12,D_S_model16F12)
plot_agedist_deaths(S_model16F12,D_S_model16F12)


most_by_date = NULL
for(i in  42-c(5,10,15,20,25,30)) {
  most_by_date = bind_rows(most_by_date,
                              summary(get(paste0("S_model16F",i)),pars=c("cfr_A_symptomatic","cfr_B_symptomatic","cfr_C_symptomatic","cfr_D_symptomatic","cfr_C_all","cfr_D_all"))[[1]] %>%
                                as.data.frame() %>%
                                rownames_to_column() %>%
                                tbl_df() %>%
                                mutate(date=i))
}
most_by_date = bind_rows(most_by_date,
                            summary(S_model16A,pars=c("cfr_A_symptomatic","cfr_B_symptomatic","cfr_C_symptomatic","cfr_D_symptomatic","cfr_C_all","cfr_D_all"))[[1]] %>%
                              as.data.frame() %>%
                              rownames_to_column() %>%
                              tbl_df() %>%
                              mutate(date=42))
filter(most_by_date,rowname %in% c("cfr_A_symptomatic","cfr_D_all")) %>%
  mutate(type=factor(rowname,levels=c("cfr_A_symptomatic","cfr_D_all"),
                     labels=c("Crude","Adjusted among infected"))) %>%
  mutate(date2=day_start+date) %>%
  ggplot() +
  geom_ribbon(aes(x=date2,ymin=`2.5%`,ymax=`97.5%`,fill=type),alpha=.5) +
  geom_line(aes(x=date2,y=`50%`,colour=type)) +
  scale_y_continuous(labels=scales::percent) +
  labs(x="Date",y="Proportion",fill="Mortality",colour="Mortality")
ggsave("format_output/figures/sst_mortality_by_date.pdf",width=7,height=5)




#### Model 16 G with lower participation of presymptomatics -----------------------------------------------


data_list_model16G = data_list_model16A

data_list_model16G$p_q_P=0.30

# Create data and bash files ----
bashfile_rdump("model16",id="G",data_list_model16G,warmup=500,iter=500,adapt_delta=0.8,max_depth=10,init=0.5,timelimit=12,chains=4)

# Copy on cluster
# system("scp /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/run_models/sb_model16G.sh /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/run_models/data_SIM_model16G* /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/run_models/data_S_model16G* UBELIX:projects/COVID_age/model/.")
# module load Boost/1.66.0-foss-2018a

# Copy back posterior samples
system("scp UBELIX:projects/COVID_age/model/S_model16G_2020-05-04-17-01-43* /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/posterior_samples/.")
system("scp UBELIX:projects/COVID_age/model/data_S_model16G_2020-05-04-17-01-43* /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/posterior_samples/.")

# Load posterior samples 
D_S_model16G = read_rdump("posterior_samples/data_S_model16G_2020-05-04-17-01-43.R")
S_model16G = read_stan_csv(paste0("posterior_samples/",dir("posterior_samples",pattern = 'S_model16G_2020-05-04-17-01-43_[[:digit:]]+.csv')))

# Checks
check_hmc_diagnostics(S_model16G)
print(S_model16G,pars=c("beta","eta","epsilon","rho","pi","xi","nu","psi","phi","mu","kappa"),digits_summary=4)
print(S_model16G,pars=c("cfr_A_symptomatic","cfr_B_symptomatic","cfr_C_symptomatic","cfr_D_symptomatic","cfr_C_all","cfr_D_all"),digits_summary=5)

# Plots
source('format_output/functions_model16.R')
plot_incidence_cases(S_model16G,D_S_model16G,start_date = day_start,end_date = day_max)
plot_total_cases(S_model16G,D_S_model16G)
plot_agedist_cases(S_model16G,D_S_model16G)
plot_incidence_deaths(S_model16G,D_S_model16G,start_date = day_start,end_date = day_max+50)
plot_total_deaths(S_model16G,D_S_model16G)
plot_agedist_deaths(S_model16G,D_S_model16G)


#### Model 16 H with no participation of pre and asymptomatics  -----------------------------------------------

data_list_model16H = data_list_model16A

data_list_model16H$p_q_P=0

# Create data and bash files ----
bashfile_rdump("model16",id="H",data_list_model16H,warmup=500,iter=500,adapt_delta=0.8,max_depth=10,init=0.5,timelimit=12,chains=4)

# Copy on cluster
# system("scp /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/run_models/sb_model16H.sh /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/run_models/data_S_model16H_* UBELIX:projects/COVID_age/model/.")
# module load Boost/1.66.0-foss-2018a

# Copy back posterior samples
system("scp UBELIX:projects/COVID_age/model/S_model16H_2020-05-04-17-02-52* /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/posterior_samples/.")
system("scp UBELIX:projects/COVID_age/model/data_S_model16H_2020-05-04-17-02-52* /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/posterior_samples/.")

# Load posterior samples 
D_S_model16H = read_rdump("posterior_samples/data_S_model16H_2020-05-04-17-02-52.R")
S_model16H = read_stan_csv(paste0("posterior_samples/",dir("posterior_samples",pattern = 'S_model16H_2020-05-04-17-02-52_[[:digit:]]+.csv')))

# Checks
check_hmc_diagnostics(S_model16H)
print(S_model16H,pars=c("beta","eta","epsilon","rho","pi","xi","nu","psi","phi"),digits_summary=4)
print(S_model16H,pars=c("cfr_A_symptomatic","cfr_B_symptomatic","cfr_C_symptomatic","cfr_D_symptomatic","cfr_C_all","cfr_D_all"),digits_summary=5)

# Plots
source('format_output/functions_model16.R')
plot_incidence_cases(S_model16H,D_S_model16H,start_date = day_start,end_date = day_max)
plot_total_cases(S_model16H,D_S_model16H)
plot_agedist_cases(S_model16H,D_S_model16H)
plot_incidence_deaths(S_model16H,D_S_model16H,start_date = day_start,end_date = day_max+50)
plot_total_deaths(S_model16H,D_S_model16H)
plot_agedist_deaths(S_model16H,D_S_model16H)





#### Model 16 I:  model16A + varying onset to mortality -----------------------------------------------

data_list_model16I = data_list_model16A

# Fom Linton et al
linton_mean = 15
linton_sd = 10
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
plot(gamma)
data_list_model16I$p_gamma = gamma

# Create data and bash files ----
bashfile_rdump("model16",id="I",data_list_model16I,warmup=500,iter=500,adapt_delta=0.8,max_depth=10,init=0.5,timelimit=12,chains=4)

# Copy on cluster
# system("scp /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/run_models/sb_model16I.sh /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/run_models/data_S_model16I_* UBELIX:projects/COVID_age/model/.")
# module load Ioost/1.66.0-foss-2018a

# Copy back posterior samples
system("scp  UBELIX:projects/COVID_age/model/S_model16I_2020-05-05-16-31-49* /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/posterior_samples/.")
system("scp  UBELIX:projects/COVID_age/model/data_S_model16I_2020-05-05-16-31-49* /home/julien/Dropbox/Unibe/covid-19/covid_adjusted_cfr/posterior_samples/.")

# Load posterior samples 
D_S_model16I = read_rdump("posterior_samples/data_S_model16I_2020-05-05-16-31-49.R")
S_model16I = read_stan_csv(paste0("posterior_samples/",dir("posterior_samples",pattern = 'S_model16I_2020-05-05-16-31-49_[[:digit:]]+.csv')))

# Checks
check_hmc_diagnostics(S_model16I)
print(S_model16I,pars=c("beta","eta","epsilon","rho","pi","xi","nu","psi","mu","kappa"),digits_summary=4)
print(S_model16I,pars=c("cfr_A_symptomatic","cfr_B_symptomatic","cfr_C_symptomatic","cfr_D_symptomatic","cfr_C_all","cfr_D_all"),digits_summary=5)

# Plots
source('format_output/functions_model16.R')
plot_incidence_cases(S_model16I,D_S_model16I,start_date = day_start,end_date = day_max)
plot_total_cases(S_model16I,D_S_model16I)
plot_agedist_cases(S_model16I,D_S_model16I)
plot_incidence_deaths(S_model16I,D_S_model16I,start_date = day_start,end_date = day_max+50)
plot_total_deaths(S_model16I,D_S_model16I)
plot_agedist_deaths(S_model16I,D_S_model16I)

# Save
save(S_model16B,D_S_model16B,
     S_model16C,D_S_model16C,
     S_model16D,D_S_model16D,
     S_model16E10,D_S_model16E10,
     S_model16E20,D_S_model16E20,
     S_model16E30,D_S_model16E30,
     S_model16E40,D_S_model16E40,
     S_model16E50,D_S_model16E50,
     S_model16E60,D_S_model16E60,
     S_model16E70,D_S_model16E70,
     S_model16E80,D_S_model16E80,
     S_model16E90,D_S_model16E90,
     S_model16F12,D_S_model16F12,
     S_model16F17,D_S_model16F17,
     S_model16F22,D_S_model16F22,
     S_model16F27,D_S_model16F27,
     S_model16F32,D_S_model16F32,
     S_model16F37,D_S_model16F37,
     S_model16G,D_S_model16G,
     S_model16H,D_S_model16H,
     S_model16I,D_S_model16I,
     file="posterior_samples/posterior_samples_supp_HU_2020-05-04.Rdata")

