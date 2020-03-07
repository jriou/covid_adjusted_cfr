source("data_management.R")
source("setup.R")


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

# Copy on cluster, get cmdstan chains
S_model10 = read_stan_csv(paste0("posterior_samples/",dir("posterior_samples",pattern = 'S_model10_2020-03-02-17-48-36_[[:digit:]]+.csv')))
D_S_model10 = read_rdump("posterior_samples/data_S_model10.R")
check_hmc_diagnostics(S_model10)
print(S_model10,pars=c("beta","eta","epsilon","rho_K","pi"),digits_summary=4)
save(S_model10,D_S_model10,
  file="model10_2020-03-02.Rdata")

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
data_list_model12 = list(
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
  p_psi=c(302,319), # Diamond princess
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
M_model12 = stan_model("model/model12.stan")
T_model12 = sampling(M_model12,
                     data = data_list_model12,
                     iter = 5,
                     chains = 1,
                     init=0,
                     control=list(max_treedepth=10,adapt_delta=0.8))
print(T_model12,pars=c("beta","eta","epsilon","rho_K","pi"))

# Create data and bash files
bashfile_rdump("model12",data_list_model12,warmup=500,iter=500,adapt_delta=0.8,max_depth=10,init=0.5,timelimit=96,chains=4,priorpredcheck=FALSE)

# Copy on cluster
system("scp /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/model12.stan /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/sb_model12.sh /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/data_SIM_model12.R /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/data_S_model12.R UBELIX:projects/COVID_age/model/.")

# Run on cluster
# >> R
# >> library(rstan)
SIM_model12 = read_stan_csv(dir(".",pattern = 'SIM_model12_2020-03-02-18-45-16_[[:digit:]]+.csv'))
D_SIM_model12 = read_rdump("data_SIM_model12.R")
check_hmc_diagnostics(SIM_model12)
print(SIM_model12,pars=c("beta","eta","epsilon","rho_K","pi"),digits_summary=4)
S_model12 = read_stan_csv(dir(".",pattern = 'S_model12_2020-03-02-18-45-16_[[:digit:]]+.csv'))
D_S_model12 = read_rdump("data_S_model12.R")
check_hmc_diagnostics(S_model12)
print(S_model12,pars=c("beta","eta","epsilon","rho_K","pi"),digits_summary=4)
save(
  # SIM_model12,D_SIM_model12,
  S_model12,D_S_model12,
  file="model12_2020-03-02.Rdata")

# Retrieve from cluster
system("scp  UBELIX:projects/COVID_age/model/model12_2020-03-02.Rdata /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/.")

# Load
load("model/model12_2020-03-02.Rdata")

# Prior predictive check
check_hmc_diagnostics(SIM_model12)
print(SIM_model12,pars=c("beta","eta","epsilon","rho_K","pi","phi"),digits_summary=4)
plot_incidence_cases(SIM_model12,D_SIM_model12)
plot_incidence_deaths(SIM_model12,D_SIM_model12)
pairs(SIM_model12,pars=c("beta","eta","output_incidence_cases[50]"))

# Posterior predictive check
check_hmc_diagnostics(S_model12)
print(S_model12,pars=c("beta","eta","epsilon","rho_K","pi","phi"),digits_summary=4)
plot_incidence_cases(S_model12,D_S_model12)
plot_incidence_deaths(S_model12,D_S_model12,show.after.tmax = TRUE)

print(S_model12,pars=c("cfr_A_by_age","cfr_B_by_age","cfr_C_by_age","cfr_D_by_age"),digits_summary = 4)
print(S_model12,pars=c("cfr_A","cfr_B","cfr_C","cfr_D"),digits_summary = 4)

print(S_model12,pars=c("predicted_total_reported_cases","predicted_total_overall_cases"),digits_summary = 4)

D_S_model12$contact

g1 = plot_incidence_cases(S_model12,D_S_model12,col1=col.cases,col2=col.cases.adj,show.asympto=TRUE) +
  # annotate("label",x=c(day_tswitch,day_tmax),y=c(13000,13000),label=c("a","b")) +
  coord_cartesian(xlim=c(as.Date("2019-12-30"),as.Date("2020-03-31")),ylim=c(0,13000*1.03)) +
  scale_y_continuous(breaks=c(0,5000,10000),labels=c("0","5K","10K"))
g1
g3 = plot_incidence_deaths(S_model12,D_S_model12,col1=col.deaths,col2=col.deaths.adj,show.after.tmax = TRUE)+
  # annotate("label",x=c(day_tswitch,day_tmax),y=c(125,125),label=c("a","b")) +
  coord_cartesian(xlim=c(as.Date("2019-12-30"),as.Date("2020-03-31")),ylim=c(0,125*1.03))
g3
g4 = plot_agedist_cases(S_model12,D_S_model12,col1=col.cases,col2=col.cases.adj,type=2)
g4
g5 = plot_agedist_deaths(S_model12,D_S_model12,col1=col.deaths,col2=col.deaths.adj)
g5






# Model 12  ----------------------------------------------------------

# with data from Bi et al https://www.medrxiv.org/content/10.1101/2020.03.03.20028423v1.full.pdf
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
data_list_model12 = list(
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
  p_incubation=5.95, # Bi et al
  p_infectious=2.4, # Bi et al
  p_psi=c(71,18), # Bi et al
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
M_model12 = stan_model("model/model12.stan")
T_model12 = sampling(M_model12,
                     data = data_list_model12,
                     iter = 5,
                     chains = 1,
                     init=0,
                     control=list(max_treedepth=10,adapt_delta=0.8))
print(T_model12,pars=c("beta","eta","epsilon","rho_K","pi"))

# Create data and bash files
bashfile_rdump("model12",data_list_model12,warmup=500,iter=500,adapt_delta=0.8,max_depth=10,init=0.5,timelimit=96,chains=4,priorpredcheck=FALSE)

# Copy on cluster
system("scp /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/model12.stan /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/sb_model12.sh /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/data_SIM_model12.R /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/data_S_model12.R UBELIX:projects/COVID_age/model/.")

# Run on cluster
# >> R
# >> library(rstan)
SIM_model12 = read_stan_csv(dir(".",pattern = 'SIM_model12_2020-03-02-18-45-16_[[:digit:]]+.csv'))
D_SIM_model12 = read_rdump("data_SIM_model12.R")
check_hmc_diagnostics(SIM_model12)
print(SIM_model12,pars=c("beta","eta","epsilon","rho_K","pi"),digits_summary=4)
S_model12 = read_stan_csv(dir(".",pattern = 'S_model12_2020-03-02-18-45-16_[[:digit:]]+.csv'))
D_S_model12 = read_rdump("data_S_model12.R")
check_hmc_diagnostics(S_model12)
print(S_model12,pars=c("beta","eta","epsilon","rho_K","pi"),digits_summary=4)
save(
  # SIM_model12,D_SIM_model12,
  S_model12,D_S_model12,
  file="model12_2020-03-02.Rdata")

# Retrieve from cluster
system("scp  UBELIX:projects/COVID_age/model/model12_2020-03-02.Rdata /home/julien/Dropbox/Unibe/covid-19/COVID_age/model/.")

# Load
load("model/model12_2020-03-02.Rdata")

# Prior predictive check
check_hmc_diagnostics(SIM_model12)
print(SIM_model12,pars=c("beta","eta","epsilon","rho_K","pi","phi"),digits_summary=4)
plot_incidence_cases(SIM_model12,D_SIM_model12)
plot_incidence_deaths(SIM_model12,D_SIM_model12)
pairs(SIM_model12,pars=c("beta","eta","output_incidence_cases[50]"))

# Posterior predictive check
check_hmc_diagnostics(S_model12)
print(S_model12,pars=c("beta","eta","epsilon","rho_K","pi","phi"),digits_summary=4)
plot_incidence_cases(S_model12,D_S_model12)
plot_incidence_deaths(S_model12,D_S_model12,show.after.tmax = TRUE)

print(S_model12,pars=c("cfr_A_by_age","cfr_B_by_age","cfr_C_by_age","cfr_D_by_age"),digits_summary = 4)
print(S_model12,pars=c("cfr_A","cfr_B","cfr_C","cfr_D"),digits_summary = 4)

print(S_model12,pars=c("predicted_total_reported_cases","predicted_total_overall_cases"),digits_summary = 4)

D_S_model12$contact

g1 = plot_incidence_cases(S_model12,D_S_model12,col1=col.cases,col2=col.cases.adj,show.asympto=TRUE) +
  # annotate("label",x=c(day_tswitch,day_tmax),y=c(13000,13000),label=c("a","b")) +
  coord_cartesian(xlim=c(as.Date("2019-12-30"),as.Date("2020-03-31")),ylim=c(0,13000*1.03)) +
  scale_y_continuous(breaks=c(0,5000,10000),labels=c("0","5K","10K"))
g1
g3 = plot_incidence_deaths(S_model12,D_S_model12,col1=col.deaths,col2=col.deaths.adj,show.after.tmax = TRUE)+
  # annotate("label",x=c(day_tswitch,day_tmax),y=c(125,125),label=c("a","b")) +
  coord_cartesian(xlim=c(as.Date("2019-12-30"),as.Date("2020-03-31")),ylim=c(0,125*1.03))
g3
g4 = plot_agedist_cases(S_model12,D_S_model12,col1=col.cases,col2=col.cases.adj,type=2)
g4
g5 = plot_agedist_deaths(S_model12,D_S_model12,col1=col.deaths,col2=col.deaths.adj)
g5

