functions {
  real switch_zero(real t, real t1) {
    return(1/(1+exp(5*(t-t1))));
  }
  real switch_eta(real t, real t1, real eta) {
    return(eta+(1-eta)/(1+exp(3*(t-t1))));
  } 
  real[] SEIR(real t,
  real[] y,
  real[] theta,
  real[] x_r,
  int[] x_i
  ) {
    int K = x_i[1];
    real tmax2 = x_r[1];
    real tswitch = x_r[2];
    real dydt[(7*K)]; //SEIR, then C (cumulative incidence) and M (mortality)
    real nI; // total infectious
    real ntot;
    
    real beta; // transmission rate
    real eta; // reduction in transmission rate after quarantine
    real epsilon[K]; // mortality per age group
    real tau; // incubation period
    real mu; // infectious period
    real delta; // time from isolation to death
    real p_tmax; // tau is 0 after tmax (stop "recruiting")
    real p_tswitch;
    
    // estimated parameters
    beta = theta[1];
    eta = theta[2];
    delta = theta[3];
    epsilon = theta[4:(4+K-1)];
    
    // fixed parameters
    tau = 1.0/x_r[3];
    mu = 1.0/x_r[4];
    
    //Total number of infectious people
    nI = sum(y[(2*K+1):(3*K)]);
    ntot = sum(y[1:(4*K)]);
    p_tmax = switch_zero(t,tmax2);
    p_tswitch = switch_eta(t,tswitch,eta);
    
    for (k in 1:K) {
      // S
      dydt[k] = -beta * p_tswitch * y[k] * nI/ntot; 
      // E
      dydt[K+k] = beta * p_tswitch * y[k] * nI/ntot - tau * p_tmax * y[K+k];
      // I
      dydt[2*K+k] = tau* p_tmax * y[K+k] - mu * y[2*K+k];
      // J
      dydt[3*K+k] = mu * y[2*K+k] - delta * y[3*K+k];
      // R
      dydt[4*K+k] = delta * (1.0 - epsilon[k]) * y[3*K+k];
      // C
      dydt[5*K+k] = tau* p_tmax * y[K+k];
      // M
      dydt[6*K+k] = delta * epsilon[k] * y[3*K+k];
    }
    return(dydt);
  }
}

data {
  int K; //number of age classes
  vector[K] age_dist;
  int pop_t; //total population
  int tmax; //number of days between date of first "infection" (late November) and 11 February (when data is collected)
  real tmax2;
  real tswitch;
  
  //Data to fit
  int D; //number of days with reported incidence
  int incidence_cases[D]; //overal incidence for W weeks
  int incidence_deaths[D]; //overal incidence for W weeks
  int agedistr_cases[K]; //number of cases at tmax for the K age classes
  int agedistr_deaths[K]; //mortality at tmax for the K age classes
  //Parameters in priors
  real p_beta;
  real p_eta[2];
  real p_delta;
  real p_pi[2];
  real p_epsilon[2];
  real p_rho[2];
  real p_phi;
  // fixed parameters
  real p_incubation;
  real p_infectious;
  real p_asymptomatic_age80plus;
  //Simulation
  real t0; //starting time
  int t_data; //time of first data
  int S;
  real ts[S]; // time bins
  
  int inference;
  int doprint;
}

transformed data {
  real x_r[4] = {tmax2,tswitch,p_incubation,p_infectious};
  int x_i[1] = {K};
}

parameters{
  real<lower=0> beta; // base transmission rate
  real<lower=0,upper=1> eta; // reduction in transmission rate after quarantine measures
  real<lower=0> delta; // time from isolation to death
  vector<lower=0,upper=1> [K] epsilon; // age-dependent mortality probability
  vector<lower=0,upper=1> [K-1] rho; // age-dependent reporting probability
  real<lower=0, upper=1> pi; // number of cases at t0
  real<lower=0> phi[2]; // variance parameters
}
transformed parameters {
  // transformed parameters
  vector[K] rho_K;
  // change of format for integrate_ode_rk45
  real theta[4+K-1]; // vector of parameters
  real init[K*7]; // initial values
  real y[S,K*7]; // raw ODE output
  vector[K] comp_S[S];
  vector[K] comp_E[S];
  vector[K] comp_I[S];
  vector[K] comp_J[S];
  vector[K] comp_R[S];
  vector[K] comp_C[S];
  vector[K] comp_M[S];
  vector[K] comp_diffC[S];
  vector[K] comp_diffM[S];
  // outcomes
  vector[D] output_incidence_cases; // overall case incidence by day
  vector[D] output_incidence_deaths; // overal mortality incidence by day 
  simplex[K] output_agedistr_cases; // final age distribution of cases
  simplex[K] output_agedistr_deaths; // final age distribution of deaths
  
  // transformed parameters
  for(i in 1:(K-1)){
    rho_K[i]=rho[i]*p_asymptomatic_age80plus;
  }
  rho_K[K]=p_asymptomatic_age80plus;
  // change of format for integrate_ode_rk45
  theta[1:3] = {beta,eta,delta};
  theta[4:(4+K-1)] = to_array_1d(epsilon);
  for(k in 1:K){
    init[k] = age_dist[k] * (1-pi);
    init[K+k] = age_dist[k] * pi;
    init[2*K+k] = 0.0;
    init[3*K+k] = 0.0;
    init[4*K+k] = 0.0;
    init[5*K+k] = 0.0;
    init[6*K+k] = 0.0;
  }
  // run ODE solver
  y = integrate_ode_rk45(
    SEIR, // ODE function
    init, // initial states
    t0, // t0
    ts, // evaluation dates (ts)
    theta, // parameters
    x_r, // real data
    x_i, // integer data
    1.0E-10, 1.0E-10, 1.0E3); // tolerances and maximum steps
    // extract and format ODE results (1.0E-9 correction to avoid negative values due to unprecise estimates of zeros as tolerance is 1.0E-10)
    for(i in 1:S) {
      comp_S[i] = (to_vector(y[i,1:K]) + 1.0E-9) * pop_t;
      comp_E[i] = (to_vector(y[i,(K+1):(2*K)]) + 1.0E-9) * pop_t;
      comp_I[i] = (to_vector(y[i,(2*K+1):(3*K)]) + 1.0E-9) * pop_t;
      comp_J[i] = (to_vector(y[i,(3*K+1):(4*K)]) + 1.0E-9) * pop_t;
      comp_R[i] = (to_vector(y[i,(4*K+1):(5*K)]) + 1.0E-9) * pop_t;
      comp_C[i] = (to_vector(y[i,(5*K+1):(6*K)]) + 1.0E-9) * pop_t;
      comp_M[i] = (to_vector(y[i,(6*K+1):(7*K)]) + 1.0E-9) * pop_t;
      comp_diffC[i] = i==1 ? comp_C[i,] : 1.0E-9*pop_t + comp_C[i,] - comp_C[i-1,]; // lagged difference of cumulative incidence
      comp_diffM[i] = i==1 ? comp_M[i,] : 1.0E-9*pop_t + comp_M[i,] - comp_M[i-1,]; // lagged difference of cumulative mortality
    }
    // compute outcomes (again, 1.0E-9 correction to avoid negative values in the lagged differences)
    for(i in t_data:tmax){
      output_incidence_cases[i-t_data+1] = sum(comp_diffC[i].*rho_K);
      output_incidence_deaths[i-t_data+1] = sum(comp_diffM[i]);
    }
    output_agedistr_cases = (comp_C[tmax,].*rho_K) ./ sum(comp_C[tmax,].*rho_K);
    output_agedistr_deaths = (comp_M[tmax,]) ./ sum(comp_M[tmax,]);
}
model {
  // priors
  beta ~ exponential(p_beta);
  eta ~ beta(p_eta[1],p_eta[2]);
  delta ~ exponential(p_delta);
  for(k in 1:K) epsilon[k] ~ beta(p_epsilon[1],p_epsilon[2]);
  for(k in 1:(K-1)) rho[k] ~ beta(p_rho[1],p_rho[2]);
  pi ~ beta(p_pi[1],p_pi[2]);
  phi ~ exponential(p_phi);
  // debug
  if(doprint==1) {
    print("beta: ",beta);
    print("eta: ",beta);
    print("epsilon: ",epsilon);
    print("rho: ",rho);
    print("pi: ",pi);
    print("y[5,]: ",y[5,]);
    print("comp_C[5,]: ",comp_C[5,]);
    print("comp_diffC[5,]: ",comp_diffC[5,]);
    print("comp_M[5,]: ",comp_M[5,]);
    print("comp_diffM[5,]: ",comp_diffM[5,]);
  }
  // likelihood
  if (inference==1) {
    for(i in 1:D) {
      target += neg_binomial_2_lpmf( incidence_cases[i] | output_incidence_cases[i], output_incidence_cases[i]/phi[1]);
      target += neg_binomial_2_lpmf( incidence_deaths[i] | output_incidence_deaths[i],output_incidence_deaths[i]/phi[2]);
    }
    target += multinomial_lpmf(agedistr_cases | output_agedistr_cases);
    target += multinomial_lpmf(agedistr_deaths | output_agedistr_deaths);
  }
}

generated quantities{
  real avg_rho = sum(age_dist .* rho_K);
  real beta2 = beta*eta;
  
  int predicted_reported_incidence_cases[S]; 
  real predicted_overall_incidence_cases[S]; 
  int predicted_overall_incidence_deaths[S];
  
  int pred_comp_reported_diffC[S,K];
  vector[K] pred_comp_overall_diffC[S];
  int pred_comp_diffM[S,K];
  
  vector[K] predicted_total_reported_cases_by_age;
  vector[K] predicted_total_overall_cases_by_age;
  vector[K] predicted_total_overall_deaths_tmax_by_age;
  vector[K] predicted_total_overall_deaths_delay_by_age;
  
  real predicted_total_reported_cases;
  real predicted_total_overall_cases;
  real predicted_total_overall_deaths_tmax;
  real predicted_total_overall_deaths_delay;
  
  vector[K] cfr_A_by_age; //cfr by age classes, no correction of underreporting, no correction of time lag
  vector[K] cfr_B_by_age; //cfr by age classes, no correction of underreporting, correction of time lag
  vector[K] cfr_C_by_age; //cfr by age classes, correction of underreporting, no correction of time lag
  vector[K] cfr_D_by_age; //cfr by age classes, correction of underreporting, correction of time lag
  real cfr_A; //cfr by age classes, no correction of underreporting, no correction of time lag
  real cfr_B; //cfr by age classes, no correction of underreporting, correction of time lag
  real cfr_C; //cfr by age classes, correction of underreporting, no correction of time lag
  real cfr_D; //cfr by age classes, correction of underreporting, correction of time lag
  
  for(i in 1:S){
    predicted_reported_incidence_cases[i] = neg_binomial_2_rng(sum(comp_diffC[i].*rho_K), sum(comp_diffC[i].*rho_K)/phi[1]);
    predicted_overall_incidence_cases[i] = predicted_reported_incidence_cases[i] / avg_rho;
    predicted_overall_incidence_deaths[i] = neg_binomial_2_rng(sum(comp_diffM[i]),sum(comp_diffM[i])/phi[2]);
  }
  for(i in 1:S) {
    pred_comp_reported_diffC[i] = predicted_reported_incidence_cases[i] == 0 ? rep_array(0,K) : multinomial_rng(output_agedistr_cases,predicted_reported_incidence_cases[i]);
    pred_comp_overall_diffC[i] = to_vector(pred_comp_reported_diffC[i]) ./ rho_K;
    pred_comp_diffM[i] = predicted_overall_incidence_deaths[i] == 0 ? rep_array(0,K) : multinomial_rng(output_agedistr_deaths,predicted_overall_incidence_deaths[i]);
  }
  for(i in 1:K) {
    predicted_total_reported_cases_by_age[i] = sum(pred_comp_reported_diffC[1:tmax,i]);
    predicted_total_overall_cases_by_age[i] = sum(pred_comp_overall_diffC[1:tmax,i]);
    predicted_total_overall_deaths_tmax_by_age[i] = sum(pred_comp_diffM[1:tmax,i]);
    predicted_total_overall_deaths_delay_by_age[i] = sum(pred_comp_diffM[1:S,i]);
  }
  predicted_total_reported_cases = sum(predicted_total_reported_cases_by_age);
  predicted_total_overall_cases = sum(predicted_total_overall_cases_by_age);
  predicted_total_overall_deaths_tmax = sum(predicted_total_overall_deaths_tmax_by_age);
  predicted_total_overall_deaths_delay = sum(predicted_total_overall_deaths_delay_by_age);
  
  cfr_A_by_age = predicted_total_overall_deaths_tmax_by_age ./ predicted_total_reported_cases_by_age;
  cfr_B_by_age = predicted_total_overall_deaths_delay_by_age ./ predicted_total_reported_cases_by_age;
  cfr_C_by_age = predicted_total_overall_deaths_tmax_by_age ./ predicted_total_overall_cases_by_age;
  cfr_D_by_age = predicted_total_overall_deaths_delay_by_age ./ predicted_total_overall_cases_by_age;
  
  cfr_A = predicted_total_overall_deaths_tmax / predicted_total_reported_cases;
  cfr_B = predicted_total_overall_deaths_delay / predicted_total_reported_cases;
  cfr_C = predicted_total_overall_deaths_tmax / predicted_total_overall_cases;
  cfr_D = predicted_total_overall_deaths_delay / predicted_total_overall_cases;
}
