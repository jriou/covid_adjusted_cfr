functions {
  real switch_eta(real t, real t1, real eta, real nu, real xi) {
    return(eta+(1-eta)/(1+exp(xi*(t-t1-nu))));
  }
  real[] SEIR(real t,
  real[] y,
  real[] theta,
  real[] x_r,
  int[] x_i
  ) {
    int K = x_i[1];
    
    real tswitch = x_r[1];
    real dydt[(4*K)]; //SEIAR, then C and D
    real nI; // total infectious
    real ntot;
    
    real beta; // transmission rate
    real eta; // reduction in transmission rate after quarantine
    real xi; // slope of quarantine implementation
    real nu; // shift of quarantine implementation
    real tau; // incubation period
    real mu; // infectious period
    real psi; // probability of symptoms
    real p_tswitch;
    real contact[K*K];//contact matrix, first K values, corresponds to number of contact between age class 1 and other classes, etc
    real n_by_age[K];
    real f_inf[K];
    
    real init[K*4];
    real age_dist[K];
    real pi; // number of cases at t0
    
    // estimated parameters
    beta = theta[1];
    eta = theta[2];
    xi = theta[3];
    nu = theta[4];
    pi = theta[5];
    psi = theta[6];
    
    // Initial conditions
    for(k in 1:K){
      age_dist[k] = x_r[3+K*K + k];
      init[k] = age_dist[k] * (1-pi);
      init[K+k] = age_dist[k] * pi;
      init[2*K+k] = 0.0;
      init[3*K+k] = 0.0;
    }
    
    // fixed parameters
    tau = 1.0/x_r[2];
    mu = 1.0/x_r[3];
    contact = x_r[4:(3+K*K)];
    
    //Total number of infectious people
    p_tswitch = switch_eta(t,tswitch,eta,nu,xi);
    
    //Force of infection by age classes: beta * p_tswitch * sum((number of infected people by age) / (total number of people by age) * (number of contact by age))
    for(k in 1:K){
      f_inf[k] = beta * p_tswitch * sum(to_vector(y[(2*K+1):(3*K)]) ./ to_vector(age_dist) .* to_vector(contact[(K*(k-1)+1):(k*K)])); //
    }
    //print("Force of infection",f_inf);
    
    for (k in 1:K) {
      // S
      dydt[k] = - f_inf[k] * (y[k]+init[k]); 
      // E
      dydt[K+k] = f_inf[k] * (y[k]+init[k])- tau * (y[K+k]+init[K+k]);
      // I
      dydt[2*K+k] = psi * tau * (y[K+k]+init[K+k]) - mu * (y[2*K+k]+init[2*K+k]);
      // C
      dydt[3*K+k] = psi * tau * (y[K+k]+init[K+k]);
    }
    return(dydt);
  }
}

data {
  int K; //number of age classes
  vector[K] age_dist;
  int pop_t; //total population
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
  real p_pi[2];
  real p_epsilon[2];
  real p_rho[2];
  real p_phi;
  real p_xi;
  real p_nu;
  real p_psi[2];
  // fixed parameters
  real p_incubation;
  real p_infectious;
  int G;
  real p_gamma[G];
  //Simulation
  real t0; //starting time
  int t_data; //time of first data
  int S;
  real ts[S]; // time bins
  real contact[K*K];
  
  int inference;
  int doprint;
}

transformed data {
  real x_r[3+K*K+K]; //4 parameters + K*K contact matrix parameters + K age_dist parameters
  int x_i[1] = {K};
  real init[K*4] = rep_array(0.0, K * 4); // initial values
  x_r[1] = tswitch;
  x_r[2] = p_incubation;
  x_r[3] = p_infectious;
  x_r[4:(3+K*K)] = contact;
  for(k in 1:K) {
    x_r[3 + K*K + k] = age_dist[k];
  }
}

parameters{
  real<lower=0,upper=1> beta; // base transmission rate
  real<lower=0,upper=1> eta; // reduction in transmission rate after quarantine measures
  vector<lower=0,upper=1> [K] epsilon; // age-dependent mortality probability
  vector<lower=0,upper=1> [K] rho; // age-dependent reporting probability
  real<lower=0, upper=1> pi; // number of cases at t0
  real<lower=0> phi[2]; // variance parameters
  real<lower=0,upper=1> xi_raw; // slope of quarantine implementation
  real<lower=0> nu; // shift of quarantine implementation
  real<lower=0,upper=1> psi; // proportion of symptomatics
}
transformed parameters {
  // transformed parameters
  real xi = xi_raw+0.5;
  // change of format for integrate_ode_rk45
  real theta[6]; // vector of parameters
  real y[S,K*4]; // raw ODE output
  vector[K] comp_S[S];
  vector[K] comp_E[S];
  vector[K] comp_I[S];
  vector[K] comp_D[S];
  vector[K] comp_diffD[S];
  vector[K] comp_C[S+G];
  vector[K] comp_diffC[S+G];
  vector[K] comp_diffM[S+G];
  vector[K] comp_M[S+G];
  // outcomes
  vector[D] output_incidence_cases; // overall case incidence by day
  vector[D] output_incidence_deaths; // overal mortality incidence by day 
  simplex[K] output_agedistr_cases; // final age distribution of cases
  simplex[K] output_agedistr_deaths; // final age distribution of deaths
  
  // change of format for integrate_ode_rk45
  theta[1:6] = {beta,eta,xi,nu,pi,psi};
  // run ODE solver
  y = integrate_ode_bdf(
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
      comp_S[i] = (to_vector(y[i,1:K]) + to_vector(age_dist) * (1-pi) + 1.0E-9) * pop_t;
      comp_E[i] = (to_vector(y[i,(K+1):(2*K)]) + to_vector(age_dist) * pi + 1.0E-9) * pop_t;
      comp_I[i] = (to_vector(y[i,(2*K+1):(3*K)]) + 1.0E-9) * pop_t;

      comp_C[i] = (to_vector(y[i,(3*K+1):(4*K)]) + 1.0E-9) * pop_t;
      comp_diffC[i] = i==1 ? comp_C[i,] : 1.0E-9*pop_t + comp_C[i,] - comp_C[i-1,]; // lagged difference of cumulative incidence of symptomatics
    }
    //Incidence and cumulative incidence after S
    for(g in 1:G){
      comp_C[S+g] = comp_C[S];
      comp_diffC[S+g] = rep_vector(1.0E-9,K);
    }
    //Mortality
    //set diffM and M to 0
    for(i in 1:(G+S)){
      comp_diffM[i] = rep_vector(1.0E-9,K);
    }
    //compute mortality
    for(i in 1:S) {
      for(g in 1:G) {
        comp_diffM[i+g] += comp_diffC[i] .* epsilon * p_gamma[g] ;
      }
    }
    // cumulative sum
    for(i in 1:(S+G)) {
      for(k in 1:K) {
        comp_M[i,k] = sum(comp_diffM[1:i,k]);
      }
    }
    //compute D and diffD
    for(i in 1:S){
      comp_D[i] = (1.0-psi)/psi * comp_C[i];
      comp_diffD[i] = (1.0-psi)/psi * comp_diffC[i];
    }
    // compute outcomes (again, 1.0E-9 correction to avoid negative values in the lagged differences)
    for(i in t_data:S){
      output_incidence_cases[i-t_data+1] = sum(comp_diffC[i].*rho);
      output_incidence_deaths[i-t_data+1] = sum(comp_diffM[i]);
    }
    output_agedistr_cases = (comp_C[S,].*rho) ./ sum(comp_C[S,].*rho);
    output_agedistr_deaths = (comp_M[S,]) ./ sum(comp_M[S,]);
}
model {
  // priors
  beta ~ beta(p_beta,p_beta);
  eta ~ beta(p_eta[1],p_eta[2]);
  for(k in 1:K) epsilon[k] ~ beta(p_epsilon[1],p_epsilon[2]);
  for(k in 1:K) rho[k] ~ beta(p_rho[1],p_rho[2]);
  pi ~ beta(p_pi[1],p_pi[2]);
  phi ~ exponential(p_phi);
  xi_raw ~ beta(p_xi,p_xi); 
  nu ~ exponential(p_nu);
  psi ~ beta(p_psi[1],p_psi[2]);
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
  real avg_rho = sum(age_dist .* rho);
  real beta2 = beta*eta;
  
  int predicted_reported_incidence_symptomatic_cases[S]; 
  real predicted_overall_incidence_symptomatic_cases[S]; 
  real predicted_overall_incidence_all_cases[S]; 
  int predicted_overall_incidence_deaths[S+G];
  
  int predicted_comp_reported_diffC[S,K];
  vector[K] predicted_comp_overall_diffC[S];
  vector[K] predicted_comp_overall_diffA[S];
  int predicted_comp_diffM[S+G,K];
  
  vector[K] predicted_total_reported_symptomatic_cases_by_age;
  vector[K] predicted_total_overall_symptomatic_cases_by_age;
  vector[K] predicted_total_overall_all_cases_by_age;
  vector[K] predicted_total_overall_deaths_tmax_by_age;
  vector[K] predicted_total_overall_deaths_delay_by_age;
  
  real predicted_total_reported_symptomatic_cases;
  real predicted_total_overall_symptomatic_cases;
  real predicted_total_overall_all_cases;
  real predicted_total_overall_deaths_tmax;
  real predicted_total_overall_deaths_delay;
  
  vector[K] cfr_A_symptomatic_by_age; //cfr by age classes, no correction of underreporting, no correction of time lag
  vector[K] cfr_B_symptomatic_by_age; //cfr by age classes, no correction of underreporting, correction of time lag
  vector[K] cfr_C_symptomatic_by_age; //cfr by age classes, correction of underreporting, no correction of time lag
  vector[K] cfr_D_symptomatic_by_age; //cfr by age classes, correction of underreporting, correction of time lag
  vector[K] cfr_C_all_by_age; //cfr by age classes, correction of underreporting and asymptomatics, no correction of time lag
  vector[K] cfr_D_all_by_age; //cfr by age classes, correction of underreporting and asymptomatics, correction of time lag
  real cfr_A_symptomatic; //cfr by age classes, no correction of underreporting, no correction of time lag
  real cfr_B_symptomatic; //cfr by age classes, no correction of underreporting, correction of time lag
  real cfr_C_symptomatic; //cfr by age classes, correction of underreporting, no correction of time lag
  real cfr_D_symptomatic; //cfr by age classes, correction of underreporting, correction of time lag
  real cfr_C_all; //cfr by age classes, correction of underreporting and asymptomatics, no correction of time lag
  real cfr_D_all; //cfr by age classes, correction of underreporting and asymptomatics, correction of time lag
  
  for(i in 1:S){
    predicted_reported_incidence_symptomatic_cases[i] = neg_binomial_2_rng(sum(comp_diffC[i].*rho), sum(comp_diffC[i].*rho)/phi[1]);
    predicted_overall_incidence_symptomatic_cases[i] = predicted_reported_incidence_symptomatic_cases[i] / avg_rho;
    predicted_overall_incidence_all_cases[i] = predicted_reported_incidence_symptomatic_cases[i] / avg_rho / psi;
  }
  for(i in 1:(S+G)) predicted_overall_incidence_deaths[i] = neg_binomial_2_rng(sum(comp_diffM[i]),sum(comp_diffM[i])/phi[2]);
  for(i in 1:S) {
    predicted_comp_reported_diffC[i] = predicted_reported_incidence_symptomatic_cases[i] == 0 ? rep_array(0,K) : multinomial_rng(output_agedistr_cases,predicted_reported_incidence_symptomatic_cases[i]);
    predicted_comp_overall_diffC[i] = to_vector(predicted_comp_reported_diffC[i]) ./ rho;
    predicted_comp_overall_diffA[i] = to_vector(predicted_comp_reported_diffC[i]) ./ rho * (1-psi) / psi;
  }
  for(i in 1:(S+G)) predicted_comp_diffM[i] = predicted_overall_incidence_deaths[i] == 0 ? rep_array(0,K) : multinomial_rng(output_agedistr_deaths,predicted_overall_incidence_deaths[i]);
  for(i in 1:K) {
    predicted_total_reported_symptomatic_cases_by_age[i] = sum(predicted_comp_reported_diffC[1:S,i]);
    predicted_total_overall_symptomatic_cases_by_age[i] = sum(predicted_comp_overall_diffC[1:S,i]);
    predicted_total_overall_all_cases_by_age[i] = sum(predicted_comp_overall_diffC[1:S,i]) + sum(predicted_comp_overall_diffA[1:S,i]);
    predicted_total_overall_deaths_tmax_by_age[i] = sum(predicted_comp_diffM[1:S,i]);
    predicted_total_overall_deaths_delay_by_age[i] = sum(predicted_comp_diffM[1:(S+G),i]);
  }
  predicted_total_reported_symptomatic_cases = sum(predicted_total_reported_symptomatic_cases_by_age);
  predicted_total_overall_symptomatic_cases = sum(predicted_total_overall_symptomatic_cases_by_age);
  predicted_total_overall_all_cases = sum(predicted_total_overall_all_cases_by_age);
  predicted_total_overall_deaths_tmax = sum(predicted_total_overall_deaths_tmax_by_age);
  predicted_total_overall_deaths_delay = sum(predicted_total_overall_deaths_delay_by_age);
  
  cfr_A_symptomatic_by_age = predicted_total_overall_deaths_tmax_by_age ./ predicted_total_reported_symptomatic_cases_by_age;
  cfr_B_symptomatic_by_age = predicted_total_overall_deaths_delay_by_age ./ predicted_total_reported_symptomatic_cases_by_age;
  cfr_C_symptomatic_by_age = predicted_total_overall_deaths_tmax_by_age ./ predicted_total_overall_symptomatic_cases_by_age;
  cfr_D_symptomatic_by_age = predicted_total_overall_deaths_delay_by_age ./ predicted_total_overall_symptomatic_cases_by_age;
  
  cfr_A_symptomatic = predicted_total_overall_deaths_tmax / predicted_total_reported_symptomatic_cases;
  cfr_B_symptomatic = predicted_total_overall_deaths_delay / predicted_total_reported_symptomatic_cases;
  cfr_C_symptomatic = predicted_total_overall_deaths_tmax / predicted_total_overall_symptomatic_cases;
  cfr_D_symptomatic = predicted_total_overall_deaths_delay / predicted_total_overall_symptomatic_cases;
  
  cfr_C_all_by_age = predicted_total_overall_deaths_tmax_by_age ./ predicted_total_overall_all_cases_by_age;
  cfr_D_all_by_age = predicted_total_overall_deaths_delay_by_age ./ predicted_total_overall_all_cases_by_age;
  
  cfr_C_all = predicted_total_overall_deaths_tmax / predicted_total_overall_all_cases;
  cfr_D_all = predicted_total_overall_deaths_delay / predicted_total_overall_all_cases;
}
