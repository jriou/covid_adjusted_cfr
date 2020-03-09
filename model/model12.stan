functions {
  real switch_zero(real t, real t1) {
    return(1/(1+exp(5*(t-t1))));
  }
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
    real tmax2 = x_r[1];
    real tswitch = x_r[2];
    real dydt[(7*K)]; //SEIAR, then C and D
    real nI; // total infectious
    real ntot;
    
    real beta; // transmission rate
    real eta; // reduction in transmission rate after quarantine
    real xi; // slope of quarantine implementation
    real nu; // shift of quarantine implementation
    real tau; // incubation rate
    real mu; // infectious rate
    real psi; // probability of symptoms
    real p_tmax; // tau is 0 after tmax (stop "recruiting")
    real p_tswitch;
    real contact[K*K];//contact matrix, first K values, corresponds to number of contact between age class 1 and other classes, etc
    real n_by_age[K];
    real f_inf[K];
    
    // estimated parameters
    beta = theta[1];
    eta = theta[2];
    xi = theta[3];
    nu = theta[4];
    psi = theta[5];
    
    // fixed parameters
    tau = 1.0/x_r[3];
    mu = 1.0/x_r[4];
    contact = x_r[5:(4+K*K)];
    
    //Total number of infectious people
    nI = sum(y[(2*K+1):(3*K)]);
    ntot = sum(y[1:(4*K)]);
    p_tmax = switch_zero(t,tmax2);
    p_tswitch = switch_eta(t,tswitch,eta,nu,xi);
    
    //Total number of people (S+E+I+J+R) by age classes
    for(k in 1:K){
      n_by_age[k] = sum(y[{k,k+K,k+2*K,k+3*K,k+4*K}]);
    }
    //Force of infection by age classes: beta * p_tswitch * sum((number of infected people by age) / (total number of people by age) * (number of contact by age))
    for(k in 1:K){
      f_inf[k] = beta * p_tswitch * sum(to_vector(y[(2*K+1):(3*K)]) ./ to_vector(n_by_age) .* to_vector(contact[(K*(k-1)+1):(k*K)])); //
    }
    //print("Force of infection",f_inf);
    
    for (k in 1:K) {
     // S
      dydt[k] = - f_inf[k] * y[k]; 
      // E
      dydt[K+k] = f_inf[k] * y[k]- tau * p_tmax * y[K+k];
      // I
      dydt[2*K+k] = psi * tau * p_tmax * y[K+k] - mu * y[2*K+k];
      // A
      dydt[3*K+k] = (1-psi) * tau * p_tmax * y[K+k] - mu * y[3*K+k];
      // R
      dydt[4*K+k] =  mu * y[2*K+k] +  mu * y[3*K+k];
      // C
      dydt[5*K+k] = psi * tau * p_tmax * y[K+k];
      // D
      dydt[6*K+k] = (1-psi) * tau * p_tmax * y[K+k];
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
  real x_r[4+K*K];
  int x_i[1] = {K};
  x_r[1:4] = {tmax2,tswitch,p_incubation,p_infectious};
  x_r[5:(4+K*K)] = contact;
}

parameters{
  real<lower=0> beta; // base transmission rate
  real<lower=0,upper=1> eta; // reduction in transmission rate after quarantine measures
  vector<lower=0,upper=1> [K] epsilon; // age-dependent mortality probability
  vector<lower=0,upper=1> [K-1] rho; // age-dependent reporting probability
  real<lower=0, upper=1> pi; // number of cases at t0
  real<lower=0> phi[2]; // variance parameters
  real<lower=0,upper=1> xi_raw; // slope of quarantine implementation
  real<lower=0> nu; // shift of quarantine implementation
  real<lower=0,upper=1> psi; // proportion of symptomatics
}
transformed parameters {
  // transformed parameters
  vector[K] rho_K;
  real xi = xi_raw+0.5;
  // change of format for integrate_ode_rk45
  real theta[5]; // vector of parameters
  real init[K*7]; // initial values
  real y[S,K*7]; // raw ODE output
  vector[K] comp_S[S];
  vector[K] comp_E[S];
  vector[K] comp_I[S];
  vector[K] comp_A[S];
  vector[K] comp_R[S];
  vector[K] comp_C[S];
  vector[K] comp_D[S];
  vector[K] comp_diffC[S];
  vector[K] comp_diffD[S];
  vector[K] comp_diffM[S];
  vector[K] comp_M[S];
  // outcomes
  vector[D] output_incidence_cases; // overall case incidence by day
  vector[D] output_incidence_deaths; // overal mortality incidence by day 
  simplex[K] output_agedistr_cases; // final age distribution of cases
  simplex[K] output_agedistr_deaths; // final age distribution of deaths
  
  // transformed parameters
  for(i in 1:(K-1)){
    rho_K[i]=rho[i];
  }
  rho_K[K]=1.0;
  // change of format for integrate_ode_rk45
  theta[1:5] = {beta,eta,xi,nu,psi};
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
      comp_A[i] = (to_vector(y[i,(3*K+1):(4*K)]) + 1.0E-9) * pop_t;
      comp_R[i] = (to_vector(y[i,(4*K+1):(5*K)]) + 1.0E-9) * pop_t;
      comp_C[i] = (to_vector(y[i,(5*K+1):(6*K)]) + 1.0E-9) * pop_t;
      comp_D[i] = (to_vector(y[i,(6*K+1):(7*K)]) + 1.0E-9) * pop_t;
      comp_diffC[i] = i==1 ? comp_C[i,] : 1.0E-9*pop_t + comp_C[i,] - comp_C[i-1,]; // lagged difference of cumulative incidence of symptomatics
      comp_diffD[i] = i==1 ? comp_D[i,] : 1.0E-9*pop_t + comp_D[i,] - comp_D[i-1,]; // lagged difference of cumulative incidence of asymptomatics
    }
    // compute mortality
    for(i in 1:S) {
      comp_diffM[i] = rep_vector(1.0E-9,K);
      comp_M[i] = rep_vector(0.0,K);
      for(g in 1:G) {
        if (i>g) comp_diffM[i] += comp_diffC[i-g] .* epsilon * p_gamma[g] ;
      }
      comp_M[i] += comp_diffM[i];
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
  for(k in 1:K) epsilon[k] ~ beta(p_epsilon[1],p_epsilon[2]);
  for(k in 1:(K-1)) rho[k] ~ beta(p_rho[1],p_rho[2]);
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
  real avg_rho = sum(age_dist .* rho_K);
  real beta2 = beta*eta;
  
  int predicted_reported_incidence_symptomatic_cases[S]; 
  real predicted_overall_incidence_symptomatic_cases[S]; 
  real predicted_overall_incidence_all_cases[S]; 
  int predicted_overall_incidence_deaths[S];
  
  int predicted_comp_reported_diffC[S,K];
  vector[K] predicted_comp_overall_diffC[S];
  vector[K] predicted_comp_overall_diffA[S];
  int predicted_comp_diffM[S,K];
  
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
    predicted_reported_incidence_symptomatic_cases[i] = neg_binomial_2_rng(sum(comp_diffC[i].*rho_K), sum(comp_diffC[i].*rho_K)/phi[1]);
    predicted_overall_incidence_symptomatic_cases[i] = predicted_reported_incidence_symptomatic_cases[i] / avg_rho;
    predicted_overall_incidence_all_cases[i] = predicted_reported_incidence_symptomatic_cases[i] / avg_rho / psi;
    predicted_overall_incidence_deaths[i] = neg_binomial_2_rng(sum(comp_diffM[i]),sum(comp_diffM[i])/phi[2]);
  }
  for(i in 1:S) {
    predicted_comp_reported_diffC[i] = predicted_reported_incidence_symptomatic_cases[i] == 0 ? rep_array(0,K) : multinomial_rng(output_agedistr_cases,predicted_reported_incidence_symptomatic_cases[i]);
    predicted_comp_overall_diffC[i] = to_vector(predicted_comp_reported_diffC[i]) ./ rho_K;
    predicted_comp_overall_diffA[i] = to_vector(predicted_comp_reported_diffC[i]) ./ rho_K / psi;
    predicted_comp_diffM[i] = predicted_overall_incidence_deaths[i] == 0 ? rep_array(0,K) : multinomial_rng(output_agedistr_deaths,predicted_overall_incidence_deaths[i]);
  }
  for(i in 1:K) {
    predicted_total_reported_symptomatic_cases_by_age[i] = sum(predicted_comp_reported_diffC[1:tmax,i]);
    predicted_total_overall_symptomatic_cases_by_age[i] = sum(predicted_comp_overall_diffC[1:tmax,i]);
    predicted_total_overall_all_cases_by_age[i] = sum(predicted_comp_overall_diffC[1:tmax,i]) + sum(predicted_comp_overall_diffA[1:tmax,i]);
    predicted_total_overall_deaths_tmax_by_age[i] = sum(predicted_comp_diffM[1:tmax,i]);
    predicted_total_overall_deaths_delay_by_age[i] = sum(predicted_comp_diffM[1:S,i]);
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
