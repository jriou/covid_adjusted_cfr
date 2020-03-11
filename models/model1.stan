functions {
  real switch_zero(real t, real t1) {
    return(1/(1+exp(5*(t-t1))));
  }
  real switch_eta(real t, real t1, real eta) {
    return(eta+(1-eta)/(1+exp(5*(t-t1))));
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
  real dydt[(6*K)]; //SEIR, then C (cumulative incidence) and M (mortality)
  real nI; // total infectious
  real ntot;
  
  real beta; //transmission rate
  real eta; // reduction in transmission rate after quarantine
  real epsilon[K]; //case fatality rate per age group
  real tau; //incubation period
  real mu; //time to recovery (or death)
  real p_tmax; //tau is 0 after tmax (stop "recruiting")
  real p_tswitch;
  
  beta = theta[1];
  eta = theta[2];
  tau = theta[3];
  mu = theta[4];
  epsilon = theta[5:(5+K-1)];

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
    // R
    dydt[3*K+k] = mu * (1.0 - epsilon[k]) * y[2*K+k];
    // Cumulative Incidence
    dydt[4*K+k] = tau* p_tmax * y[K+k];
    //Cumulative mortality
    dydt[5*K+k] = mu * epsilon[k] * y[2*K+k];
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
    real p_onset_to_death[2];
    real p_pi[2];
    real p_incubation[2];
    real p_epsilon[2];
    real p_rho[2];
    real p_phi;
    //Simulation
    real t0; //starting time
    int t_data; //time of first data
    int S;
    real ts[S]; // time bins
    
    int inference;
    int doprint;
}

transformed data {
    real x_r[2] = {tmax2,tswitch};
    int x_i[1] = {K};
}

parameters{
    real<lower=0> beta; // base transmission rate
    real<lower=0,upper=1> eta; // reduction in transmission rate after quarantine measures
    real<lower=0> incubation; // incubation period (Lindon 2020)
    real<lower=0> onset_to_death; // duration from disease onset to death (Lindon 2020)
    vector<lower=0,upper=1> [K] epsilon; // age-dependent mortality probability
    vector<lower=0,upper=1> [K-1] rho; // age-dependent reporting probability
    real<lower=0, upper=1> pi; // number of cases at t0
    real<lower=0> phi[2]; // variance parameters
}
transformed parameters {
  // transformed parameters
  vector[K] rho_K;
  real<lower=0> tau;
  real<lower=0> mu;
  // change of format for integrate_ode_rk45
  real theta[5+K-1]; // vector of parameters
  real init[K*6]; // initial values
  real y[S,K*6]; // raw ODE output
  real y_pop[S,K*6]; // transformed ODE output 
  // outcomes
  vector[D] output_incidence_cases; // overall case incidence by day
  vector[D] output_incidence_deaths; // overal mortality incidence by day 
  simplex[K] output_agedistr_cases; // final age distribution of cases
  simplex[K] output_agedistr_deaths; // final age distribution of deaths
  
  // transformed parameters
  tau = 1.0/incubation;
  mu = 1.0/onset_to_death;
  for(i in 1:(K-1)){
    rho_K[i]=rho[i];
  }
  rho_K[K]=1.0;
  // change of format for integrate_ode_rk45
  theta[1:4] = {beta,eta,tau,mu};
  theta[5:(5+K-1)] = to_array_1d(epsilon);
  for(k in 1:K){
    init[k] = age_dist[k] * (1-pi);
    init[K+k] = age_dist[k] * pi;
    init[2*K+k] = 0.0;
    init[3*K+k] = 0.0;
    init[4*K+k] = 0.0;
    init[5*K+k] = 0.0;
  }
  // run ODE solver
  y = integrate_ode_rk45(
    SEIR, 
    init, // initial states
    t0, // t0
    ts, // evaluation dates (ts)
    theta, // parameters
    x_r, // real data
    x_i, // integer data
    1.0E-10, 1.0E-10, 1.0E3);
  // correct for imprecision of 10e-10 leading to some negative values
  for(k in 1:(K*6)) {
    for(j in 1:S){
      y_pop[j,k] = (y[j,k] + 1.0E-9)*pop_t;
    }
  }
  // compute outcomes
  for(i in t_data:tmax){
    output_incidence_cases[i-t_data+1] = sum((to_vector(y_pop[i,(4*K+1):(5*K)]) + 1.0E-9*pop_t - to_vector(y_pop[i-1,(4*K+1):(5*K)])).*rho_K); 
    output_incidence_deaths[i-t_data+1] =  sum((to_vector(y_pop[i,(5*K+1):(6*K)]) + 1.0E-9*pop_t - to_vector(y_pop[i-1,(5*K+1):(6*K)])));
  }
  output_agedistr_cases = (to_vector(y_pop[tmax,(4*K+1):(5*K)]).*rho_K)./sum(to_vector(y_pop[tmax,(4*K+1):(5*K)]).*rho_K);
  output_agedistr_deaths = (to_vector(y_pop[tmax,(5*K+1):(6*K)]))./sum(to_vector(y_pop[tmax,(5*K+1):(6*K)]));
}
model {
    beta ~ exponential(p_beta);
    eta ~ beta(p_eta[1],p_eta[2]);
    incubation ~ gamma(p_incubation[1]^2/p_incubation[2]^2,p_incubation[1]/p_incubation[2]^2);
    onset_to_death ~ gamma(p_onset_to_death[1]^2/p_onset_to_death[2]^2,p_onset_to_death[1]/p_onset_to_death[2]^2);
    for(k in 1:K) epsilon[k] ~ beta(p_epsilon[1],p_epsilon[2]);
    for(k in 1:(K-1)) rho[k] ~ beta(p_rho[1],p_rho[2]);
    pi ~ beta(p_pi[1],p_pi[2]);
    phi ~ exponential(p_phi);
    
    if(doprint==1) {
      print("beta: ",beta);
      print("eta: ",beta);
      print("tau: ",tau);
      print("incubation: ",incubation);
      print("mu: ",mu);
      print("onset_to_death: ",onset_to_death);
      print("epsilon: ",epsilon);
      print("rho: ",rho);
      print("pi: ",pi);
  
      print("y[1,]: ",y[1,]);
      print("output_incidence_cases: ", output_incidence_cases);
      print("output_incidence_deaths: ", output_incidence_deaths);
    }

    if (inference==1){
      for(i in 1:D){
        target += neg_binomial_2_lpmf( incidence_cases[i] | output_incidence_cases[i], output_incidence_cases[i]/phi[1]);
        target += neg_binomial_2_lpmf( incidence_deaths[i] | output_incidence_deaths[i],output_incidence_deaths[i]/phi[2]);
      }
      target +=multinomial_lpmf(agedistr_cases | output_agedistr_cases);
      target +=multinomial_lpmf(agedistr_deaths | output_agedistr_deaths);
    }
  }

generated quantities{
  real beta2 = beta*eta;
  vector[K] predicted_reported_incidence_cases[S-1]; 
  vector[K] predicted_overall_incidence_cases[S-1]; 
  vector[K] predicted_overall_incidence_deaths[S-1]; 

  vector [K] cfr_A; //cfr by age classes, no correction of underreporting, no correction of time lag
  vector [K] cfr_B; //cfr by age classes, no correction of underreporting, correction of time lag
  vector [K] cfr_C; //cfr by age classes, correction of underreporting, no correction of time lag
  vector [K] cfr_D; //cfr by age classes, correction of underreporting, correction of time lag
  real cfr_overall[4]; //cfr overall, for method A, B, C, D respectively 
  
  for(i in 2:S){
      predicted_reported_incidence_cases[i-1] = to_vector(y[i,(4*K+1):(5*K)]).*rho_K - to_vector(y[i-1,(4*K+1):(5*K)]).*rho_K;
      predicted_overall_incidence_cases[i-1] = to_vector(y[i,(4*K+1):(5*K)]) - to_vector(y[i-1,(4*K+1):(5*K)]);
      predicted_overall_incidence_deaths[i-1] = to_vector(y[i,(5*K+1):(6*K)]) - to_vector(y[i-1,(5*K+1):(6*K)]);
  }
  
  cfr_A=(to_vector(y[tmax,(5*K+1):(6*K)]))./(to_vector(y[tmax,(4*K+1):(5*K)]).*rho_K);
  cfr_B=(to_vector(y[S,(5*K+1):(6*K)]))./(to_vector(y[S,(4*K+1):(5*K)]).*rho_K);
  cfr_C=(to_vector(y[tmax,(5*K+1):(6*K)]))./(to_vector(y[tmax,(4*K+1):(5*K)]));
  cfr_D=(to_vector(y[S,(5*K+1):(6*K)]))./(to_vector(y[S,(4*K+1):(5*K)]));
  
  cfr_overall = {
    sum(to_vector(y[tmax,(5*K+1):(6*K)]))/sum(to_vector(y[tmax,(4*K+1):(5*K)]).*rho_K),
    sum(to_vector(y[S,(5*K+1):(6*K)]))/sum(to_vector(y[S,(4*K+1):(5*K)]).*rho_K),
    sum(to_vector(y[tmax,(5*K+1):(6*K)]))/sum(to_vector(y[tmax,(4*K+1):(5*K)])),
    sum(to_vector(y[S,(5*K+1):(6*K)]))/sum(to_vector(y[S,(4*K+1):(5*K)]))
    };
 }
