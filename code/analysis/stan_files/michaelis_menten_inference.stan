functions {
    // Define Michaelis-Menten model
    vector enzyme_kinetics(
        real t,  // time
        vector s,  // substrate concentration
        real[] theta // parameters (V_max, K_M)
    ) {
        // Initialize dS/dt
        vector[1] ds_dt;
        // Compute dynamics
        ds_dt = - theta[1] * s / (theta[2] + s[1]);
          
    return ds_dt;
    }
}

data {
    // Inference information
    int<lower=1> n_sample;  // number of time samples
    real t_sample[n_sample];  // time points observed
    real s_[n_sample];  // data, substrate concentration
    vector[1] s0;  // Initial substrate concentration
    real t0;  // Initial time

    //  Simulation information
    int<lower=1> n_sim;  // number of time points for numerical integration
    real t_sim[n_sim];  // time points for numerical integration

    // Parameters for prior distributions
    real vmax_param[2];  // parameters for prior dis. of V_max
    real km_param[2];  // parameters for prior dis. of K_M
    real sigma_param[2];  // parameters for prior dis. of sigma^2 
}

transformed data {

}

parameters {
    // Initialize parameters
    real<lower=0> vmax_;  // maximum catalytic rate
    real<lower=0> km_;  // Michaelis-Menten constant
    real<lower=0> sigma_;  // Observational error
}

transformed parameters{
    // Pack parameters together for ODE solver
    real<lower=0> theta[2];  
    theta[1] = vmax_;
    theta[2] = km_;

    // Numerically integrate the SIR ODEs
    vector[1] s_hat[n_sample];   // solution from the ODE solver
    s_hat = ode_rk45(enzyme_kinetics, s0, t0, t_sample, theta);
}

model {
    // Priors
    vmax_ ~ normal(vmax_param[1], vmax_param[2]);
    km_ ~ normal(km_param[1], km_param[2]);
    sigma_ ~ normal(sigma_param[1], sigma_param[2]); 
    
    // Likelihood
    // Unpack predictions into array
    real s_int[n_sample];
    for (i in 1:n_sample){
        s_int[i] = s_hat[i][1];
    }
    s_ ~ normal(s_int, sigma_);
}

generated quantities {
    // Run numerical integration with a finer grid 
    vector[1] s_sim[n_sim];
    s_sim = ode_rk45(enzyme_kinetics, s0, t0, t_sim, theta);
    
    // Simulate observational model
    real s_tilde[n_sim];
    for (i in 1:n_sim) {
        s_tilde[i] = normal_rng(s_sim[i][1], sigma_);
    }
}