functions {
    // Define Michaelis-Menten model
    real[] enzyme_kinetics(
        real t,  // time
        real[] s,  // substrate concentration
        real[] theta, // parameters (V_max, K_M)
        data real[] x_r,
        data int[] x_i
    ) {
        // Initialize dS/dt
        real ds_dt[2];
        // Compute dynamics
        ds_dt[1] = - (theta[1] * s[2]) * s[1] / (theta[2] + s[1]);
        ds_dt[2] = - theta[3] * s[2];
          
    return ds_dt;
    }
}

data {
    // Inference information
    int<lower=1> n_replicate;  // number of replicates
    int<lower=1> n_sample;  // number of time samples
    real t_sample[n_replicate, n_sample];  // time points observed
    real s_[n_replicate, n_sample];  // data, substrate concentration
    real t0;  // Initial time

    //  Simulation information
    int<lower=1> n_sim;  // number of time points for numerical integration
    real t_sim[n_sim];  // time points for numerical integration

    // Parameters for prior distributions
    real vmax_param[2];  // parameters for prior dis. of V_max
    real km_param[2];  // parameters for prior dis. of K_M
    real s0_param[2];  // parameters for prior dis. of s0
    real deg_param[2];  // parameters for prior dis. of enzyme degradation rate
    real sigma_param[2];  // parameters for prior dis. of sigma^2 
}

transformed data {
    real x_r[0];
    int x_i[0];
}

parameters {
    // Initialize parameters
    real<lower=0> vmax_;  // maximum catalytic rate
    real<lower=0> km_;  // Michaelis-Menten constant
    real<lower=0> s0_;  // Michaelis-Menten constant
    real<lower=0> deg_;  // enzyme degradation rate
    real<lower=0> sigma_;  // observational error
}

transformed parameters{
    // Pack parameters together for ODE solver
    real<lower=0> theta[3];  
    theta[1] = vmax_;
    theta[2] = km_;
    theta[3] = deg_;

    // Convert s0_ to the right format
    real s0[2];
    s0[1] = s0_;  // initial concentration
    s0[2] = 1;  // initial total enzyme
    
}

model {
    // Priors
    vmax_ ~ normal(vmax_param[1], vmax_param[2]);
    km_ ~ normal(km_param[1], km_param[2]);
    s0_ ~ normal(s0_param[1], s0_param[2]);
    deg_ ~ normal(deg_param[1], deg_param[2]);
    sigma_ ~ normal(sigma_param[1], sigma_param[2]); 
    
    // Likelihood
    
    // Loop through replicates
    for (j in 1:n_replicate) {
        // Numerically integrate the SIR ODEs
        real s_hat[n_sample, 2];   // solution from the ODE solver
        s_hat = integrate_ode_rk45(
            enzyme_kinetics, s0, t0, t_sample[j, :], theta, 
            x_r, x_i
        );

        s_[j,:] ~ normal(s_hat[:, 1], sigma_);
    }
}

generated quantities {
    // Run numerical integration with a finer grid 
    real s_sim[n_sim, 2];
    s_sim = integrate_ode_rk45(enzyme_kinetics, s0, t0, t_sim, theta, x_r, x_i);
    
    // Simulate observational model
    real s_tilde[n_sim];
    for (i in 1:n_sim) {
        s_tilde[i] = normal_rng(s_sim[i, 1], sigma_);
    }
}