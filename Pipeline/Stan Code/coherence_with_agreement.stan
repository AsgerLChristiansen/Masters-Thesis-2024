

data {
  int N_channels; // How many channels?
  int N_visits; // How many repeat measures per channel? That is to say, visits.
  int N_pairs;  // How many pairs?
  int N_pairs_real; // how many of them are real?
  int N_pairs_surrogate; // and how many are surrogate=?
  int N_conditions;
  // matrix of coherence estimates per channel per visit.
  array[N_visits, N_conditions, N_channels, N_pairs] real<lower=0,upper=1> coherence;
  array[N_pairs] vector<lower=0,upper=1>[N_visits]  SVS_agreement;
  // 2 visits. Coded 0 and 1.
  vector[N_visits] visit;
}




parameters { 
  
  //Group level parameters for group 1 and 2
  real mu_beta_1_real;
  real mu_beta_1_surrogate;
  real mu_beta_2_real;
  real mu_beta_2_surrogate;
  real mu_beta_3_real;
  real mu_beta_3_surrogate;
  
  real<lower=0.1,upper=2> sigma_beta_1_real;
  real<lower=0.1,upper=2> sigma_beta_1_surrogate;
  real<lower=0.1,upper=2> sigma_beta_2_real;
  real<lower=0.1,upper=2> sigma_beta_2_surrogate;
  real<lower=0.1,upper=2> sigma_beta_3_real;
  real<lower=0.1,upper=2> sigma_beta_3_surrogate;
  
  // Pair level parameters for beta 1 and beta 2
  vector[N_pairs] mu_beta_1_pair;
  vector[N_pairs] mu_beta_2_pair;
  vector[N_pairs] mu_beta_3_pair;
  
  
  vector<lower=0.1,upper=2>[N_pairs] sigma_beta_1_pair;
  vector<lower=0.1,upper=2>[N_pairs] sigma_beta_2_pair;
  vector<lower=0.1,upper=2>[N_pairs] sigma_beta_3_pair;
  // Channel level parameters for beta_1 and beta_2 (Visit4)
  matrix[N_channels, N_pairs] mu_beta_1_channel;
  matrix[N_channels, N_pairs] mu_beta_2_channel;
  matrix[N_channels, N_pairs] mu_beta_3_channel;
  matrix<lower=0.1,upper=2>[N_channels, N_pairs] sigma_beta_1_channel;
  matrix<lower=0.1,upper=2>[N_channels, N_pairs] sigma_beta_2_channel;
  matrix<lower=0.1,upper=2>[N_channels, N_pairs] sigma_beta_3_channel;


  // Condition level parameters  (The lowest level)
  // Note that matrices are no longer possible.
  array[N_conditions, N_channels, N_pairs] real beta_1;
  
  array[N_conditions, N_channels, N_pairs] real beta_2;
  array[N_conditions, N_channels, N_pairs] real beta_3;
  
  
  vector<lower=0.1,upper=2>[N_conditions] sigma_condition;

}


transformed parameters {
  
  // Mu is dependent on all of the above.
  array[N_conditions, N_channels, N_pairs] vector[N_visits] mu;
  
  
  
  
  
  // Tried and failed to avoid this many for-loops.
  for (p in 1:N_pairs){
  for (i in 1:N_channels){
  for (c in 1:N_conditions){ 
    if (ceil(coherence[1, c, i, p]) !=0 && ceil(coherence[2, c, i, p])!=0 ){
    
  //mu[c,i,p] = beta_1[c,i,p] - beta_1[c,i,p] * visit +  beta_2[c,i,p] * visit; // For visit in {0,1}
  mu[c,i,p] = beta_1[c,i,p] - beta_1[c,i,p] * visit +  beta_2[c,i,p] * visit + beta_3[c,i,p]*SVS_agreement[p]; // For visit in {0,1}
  
  
  
  }
  
  else if (ceil(coherence[1, c, i, p])!=0 ){
    mu[c,i,p,1] = beta_1[c,i,p] - beta_1[c,i,p] * visit[1] +  beta_2[c,i,p] * visit[1]+ beta_3[c,i,p]*SVS_agreement[p,1];
    
  }
  else if (ceil(coherence[2, c, i, p])!=0 ){
    mu[c,i,p,2] = beta_1[c,i,p] - beta_1[c,i,p] * visit[2] +  beta_2[c,i,p] * visit[2]+ beta_3[c,i,p]*SVS_agreement[p,2];
    
  }
  
  }
  }
}}

model {
  
  // Start at group level - real pairs versus surrogate pairs.
  
  // expected baseline coherence at visit 1 (In logit. -8 ~= 0.3)
  mu_beta_1_real ~ normal(-0.80, 0.5);
  // expected baseline coherence at visit 2 (Slightly larger)
  mu_beta_2_real ~  normal(-0.6, 0.5);
  // expected effect of hypothetical 100% agreement (roughly one SD in logit scale).
  // Generally, it should be very surprising (though not impossible) if it was negative, but its perfectly
  // plausible that the effect is close to 0.
  mu_beta_3_real ~  normal(0.4, 0.2);
 
  sigma_beta_1_real ~ gamma(1, 3);
  
  sigma_beta_2_real ~ gamma(1, 3);
  sigma_beta_3_real ~ gamma(1, 3);
  
  // We expect no increase in coherence from visit 1 to visit 4.
  mu_beta_1_surrogate ~ normal(-0.80,0.5);
  mu_beta_2_surrogate ~ normal(-0.80,0.5);
  
  // In the surrogate case, we generally expect great variation due to
  // random pairings. We expect that the result is close to zero, or even
  // negative, since the relation between surrogate agreement and synchrony
  // should be unsystematic in nature.
  mu_beta_3_surrogate ~ normal(0,0.3);
  
  sigma_beta_1_surrogate ~  gamma(1, 3);
  sigma_beta_2_surrogate ~  gamma(1, 3);
  sigma_beta_3_surrogate ~  gamma(1, 3);
  
  //Moving into pair level parameters
  target += normal_lpdf(mu_beta_1_pair[1:N_pairs_real] | mu_beta_1_real, sigma_beta_1_real);
  
  target += normal_lpdf(  mu_beta_1_pair[(N_pairs_real+1):(N_pairs_real + N_pairs_surrogate)] |mu_beta_1_surrogate ,  sigma_beta_1_surrogate);
  
  target += normal_lpdf(mu_beta_2_pair[1:N_pairs_real] | mu_beta_2_real, sigma_beta_2_real);
  
  target += normal_lpdf(  mu_beta_2_pair[(N_pairs_real+1):(N_pairs_real + N_pairs_surrogate)] |mu_beta_2_surrogate ,  sigma_beta_2_surrogate);

  target += normal_lpdf(mu_beta_3_pair[1:N_pairs_real] | mu_beta_3_real, sigma_beta_3_real);
  
  target += normal_lpdf(  mu_beta_3_pair[(N_pairs_real+1):(N_pairs_real + N_pairs_surrogate)] |mu_beta_3_surrogate ,  sigma_beta_3_surrogate);

  
  sigma_beta_1_pair  ~ gamma(1, 3);
  sigma_beta_2_pair  ~ gamma(1, 3);
  sigma_beta_3_pair  ~ gamma(1, 3);
  
  sigma_condition  ~ gamma(1, 3);
  
  
  // Channel-level parameters
  for (p in 1:N_pairs){
    
  target += normal_lpdf( mu_beta_1_channel[,p]  |mu_beta_1_pair[p], sigma_beta_1_pair[p]);
  
  target += normal_lpdf(  mu_beta_2_channel[,p] |mu_beta_2_pair[p], sigma_beta_2_pair[p]);
  target += normal_lpdf(  mu_beta_3_channel[,p] |mu_beta_3_pair[p], sigma_beta_3_pair[p]);
 
  sigma_beta_1_channel[,p]  ~ gamma(1, 3);
  sigma_beta_2_channel[,p]  ~ gamma(1, 3);
  sigma_beta_3_channel[,p]  ~ gamma(1, 3);
 
  // Condition-level parameters
  for (i in 1:N_channels){
    
    target += normal_lpdf(  beta_1[,i,p]   |mu_beta_1_channel[i,p], sigma_beta_1_channel[i,p]);
    target += normal_lpdf(  beta_2[,i,p]   |mu_beta_2_channel[i,p], sigma_beta_2_channel[i,p]);
    target += normal_lpdf(  beta_3[,i,p]   |mu_beta_3_channel[i,p], sigma_beta_3_channel[i,p]);
    
  // And finally, for every condition, predict logit-transformed coherence
  // for both visits using the appropriate mu-vector.
  for (c in 1:N_conditions){
    
    
    if (ceil(coherence[1, c, i, p])!=0  && ceil(coherence[2, c, i, p])!=0 ){
  target += normal_lpdf( logit(coherence[,c,i,p]) |mu[c,i,p], sigma_condition[c]);
    }
    
    else if (ceil(coherence[1, c, i, p]) !=0 ){
  target += normal_lpdf( logit(coherence[1,c,i,p]) |mu[c,i,p,1], sigma_condition[c]);
      
    }
    else if (ceil(coherence[2, c, i, p]) !=0 ){
  target += normal_lpdf( logit(coherence[2,c,i,p]) |mu[c,i,p,2], sigma_condition[c]);
      
    }
}}

}}



generated quantities {
  real mu_beta_1_real_prior;
  real mu_beta_1_surrogate_prior;
  real mu_beta_2_real_prior;
  real mu_beta_2_surrogate_prior;
  real mu_beta_3_real_prior;
  real mu_beta_3_surrogate_prior;
   mu_beta_1_real_prior =  normal_rng(-0.80, 0.5);
   mu_beta_1_surrogate_prior = normal_rng(-0.6, 0.5);
   mu_beta_2_real_prior = normal_rng(-0.80, 0.5);
   mu_beta_2_surrogate_prior = normal_rng(-0.80, 0.5);
   mu_beta_3_real_prior = normal_rng(0.4, 0.2);
   mu_beta_3_surrogate_prior = normal_rng(0, 0.3);

  
  real sigma_beta_1_real_prior;
  real sigma_beta_1_surrogate_prior;
  real sigma_beta_2_real_prior;
  real sigma_beta_2_surrogate_prior;
  real sigma_beta_3_real_prior;
  real sigma_beta_3_surrogate_prior;
   sigma_beta_1_real_prior =  gamma_rng(1, 3);
   sigma_beta_1_surrogate_prior = gamma_rng(1, 3);
   sigma_beta_2_real_prior = gamma_rng(1, 3);
   sigma_beta_2_surrogate_prior = gamma_rng(1, 3);
   sigma_beta_3_real_prior = gamma_rng(1, 3);
   sigma_beta_3_surrogate_prior = gamma_rng(1, 3);
  
}
