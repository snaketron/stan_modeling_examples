data {
  int <lower = 0> N_sample; // number of samples (repertoires)
  int <lower = 0> N_gene; // number of genes
  int Y_1 [N_gene, N_sample]; // number of successes (cells) in samples x gene
  int Y_2 [N_gene, N_sample]; // number of successes (cells) in samples x gene
  int N [N_sample, 2]; // number of total tries (repertoire size)
  // test data
  int <lower = 0> N_sample_test; // number of samples (repertoires)
  int Y_1_test [N_gene, N_sample_test]; // number of successes (cells) in samples x gene
  int Y_2_test [N_gene, N_sample_test]; // number of successes (cells) in samples x gene
  int N_test [N_sample_test, 2]; // number of total tries (repertoire size)
}

transformed data {
  real N_real [N_sample, 2];
  N_real = N;
}

parameters {
  vector [N_gene] alpha_mu_gene;
  real <lower=0> beta_sigma_gene;
  real <lower=0> alpha_sigma_gene;
  real <lower=0> beta_sigma_pop;
  vector [N_gene] beta_z [N_sample];
  vector [N_gene] alpha_z [N_sample];
  vector [N_gene] beta_z_gene;
}

transformed parameters {
  vector [N_gene] alpha [N_sample];
  vector [N_gene] beta [N_sample];
  vector [N_gene] beta_mu_gene;
  
  beta_mu_gene = 0+beta_sigma_pop*beta_z_gene;
  for(i in 1:N_sample) {
    beta[i]=beta_mu_gene+beta_sigma_gene*beta_z[i];
    alpha[i]=alpha_mu_gene+alpha_sigma_gene*alpha_z[i];
  }
}

model {
  target += cauchy_lpdf(beta_sigma_pop | 0, 1);
  target += cauchy_lpdf(alpha_sigma_gene | 0, 1);
  target += cauchy_lpdf(beta_sigma_gene | 0, 1);
  for(i in 1:N_sample) {
    target += normal_lpdf(alpha_z[i] | 0, 1);
    target += normal_lpdf(beta_z[i] | 0, 1);
  }
  target += normal_lpdf(beta_z_gene | 0, 1);
  target += normal_lpdf(alpha_mu_gene | 0, 5);
  
  // likelihood
  for(i in 1:N_sample) {
    target += multinomial_lpmf(Y_1[,i] | softmax(alpha[i]-beta[i]));
    target += multinomial_lpmf(Y_2[,i] | softmax(alpha[i]+beta[i]));
  }
}

generated quantities {
  int Y_hat_1 [N_gene, N_sample];
  int Y_hat_2 [N_gene, N_sample];
  int Y_hat_group_1 [N_gene, N_sample];
  int Y_hat_group_2 [N_gene, N_sample];
  real log_lik [N_sample, 2];
  real log_lik_test [N_sample_test, 2];
  real log_lik_train [N_sample, 2];
  
  vector [N_gene] mu [2];
  
  // PPC
  for(i in 1:N_sample) {
    mu[1] = softmax(alpha[i] - beta[i]);
    mu[2] = softmax(alpha[i] + beta[i]);
    
    Y_hat_1[,i] = multinomial_rng(mu[1], sum(Y_1[,i]));
    Y_hat_2[,i] = multinomial_rng(mu[2], sum(Y_2[,i]));
    log_lik[i,1] = multinomial_lpmf(Y_1[,i] | mu[1]);
    log_lik[i,2] = multinomial_lpmf(Y_2[,i] | mu[2]);
  }
  
  // PPC - test (condition-specific)
  mu[1]=softmax((alpha_mu_gene+alpha_sigma_gene)-(beta_mu_gene+beta_sigma_gene));
  mu[2]=softmax((alpha_mu_gene+alpha_sigma_gene)+(beta_mu_gene+beta_sigma_gene));
  for(i in 1:N_sample_test) {
    log_lik_test[i,1] = multinomial_lpmf(Y_1_test[,i]|mu[1]);
    log_lik_test[i,2] = multinomial_lpmf(Y_2_test[,i]|mu[2]);
  }
  for(i in 1:N_sample) {
    log_lik_train[i,1] = multinomial_lpmf(Y_1[,i]|mu[1]);
    log_lik_train[i,2] = multinomial_lpmf(Y_2[,i]|mu[2]);
    Y_hat_group_1[,i] = multinomial_rng(mu[1], N[i,1]);
    Y_hat_group_2[,i] = multinomial_rng(mu[2], N[i,2]);
  }
}