data {
  int<lower=0> K; // categories (genes)
  vector[K] alpha; // intercept
  vector[K] beta; // effect
  real xi; // precision
  int n; // total counts (tries)
}

generated quantities {
  // counts
  int y_a [K];
  int y_b [K];
  // intermediate dirichlet simulations
  vector [K] t [2];
  
  // simulate from dirichlet dist.
  t[1] = dirichlet_rng(xi*softmax(alpha + beta));
  t[2] = dirichlet_rng(xi*softmax(alpha - beta));
  
  // simulate from multinomial dist.
  y_a = multinomial_rng(t[1], n);
  y_b = multinomial_rng(t[2], n);
}
