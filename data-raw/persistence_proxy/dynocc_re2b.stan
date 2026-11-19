// Dynamic occupancy, v2: detection per visit with fixed effects of survey window (run 1-4),
// route and stop random effects; colonisation and extinction with route and year random
// effects. Latent states marginalised with the forward algorithm.
data {
  int<lower=1> N; int<lower=1> T; int<lower=1> R; int<lower=1> W;
  array[N] int<lower=1, upper=R> route;
  array[N] int<lower=1, upper=T> first;
  array[N, T] int<lower=0> cell;          // 0 = no visit, else index into cell ranges
  int<lower=1> K; array[K] int<lower=1> vstart; array[K] int<lower=0> vend;
  int<lower=1> V; array[V] int<lower=0, upper=1> y; array[V] int<lower=1, upper=W> win;
}
parameters {
  real a_psi; real a_col; real a_ext; vector[W] b_p;
  real<lower=0> s_col_r; real<lower=0> s_ext_r; real<lower=0> s_col_t; real<lower=0> s_ext_t;
  real<lower=0> s_p_r; real<lower=0> s_p_i;
  vector[R] z_col_r; vector[R] z_ext_r; vector[T] z_col_t; vector[T] z_ext_t;
  vector[R] z_p_r; vector[N] z_p_i;
}
model {
  a_psi ~ normal(0, 1.5); a_col ~ normal(-1.5, 1.5); a_ext ~ normal(-1.5, 1.5); b_p ~ normal(0, 1.5);
  s_col_r ~ normal(0, 1); s_ext_r ~ normal(0, 1); s_col_t ~ normal(0, 1); s_ext_t ~ normal(0, 1);
  s_p_r ~ normal(0, 1); s_p_i ~ normal(0, 1);
  z_col_r ~ std_normal(); z_ext_r ~ std_normal(); z_col_t ~ std_normal(); z_ext_t ~ std_normal();
  z_p_r ~ std_normal(); z_p_i ~ std_normal();
  for (i in 1:N) {
    int r = route[i];
    real off_p = s_p_r * z_p_r[r] + s_p_i * z_p_i[i];
    vector[2] la;
    la[1] = log1m_inv_logit(a_psi); la[2] = log_inv_logit(a_psi);
    for (t in first[i]:T) {
      if (t > first[i]) {
        real lc = a_col + s_col_r * z_col_r[r] + s_col_t * z_col_t[t];
        real le = a_ext + s_ext_r * z_ext_r[r] + s_ext_t * z_ext_t[t];
        vector[2] nw;
        nw[1] = log_sum_exp(la[1] + log1m_inv_logit(lc), la[2] + log_inv_logit(le));
        nw[2] = log_sum_exp(la[1] + log_inv_logit(lc), la[2] + log1m_inv_logit(le));
        la = nw;
      }
      if (cell[i, t] > 0) {
        int k = cell[i, t]; int any = 0; real ll = 0;
        for (v in vstart[k]:vend[k]) {
          ll += bernoulli_logit_lpmf(y[v] | b_p[win[v]] + off_p);
          any += y[v];
        }
        if (any > 0) la[1] = negative_infinity();
        la[2] += ll;
      }
    }
    target += log_sum_exp(la);
  }
}
generated quantities {
  real col_m = inv_logit(a_col);
  real ext_m = inv_logit(a_ext);
  real eps = logit(col_m / (col_m + ext_m));
  real tau = log(-log1m(col_m + ext_m));
  real col_bar; real ext_bar; real eps_bar; real tau_bar;
  {
    real sc = sqrt(square(s_col_r) + square(s_col_t));
    real se = sqrt(square(s_ext_r) + square(s_ext_t));
    real ac = 0; real ae = 0;
    for (m in 1:500) {
      ac += inv_logit(a_col + sc * std_normal_rng());
      ae += inv_logit(a_ext + se * std_normal_rng());
    }
    col_bar = ac / 500; ext_bar = ae / 500;
    eps_bar = logit(col_bar / (col_bar + ext_bar));
    tau_bar = log(-log1m(col_bar + ext_bar));
  }
}
