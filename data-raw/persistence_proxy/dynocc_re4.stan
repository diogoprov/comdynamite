// Dynamic occupancy, v4. Detection per visit: fixed effects of survey window (run 1-4), random
// state x window and year x window effects (regional and interannual phenology), plus route,
// stop and observer random effects. Initial occupancy, colonisation and extinction:
// route and stop random effects; colonisation and extinction also year random effects.
// Latent states marginalised with the forward algorithm. Species-level colonisation and
// extinction are averaged over the random-effect distribution (Monte Carlo in generated
// quantities), in addition to the median-unit values.
data {
  int<lower=1> N; int<lower=1> T; int<lower=1> R; int<lower=1> W; int<lower=1> O;
  array[N] int<lower=1, upper=R> route;
  array[N] int<lower=1, upper=T> first;
  array[N, T] int<lower=0> cell;
  int<lower=1> K; array[K] int<lower=1> vstart; array[K] int<lower=0> vend;
  int<lower=1> V; array[V] int<lower=0, upper=1> y; array[V] int<lower=1, upper=W> win;
  array[V] int<lower=1, upper=O> obs;
  int<lower=1> G; array[V] int<lower=1, upper=G> gst; array[V] int<lower=1, upper=T> vyr;
  int<lower=1> M;                          // Monte Carlo draws for the marginal rates
}
parameters {
  real a_psi; real a_col; real a_ext; vector[W] b_p;
  real<lower=0> s_psi_i;
  real<lower=0> s_col_r; real<lower=0> s_ext_r; real<lower=0> s_col_i; real<lower=0> s_ext_i;
  real<lower=0> s_col_t; real<lower=0> s_ext_t;
  real<lower=0> s_p_r; real<lower=0> s_p_i; real<lower=0> s_p_o;
  vector<lower=0>[W] s_sw; real<lower=0> s_yw;
  matrix[G, W] z_sw; matrix[T, W] z_yw;
  vector[N] z_psi_i;
  vector[R] z_col_r; vector[R] z_ext_r; vector[N] z_col_i; vector[N] z_ext_i;
  vector[T] z_col_t; vector[T] z_ext_t;
  vector[R] z_p_r; vector[N] z_p_i; vector[O] z_p_o;
}
model {
  vector[V] eta_p;
  for (v in 1:V) eta_p[v] = b_p[win[v]] + s_sw[win[v]] * z_sw[gst[v], win[v]] + s_yw * z_yw[vyr[v], win[v]];
  eta_p += s_p_o * z_p_o[obs];
  a_psi ~ normal(0, 1.5); a_col ~ normal(-1.5, 1.5); a_ext ~ normal(-1.5, 1.5); b_p ~ normal(0, 1.5);
  s_psi_i ~ normal(0, 1);
  s_col_r ~ normal(0, 1); s_ext_r ~ normal(0, 1); s_col_i ~ normal(0, 1); s_ext_i ~ normal(0, 1);
  s_col_t ~ normal(0, 1); s_ext_t ~ normal(0, 1);
  s_p_r ~ normal(0, 1); s_p_i ~ normal(0, 1); s_p_o ~ normal(0, 1);
  s_sw ~ normal(0, 1.5); s_yw ~ normal(0, 1);
  to_vector(z_sw) ~ std_normal(); to_vector(z_yw) ~ std_normal();
  z_psi_i ~ std_normal();
  z_col_r ~ std_normal(); z_ext_r ~ std_normal(); z_col_i ~ std_normal(); z_ext_i ~ std_normal();
  z_col_t ~ std_normal(); z_ext_t ~ std_normal();
  z_p_r ~ std_normal(); z_p_i ~ std_normal(); z_p_o ~ std_normal();
  for (i in 1:N) {
    int r = route[i];
    real off_p = s_p_r * z_p_r[r] + s_p_i * z_p_i[i];
    real oc = a_col + s_col_r * z_col_r[r] + s_col_i * z_col_i[i];
    real oe = a_ext + s_ext_r * z_ext_r[r] + s_ext_i * z_ext_i[i];
    real lp = a_psi + s_psi_i * z_psi_i[i];
    vector[2] la;
    la[1] = log1m_inv_logit(lp); la[2] = log_inv_logit(lp);
    for (t in first[i]:T) {
      if (t > first[i]) {
        real lc = oc + s_col_t * z_col_t[t];
        real le = oe + s_ext_t * z_ext_t[t];
        vector[2] nw;
        nw[1] = log_sum_exp(la[1] + log1m_inv_logit(lc), la[2] + log_inv_logit(le));
        nw[2] = log_sum_exp(la[1] + log_inv_logit(lc), la[2] + log1m_inv_logit(le));
        la = nw;
      }
      if (cell[i, t] > 0) {
        int k = cell[i, t];
        int any = sum(y[vstart[k]:vend[k]]);
        la[2] += bernoulli_logit_lpmf(y[vstart[k]:vend[k]] | eta_p[vstart[k]:vend[k]] + off_p);
        if (any > 0) la[1] = negative_infinity();
      }
    }
    target += log_sum_exp(la);
  }
}
generated quantities {
  // median unit (all random effects at zero)
  real col_m = inv_logit(a_col);
  real ext_m = inv_logit(a_ext);
  real eps_m = logit(col_m / (col_m + ext_m));
  real tau_m = log(-log1m(col_m + ext_m));
  // averaged over route, stop and year random effects
  real col_bar; real ext_bar; real eps_bar; real tau_bar;
  {
    real sc = sqrt(square(s_col_r) + square(s_col_i) + square(s_col_t));
    real se = sqrt(square(s_ext_r) + square(s_ext_i) + square(s_ext_t));
    real ac = 0; real ae = 0;
    for (m in 1:M) {
      ac += inv_logit(a_col + sc * std_normal_rng());
      ae += inv_logit(a_ext + se * std_normal_rng());
    }
    col_bar = ac / M; ext_bar = ae / M;
    eps_bar = logit(col_bar / (col_bar + ext_bar));
    tau_bar = log(-log1m(col_bar + ext_bar));
  }
}
