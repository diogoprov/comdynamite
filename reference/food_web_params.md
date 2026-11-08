# Parameters of the three-trophic-level eco-evolutionary model

Collects the parameters of
[`food_web_dynamics()`](https://provetelab.org/comdynamite/reference/food_web_dynamics.md),
named as in chapter 3 of McPeek (2017), with the values of its figure
3.2 for the resource and the consumer and of its figure 3.11
(bidirectional-dependent traits) for the predator. Greek letters are
spelt out and the optima are prefixed with `z`: `zc_R` is the trait
value at which the intrinsic birth rate of the resource is highest,
`zf_N` and `zx_P` the values at which the intrinsic death rates of the
consumer and the predator are lowest.

## Usage

``` r
food_web_params(
  c0 = 1,
  gamma = 0.002,
  zc_R = 10,
  d = 0.02,
  a0 = 0.1,
  b = 0.1,
  h = 0.1,
  f0 = 0.1,
  theta = 0.01,
  zf_N = 4,
  g = 0,
  m0 = 0.05,
  n = 0.1,
  l = 0.1,
  v0 = 0.001,
  w = 0.1,
  u = 0.1,
  x0 = 0.01,
  delta = 0.01,
  zx_P = 1,
  y = 0,
  beta = 5,
  phi = 5,
  psi = 5,
  alpha = 0.1,
  rho = 0.1,
  tau = 0.1,
  eps_R = 20,
  eps_N = 20,
  eta_N = 20,
  eta_P = 20,
  kappa_R = 20,
  kappa_P = 20
)
```

## Arguments

- c0, gamma, zc_R:

  Maximum intrinsic birth rate of the resource, the strength of
  stabilizing selection on it and its optimum.

- d:

  Density-dependent death rate of the resource.

- a0, b, h:

  Maximum attack coefficient, conversion efficiency and handling time of
  the consumer feeding on the resource.

- f0, theta, zf_N:

  Minimum intrinsic death rate of the consumer, the strength of
  stabilizing selection on it and its optimum.

- g:

  Density-dependent death rate of the consumer.

- m0, n, l:

  Maximum attack coefficient, conversion efficiency and handling time of
  the predator feeding on the consumer.

- v0, w, u:

  The same for the predator feeding on the resource, used only with
  `omnivory = TRUE`.

- x0, delta, zx_P:

  Minimum intrinsic death rate of the predator, the strength of
  stabilizing selection on it and its optimum.

- y:

  Density-dependent death rate of the predator.

- beta, phi, psi:

  Widths of the Gaussian attack coefficients of the consumer on the
  resource, the predator on the consumer and the predator on the
  resource.

- alpha, rho, tau:

  Steepnesses of the logistic attack coefficients, in the same order.

- eps_R, eps_N, eta_N, eta_P, kappa_R, kappa_P:

  Half-saturation constants of the Michaelis-Menten attack coefficients:
  of the resource and consumer traits for the consumer on the resource,
  of the consumer and predator traits for the predator on the consumer,
  and of the resource and predator traits for the predator on the
  resource.

## Value

A named list of class `food_web_params`.

## Details

The intrinsic rates are quadratic in the traits (equations 3.12, 3.13
and 3.21): \\c(z_R) = c_0 (1 - \gamma (z_R - zc_R)^2)\\, \\f(z_N) = f_0
(1 + \theta (z_N - zf_N)^2)\\ and \\x(z_P) = x_0 (1 + \delta (z_P -
zx_P)^2)\\. Each feeding link has an attack coefficient with a maximum
(`a0` for the consumer on the resource, `m0` for the predator on the
consumer, `v0` for the predator on the resource) and a parameter that
sets the strength of selection on it, which one depends on the `attack`
form chosen in
[`food_web_dynamics()`](https://provetelab.org/comdynamite/reference/food_web_dynamics.md):

- `"gaussian"` (bidirectional-dependent traits, equation 3.16 and box
  3.1): the attack coefficient is highest when the traits match and
  declines as \\a_0 \exp(-((z_N - z_R) / \beta)^2)\\; the widths are
  `beta`, `phi` and `psi`.

- `"logistic"` (unidirectional-dependent traits, equation 3.15): the
  attack coefficient rises with the difference between the traits of the
  consumer and its prey, as \\a_0 / (1 + \exp(-\alpha (z_N - z_R)))\\;
  the steepnesses are `alpha`, `rho` and `tau`.

- `"michaelis"` (unidirectional-independent traits, equation 3.14): each
  trait raises the attack coefficient on its own, as \\a_0 \\ z_R /
  (\epsilon_R + z_R) \cdot z_N / (\epsilon_N + z_N)\\; the
  half-saturation constants are `eps_R` and `eps_N`, `eta_N` and
  `eta_P`, `kappa_R` and `kappa_P`. Traits must then be positive.

## References

McPeek, M. A. (2017). *Evolutionary Community Ecology*. Princeton
University Press.
[doi:10.23943/princeton/9780691088778.001.0001](https://doi.org/10.23943/princeton/9780691088778.001.0001)

## See also

[`food_web_dynamics()`](https://provetelab.org/comdynamite/reference/food_web_dynamics.md)

## Examples

``` r
p <- food_web_params(c0 = 10, d = 0.2)
p$c0
#> [1] 10
```
