library(brms)


# Your data
mydata <- data.frame(
  positive = c(22, 5, ...),
  tested = c(100, 80, ...),
  predictor1 = c(...),
  predictor2 = c(...)
)

# Fit binomial model
fit_binom <- brm(
  positive | trials(tested) ~ predictor1 + predictor2,
  family = binomial(),
  data = mydata
)

# Fit beta-binomial model
fit_bb <- brm(
  positive | trials(tested) ~ predictor1 + predictor2,
  family = beta_binomial(),
  data = mydata
)

# Model comparison
loo(fit_binom, fit_bb)

# Posterior predictive check
pp_check(fit_bb)

# Explore effect of parameters
conditional_effects(fit_bb)