#-----------------------------------
# Estimating Chi-Square distribution
# parameters using R
#-----------------------------------

# load libraries
library(tidyverse)
library(bbmle)

# dfine parameters and grid

df = 1:10
ncp = 1:10
n = runif(10, 250, 500) |> trunc()
param_grid = expand_grid(n = n, df = df, ncp = ncp) # Create a tibble from all combinations of inputs
head(param_grid)

# functions to estimate the parameters of a chisq distribution
# dof
mean_x <- function(x) mean(x)
mean_minus_1 <- function(x) mean(x) - 1
var_div_2 <- function(x) var(x) / 2
length_minus_1 <- function(x) length(x) - 1
# ncp
mean_minus_mean_minus_1 <- function(x) mean(x) - (mean(x) - 1)
ie_mean_minus_var_div_2 <- function(x) ifelse((mean(x) - (var(x) / 2)) < 0, 0, mean(x) - var(x)/2)
ie_optim <- function(x) optim(par = 0,
                              fn = function(ncp) {
                                -sum(dchisq(x, df = var(x)/2, ncp = ncp, log = TRUE))
                              },
                              method = "Brent",
                              lower = 0,
                              upper = 10 * var(x)/2)$par
# both
estimate_chisq_params <- function(data) {
  # Negative log-likelihood function
  negLogLik <- function(df, ncp) {
    -sum(dchisq(data, df = df, ncp = ncp, log = TRUE))
  }
  
  # Initial values (adjust based on your data if necessary)
  start_vals <- list(df = trunc(var(data)/2), ncp = trunc(mean(data)))
  
  # MLE using bbmle
  mle_fit <- bbmle::mle2(negLogLik, start = start_vals)
  # Return estimated parameters as a named vector
  df <- dplyr::tibble(
    est_df = coef(mle_fit)[1],
    est_ncp = coef(mle_fit)[2]
  )
  return(df)
}
safe_estimates <- {
  purrr::possibly(
    estimate_chisq_params,
    otherwise = NA_real_,
    quiet = TRUE
  )
}

# estimation

# Simulate data ----
set.seed(2024)
dff <- param_grid |>
  mutate(x = pmap(pick(everything()), match.fun("rchisq"))) |>
  mutate(
    safe_est_parms = map(x, safe_estimates),
    dfa = map_dbl(x, mean_minus_1),
    dfb = map_dbl(x, var_div_2),
    dfc = map_dbl(x, length_minus_1),
    ncpa = map_dbl(x, mean_minus_mean_minus_1),
    ncpb = map_dbl(x, ie_mean_minus_var_div_2),
    ncpc = map_dbl(x, ie_optim)
  ) |>
  select(-x) |>
  filter(map_lgl(safe_est_parms, ~ any(is.na(.x))) == FALSE) |>
  unnest(cols = safe_est_parms) |>
  mutate(
    dfa_resid = dfa - df,
    dfb_resid = dfb - df,
    dfc_resid = dfc - df,
    dfd_resid = est_df - df,
    ncpa_resid = ncpa - ncp,
    ncpb_resid = ncpb - ncp,
    ncpc_resid = ncpc - ncp,
    ncpd_resid = est_ncp - ncp
  )
glimpse(dff)

#glimpse() is like a transposed version of print(): columns run down the page,
# and data runs across

# visual insights

par(mfrow = c(1, 2))
boxplot(dff$dfa ~ dff$df, main = "mean(x) -1 ~ df")
boxplot(dff$dfa_resid ~ dff$df, main = "mean(x) -1 ~ df Residuals")

par(mfrow = c(1, 1))
par(mfrow = c(1, 2))
boxplot(dff$dfb ~ dff$df, main = "var(x) / 2 ~ df")
boxplot(dff$dfb_resid ~ dff$df, main = "var(x) / 2 ~ df Residuals")

par(mfrow = c(1, 1))
par(mfrow = c(1, 2))
boxplot(dff$dfc ~ dff$df, main = "length(x) - 1 ~ df")
boxplot(dff$dfc_resid ~ dff$df, main = "length(x) - 1 ~ df Residuals")

par(mfrow = c(1, 2))
boxplot(dff$est_df ~ dff$df, main = "negloglik ~ df - Looks Good")
boxplot(dff$dfd_resid ~ dff$df, main = "negloglik ~ df Residuals")

par(mfrow = c(1, 1))
par(mfrow = c(1, 2))
boxplot(dff$ncpb ~ dff$ncp, main = "mean(x) - var(x)/2 ~ nc")
boxplot(dff$ncpb_resid ~ dff$ncp, main = "mean(x) - var(x)/2 ~ ncp Residuals")

par(mfrow = c(1, 1))
par(mfrow = c(1, 2))
boxplot(dff$ncpc ~ dff$ncp, main = "optim ~ ncp")
boxplot(dff$ncpc_resid ~ dff$ncp, main = "optim ~ ncp Residuals")

par(mfrow = c(1, 1))
par(mfrow = c(1, 2))
boxplot(dff$est_ncp ~ dff$ncp, main = "negloglik ~ ncp - Looks Good")
boxplot(dff$ncpd_resid ~ dff$ncp, main = "negloglik ~ ncp Residuals")

# link to R-bloggers article:

https://www.r-bloggers.com/2024/04/estimating-chi-square-distribution-parameters-using-r/?utm_source=phpList&utm_medium=email&utm_campaign=R-bloggers-daily&utm_content=HTML&fbclid=IwZXh0bgNhZW0CMTAAAR3d7GuKRRNYewxzDxBddefFmCLeHY_vv7Ed8W1bv7L-hX9TXeaxy1wb4Jc_aem_ASlpTgGI-gxGU3O7aCzR8hRGVAwvc1NbAlugknjYkmCHGC7tbZIGJv-7g7UKVVYRSsgR-YCG4scfw8gHNbdl6CdQ

#----
# end
#----