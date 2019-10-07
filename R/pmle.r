library(dplyr)
library(dtplyr)
library(tidyr)
library(data.table)

#' @export
pmle <- function(cliques, directional=F, return_model=F){
  params <- all_param_names(drop_alpha=T, drop_beta=F, directional=directional)
  formula <- reformulate(params, 'event')
  log_reg <- glm(
    formula,
    data = cliques,
    family = binomial(link='logit'),
    control = glm.control(maxit = 100)
  )
  
  if (return_model) return(log_reg)
  ple <- log_reg$coefficients
  
  # Get standard error, if the coef blew up then set to 0
  se <- sqrt(diag(vcov(log_reg)))
  ple[se > 100] <- 0
  return(ple)
}

#' @export
get_lik <- function(glm, df){
  p <- predict(glm, df, type='response')
  log_lik <- sum(log(p*df$event + (1-p)*(1-df$event)))
  return(log_lik)
}

#' @export
pseudoliklihood <- function(grid_i, prior_grid, grid_size){
  cliques_full <- grid_cliques(
    grid_i, 
    prior_grid, 
    grid_size,
    directional = T
  )

  cliques_reduced <- cliques_full
  cliques_reduced$gamma <- cliques_reduced$gamma_d + cliques_reduced$gamma_u
  cliques_reduced$lambda <- cliques_reduced$lambda_l + cliques_reduced$lambda_r
  cliques_reduced$kappa <- cliques_reduced$kappa_dr + cliques_reduced$kappa_ul
  cliques_reduced$delta <- cliques_reduced$delta_ur + cliques_reduced$delta_dl

  full_model <- pmle(cliques_full, directional=T, return_model=T)
  reduced_model <- pmle(cliques_reduced, directional=F, return_model=T)

  log_lik_full <- get_lik(full_model, cliques_full)
  log_lik_red <- get_lik(reduced_model, cliques_reduced)

  lr <- -2*(log_lik_red - log_lik_full)
  df <- length(full_model$coefficients) - length(reduced_model$coefficients)
  p <- 1-pchisq(lr, df)

  results <- c(
    lr = lr,
    p = p
  )

  return(results)
}

#' @export
ple_df <- function(directional=F){
  if (directional){
    df <- data.frame(
      alpha=0,
      beta=0,
      gamma_d=0,
      gamma_u=0,
      lambda_l=0,
      lambda_r=0,
      kappa_dr=0,
      kappa_ul=0,
      delta_ur=0,
      delta_dl=0,
      t=1
    )
  } else {
    df <- data.frame(
      alpha=0,
      beta=0,
      delta=0,
      gamma=0,
      kappa=0,
      lambda=0,
      t=1
    )
  }
  
  return(df)
}