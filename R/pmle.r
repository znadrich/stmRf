

#' @export
pmle <- function(cliques, params, return_model = F){
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
  # Not alpha though
  alpha <- ple[1]
  se <- sqrt(diag(vcov(log_reg)))
  ple[se > 100] <- 0
  ple[1] <- alpha
  return(ple)
}

#' @export
pmle_lasso <- function(cliques, params, return_model = F){
  # glmnet throws error due to CV not having enough points in class
  if(sum(cliques$event) < 10){
    ple <- pmle(cliques, params, return_model)
  } else {
    formula <- reformulate(params)
    mm <- model.matrix(formula, cliques)
    cv.lasso <- glmnet(
      x=mm[, -1],
      y=as.factor(cliques$event),
      family='binomial',
      maxit=1000,
      standardize=T
    )
    s <- cv.lasso$lambda[which.min(deviance(cv.lasso))]
    coef_choose <- as.numeric(coef(cv.lasso, s=s))[-1] # drop alpha
    param_choose <- params[which(coef_choose != 0)]
    if(length(param_choose) != 0){
      ple_chosen <- pmle(cliques, param_choose, return_model)
      ple <- rep(0, length(params)+1)
      names(ple) <- c('alpha', params)
      ple[1] <- ple_chosen[1]
      
      for(p in params){
        if(p %in% names(ple_chosen)){
          ple[p] <- ple_chosen[p]
        }
      }
    } else {
      ple <- pmle(cliques, params, return_model)
    }
    
  }
  return(ple)
}

#' @export
pmle_ridge <- function(cliques, params, return_model = F){
  # glmnet throws error due to CV not having enough points in class
  if(sum(cliques$event) < 10){
    ple <- pmle(cliques, params, return_model)
  } else {
    formula <- reformulate(params)
    mm <- model.matrix(formula, cliques)
    cv.lasso <- glmnet(
      x=mm[, -1],
      y=as.factor(cliques$event),
      family='binomial',
      alpha=0,
      maxit=1000,
      standardize=T
    )
    s <- cv.lasso$lambda[which.min(deviance(cv.lasso))]
    ple <- as.numeric(coef(cv.lasso, s=s))
    
  }
  return(ple)
}

#' @export
pseudoliklihood <- function(grid_i, prior_grid, grid_size, full_params, reduced_params){
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

  full_model <- pmle(cliques_full, full_params, return_model = TRUE)
  reduced_model <- pmle(cliques_reduced, reduced_params, return_model = TRUE)

  lr <- reduced_model$deviance - full_model$deviance

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
      gamma_dd=0,
      gamma_uu=0,
      lambda_ll=0,
      lambda_rr=0,
      kappa_drdr=0,
      kappa_ulul=0,
      delta_urur=0,
      delta_dldl=0,
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