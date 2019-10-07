library(dplyr)
library(tidyr)
library(ggplot2)
library(gganimate)
library(dtplyr)
library(data.table)

generate_params <- function(p_interact = .05, prev_interact = F){
  if(prev_interact) {
    p_interact <- .8
  }
  
  if(runif(1) <= p_interact){
    interact <- T
    params <- rgamma(5, 1, .5)
  } else {    
    interact <- F
    params <- rgamma(5, 1, 20)
  }

  neighborhood_params <- list(
    beta = params[1]*0,
    gamma = params[2],
    lambda = params[3],
    kappa = params[4],
    delta = params[5]
  )
  
  return_list <- list(params=neighborhood_params, interact=interact)
  return(return_list)
}

get_param_maginitude <- function(params){
  magnitude <- sum(as.numeric(params))
  
  return(magnitude)
}

generate_pixels <- function(alpha, size){
  p <- runif(size)
  x <- as.integer(p < boot::inv.logit(alpha))
  return(x)
}

generate_pixels_i_j <- function(theta_i_j, size){
  p <- runif(size)
  x <- as.integer(p < boot::inv.logit(theta_i_j))
  return(x)
}

empty_grid <- function(grid_size){
  grid <- expand.grid(
    latitude = round(seq(0, 1, 1/grid_size), 2),
    longitude = round(seq(0, 1, 1/grid_size), 2)
  )
  
  return(grid)
}

generate_grid_init <- function(alpha, t, grid_size = 100){
  grid <- empty_grid(grid_size)
  
  pixels <- generate_pixels(alpha, grid_size**2)
  has_event <- which(pixels == 1)
  
  grid <- grid[has_event, ]
  
  if(nrow(grid) == 0){    
    grid <- empty_grid(grid_size)
    grid <- grid[sample(nrow(grid), 1), ]
  }
  
  grid$t <- t
  return(grid)
}

map_params <- function(loc, params, directional=F){
  if (directional){
    if(loc == 'c') as.numeric(params['beta'])
    else if (loc == 'd') as.numeric(params['gamma_d'])
    else if (loc == 'u') as.numeric(params['gamma_u'])
    else if (loc == 'l') as.numeric(params['lambda_l'])
    else if (loc == 'r') as.numeric(params['lambda_r'])
    else if (loc == 'dr') as.numeric(params['kappa_dr'])
    else if (loc =='ul') as.numeric(params['kappa_ul'])
    else if (loc == 'ur') as.numeric(params['delta_ur'])
    else if (loc == 'dl') as.numeric(params['delta_dl'])
  } else {
    if(loc == 'c') as.numeric(params['beta'])
    else if (loc %in% c('u', 'd')) as.numeric(params['gamma'])
    else if (loc %in% c('r', 'l')) as.numeric(params['lambda'])
    else if (loc %in% c('dr', 'ul')) as.numeric(params['kappa'])
    else if (loc %in% c('dl', 'ur')) as.numeric(params['delta'])
  }
  
}

map_param_names <- function(loc, directional=F){
  if (directional){
    inner_func <- function(loc){
      if(loc == 'c') 'beta'
      else if (loc == 'd') 'gamma_d'
      else if (loc == 'u') 'gamma_u'
      else if (loc == 'l') 'lambda_l'
      else if (loc == 'r') 'lambda_r'
      else if (loc == 'dr') 'kappa_dr'
      else if (loc =='ul') 'kappa_ul'
      else if (loc == 'ur') 'delta_ur'
      else if (loc == 'dl') 'delta_dl'
    }
  } else {
    inner_func <- function(loc){
      if(loc == 'c') 'beta'
      else if (loc %in% c('u', 'd')) 'gamma'
      else if (loc %in% c('r', 'l')) 'lambda'
      else if (loc %in% c('dr', 'ul')) 'kappa'
      else if (loc %in% c('dl', 'ur')) 'delta'
    }
  } 
  
  nm <- sapply(loc, inner_func)
  return(nm)
}

all_param_names <- function(drop_alpha=F, drop_beta=F, directional=F){
  if (directional){
    directional_params <- c(
      'gamma_d', 'gamma_u', 
      'lambda_l', 'lambda_r', 
      'kappa_dr', 'kappa_ul', 
      'delta_ur', 'delta_dl'  
    )

    v <- c('alpha', 'beta', directional_params)
  } else {
    v <- c('alpha', 'beta', 'gamma', 'lambda', 'kappa', 'delta')
  }
  
  if (drop_alpha) v <- v[v != 'alpha']
  if (drop_beta) v <- v[v != 'beta']

  return(v)
}

decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

eval_neighbors <- function(prior_grid, grid_size){
  g <- mapply(
    FUN = gen_neighbors,
    i = prior_grid$longitude, 
    j = prior_grid$latitude,
    grid_size = grid_size
  ) %>%
    do.call(rbind, .) %>%
    as.data.frame(stringsAsFactors = F) %>%
    lapply(unlist) %>%
    as.data.frame(stringsAsFactors = F)
  
  names(g) <- c('longitude', 'latitude', 'eta')
  
  g$latitude <- round(g$latitude, decimalplaces(1/grid_size))
  g$longitude <- round(g$longitude, decimalplaces(1/grid_size))

  return(g)
}

gen_neighbors <- function(i, j, grid_size){
  size <- 1/grid_size
  eta <- list(
    ul = c(-1, 1), u = c(0, 1), ur = c(1, 1),
    l = c(-1, 0), c = c(0, 0), r = c(1, 0),
    dl = c(-1, -1), d = c(0, -1), dr = c(1, -1)
  ) %>% lapply(function(x) x*size)
  
  f <- function(t, size, eta, i, j){
    x <- eta[[t]]
    list(as.numeric(x[1] + i), as.numeric(x[2] + j), names(eta)[t])
  }
  
  neighbors <- lapply(
    seq_along(eta), 
    f,
    size = size,
    eta = eta,
    i = i,
    j = j
  )
  
  return(neighbors)
}

has_neighbors <- function(grid, prior_grid){
  prior_grid$c <- 1
  grid_pairs <- grid %>%
    left_join(prior_grid, by = c('latitude', 'longitude'))
  x <- which(grid_pairs$c == 1)
  return(x)
}

which_neighbors <- function(grid, neighbors_prior){
  index_has_neighbors <- has_neighbors(grid, neighbors_prior)
  with_n_type <- grid[index_has_neighbors, ] %>%
    inner_join(neighbors_prior, by = c('latitude', 'longitude'))
  return(with_n_type)
}

get_params <- function(i, neighbors, params){
  x <- neighbors[i, ]
  x$param_val <- map_params(x$eta, params)
  return(x)
}

neighbor_hood_calculations <- function(grid, prior_grid, grid_size, neighborhood_params){
  neighbors_prior <- eval_neighbors(prior_grid, grid_size)
  neighbors <- which_neighbors(grid, neighbors_prior) %>%
    lapply(1:nrow(.), get_params, neighbors = ., params = neighborhood_params) %>%
    do.call(rbind, .)
  
  neighbors <- neighbors %>%
    group_by(latitude, longitude) %>%
    summarize(eta_i_j = sum(param_val)) %>%
    ungroup
  return(neighbors)
}

generate_grid_main <- function(alpha, prior_grid, t, neighborhood_params, grid_size = 100){
  grid <- empty_grid(grid_size)
  
  if(nrow(prior_grid) > 0){
    neighbors <- neighbor_hood_calculations(grid, prior_grid, grid_size, neighborhood_params)
    theta_i_j <- grid %>%
      left_join(neighbors, by = c('latitude', 'longitude')) %>%
      mutate(theta = coalesce(eta_i_j, 0) + alpha) %>%
      select(theta) %>%
      unlist
    
    pixels <- generate_pixels_i_j(theta_i_j, nrow(grid))
  } else {
    pixels <- generate_pixels(alpha, nrow(grid))
  }
  
  has_event <- which(pixels == 1)
  
  grid <- grid[has_event, ]
  
  if(nrow(grid) == 0){    
    grid <- empty_grid(grid_size)
    grid <- grid[sample(nrow(grid), 1), ]
  }
  
  grid$t <- t
  return(grid)
}

get_prior_grid <- function(x, t){
  prior_grid <- x[x$t == t-1, 1:2]
  if(is.null(prior_grid)){
    prior_grid <- data.frame()
  }
  
  return(prior_grid)
}

update_data <- function(x_i, x, params, interact){
  if(!is.null(x_i)){
    x_i$magnitude <- ifelse(interact, "high", "low")
    x <- rbind(x, x_i)
  }
  
  return(x)
}

eval_cliques <- function(grid_i, prior_grid, grid_size, directional=F){
  neighbors_prior <- eval_neighbors(prior_grid, grid_size) %>%
    mutate(eta = map_param_names(eta, directional)) %>%
    lazy_dt()
  
  if(!is.null(grid_i)){
    grid_i <- grid_i %>%
      select(-t) %>%
      mutate(eta = 'alpha')
    neighbors_prior <- neighbors_prior %>%
      as.data.frame() %>%
      rbind(grid_i) %>%
      lazy_dt()
  }
  
  grid <- empty_grid(grid_size) %>%
    lazy_dt() %>%
    left_join(neighbors_prior, by = c('latitude', 'longitude'))

  return(grid)
}

n_clique <- function(cliques){
  cliques %>% 
    filter(!is.na(eta)) %>%
    group_by(eta) %>% 
    summarize(n=n()) %>%
    return()
}

clique_delta <- function(param, event){
  d <- (1-event)*param - event*param
  return(d)
}

grid_cliques <- function(grid_i, prior_grid, grid_size, directional=F){
  cliques <- eval_cliques(grid_i, prior_grid, grid_size, directional)
  
  grid <- cliques %>% 
    mutate(eta=coalesce(eta, 'none')) %>% 
    group_by(latitude, longitude, eta) %>% 
    summarise(n=n()) %>%
    as.data.frame() %>%
    spread(eta, n)
    
  if(!('alpha' %in% colnames(grid))){
    grid$alpha <- 0
  } 
  
  grid$event <- ifelse(grid$alpha == 0, 0, 1)
  
  grid[is.na(grid)] <- 0

  missing_params <- all_param_names(directional=directional)
  missing_params <- missing_params[!(missing_params %in% colnames(grid))]
  for (p in missing_params){
    grid[, p] <- 0
  }

  if (directional){
    grid$alpha <- 
        grid$beta +
        grid$gamma_d +
        grid$gamma_u +
        grid$lambda_l +
        grid$lambda_r +
        grid$kappa_dr +
        grid$kappa_ul +
        grid$delta_ur +
        grid$delta_dl
  } else {
    grid$alpha <- grid$beta+grid$delta+grid$gamma+grid$kappa+grid$lambda
  }
  
  return(grid)
}

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

get_lik <- function(glm, df){
  p <- predict(glm, df, type='response')
  log_lik <- sum(log(p*df$event + (1-p)*(1-df$event)))
  return(log_lik)
}

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

plot_param <- function(estimate, param, dates){
  colors = c('red', 'blue', 'green', 'cyan')
  dates <- as.factor(dates)
  ix <- estimate$t
  v_est <- estimate[, param]
  lvls <- levels(dates)
  plot(ix[dates == lvls[1]], v_est[dates == lvls[1]],
       col=colors[1], type='l',
       xlim = range(ix), ylim = range(v_est),
       main=param, xlab='index', ylab='estimate')
  for(i in 2:length(lvls)){
    lines(ix[dates == lvls[i]], v_est[dates == lvls[i]], col=colors[i])
  }
  if(param == 'alpha') legend('topleft', legend = lvls, col = colors, lty=rep(1, length(lvls)))
}
