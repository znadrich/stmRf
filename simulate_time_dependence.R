library(dplyr)
library(tidyr)
library(ggplot2)
library(gganimate)

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
    kappa = params[4]*0,
    delta = params[5]*0
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

map_params <- function(loc, params){
  if(loc == 'c') as.numeric(params['beta'])
  else if (loc %in% c('u', 'd')) as.numeric(params['gamma'])
  else if (loc %in% c('r', 'l')) as.numeric(params['lambda'])
  else if (loc %in% c('dr', 'ul')) as.numeric(params['kappa'])
  else if (loc %in% c('dl', 'ur')) as.numeric(params['delta'])
}

map_param_names <- function(loc){
  inner_func <- function(loc){
    if(loc == 'c') 'beta'
    else if (loc %in% c('u', 'd')) 'gamma'
    else if (loc %in% c('r', 'l')) 'lambda'
    else if (loc %in% c('dr', 'ul')) 'kappa'
    else if (loc %in% c('dl', 'ur')) 'delta'
  }
  
  nm <- sapply(loc, inner_func)
  return(nm)
}

all_param_names <- function(drop_beta=F){
  if(drop_beta) v <- c('alpha', 'gamma', 'lambda', 'kappa', 'delta')
  else v <- c('alpha', 'beta', 'gamma', 'lambda', 'kappa', 'delta')
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
    i = prior_grid$latitude, 
    j = prior_grid$longitude,
    grid_size = grid_size
  ) %>%
    do.call(rbind, .) %>%
    as.data.frame(stringsAsFactors = F) %>%
    lapply(unlist) %>%
    as.data.frame(stringsAsFactors = F)
  
  names(g) <- c('latitude', 'longitude', 'eta')
  
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

eval_cliques <- function(grid_i, prior_grid, grid_size){
  neighbors_prior <- eval_neighbors(prior_grid, grid_size) %>%
    mutate(eta = map_param_names(eta))
  
  if(!is.null(grid_i)){
    grid_i <- grid_i %>%
      select(-t) %>%
      mutate(eta = 'alpha')
    neighbors_prior <- rbind(neighbors_prior, grid_i)
  }
  
  grid <- empty_grid(grid_size) %>%
    left_join(neighbors_prior, by = c('latitude', 'longitude'))
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

grid_cliques <- function(grid_i, prior_grid, grid_size){
  cliques <- eval_cliques(grid_i, prior_grid, grid_size)
  
  grid <- cliques %>% 
    mutate(eta=coalesce(eta, 'none')) %>% 
    group_by(latitude, longitude, eta) %>% 
    summarise(n=n()) %>%
    spread(eta, n)
    
    if(!('alpha' %in% colnames(grid))){
      grid$alpha <- 0
    } 
    
    grid <- grid %>%
      mutate(event=ifelse(alpha == 0, 0, 1)) 
  
  grid[is.na(grid)] <- 0
  grid$alpha <- grid$beta+grid$delta+grid$gamma+grid$kappa+grid$lambda
  # grid <- grid %>%
  #   mutate(
  #     alpha=clique_delta(beta+delta+gamma+kappa+lambda, event),
  #     beta=clique_delta(beta, event),
  #     delta=clique_delta(delta, event),
  #     gamma=clique_delta(gamma, event),
  #     kappa=clique_delta(kappa, event),
  #     lambda=clique_delta(lambda, event),
  #   )
  
  return(grid)
}

pmle <- function(grid_i, prior_grid, grid_size){
  grid <- grid_cliques(grid_i, prior_grid, grid_size)
  log_reg <- glm(
    event ~ beta+delta+gamma+kappa+lambda,
    data = grid,
    family = binomial(link='logit'),
    control = glm.control(maxit = 100)
  )
  
  ple <- log_reg$coefficients
  
  # Get standard error, if the coef blew up then set to 0
  se <- sqrt(diag(vcov(log_reg)))
  ple[se > 1000] <- 0
  return(ple)
}

plot_param <- function(real, estimate, param){
  v_real <- real[, param]
  v_est <- estimate[, param]  
  plot(v_real, type='l', col='red', main=param)
  lines(v_est, col='blue')
  legend(
    'topleft', 
    c('real', 'estimate'), 
    lty=c(1, 1), 
    col=c('red', 'blue'),
    cex=.8
  )
}
