library(dplyr)
library(ggplot2)
library(gganimate)

p <- .0003
alpha <- boot::logit(p)
params <- rgamma(5, 10, 3)
neighborhood_params = list(
  beta = params[1],
  gamma = params[2],
  lambda = params[3],
  kappa = params[4],
  delta = params[5]
)

generate_params <- function(p_interact = .05, prev_param_mag = 0){
  interact <- runif(1)
  if(prev_param_mag > 5) {
    p_interact <- .8
  }
  
  if(interact <= p_interact){
    params <- rgamma(5, 10, 3)
  } else {
    params <- rgamma(5, 1, 5)
  }
  
  return(params)
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
    latitude = seq(0, 1, 1/grid_size),
    longitude = seq(0, 1, 1/grid_size)
  )
  
  return(grid)
}

generate_grid_init <- function(alpha, t, grid_size = 100){
  grid <- empty_grid(grid_size)
  
  pixels <- generate_pixels(alpha, grid_size**2)
  has_event <- which(pixels == 1)
  
  grid <- grid[has_event, ]
  if(nrow(grid) > 0){
    grid$t <- t
    return(grid)
  }
}

map_params <- function(loc, params){
  if(loc == 'c') as.numeric(params['beta'])
  else if (loc %in% c('u', 'd')) as.numeric(params['gamma'])
  else if (loc %in% c('r', 'l')) as.numeric(params['lambda'])
  else if (loc %in% c('dr', 'ul')) as.numeric(params['kappa'])
  else if (loc %in% c('dl', 'ur')) as.numeric(params['delta'])
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
  x <- which(grid$latitude %in% prior_grid$latitude & grid$longitude %in% prior_grid$longitude)
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
  if(nrow(grid) > 0){
    grid$t <- t
    return(grid)
  }
}

t_v <- 1:100

params <- generate_params(p_interact = 0.1)
param_magnitude <- sum(params)
neighborhood_params = list(
  beta = params[1],
  gamma = params[2],
  lambda = params[3],
  kappa = params[4],
  delta = params[5]
)

grid_size <- 100
x <- generate_grid_init(alpha, t = 1, grid_size = grid_size)
x$magnitude <- ifelse(sum(params) > 5, "high", "low")
for(t in t_v[-1]){
  prior_grid <- x[x$t == t-1, 1:2]
  if(is.null(prior_grid)){
    prior_grid <- data.frame()
  }
  
  params <- generate_params(
    p_interact = 0.1, 
    prev_param_mag = tail(param_magnitude, 1)
  )
  
  neighborhood_params = list(
    beta = params[1]*0,
    gamma = params[2],
    lambda = params[3],
    kappa = params[4],
    delta = params[5]
  )
  param_magnitude <- c(param_magnitude, sum(params))
  x_i <- generate_grid_main(alpha, prior_grid, t, neighborhood_params, grid_size = grid_size)
  if(!is.null(x_i)){
    x_i$magnitude <- ifelse(sum(params) > 5, "high", "low")
    x <- rbind(x, x_i)
  }
}

plot(param_magnitude, type = 'l')

anim <- x %>%
  mutate(begin = as.integer(t),
         length = as.integer(1*2),
         exit = as.integer(1)) %>%
  ggplot(aes(x = latitude, y = longitude)) +
  geom_point(aes(color = magnitude)) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
    transition_events(start = begin,
                      end = begin + length,
                      enter_length = as.integer(1),
                      exit_length = as.integer(1)) +
    ggtitle('{frame_time}')


fps <- 2
animate(anim, fps=fps)
