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