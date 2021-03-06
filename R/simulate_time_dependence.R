#' @export
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

#' @export
generate_pixels <- function(alpha, size){
  p <- runif(size)
  x <- as.integer(p < boot::inv.logit(alpha))
  return(x)
}

#' @export
generate_pixels_i_j <- function(theta_i_j){
  p <- boot::inv.logit(theta_i_j)
  x <- sapply(1:length(p), function(x) rbinom(1, 1, p[x]))
  return(x)
}

#' @export
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

#' @export
generate_grid_main <- function(alpha, prior_grid, neighborhood_params, grid_size = 100, directional=F){
  x_i_tmp <- empty_grid(grid_size)
  x_i_tmp$t <- 0
  grid <- grid_cliques(
    grid_i = x_i_tmp,
    prior_grid = prior_grid,
    grid_size = grid_size,
    directional = directional
  )
  grid <- grid[grid$alpha > 0, ] 
  nms <- all_param_names(drop_alpha = T, drop_beta=F, directional = directional)
  params_v <- unlist(neighborhood_params[nms])
  n_v <- as.matrix(grid[, nms])
  theta_i_j <- (n_v %*% params_v) + alpha
  grid$theta_i_j <- theta_i_j[, 1]
  pixels <- generate_pixels_i_j(theta_i_j = grid$theta_i_j)
  
  has_event <- which(pixels == 1)
  
  grid <- grid[has_event, ]
  
  if(nrow(grid) == 0){
    grid <- empty_grid(grid_size)
    grid <- grid[sample(nrow(grid), 1), ]
  }
  
  return(grid)
}

#' @export
update_data <- function(x_i, x, params, interact){
  if(!is.null(x_i)){
    x_i$magnitude <- ifelse(interact, "high", "low")
    x <- rbind(x, x_i)
  }
  
  return(x)
}