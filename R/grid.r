

#' @export
empty_grid <- function(grid_size){
  grid <- expand.grid(
    latitude = round(seq(0, 1, 1/grid_size), 2),
    longitude = round(seq(0, 1, 1/grid_size), 2)
  )
  
  return(grid)
}

#' @export
get_prior_grid <- function(x, t){
  
  prior_grid <- x[x$t == t-1, c('latitude', 'longitude')]
  if(is.null(prior_grid)){
    prior_grid <- data.frame()
  }
  
  return(prior_grid)
}