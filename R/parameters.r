

#' @export
get_params <- function(i, neighbors, params, directional=F){
  x <- neighbors[i, ]
  x$param_val <- map_params(x$eta, params, directional=directional)
  return(x)
}

#' @export
all_param_names <- function(drop_alpha=F, drop_beta=F, directional=F){
  if (directional){
    directional_params <- c(
      'gamma_d', 'gamma_u', 'gamma_dd', 'gamma_uu',
      'lambda_l', 'lambda_r', 'lambda_ll', 'lambda_rr',
      'kappa_dr', 'kappa_ul', 'kappa_drdr', 'kappa_ulul',
      'delta_ur', 'delta_dl', 'delta_urur', 'delta_dldl'
    )

    v <- c('alpha', 'beta', directional_params)
  } else {
    v <- c('alpha', 'beta', 'gamma', 'lambda', 'kappa', 'delta')
  }
  
  if (drop_alpha) v <- v[v != 'alpha']
  if (drop_beta) v <- v[v != 'beta']

  return(v)
}

#' @export
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
    else if (loc == 'dd') as.numeric(params['gamma_dd'])
    else if (loc == 'uu') as.numeric(params['gamma_uu'])
    else if (loc == 'll') as.numeric(params['lambda_ll'])
    else if (loc == 'rr') as.numeric(params['lambda_rr'])
    else if (loc == 'drdr') as.numeric(params['kappa_drdr'])
    else if (loc =='ulul') as.numeric(params['kappa_ulul'])
    else if (loc == 'urur') as.numeric(params['delta_urur'])
    else if (loc == 'dldl') as.numeric(params['delta_dldl'])
  } else {
    if(loc == 'c') as.numeric(params['beta'])
    else if (loc %in% c('u', 'd')) as.numeric(params['gamma'])
    else if (loc %in% c('r', 'l')) as.numeric(params['lambda'])
    else if (loc %in% c('dr', 'ul')) as.numeric(params['kappa'])
    else if (loc %in% c('dl', 'ur')) as.numeric(params['delta'])
  }
  
}

#' @export
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
      else if (loc == 'dd') 'gamma_dd'
      else if (loc == 'uu') 'gamma_uu'
      else if (loc == 'll') 'lambda_ll'
      else if (loc == 'rr') 'lambda_rr'
      else if (loc == 'drdr') 'kappa_drdr'
      else if (loc =='ulul') 'kappa_ulul'
      else if (loc == 'urur') 'delta_urur'
      else if (loc == 'dldl') 'delta_dldl'
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

#' @export
get_param_maginitude <- function(params){
  magnitude <- sum(as.numeric(params))
  
  return(magnitude)
}