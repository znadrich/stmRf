decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
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

latlong2state <- function(pointsDF) {
  # Prepare SpatialPolygons object with one SpatialPolygon
  # per state (plus DC, minus HI & AK)
  states <- map('state', fill=TRUE, col="transparent", plot=FALSE)
  IDs <- sapply(strsplit(states$names, ":"), function(x) x[1])
  states_sp <- map2SpatialPolygons(states, IDs=IDs,
                                   proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  # Convert pointsDF to a SpatialPoints object 
  pointsSP <- SpatialPoints(pointsDF, 
                            proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  # Use 'over' to get _indices_ of the Polygons object containing each point 
  indices <- over(pointsSP, states_sp)
  
  # Return the state names of the Polygons object containing each point
  stateNames <- sapply(states_sp@polygons, function(x) x@ID)
  stateNames[indices]
}

rescale_round <- function(x){
  x <- round(scales::rescale(x, c(0, 1)), 2)
  x <- round(x/.01)*.01
  return(x)
}
