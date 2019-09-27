setwd('/Users/Zack/Dropbox/Thesis/R')
source('simulate_time_dependence.R')
library(sf)
library(stringr)
library(maps)
library(maptools)

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

get_prior_grid <- function(x, t){
  
  prior_grid <- x[x$t == t-1, c('latitude', 'longitude')]
  if(is.null(prior_grid)){
    prior_grid <- data.frame()
  }
  
  return(prior_grid)
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
       main=param)
  for(i in 2:length(lvls)){
    lines(ix[dates == lvls[i]], v_est[dates == lvls[i]], col=colors[i])
  }
  if(param == 'alpha') legend('topleft', legend = lvls, col = colors, lty=rep(1, length(lvls)))
}

rescale_round <- function(x){
  x <- round(scales::rescale(x, c(0, 1)), 2)
  x <- round(x/.01)*.01
  return(x)
}

df <- as.data.frame(dbf)
df$latitude <- df$LAT
df$longitude <- df$LONG


# Test the function using points in Wisconsin and Oregon.
testPoints <- data.frame(x = df$LONG, y = df$LAT)

df$state <- latlong2state(testPoints)
df <- df[!is.na(df$state), ]
df <- df[df$state == 'california', ]

lat_boundaries <- range(df$latitude)
long_boundaries <- range(df$longitude)

df$og_latitude <- df$latitude
df$og_longitude <- df$longitude
df$latitude <- rescale_round(df$latitude)
df$longitude <- rescale_round(df$longitude)
df$hr <- str_pad(round(df$GMT/100)*100, width = 4, side = 'left', pad = '0')
df <- df[!(df$hr %in% c('0800', '1100', '1900', '2200')), ]

filtered <- df %>% 
  group_by(JULIAN) %>% 
  summarize(n=n()) %>%
  filter(n < 10) %>%
  select(JULIAN) %>% 
  unlist

df$index <- as.numeric(paste(df$JULIAN, df$hr, sep = ''))
df <- df[!(df$JULIAN %in% filtered), ]
df$t <- dense_rank(df$JULIAN)
dates <- unique(df[, c('t', 'DATE')]) %>%
  arrange(t)
dates$notable <- ifelse(month(dates$DATE) %in% c(7,8), "Jul-Aug", ifelse(month(dates$DATE) == 11, 'Nov', ifelse(month(dates$DATE) < 7, 'pre-Jul', 'Sep-Oct')))
dates$notable <- factor(dates$notable, levels = c('pre-Jul', 'Jul-Aug', 'Sep-Oct', 'Nov'))
t_v <- 1:max(df$t)

ple <- data.frame(
  alpha=0,
  beta=0,
  delta=0,
  gamma=0,
  kappa=0,
  lambda=0,
  t=1
)

ple_names <- names(ple)

grid_size <- 100

for(t in t_v[-1]){
  print(t)
  x_i <- unique(df[df$t == t, c('latitude', 'longitude', 't')])
    
  prior_grid <- unique(get_prior_grid(df, t))
  
  cliques <- grid_cliques(grid_i = x_i, prior_grid = prior_grid, grid_size)
  
  ple_i <- pmle(cliques)
  ple_i <- c(ple_i, t)
  
  ple <- rbind(ple, ple_i)
}

names(ple) <- ple_names
par(mfrow=c(2,3))
for(p in all_param_names(drop_beta=F)){
  plot_param(ple, p, dates=dates$notable)
}
par(mfrow=c(1,1))

which(ple$delta > 2)
index_look <- 300
ple[index_look, ]
table(df$DATE[df$t == index_look])
table(df$DATE[df$t == index_look-1])
look <- c(ple$t[index_look], ple$t[index_look-1])
df %>%
  filter(t %in% look) %>%
  mutate(t = as.factor(t)) %>%
  ggplot(aes(x = longitude, y = latitude, color = t)) +
  geom_point() +
  theme_classic() 


df %>%
  filter(t %in% look) %>%
  mutate(t = as.factor(t)) %>%
  ggplot(aes(x = og_longitude, y = og_latitude, color = t)) +
  geom_point() +
  theme_classic() +
  scale_x_continuous(limits = long_boundaries) +
  scale_y_continuous(limits = lat_boundaries)
