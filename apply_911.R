setwd('/Users/Zack/Dropbox/Thesis/R')
source('simulate_time_dependence.R')
get_prior_grid <- function(x, t){
  
  prior_grid <- x[x$t == t-1, c('latitude', 'longitude')]
  if(is.null(prior_grid)){
    prior_grid <- data.frame()
  }
  
  return(prior_grid)
}

plot_param <- function(estimate, param){
  v_est <- estimate[, param]  
  plot(v_est, type='l', col='red', main=param)
}

rescale_round <- function(x){
  x <- round(scales::rescale(x, c(0, 1)), 2)
  x <- round(x/.02)*.02
  return(x)
}

ple <- data.frame(
  alpha=numeric(0),
  beta=numeric(0),
  delta=numeric(0),
  gamma=numeric(0),
  kappa=numeric(0),
  lambda=numeric(0)
)

ple_names <- names(ple)

filtered_df <- filtered_df[filtered_df$latitude >= lat_boundaries[1] & filtered_df$latitude <= lat_boundaries[2], ]
filtered_df <- filtered_df[filtered_df$longitude >= long_boundaries[1] & filtered_df$longitude <= long_boundaries[2], ]
filtered_df$og_latitude <- filtered_df$latitude
filtered_df$og_longitude <- filtered_df$longitude
filtered_df$latitude <- rescale_round(filtered_df$latitude)
filtered_df$longitude <- rescale_round(filtered_df$longitude)
filtered_df$t <- dense_rank(filtered_df$nearest_5_min)
t_v <- 1:max(filtered_df$t)

grid_size <- 50

for(t in t_v[-1]){
  print(t)
  x_i <- filtered_df[filtered_df$t == t, c('latitude', 'longitude', 't')]
  prior_grid <- get_prior_grid(filtered_df, t)
  ple <- rbind(ple, pmle(x_i, prior_grid, grid_size))
}

names(ple) <- ple_names

par(mfrow=c(2,3))
for(p in all_param_names(drop_beta=F)){
  plot_param(ple, p)
}
par(mfrow=c(1,1))

which(ple$delta > 0)+1
look <- c(250, 249)
filtered_df %>%
  filter(t %in% look) %>%
  mutate(t = as.factor(t)) %>%
  ggplot(aes(x = og_latitude, y = og_longitude, color = t)) +
  geom_point() +
  scale_x_continuous(limits = lat_boundaries) +
  scale_y_continuous(limits = long_boundaries) +
  theme_classic() 
