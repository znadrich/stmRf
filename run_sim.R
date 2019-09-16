setwd('/Users/Zack/Dropbox/Thesis/R')
source('simulate_time_dependence.R')

p <- .005
alpha <- boot::logit(p)
params <- rgamma(5, 10, 2)
neighborhood_params = list(
  beta = params[1],
  gamma = params[2],
  lambda = params[3],
  kappa = params[4],
  delta = params[5]
)

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

grid_size <- 50
x <- generate_grid_init(alpha, t = 1, grid_size = grid_size)

x$magnitude <- ifelse(sum(params) > 5, "high", "low")

ple <- data.frame(
  int=numeric(0),
  beta=numeric(0),
  delta=numeric(0),
  gamma=numeric(0),
  kappa=numeric(0),
  lambda=numeric(0)
)
ple_names <- names(ple)

real_params <- data.frame(
  beta=numeric(0),
  delta=numeric(0),
  gamma=numeric(0),
  kappa=numeric(0),
  lambda=numeric(0)
)

for(t in t_v[-1]){
  print(t)
  prior_grid <- get_prior_grid(x, t)
  
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
  
  ple <- rbind(ple, pmle(x_i, prior_grid, grid_size))
  real_params <- rbind(real_params, as.data.frame(neighborhood_params))
  x <- update_data(x_i, x, params, 2)
}
names(ple) <- ple_names
real_params$alpha <- alpha
plot(param_magnitude, type='l')

par(mfrow=c(2,2))
for(p in all_param_names()[-1]){
  plot_param(real_params, ple, p)
}
par(mfrow=c(1,1))

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
                    enter_length = as.integer(0),
                    exit_length = as.integer(0)) +
  ggtitle('{frame_time}')


fps <- 2
animate(anim, fps=fps)
