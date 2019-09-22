library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(lubridate)
library(gganimate)
library(purrr)
library(magick)

extract_coord <- function(location){
  paren_regex <- '\\((.*)\\, (.*)\\)'
  coords <- str_extract(location, paren_regex) %>%
    str_replace("\\(", "") %>%
    str_replace("\\)", "") %>%
    str_split("\\,\\ ")
  return(coords)
}

get_lat_long <- function(location){
  coords <- extract_coord(location)
  lat <- sapply(coords, function(x) as.numeric(x[1]))
  long <- sapply(coords, function(x) as.numeric(x[2]))
  
  coord_df <- data.frame(latitude = lat, longitude = long)
  
  return(coord_df)
}

extract_time <- function(dt){
  hour <- hour(dt)
  min <- str_pad(string = minute(dt), width = 2, side = 'left', "0")
  str <- paste(hour, min, sep=':')
  return(str)
}

plot_title <- function(df){
  dt <- fast_strptime(df$CallDateTime, format="%m/%d/%Y %I:%M:%S %p")
  min_time <- extract_time(min(dt))
  max_time <- extract_time(max(dt))
  date <- unique(df$date)
  title <- paste(date, min_time, max_time)
  
  return(title)
}

include_miss <- function(df){
  rng <- range(df$minute)
  all_mins <- seq(rng[1], rng[2])
  miss_mins <- all_mins[!(all_mins %in% unique(df$minute))]
  miss_df <- data.frame(latitude=0, longitude=0, minute=miss_mins)
  new_df <- rbind(df, miss_df)
  return(new_df)
}

frame_ts <- function(frame_time){
  hours_mins <- (3:26)*60
  hour <- which.min(frame_time >= hours_mins)
  minute <- frame_time-(hour+1)*60 
  hour <- hour-1 %>%
    str_pad(2, 'left', '0')
  minute <- minute %>%
    round() %>%
    str_pad(2, 'left', '0')
    
  ts <- paste(hour, minute, sep=':')
}

frame_title <- function(frame_time){
  begin <- frame_time-5
  end <- frame_time+5
  begin_ts <- frame_ts(begin)
  end_ts <- frame_ts(end)
  title <- paste(begin_ts, "-", end_ts)
  return(title)
}

df <- data.table::fread("C:/Users/Zack/Downloads/911_Police_Calls_for_Service.csv", nThread = 4 , data.table = F)


coords <- get_lat_long(df$Location)
df <- cbind(df, coords)
df <- df[!is.na(df$latitude), ]

df$date <- str_sub(df$CallDateTime, 1, 10)
df$datetime <- fast_strptime(df$CallDateTime, format="%m/%d/%Y %I:%M:%S %p")
df$hour <- hour(df$datetime)
df$minute <- minute(df$datetime)+(df$hour+1)*60

lat_boundaries <- c(
  quantile(df$latitude[df$latitude != 0], .05),
  quantile(df$latitude[df$latitude != 0], .95)
)

long_boundaries <- c(
    quantile(df$longitude[df$longitude != 0], .05),
    quantile(df$longitude[df$longitude != 0], .95)
  )

dates <- '04/18/2015' 

hours <- 0:23

filtered_df <- df %>%
  select(-datetime) %>%
  filter(date %in% dates,
         hour %in% hours,
         latitude != 0) %>%
  mutate(nearest_5_min = 5*floor(minute/5),
         nearest_10_min = 10*floor(minute/10)) 

ts <- 5

plot_df <- filtered_df %>%
  select(latitude, longitude, minute) %>%
  include_miss %>%
  mutate(begin = minute,
         length = ts*2,
         exit = 1)

anim <- plot_df %>%
  ggplot(aes(x=latitude, y=longitude)) +
  geom_point() +
  scale_x_continuous(limits = lat_boundaries) +
  scale_y_continuous(limits = long_boundaries) +
  theme_classic() +
  transition_events(start = begin,
                    end = begin + length,
                    enter_length = 0,
                    exit_length = 0) +
  ggtitle(paste(dates, '{frame_title(frame_time)}'))

fps <- 5

animate(anim, fps=fps)
frame_vars()
