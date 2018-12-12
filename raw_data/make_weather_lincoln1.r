
library(tidyverse)
library(lubridate)
library(readxl)
# "lubridate" %in% (.packages()) # how to check if package loaded

options(warn=1)

# read our weather data
head <- read_tsv('raw_data/Lincoln.txt', col_names=FALSE)[1:6,] # header rows
data <- read_tsv('raw_data/Lincoln.txt', skip=6, col_names=FALSE) # data rows
names(data) <- c('date', 
                 'tmax_C', 'tmin_C', 
                 'radn_MJ_m2', 'sun_hrs', 
                 'windrun_km', 
                 'rain_mm', 'rain_hrs', 
                 'penman_mm', 'priestley_mm', 'open_pan_mm', 
                 'blank1', 'blank2')

# erroneous dates?
data$date <- dmy(data$date)
plot(data$date)

# add extra variables
data <- data %>%
  mutate(station=1, # 
         newdate=date, # parse date
         year=format(year(newdate), scientific=FALSE),
         doy=yday(newdate))

# test for missing rows
# temp1 <- data$newdate[2:nrow(data)] 
# temp2 <- data$newdate[1:(nrow(data)-1)]
temp3 <- as.tibble(c(difftime(data$newdate[2:nrow(data)], data$newdate[1:(nrow(data)-1)], units='days'), 1))
missing <- which(temp3 != 1)
if (length(missing)>0) {
  print('Error: Illegal step sizes in time series (should all be 1)')
  print(missing + 6L)
  print(temp3$value[missing])
  # print('Truncating time series')
  # data <- data[1:missing[1], ]
} else {  
  print('Confirmed: No missing steps in time series')
}

# get variables for BASGRA weather file
data2 <- data %>% 
  dplyr::select(station, year, doy, tmin_C, tmax_C, rain_mm, radn_MJ_m2, priestley_mm)

# simple handling of missing data
missing <- which(data2$tmin_C == -999)
print('Warning: Imputing data2$tmin_C == -999')
print(missing)
for (i in missing) {
  data2$tmin_C[i] <- data2$tmin_C[i-1]
}
missing <- which(data2$tmax_C == -999)
print('Warning: Imputing data2$tmax_C == -999')
print(missing)
for (i in missing) {
  data2$tmax_C[i] <- data2$tmax_C[i-1]
}
missing <- which(data2$rain_mm < 0)
print('Warning: Imputing data2$rain_mm < 0')
print(missing)
for (i in missing) {
  data2$rain_mm[i] <- data2$rain_mm[i-1]
}
missing <- which(data2$radn_MJ_m2 < 0)
print('Warning: Imputing data2$radn_MJ_m2 < 0')
print(missing)
for (i in missing) {
  data2$radn_MJ_m2[i] <- data2$radn_MJ_m2[i-1]
}
missing <- which(data2$priestley_mm < 0)
print('Warning: Imputing data2$priestley_mm < 0')
print(missing)
for (i in missing) {
  data2$priestley_mm[i] <- data2$priestley_mm[i-1]
}

# get irrigation data and add to rain
print('Reading and adding irrigation')
irrig <- read_excel('raw_data/Canterbury irrigation 2011-2017 Simon.xlsx', col_names=TRUE)[,1:4]
names(irrig) <- c("year", "month", "date", "irrig_mm")
irrig <- irrig %>% 
  mutate(year=format(year(date), scientific=FALSE),
         doy = yday(date))
data3 <- data2 %>% 
  left_join(irrig, by=c("year", "doy")) %>%
  select(-month, -date)
data3 <- data3 %>% 
  mutate(rain_mm = if_else(is.na(irrig_mm), rain_mm, rain_mm + irrig_mm)) %>% 
  select(-irrig_mm)

# write BASGRA weather file
write_tsv(data3, 'raw_data/weather_Lincoln.txt')
print('Remember to copy weather to scenario folder!')

