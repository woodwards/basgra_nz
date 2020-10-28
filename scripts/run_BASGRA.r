# run_BASGRA
cat(file = stderr(), "Starting run_BASGRA.r", "\n")

# load packages
library(tidyverse)

# Load DLL here (parallel version uses BASGRA package instead)
BASGRA_DLL <- "model/BASGRA_WG.DLL" # use _WG version to use own PET
if (.Machine$sizeof.pointer == 4) {
  print("32-bit R detected")
  print("Warning: Might have memory problems")
  dyn.load(BASGRA_DLL)
} else if (.Machine$sizeof.pointer == 8) {
  print("64-bit R detected")
  # print(paste("Can't load", BASGRA_DLL))
}

#### user options ####
run_parallel <- TRUE
fixpars <- FALSE
fixpar <- c()
single_year <- 2013
base_s <- 1
scenarios <- c("run_mean", "run_mean/scenario_northland", "run_mean/scenario_scott")
parameters <- c("run_mean", "run_mean", "run_mean")
fit_mcmcs <- c(FALSE, FALSE, FALSE)
plot_varsets <- c(1, 2, 2)

# loop here
i <- 1
s <- 2

# options
scenario <- scenarios[i]
fit_mcmc <- fit_mcmcs[i]
plot_varset <- plot_varsets[i]
parameter_location <- parameters[i]
cat(file = stderr(), "\nRunning scenario", scenario, "\n")
cat(file = stderr(), "Using parameters from", parameter_location, "\n")

# initialise run
cat(file = stderr(), "Running initialise_BASGRA_general.R", "\n")
source("scripts/initialise_BASGRA_general.R")

# modify simulation length if desired
year_start <- as.integer(2011)
doy_start <- as.integer(121) # 1 May
year_stop <- as.integer(2017)
doy_stop <- as.integer(120) # 30 April
NDAYS <- as.integer(difftime(
  as.Date(doy_stop - 1, origin = paste(year_stop, 1, 1, sep = "-")),
  as.Date(doy_start - 1, origin = paste(year_start, 1, 1, sep = "-")),
  units = "days"
)) + 1L

# input files
file_weather <- paste(scenario, "/weather_Scott.txt", sep = "")
file_params <- paste(parameter_location, "/BASGRA_parModes.txt", sep = "")
parcol <- 2 + 8*(s-1)
file_harvest <- paste(scenario, "/harvest_Scott_0.txt", sep = "")
days_harvest <- as.matrix(read.table(file = file_harvest, sep = "\t"))

## 3. CREATE HARVEST CALENDAR AND WEATHER INPUT ##
days_harvest <- as.integer(days_harvest)
matrix_weather <- read_weather_WG(year_start, doy_start, NDAYS, file_weather) # function in initialise_BASGRA_general.R

## 4. CREATE VECTOR "PARAMS" ##
df_params <- read.table(file_params, header = T, sep = "\t", row.names = 1)
params <- df_params[, parcol]

## 5. CREATE EMPTY MATRIX y ##
y <- matrix(0, NDAYS, NOUT)

## 2. RUNNING ##
cat(file = stderr(), "Calling run_model()", "\n")
output <- run_model()

# prepare output for plotting
output <- as.data.frame(output)
names(output) <- outputNames
write.csv(output, paste0(scenario, "/sample_output.csv"))
x <- output$Time
plot(x, output$YIELD)
