#### Calibrate BASGRA_N using BayesianTools ####
library(BayesianTools)

#### read input data ####
year_start     <- as.integer(2011)
doy_start      <- as.integer(152) # 1 June
year_stop      <- as.integer(2017)
doy_stop       <- as.integer(151) # 31 May
NDAYS          <- as.integer(difftime(
  as.Date(doy_stop -1, origin=paste(year_stop ,1,1,sep='-')), 
  as.Date(doy_start-1, origin=paste(year_start,1,1,sep='-')),
  units="days")) + 1L
days_harvest <- matrix( as.integer(-1), nrow=100, ncol=3 ) # Simon added harv column
days_harvest   <- as.matrix(read.table(file="model_inputs/harvest_Scott.txt", sep='\t'))
days_harvest   <- as.integer(days_harvest)
file_weather   <- 'model_inputs/weather_Scott.txt'
read_weather_WG <- function(y = year_start,
                            d = doy_start,
                            n = NDAYS,
                            f = file_weather) {
  df_weather            <- read.table( f, header=TRUE, sep="\t" )
  colnames(df_weather)  <- c("ST","YR","doy","TMMNI","TMMXI","RAINI","GR","PETI")
  row_start             <- 1
  while( df_weather[row_start,]$YR  < y ) { row_start <- row_start+1 }
  while( df_weather[row_start,]$doy < d ) { row_start <- row_start+1 }
  df_weather_sim        <- df_weather[row_start:(row_start+n-1),]
  NMAXDAYS              <- as.integer(10000)
  NWEATHER              <- as.integer(7)
  matrix_weather        <- matrix( 0., nrow=NMAXDAYS, ncol=NWEATHER )
  matrix_weather[1:n,1] <- df_weather_sim$YR
  matrix_weather[1:n,2] <- df_weather_sim$doy
  matrix_weather[1:n,3] <- df_weather_sim$GR
  matrix_weather[1:n,4] <- df_weather_sim$TMMNI
  matrix_weather[1:n,5] <- df_weather_sim$TMMXI
  matrix_weather[1:n,6] <- df_weather_sim$RAINI
  matrix_weather[1:n,7] <- df_weather_sim$PETI
  return(matrix_weather)
}
matrix_weather <- read_weather_WG(year_start,doy_start,NDAYS,file_weather) # function in initialise_BASGRA_general.R

#### identify variables ####
outputNames <- c(
  "Time"      , "year"     , "doy"      , "DAVTMP"    , "CLV"      , "CLVD"     ,
  "YIELD"     , "CRES"     , "CRT"      , "CST"       , "CSTUB"    , "DRYSTOR"  ,
  "Fdepth"    , "LAI"      , "LT50"     , "PRECIP"    , "PHEN"     , "ROOTD"    ,
  "Sdepth"    , "TILG2"    , "TILG1"    , "TILV"      , "WAL"      , "WCL"      ,
  "WAPS"      , "WAS"      , "WETSTOR"  , "DM"        , "RES"      , "LERG"     , 
  "NELLVG"    , "RLEAF"    , "SLA"      , "TILTOT"    , "FRTILG"   , "FRTILG1"  ,
  "FRTILG2"   , "RDRT"     , "VERN"     ,
  "DRAIN"     , "RUNOFF"   , "EVAP"     , "TRAN"      , "LINT" )
easyNames <- c(
  "Time", "Year", "Day of Year", "Av. Temp.", "Leaf C", "Dead Leaf C",
  "Yield", "Reserve C", "Root C", "Stem C", "Stubble C", "Dry Snow",
  "Frost Depth", "LAI", "Leaf Death Temp.", "Rain", "Phen. Stage", "Root Depth",
  "Snow Depth", "Elong. Tillers", "Vernalised Tillers", "Veg. Tillers", "Soil Water", "Soil Water",
  "Surface Ice", "Soil Ice", "Wet Snow", "Herbage Mass", "Reserve C", "Leaf Elong. Rate", 
  "Elong. Leaves", "Rel. Leaf App. Rate", "Spec. Leaf Area", "Total Tillers", "Frac. Repro. Tillers", "Frac. Vern. Tillers",
  "Frac. Elong. Tillers", "Rel. Leaf Death Rate", "Vernalisation",
  "Drainage", "Runoff", "Evaporation", "Transpiration", "Light Intercep." )
outputUnits <- c(
  "(y)"       , "(y)"      , "(d)"      , "(degC)"    , "(g C m-2)", "(g C m-2)",
  "(g DM m-2)", "(g C m-2)", "(g C m-2)", "(g C m-2)" , "(g C m-2)", "(mm)"     ,
  "(m)"       , "(m2 m-2)" , "(degC)"   , "(mm d-1)"  , "(-)"      , "(m)"      ,
  "(m)"       , "(m-2)"    , "(m-2)"    , "(m-2)"     , "(mm)"     , "(%)"     ,
  "(mm)"      , "(mm)"     , "(mm)"     , "(g DM m-2)", "(g g-1)"  , "(m d-1)"  ,
  "(tiller-1)", "(d-1)"    , "(m2 g-1)" , "(m-2)"     , "(-)"      , "(-)"      ,
  "(-)"       , "(d-1)"    , "(-)"      ,
  "(mm d-1)"  , "(mm d-1)" , "(mm d-1)" , "(mm d-1)"  , "(%)"  )
NOUT <- as.integer( length(outputNames) )
chooseNames <- c( 
  "DAVTMP", "PRECIP", "EVAP", "TRAN", "RUNOFF", "DRAIN", "WAL", "WCL",
  "CLV", "CLVD", "CRES", "CRT", "CST", "CSTUB", 
  "LINT", "LAI", "DM", "SLA", "RES", 
  "TILTOT", "TILV", "TILG1", "TILG2", 
  "PHEN", "VERN", "RDRT"
)
  
#### read calibration data ####
sitesettings_filenames <- c("scripts/initialise_BASGRA_Scott.r")
sitedata_filenames     <- c("model_inputs/data_calibration_Scott.txt")
nSites                 <- length(sitedata_filenames)
s <- 1
dataset_all      <- read.table(sitedata_filenames[s],header=F,sep="")
data_year <- dataset_all[, 2]
data_doy <- dataset_all[, 3]
keeps <- which(
  (data_year == year_start & data_doy  >= doy_start) |
    (data_year  > year_start & data_year <  year_stop) |
    (data_year == year_stop  & data_doy  <= doy_stop )
)
dataset_all <- dataset_all[keeps, ]
data_name  [[s]] <-     database[[s]][,1]
data_year  [[s]] <-     database[[s]][,2]
data_doy   [[s]] <-     database[[s]][,3]
data_value [[s]] <-     database[[s]][,4]
data_sd    [[s]] <-     database[[s]][,5]
data_type  [[s]] <-     database[[s]][,6]
data_weight[[s]] <-     database[[s]][,7]
ndata         <- dim(dataset_all)[[1]]
data_index     <- sapply( 1:ndata   [s], function(i)
  which(as.character(outputNames)==data_name[[s]][i]) )

#### read all parameters ####
file_params    <- 'model_inputs/parameters_Scott.txt' # can contain multiple columns
parcol       <- 1 # which one are we going to use? (row names are ignored)
df_params      <- read.table(file_params,header=T,sep="\t",row.names=1)
params         <- df_params[,parcol]

# read calib parameter priors
file_prior    <- "model_inputs/parameters_BC_Scott.txt"

#### define likelihood function ####
BASGRA_DLL <- "Release/BASGRA_WG.DLL" # use _WG version to use own PET
dyn.load(BASGRA_DLL) 
y              <- matrix(0,NDAYS,NOUT)
run_model <- function(p = params,
                      w = matrix_weather,
                      h = days_harvest,
                      n = NDAYS) {
  .Fortran('BASGRA', p,w,h,n, NOUT,matrix(0,n,NOUT))[[6]]
}
output  <- run_model(params,matrix_weather,days_harvest,NDAYS)

#### run MCMC ####
