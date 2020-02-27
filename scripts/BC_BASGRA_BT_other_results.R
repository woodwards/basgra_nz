# other outputs

library(scales)

cat(file=stderr(), "Reading checkpoint after BASGRA results", "\n")
load(file=paste(scenario, "/checkpoint_after_results.RData", sep=""))

# source("scripts/BC_plot_parameters_traceplots.R") # very slow! but useful for convergence check
# source("scripts/BC_plot_parameters_priorbeta_histograms.R") # now using gpplot version
# outputMax[which(outputNames=="CLV")] <- 200 # it's pretty awkward to set these programmatically
# outputMax[which(outputNames=="CST")] <- 100
# outputMax[which(outputNames=="TILTOT")] <- 10000
# outputMax[which(outputNames=="WCL")] <- 60
# source("scripts/BC_plot_outputs_data.R")

# run with MAP parameters (can run shorter time window if desired to see detail detail)
s <- 1
for (s in 1:nSites){
  
  # set up
  cat(file=stderr(), paste("Running BASGRA_All with BC MAP parameters, site",s), "\n")
  # dyn.load(BASGRA_DLL) # useful for rerunning this loop
  source(sitesettings_filenames[[s]]) # site initialisation
  # modify simulation length if desired
  year_start     <- as.integer(2011)
  doy_start      <- as.integer(121) # 1 May
  year_stop      <- as.integer(2014)
  doy_stop       <- as.integer(120) # 30 April
  NDAYS_all <- NDAYS # will need to restore this for other functions
  NDAYS          <- as.integer(difftime(
    as.Date(doy_stop -1, origin=paste(year_stop ,1,1,sep="-")), 
    as.Date(doy_start-1, origin=paste(year_start,1,1,sep="-")),
    units="days")) + 1L
  y              <- matrix(0,NDAYS,NOUT)
  
  # get MAP parameters 
  file_params    <- paste(scenario, "/BASGRA_parModes.txt", sep="") 
  parcol         <- 2 + 8*(s-1)
  df_params      <- read.table(file_params,header=T,sep="\t",row.names=1)
  params         <- df_params[,parcol]
  
  # create output lists
  leg <- vector("character", 1)
  list_output <- vector("list", 1)

  # run MAP and scenarios
  k <- 1
  leg[k] <- "Current"
  list_output[[k]] <- run_model(y = matrix(0, NDAYS, NOUT))
  
  # k <- 2
  # leg[k] <- "+ Water"
  # matrix_weather <- list_matrix_weather[[s]]
  # matrix_weather[,6] <- list_matrix_weather[[3]][,6] # use Lincoln rainfall/irrigation 
  # list_output[[k]] <- run_model(y = matrix(0, NDAYS, NOUT))
  # matrix_weather <- list_matrix_weather[[s]]
  
  k <- 2
  leg[k] <- "+ Grazing"
  days_harvest <- list_days_harvest[[3]] # use Lincoln grazing pressure
  list_output[[k]] <- run_model(y = matrix(0, NDAYS, NOUT))
  days_harvest <- list_days_harvest[[s]]

  # k <- 4
  # leg[k] <- "+ Temp"
  # matrix_weather <- list_matrix_weather[[s]]
  # matrix_weather[,4] <- list_matrix_weather[[3]][,4] # use Lincoln temperature
  # matrix_weather[,5] <- list_matrix_weather[[3]][,5] # use Lincoln temperature 
  # list_output[[k]] <- run_model(y = matrix(0, NDAYS, NOUT))
  # matrix_weather <- list_matrix_weather[[s]]

  # export output using functions in initialise_BASGRA_general.R
  file_table  = paste(scenario, "/basgra_trace_table_", s, ".txt", sep="")
  file_plot   = paste(scenario, "/basgra_trace_plots_", s, ".png", sep="")
  main_title <- paste("SITE ",s," MAP (",region[s],")",sep="") 
  export_output(file_table=file_table, file_plot=file_plot, main_title=main_title,
                list_output=list_output, leg=leg, leg_title="SCENARIOS")
  
  # restore 
  NDAYS <- NDAYS_all
  y     <- matrix(0,NDAYS,NOUT)
  # dyn.unload(BASGRA_DLL) # useful for rerunning this loop 
}

# memory management
cat(file=stderr(), "Saving checkpoint after BASGRA other results", "\n")
save.image(file=paste(scenario, "/checkpoint_after_other_results.RData", sep=""))
