# other outputs

cat(file=stderr(), "Reading checkpoint after BASGRA results", "\n")
load(file=paste(scenario, "/checkpoint_after_results.RData", sep=""))

source("scripts/BC_export_parModes.R")
# source("scripts/BC_plot_parameters_traceplots.R") # very slow! but useful for convergence check
source("scripts/BC_plot_parameters_priorbeta_histograms.R")
# outputMax[which(outputNames=="CLV")] <- 200 # it's pretty awkward to set these programmatically
# outputMax[which(outputNames=="CST")] <- 100
# outputMax[which(outputNames=="TILTOT")] <- 10000
# outputMax[which(outputNames=="WCL")] <- 60
# source("scripts/BC_plot_outputs_data.R")

# run with MAP parameters (can run shorter time window)
s <- 1
for (s in 1:nSites){
  
  # set up
  cat(file=stderr(), paste("Running BASGRA_All with BC MAP parameters, site",s), "\n")
  # dyn.load(BASGRA_DLL) # useful for rerunning this loop
  source(sitesettings_filenames[[s]]) # site initialisation
  # modify simulation length if desired
  # year_start     <- as.integer(2011)
  # doy_start      <- as.integer(244) # 1 September
  # year_stop      <- as.integer(2017)
  # year_stop      <- as.integer(2014)
  # doy_stop       <- as.integer(151) # 31 May
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
  leg <- vector("character", 5)
  list_output <- vector("list", 5)
  
  # run MAP and scenarios
  leg[1] <- "Current"
  list_output[[1]] <- run_model()
  
  # leg[2] <- "Irrig+50mm"
  # matrix_weather <- list_matrix_weather[[s]]
  # ii <- matrix_weather[,2] %in% c(274,305,335,1,32,60,91,121) # irrig days
  # matrix_weather[ii,6] <- matrix_weather[ii,6] + 50.0 # add irrig
  # list_output[[2]] <- run_model()
  # matrix_weather <- list_matrix_weather[[s]]
  
  leg[2] <- "Destock-5%"
  days_harvest <- list_days_harvest[[s]]
  ii <- seq(length(days_harvest)/3*2+1, length(days_harvest))
  days_harvest[ii] <- as.integer(days_harvest[ii] * 0.95) # destock
  list_output[[2]] <- run_model()
  days_harvest <- list_days_harvest[[s]]
  
  leg[3] <- "Destock-10%"
  days_harvest <- list_days_harvest[[s]]
  ii <- seq(length(days_harvest)/3*2+1, length(days_harvest))
  days_harvest[ii] <- as.integer(days_harvest[ii] * 0.9) # destock
  list_output[[3]] <- run_model()
  days_harvest <- list_days_harvest[[s]]
  
  leg[4] <- "Destock-20%"
  days_harvest <- list_days_harvest[[s]]
  ii <- seq(length(days_harvest)/3*2+1, length(days_harvest))
  days_harvest[ii] <- as.integer(days_harvest[ii] * 0.8) # destock
  list_output[[4]] <- run_model()
  days_harvest <- list_days_harvest[[s]]
  
  leg[5] <- "Destock-30%"
  days_harvest <- list_days_harvest[[s]]
  ii <- seq(length(days_harvest)/3*2+1, length(days_harvest))
  days_harvest[ii] <- as.integer(days_harvest[ii] * 0.7) # destock
  list_output[[5]] <- run_model()
  days_harvest <- list_days_harvest[[s]]
  
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
