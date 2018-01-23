## 1. INITIALISATION ##
  cat(file=stderr(), 'Starting run_BASGRA_scott.r', "\n")
  #rm(list=ls())
  file_params    <- 'model_inputs/parameters_Scott.txt' # can contain multiple columns
  parcol       <- 1 # which one are we going to use? (row names are ignored)
  source('scripts/initialise_BASGRA_Scott_3.r')

## 2. RUNNING ##
  cat(file=stderr(), 'Calling run_model()', "\n")
   output <- run_model()
  
   # run_model <- function(p = params,
   #                       w = matrix_weather,
   #                       h = days_harvest,
   #                       n = NDAYS) {
   #   .Fortran('BASGRA', p,w,h,n, NOUT,matrix(0,n,NOUT))[[6]]
   # }
   
## 3. OUTPUT ##
   cat(file=stderr(), 'Calling export_output()', "\n")
   export_output()
   
   # export_output <- function(
   #   list_output = list(output),
   #   vars        = outputNames[-(1:3)],
   #   #  file_table  = paste( "output/output_", format(Sys.time(),"%H_%M.txt"), sep="" ),
   #   #  file_plot   = paste( "output/plot_", format(Sys.time(),"%H_%M.png"), sep="" ),
   #   file_table  = "output/output_table.txt",
   #   file_plot   = "output/output_plots.png",
   #   leg         = paste( "Run", 1:length(list_output) ),
   #   leg_title   = "LEGEND",
   #   nrow_plot   = ceiling( sqrt((length(vars)+1) * 8/11) ),
   #   ncol_plot   = ceiling( (length(vars)+1)/nrow_plot ),
   #   lty         = rep(1,length(list_output)),
   #   lwd         = rep(3,length(list_output))
   # ) 

   dyn.unload(BASGRA_DLL) 
   