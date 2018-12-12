# BASGRA Parameter check

library(matrixStats)
library(tidyverse)
library(kableExtra)
library(broom)
library(flextable)

outfile <- "check_variables"

# my parameter files
my_par <- arrange(read_tsv(paste(scenario, "/parameters_All.txt", sep="")), PARAMETER)
my_bc <- arrange(read_tsv(paste(scenario, "/parameters_BC.txt", sep=""), skip=1, col_names=FALSE), X1)
names(my_bc) <- c("PARAMETER", "MIN", "MODE", "MAX", "SHAPE", "SITES")

# check against fortran parameter list (order matters)
my_par_raw <- read_tsv(paste(scenario, "/parameters_All.txt", sep=""))
for_par <- tibble(line=readLines("model/set_params.f90")) %>%
  filter(str_detect(line, "=[:blank:]pa\\(")) %>%
  mutate(var=str_extract(line,"[:alnum:]+"),
         pa=str_extract(line,"pa\\([:digit:]+\\)"),
         index=as.numeric(str_sub(pa,4,-2)))
stopifnot(all(my_par_raw$PARAMETER==for_par$var)) # check my_par matches fortran

# check fortan output list
out_var <- read_tsv("model/output_names.tsv")
for_out<- tibble(line=readLines("model/BASGRA.f90")) %>%
  filter(str_detect(line, "y\\(day,[:space:]*[:digit:]+\\)")) %>%
  mutate(yday=str_extract(line,"y\\(day,[:space:]*[:digit:]+\\)"),
         index=as.numeric(str_sub(yday,7,-2)),
         rhs=str_sub(str_extract(line,"=\\s.+"),3,-1),
         var=str_extract(rhs,"[:alnum:]+")
  )
stopifnot(all(out_var$varname==for_out$var)) # check my_par matches fortran

# variables Simon renamed
newnames <- c("K"="KLAI", "SHAPE"="LSHAPE", "gamma"="KTSNOW")

# official parameter files (original version)
orig_par <- read_tsv("model_inputs/parameters_Gri_Saerheim.txt", skip=1, col_names=FALSE)
names(orig_par) <- c("PARAMETER", paste("SITE", 1:(ncol(orig_par)-1), sep=""))
orig_par$MINSITE <- rowMins(as.matrix(orig_par[,2:ncol(orig_par)]))
orig_par$MAXSITE <- rowMaxs(as.matrix(orig_par[,2:ncol(orig_par)]))
i <- match(names(newnames), orig_par$PARAMETER)
orig_par$PARAMETER[i[!is.na(i)]] <- newnames[!is.na(i)]
orig_par <- arrange(orig_par, PARAMETER)
orig_bc <- read_tsv("model_inputs/parameters_BC_Gri_Saerheim.txt", col_names=FALSE)
names(orig_bc) <- c("PARAMETER", "MIN", "MODE", "MAX", "SITES")
i <- match(names(newnames), orig_bc$PARAMETER)
orig_bc$PARAMETER[i[!is.na(i)]] <- newnames[!is.na(i)]
orig_bc <- arrange(orig_bc, PARAMETER)

# official parameter files (nutritive version)
new_par <- read_tsv("model_inputs/parameters_Saerheim_nutritive.txt", skip=1, col_names=FALSE)
names(new_par) <- c("PARAMETER", paste("SITE", 1:(ncol(new_par)-1), sep=""))
new_par$MINSITE <- rowMins(as.matrix(new_par[,2:ncol(new_par)]))
new_par$MAXSITE <- rowMaxs(as.matrix(new_par[,2:ncol(new_par)]))
i <- match(names(newnames), new_par$PARAMETER)
new_par$PARAMETER[i[!is.na(i)]] <- newnames[!is.na(i)]
new_par <- arrange(new_par, PARAMETER)
new_bc <- read_tsv("model_inputs/parameters_BC_Saerheim_nutritive.txt", col_names=FALSE)
names(new_bc) <- c("PARAMETER", "MIN", "MODE", "MAX", "SITES")
i <- match(names(newnames), new_bc$PARAMETER)
new_bc$PARAMETER[i[!is.na(i)]] <- newnames[!is.na(i)]
new_bc <- arrange(new_bc, PARAMETER)

# check internal consistency
print("Consistency of my_bc")
errors <- my_bc %>%
  filter(!((MIN<=MODE)&(MODE<=MAX)))
print(errors)

print("Consistency of my_par (Scott) and my_bc")
errors <- left_join(my_bc, my_par) %>%
  filter(!((MIN<=Scott)&(Scott<=MAX)))
print(errors)

print("Consistency of orig_bc")
errors <- orig_bc %>%
  filter(!((MIN<=MODE)&(MODE<=MAX)))
print(errors)

print("Consistency of orig_par and orig_bc")
errors <- left_join(orig_bc, orig_par) %>%
  dplyr::select(PARAMETER, MIN, MODE, MAX, MINSITE, MAXSITE) %>%
  filter(!((MIN<=MINSITE)&(MAXSITE<=MAX)))
print(errors)

print("Consistency of new_bc")
errors <- new_bc %>%
  filter(!((MIN<=MODE)&(MODE<=MAX)))
print(errors)

print("Consistency of new_par and new_bc")
errors <- left_join(new_bc, new_par) %>%
  dplyr::select(PARAMETER, MIN, MODE, MAX, MINSITE, MAXSITE) %>%
  filter(!((MIN<=MINSITE)&(MAXSITE<=MAX)))
print(errors)

# consistency between sets
print("Consistency of my_par (Scott) and orig_bc")
errors <- left_join(orig_bc, my_par) %>%
  filter(!((MIN<=Scott)&(Scott<=MAX)))
print(errors)

print("Consistency of my_par (Scott) and new_bc")
errors <- left_join(new_bc, my_par) %>%
  filter(!((MIN<=Scott)&(Scott<=MAX)))
print(errors)

print("Comparison of my_par (Scott) and orig_par")
errors <- left_join(dplyr::select(orig_par, PARAMETER, MINSITE, MAXSITE), my_par) %>%
  mutate(change=pmax(abs(Scott/MINSITE), 
                     abs(MINSITE/Scott), 
                     abs((MINSITE-Scott)/MINSITE), 
                     abs((MINSITE-Scott)/Scott),
                     abs(Scott/MAXSITE), 
                     abs(MAXSITE/Scott), 
                     abs((MAXSITE-Scott)/MAXSITE), 
                     abs((MAXSITE-Scott)/Scott)
         )) %>%
  arrange(desc(change))
print(errors, n=100)

print("Comparison of my_par (Scott) and new_par")
errors <- left_join(dplyr::select(new_par, PARAMETER, MINSITE, MAXSITE), my_par) %>%
  mutate(change=pmax(abs(Scott/MINSITE), 
                     abs(MINSITE/Scott), 
                     abs((MINSITE-Scott)/MINSITE), 
                     abs((MINSITE-Scott)/Scott),
                     abs(Scott/MAXSITE), 
                     abs(MAXSITE/Scott), 
                     abs((MAXSITE-Scott)/MAXSITE), 
                     abs((MAXSITE-Scott)/Scott)
  )) %>%
  arrange(desc(change))
print(errors, n=100)

# print("Comparison of my_bc and official orig_par")
# errors <- left_join(dplyr::select(orig_par, PARAMETER, MINSITE, MAXSITE), my_bc) %>%
#   mutate(change=pmax(abs(Scott/MINSITE), 
#                      abs(MINSITE/Scott), 
#                      abs((MINSITE-Scott)/MINSITE), 
#                      abs((MINSITE-Scott)/Scott),
#                      abs(Scott/MAXSITE), 
#                      abs(MAXSITE/Scott), 
#                      abs((MAXSITE-Scott)/MAXSITE), 
#                      abs((MAXSITE-Scott)/Scott)
#   )) %>%
#   arrange(desc(change))
# print(errors, n=20)

# read old parameters
# file_params    <- 'model_inputs/parameters.txt' # can contain multiple columns
# parcol       <- 1 # which one are we going to use? (row names are ignored)
# orig_params      <- read.table(file_params,header=T,sep="\t",row.names=1)

# summary
all <- dplyr::select(arrange(orig_par, PARAMETER), PARAMETER, origMINSITE=MINSITE, origMAXSITE=MAXSITE) %>%
  left_join(dplyr::select(orig_bc, PARAMETER, origMIN=MIN, origMODE=MODE, origMAX=MAX)) %>%
  full_join(dplyr::select(arrange(new_par, PARAMETER), PARAMETER, newMINSITE=MINSITE, newMAXSITE=MAXSITE)) %>%
  left_join(dplyr::select(new_bc, PARAMETER, newMIN=MIN, newMODE=MODE, newMAX=MAX)) %>%
  full_join(dplyr::select(arrange(my_par, PARAMETER), PARAMETER, ScottSITE=Scott)) %>%
  left_join(dplyr::select(my_bc, PARAMETER, ScottMIN=MIN, ScottMODE=MODE, ScottMAX=MAX, ScottSHAPE=SHAPE)) %>%
  left_join(dplyr::select(my_par, PARAMETER, Units, Description)) 
write_tsv(all, paste(scenario, "/check_variables.tsv", sep=""))  
headings <- c("PARAMETER", 
                "MIN", "MAX", "LOWER", "MODE", "UPPER",
                "MIN", "MAX", "LOWER", "MODE", "UPPER",
                "VAL", "LOWER", "MODE", "UPPER", "SHAPE",
                "UNITS", "DESCRIPTION")
                
# https://cran.r-project.org/web/packages/kableExtra/vignettes/awesome_table_in_html.html
my_format <- function(x){
  y <- format(x, digits=3, scientific=FALSE, TRIM=TRUE)
  y[is.na(x)] <- ""
  y
}

# format values by row
all_formatted <- all %>%
  `row.names<-`(.$PARAMETER) %>% # set row names for transpose
  dplyr::select(-PARAMETER, -Units, -Description) %>% # only numeric columns
  t() %>% # transpose
  tidy() %>% # convert to tibble (creates .rownames column)
  modify(my_format) %>% # apply format function to each column of values in place
  `row.names<-`(.$.rownames) %>% # set row names for transpose
  dplyr::select(-.rownames) %>% # drop rownames column
  t() %>% # transpose
  tidy() %>% # convert to tibble (creates .rownames column)
  dplyr::select(-.rownames) %>% # drop rownames
  add_column(PARAMETER=all$PARAMETER, .before=1) %>% # add back nonnumeric columns
  add_column(UNITS=if_else(!is.na(all$Units), all$Units, ""),
             DESCRIPTION=if_else(!is.na(all$Description), all$Description, ""))

this_dir <- getwd() # avoid bug in save_kable
all_formatted %>%
  kable(col.names=headings) %>%
  kable_styling(
    bootstrap_options = c("condensed", "striped", "hover"),
    full_width=FALSE, position="left", font_size=12) %>%
  add_header_above(header=c("Parameter" = 1, "Original" = 5, "Nutritive" = 5, "Scott" = 5, "Details"=2), 
                   align="left", bold=TRUE) %>%
  column_spec(18, width="20cm") %>%
  save_kable(file=paste(this_dir, "/", scenario, "/", outfile, ".html", sep=""), 
             bs_theme="Cosmo", # https://bootswatch.com/
             self_contained=TRUE) # very slow

write_tsv(all_formatted, paste(scenario, "/", outfile, ".tsv", sep=""))
