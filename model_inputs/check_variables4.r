# BASGRA Parameter check

library(matrixStats)
library(tidyverse)
library(kableExtra)

outfile <- "check_variables4"

# my parameter files
my_par <- arrange(read_tsv("model_inputs/parameters_Scott.txt"), PARAMETER)
my_bc <- arrange(read_tsv("model_inputs/parameters_BC_Scott.txt", skip=1, col_names=FALSE), X1)
names(my_bc) <- c("PARAMETER", "MIN", "MODE", "MAX", "SHAPE", "SITES")

# check against fortran parameter list (order matters)
my_par_raw <- read_tsv("model_inputs/parameters_Scott.txt")
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

print("Consistency of my_par and my_bc")
errors <- left_join(my_bc, my_par) %>%
  filter(!((MIN<=Scott)&(Scott<=MAX)))
print(errors)

print("Consistency of orig_bc")
errors <- orig_bc %>%
  filter(!((MIN<=MODE)&(MODE<=MAX)))
print(errors)

print("Consistency of orig_par and orig_bc")
errors <- left_join(orig_bc, orig_par) %>%
  select(PARAMETER, MIN, MODE, MAX, MINSITE, MAXSITE) %>%
  filter(!((MIN<=MINSITE)&(MAXSITE<=MAX)))
print(errors)

print("Consistency of new_bc")
errors <- new_bc %>%
  filter(!((MIN<=MODE)&(MODE<=MAX)))
print(errors)

print("Consistency of new_par and new_bc")
errors <- left_join(new_bc, new_par) %>%
  select(PARAMETER, MIN, MODE, MAX, MINSITE, MAXSITE) %>%
  filter(!((MIN<=MINSITE)&(MAXSITE<=MAX)))
print(errors)

# consistency between sets
print("Consistency of my_par and orig_bc")
errors <- left_join(orig_bc, my_par) %>%
  filter(!((MIN<=Scott)&(Scott<=MAX)))
print(errors)

print("Consistency of my_par and new_bc")
errors <- left_join(new_bc, my_par) %>%
  filter(!((MIN<=Scott)&(Scott<=MAX)))
print(errors)

print("Comparison of my_par and orig_par")
errors <- left_join(select(orig_par, PARAMETER, MINSITE, MAXSITE), my_par) %>%
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

print("Comparison of my_par and new_par")
errors <- left_join(select(new_par, PARAMETER, MINSITE, MAXSITE), my_par) %>%
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
# errors <- left_join(select(orig_par, PARAMETER, MINSITE, MAXSITE), my_bc) %>%
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
all <- select(orig_par, PARAMETER, origMINSITE=MINSITE, origMAXSITE=MAXSITE) %>%
  left_join(select(orig_bc, PARAMETER, origMIN=MIN, origMODE=MODE, origMAX=MAX)) %>%
  full_join(select(new_par, PARAMETER, newMINSITE=MINSITE, newMAXSITE=MAXSITE)) %>%
  left_join(select(new_bc, PARAMETER, newMIN=MIN, newMODE=MODE, newMAX=MAX)) %>%
  full_join(select(my_par, PARAMETER, ScottSITE=Scott)) %>%
  left_join(select(my_bc, PARAMETER, ScottMIN=MIN, ScottMODE=MODE, ScottMAX=MAX, ScottSHAPE=SHAPE)) %>%
  left_join(select(my_par, PARAMETER, Units, Description)) %>%
  arrange(PARAMETER)
write_tsv(all, "model_inputs/check_variables.tsv")  
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

all_formatted <- all %>% 
  mutate_at(vars(starts_with("orig")), funs(my_format)) %>%
  mutate_at(vars(starts_with("new")), funs(my_format)) %>%
  mutate_at(vars(starts_with("Scott")), funs(my_format)) 
  
# setwd("model_inputs/")
all_formatted %>%
  kable(col.names=headings) %>%
  # kable_styling(
  #   # bootstrap_options = c("condensed", "striped", "hover"),
  #   full_width=FALSE, position="left", font_size=12) %>%
  add_header_above(header=c("Parameter" = 1, "Original" = 5, "Nutritive" = 5, "Scott" = 5, "Details"=2), 
                   align="left", bold=TRUE) %>%
  save_kable(file=paste("model_inputs/", outfile, ".html", sep=""), 
             self_contained=TRUE) # very slow


