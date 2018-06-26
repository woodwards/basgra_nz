
library(tidyverse)
library(broom)
library(kableExtra)

# test data
all <- tibble(PARAMETER=c("A","B", "C"), 
              Value1=c(0.0123, 1230, NA), 
              Value2=c(0.0234, 2340, 1.23), 
              Units=c("m", "Pa", "ha"), 
              Description=c("Length", "Pressure", "Area"))

# my formatting function
my_format <- function(x){
  y <- format(x, digits=3, scientific=FALSE, TRIM=TRUE)
  y[is.na(x)] <- ""
  y
}

# format values by row
all_formatted <- all %>%
  `row.names<-`(.$PARAMETER) %>% # set row names for transpose
  select(-PARAMETER, -Units, -Description) %>% # only numeric columns
  t() %>% # transpose
  tidy() %>% # convert to tibble (creates .rownames column)
  modify(my_format) %>% # apply format function to each column of values in place
  `row.names<-`(.$.rownames) %>% # set row names for transpose
  select(-.rownames) %>% # drop rownames column
  t() %>% # transpose
  tidy() %>% # convert to tibble (creates .rownames column)
  select(-.rownames) %>% # drop rownames
  add_column(PARAMETER=all$PARAMETER, .before=1) %>% # add back nonnumeric columns
  add_column(UNITS=all$Units, 
             DESCRIPTION=all$Description)

# print formatted table
all_formatted %>%
  kable() %>%
  kable_styling(
    bootstrap_options = c("condensed", "striped", "hover"),
    full_width=FALSE, position="left", font_size=12) %>%
  save_kable(file="temp.html", self_contained=TRUE) # very slow

