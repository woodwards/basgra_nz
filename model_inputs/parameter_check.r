# BASGRA Parameter check

library(tidyverse)

# official parameter files
par <- arrange(read_tsv("model_inputs/parameters.txt", skip=1, col_names=FALSE), X1)
n <- ncol(par) - 1
names(par) <- c("PARAMETER", paste("SITE", 1:n, sep=""))
bc <- arrange(read_tsv("model_inputs/parameters_BC_Saerheim_nutritive.txt", col_names=FALSE), X1)
names(bc) <- c("PARAMETER", "MIN", "MODE", "MAX", "SITES")

# my parameter files
my_par <- arrange(read_tsv("model_inputs/parameters_Scott.txt"), PARAMETER)
my_bc <- arrange(read_tsv("model_inputs/parameters_BC_Scott.txt", col_names=FALSE), X1)
names(my_bc) <- c("PARAMETER", "MIN", "MODE", "MAX", "SITES")

print("Consistency of official bc")
errors <- bc %>%
  filter(!((MIN<=MODE)&(MODE<=MAX)))
print(errors)

print("Consistency of official par and bc")
errors <- left_join(bc, par) %>%
  filter(!((MIN<=SITE1)&(SITE1<=MAX)))
print(errors)

print("Consistency of my_bc")
errors <- my_bc %>%
  filter(!((MIN<=MODE)&(MODE<=MAX)))
print(errors)

print("Consistency of my_par and my_bc")
errors <- left_join(my_bc, my_par) %>%
  filter(!((MIN<=Scott)&(Scott<=MAX)))
print(errors)

print("Consistency of my_par and official bc")
errors <- left_join(bc, my_par) %>%
  filter(!((MIN<=Scott)&(Scott<=MAX)))
print(errors)

print("Comparison of my_par and official par")
errors <- left_join(select(par, PARAMETER, SITE1), my_par) %>%
  mutate(change=pmax(abs(Scott/SITE1), abs(SITE1/Scott), abs((SITE1-Scott)/SITE1), abs((SITE1-Scott)/Scott))) %>%
  arrange(desc(change))
print(errors, n=20)

print("Comparison of my_bc and official par")
errors <- left_join(select(par, PARAMETER, SITE1), my_bc) %>%
  mutate(change=pmax(abs(MODE/SITE1), abs(SITE1/MODE), abs((SITE1-MODE)/MODE), abs((SITE1-MODE)/MODE))) %>%
  arrange(desc(abs(change)))
print(errors, n=20)

