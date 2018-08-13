
library(reprex)

library(ggplot2)
ggplot(mtcars) +
  geom_col(aes(x=disp, y=mpg)) +
  geom_point(aes(x=disp, y=mpg)) 

# this is the reprex
# copy code to clipboard then type reprex() at console
# reprex code will be in clipboard and can be pasted to e.g. Github
# https://github.com/tidyverse/ggplot2/issues/2773

