# test facet_grid

library(tidyverse)
library(lubridate)
library(reprex)

# data_bot2 <- read_tsv("raw_data/data_bot2_Scott.tsv")
# data_bot2 <- read_tsv("raw_data/data_bot2_Lincoln.tsv")
# data_bot2 <- read_tsv("raw_data/data_bot2_Northland.tsv")
data_bot2 <- read_tsv("raw_data/data_bot2_Northland_test.txt")

data_bot2 <- data_bot2 %>%
  select(date, fraction, species, cultivar, block, seed_rate) %>% 
  filter(block < 3) %>% 
  mutate(species = factor(species, 
                          levels=c('leaf', 'stem', 'ann', 'og', 'wc', 'weed', 'dead')))

dput(data_bot2) # creates ascii string for data object

# this is the reprex
# copy code to clipboard then type reprex() at console
# reprex code will be in clipboard and can be pasted to e.g. Github
# https://github.com/tidyverse/ggplot2/issues/2773

library(tidyverse)
library(lubridate)

data_bot2 <- structure(list(date = structure(c(1316390400, 1316390400, 1316390400,  
                  1316390400, 1316390400, 1316390400, 1316390400, 1316390400, 1316390400, 
                  1316390400, 1316390400, 1316390400, 1316390400, 1316390400, 1316390400, 
                  1316390400, 1316390400, 1316390400, 1316390400, 1316390400, 1316390400, 
                  1316390400, 1316390400, 1316390400, 1316390400, 1316390400, 1316390400, 
                  1316390400, 1316390400, 1316390400, 1316390400, 1316390400, 1316390400, 
                  1316390400, 1316390400, 1320796800, 1320796800, 1320796800, 1320796800, 
                  1320796800, 1320796800, 1320796800, 1316390400, 1316390400, 1316390400, 
                  1316390400, 1316390400, 1316390400, 1316390400, 1320796800, 1320796800, 
                  1320796800, 1320796800, 1320796800, 1320796800, 1320796800, 1316390400, 
                  1316390400, 1316390400, 1316390400, 1316390400, 1316390400, 1316390400, 
                  1320796800, 1320796800, 1320796800, 1320796800, 1320796800, 1320796800, 
                  1320796800, 1316390400, 1316390400, 1316390400, 1316390400, 1316390400, 
                  1316390400, 1316390400, 1320796800, 1320796800, 1320796800, 1320796800, 
                  1320796800, 1320796800, 1320796800), class = c("POSIXct", "POSIXt"
                  ), tzone = "UTC"), fraction = c(76.85279188, 0, 2.984771574, 
                  3.736040609, 12.40609137, 0, 4.020304569, 84.01579713, 0, 0, 
                  6.110995635, 3.824568697, 0, 6.048638537, 86.14848032, 0, 1.868460389, 
                  3.836571998, 5.605381166, 0, 2.541106129, 78.49196539, 0, 0, 
                  11.42150803, 4.969097651, 0, 5.117428925, 90.9202454, 0, 0, 2.527607362, 
                  4.883435583, 0, 1.668711656, 56.36016139, 0, 28.62603525, 2.951794436, 
                  6.052240391, 0, 6.009768528, 54.99892957, 0, 40.69792336, 2.389210019, 
                  0.464568615, 0, 1.449368444, 60.97046414, 7.243319269, 10.33755274, 
                  8.509142053, 5.977496484, 0, 6.962025316, 78.43915344, 0, 6.084656085, 
                  0, 4.872134039, 0, 10.60405644, 65.43624161, 0, 26.41685309, 
                  2.125279642, 2.218493661, 0, 3.803131991, 74.36273611, 0, 15.22097928, 
                  3.300935265, 2.461947552, 0, 4.653401797, 47.54063753, 15.6217015, 
                  16.38167617, 12.64513405, 1.688832594, 0, 6.122018155), species = structure(c(1L, 
                  2L, 3L, 5L, 4L, 6L, 7L, 1L, 2L, 3L, 5L, 4L, 6L, 7L, 1L, 2L, 3L, 
                  5L, 4L, 6L, 7L, 1L, 2L, 3L, 5L, 4L, 6L, 7L, 1L, 2L, 3L, 5L, 4L, 
                  6L, 7L, 1L, 2L, 3L, 5L, 4L, 6L, 7L, 1L, 2L, 3L, 5L, 4L, 6L, 7L, 
                  1L, 2L, 3L, 5L, 4L, 6L, 7L, 1L, 2L, 3L, 5L, 4L, 6L, 7L, 1L, 2L, 
                  3L, 5L, 4L, 6L, 7L, 1L, 2L, 3L, 5L, 4L, 6L, 7L, 1L, 2L, 3L, 5L, 
                  4L, 6L, 7L), .Label = c("leaf", "stem", "ann", "og", "wc", "weed", 
                  "dead"), class = "factor"), cultivar = c("Alto", "Alto", "Alto", 
                  "Alto", "Alto", "Alto", "Alto", "Commando", "Commando", "Commando", 
                  "Commando", "Commando", "Commando", "Commando", "Alto", "Alto", 
                  "Alto", "Alto", "Alto", "Alto", "Alto", "Commando", "Commando", 
                  "Commando", "Commando", "Commando", "Commando", "Commando", "Alto", 
                  "Alto", "Alto", "Alto", "Alto", "Alto", "Alto", "Alto", "Alto", 
                  "Alto", "Alto", "Alto", "Alto", "Alto", "Commando", "Commando", 
                  "Commando", "Commando", "Commando", "Commando", "Commando", "Commando", 
                  "Commando", "Commando", "Commando", "Commando", "Commando", "Commando", 
                  "Alto", "Alto", "Alto", "Alto", "Alto", "Alto", "Alto", "Alto", 
                  "Alto", "Alto", "Alto", "Alto", "Alto", "Alto", "Commando", "Commando", 
                  "Commando", "Commando", "Commando", "Commando", "Commando", "Commando", 
                  "Commando", "Commando", "Commando", "Commando", "Commando", "Commando"
                  ), block = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                  1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                  2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
                  2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
                  2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
                  2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L), seed_rate = c("12kg", "12kg", 
                  "12kg", "12kg", "12kg", "12kg", "12kg", "12kg", "12kg", "12kg", 
                  "12kg", "12kg", "12kg", "12kg", "18kg", "18kg", "18kg", "18kg", 
                  "18kg", "18kg", "18kg", "18kg", "18kg", "18kg", "18kg", "18kg", 
                  "18kg", "18kg", "12kg", "12kg", "12kg", "12kg", "12kg", "12kg", 
                  "12kg", "12kg", "12kg", "12kg", "12kg", "12kg", "12kg", "12kg", 
                  "12kg", "12kg", "12kg", "12kg", "12kg", "12kg", "12kg", "12kg", 
                  "12kg", "12kg", "12kg", "12kg", "12kg", "12kg", "18kg", "18kg", 
                  "18kg", "18kg", "18kg", "18kg", "18kg", "18kg", "18kg", "18kg", 
                  "18kg", "18kg", "18kg", "18kg", "18kg", "18kg", "18kg", "18kg", 
                  "18kg", "18kg", "18kg", "18kg", "18kg", "18kg", "18kg", "18kg", 
                  "18kg", "18kg")), class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA, 
                  -84L))

ybreaks <- seq(0, 100, 20)
xbreaks <- seq(floor_date(min(data_bot2$date), "years"), 
               ceiling_date(max(data_bot2$date), "years"), by="1 year")

data_bot2 %>% # FIXME not all data getting plotted for some reason??? Text coding??
  split(.$block) %>%
  map(~ggplot(.) +
        labs(x='Date', y='Dry weight %', fill='Species',
             title=paste('Botanical composition (cut to 4cm), Block', unique(.$block))) +
        theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
        geom_bar(mapping=aes(x=date, y=fraction, fill=species), stat='identity') +
        scale_fill_brewer(palette='Paired') +
        facet_grid(cultivar ~ seed_rate ) +
        scale_x_datetime(breaks=xbreaks, labels=year(xbreaks), limits=range(xbreaks)) +
        scale_y_continuous(breaks=ybreaks)
  )

data_bot2 <- data_bot2 %>% 
  filter(block==2, seed_rate=="18kg", cultivar=="Alto")

library(tidyverse)
library(lubridate)

data_bot2 <- structure(list(date = structure(c(1316390400, 1316390400, 1316390400, 
  1316390400, 1316390400, 1316390400, 1316390400), class = c("POSIXct", 
  "POSIXt"), tzone = "UTC"), fraction = c(86.14848032, 0, 1.868460389, 
  3.836571998, 5.605381166, 0, 2.541106129), species = structure(c(1L, 
  2L, 3L, 5L, 4L, 6L, 7L), .Label = c("leaf", "stem", "ann", "og", 
  "wc", "weed", "dead"), class = "factor"), cultivar = c("Alto", 
  "Alto", "Alto", "Alto", "Alto", "Alto", "Alto"), block = c(1L, 
  1L, 1L, 1L, 1L, 1L, 1L), seed_rate = c("18kg", "18kg", "18kg", 
  "18kg", "18kg", "18kg", "18kg")), class = c("tbl_df", "tbl",
                                              "data.frame"), row.names = c(NA, -7L))

xbreaks <- seq(floor_date(min(data_bot2$date), "years"), 
               ceiling_date(max(data_bot2$date), "years"), by="1 year")

ggplot(data_bot2) +
  geom_bar(mapping=aes(x=date, y=fraction, fill=species), stat='identity') +
  scale_x_datetime(breaks=xbreaks, labels=year(xbreaks))   

ggplot(data_bot2) +
  geom_bar(mapping=aes(x=date, y=fraction, fill=species), stat='identity') +
  scale_x_datetime(breaks=xbreaks, labels=year(xbreaks), limits=range(xbreaks)) 
  
