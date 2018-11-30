
temp <- data_bot2 %>%
  filter(block==3) 

ggplot(temp) +
  geom_bar(mapping=aes(x=date_bot, y=fraction, fill=species), stat='identity') +
  geom_point(mapping=aes(x=date_bot, y=fraction, fill=species), alpha=0.1) +
  facet_grid(cultivar ~ seed_rate ) 
      
temp <- filter(temp, cultivar=="Alto", seed_rate=="18kg")

ggplot(temp) +
  geom_bar(mapping=aes(x=date_bot, y=fraction, fill=species), stat='identity') +
  facet_grid(cultivar ~ seed_rate ) 
