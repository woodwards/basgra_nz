
library(reprex)

df <- prior_df %>% 
  select(-xmap, -ymap) %>% 
  filter(key %in% c("ABASAL", "FGRESSI")) %>% 
  group_by(key) %>% 
  filter(y < max(y)/20)

dput(df)

library(tidyverse)
library(scales)

df <- structure(list(key = c("ABASAL", "ABASAL", "ABASAL", "ABASAL", 
"ABASAL", "ABASAL", "ABASAL", "ABASAL", "ABASAL", "ABASAL", "ABASAL", 
"ABASAL", "ABASAL", "ABASAL", "ABASAL", "ABASAL", "ABASAL", "ABASAL", 
"ABASAL", "ABASAL", "ABASAL", "ABASAL", "FGRESSI", "FGRESSI", 
"FGRESSI", "FGRESSI", "FGRESSI", "FGRESSI", "FGRESSI", "FGRESSI", 
"FGRESSI", "FGRESSI", "FGRESSI", "FGRESSI"), x = c(0.001, 0.0082, 
0.00829, 0.00838, 0.00847, 0.00856, 0.00865, 0.00874, 0.00883, 
0.00892, 0.00901, 0.0091, 0.00919, 0.00928, 0.00937, 0.00946, 
0.00955, 0.00964, 0.00973, 0.00982, 0.00991, 0.01, 0.2, 0.206, 
0.212, 0.218, 0.224, 0.23, 0.77, 0.776, 0.782, 0.788, 0.794, 
0.8), y = c(0, 11.0493287641799, 9.52417004372928, 8.13788565579873, 
6.88592012436097, 5.763323666614, 4.76476082369701, 3.88451964303298, 
3.11652152228274, 2.45433185328381, 1.89117164322408, 1.41993034474676, 
1.03318020495593, 0.723192559340674, 0.481956675150222, 0.301202035928486, 
0.172425447313178, 0.0869252360804243, 0.0358466051992872, 0.0102463214066776, 
0.00119660118123218, 0, 0, 0.00490050000000002, 0.0192080000000001, 
0.0423405000000001, 0.073728, 0.1128125, 0.1128125, 0.073728, 
0.0423405, 0.019208, 0.0049005, 0)), row.names = c(NA, -34L), class = c("grouped_df", 
"tbl_df", "tbl", "data.frame"), groups = structure(list(key = c("ABASAL", 
"FGRESSI"), .rows = list(1:22, 23:34)), row.names = c(NA, -2L
), class = c("tbl_df", "tbl", "data.frame"), .drop = TRUE))

plot1 <- ggplot(data=df) +
  geom_line(mapping=aes(x=x, y=y)) +
  scale_x_continuous(breaks=pretty_breaks(n=2, min.n=2)) +
  facet_wrap(vars(key), scales="free")
print(plot1)
