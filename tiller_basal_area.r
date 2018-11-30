# can we estimate basal area from variation in tiller count data?

library(tidyverse)

raw <- c(800,5000,8000,9000,12000)
raw <- runif(5)*12000
data <- data_frame(y=raw) %>% 
  arrange(desc(y)) %>% 
  mutate(
    xmin=seq(0.0, 0.8, 0.2),
    xmax=seq(0.2, 1.0, 0.2),
    ymean=mean(y),
    sd=sd(y),
    yes= y > ymean-2000,
    prop=sum(yes)/n(),
    ymean2=if_else(yes, ymean/prop, 0)
  )


p1 <- ggplot() +
  geom_segment(data=data, mapping=aes(x=xmin, y=y, xend=xmax, yend=y)) +
  geom_segment(data=data, mapping=aes(x=0, y=ymean, xend=1, yend=ymean), colour="blue") +
  geom_segment(data=data, mapping=aes(x=xmin, y=ymean2, xend=xmax, yend=ymean2), colour="red") 
  
print(p1)

area <- function(p, data){ # constant and proportion
  sum(case_when(
    data$xmax<=p[2] ~ (data$xmax-data$xmin)*(data$y-p[1])^2, 
    data$xmin>=p[2] ~ (data$xmax-data$xmin)*(data$y-0)^2,
    TRUE ~           (p[2]-data$xmin)*(data$y-p[1])^2 + 
      (data$xmax-p[2])*(data$y-0)^2
  ))
}

p0 <- c(mean(data$y), 0.5)
# res <- nlm(area, p0, data) 
# yc <- res$estimate[1]
# pc <- res$estimate[2]
res <- optim(p0, area, gr=NULL, data, method="Nelder-Mead")
yc <- res$par[1]
pc <- max(0,min(1,res$par[2]))

p2 <- ggplot() +
  geom_segment(data=data, mapping=aes(x=xmin, y=y, xend=xmax, yend=y)) +
  geom_segment(mapping=aes(x=0, y=yc, xend=pc, yend=yc), colour="red") +
  geom_segment(mapping=aes(x=pc, y=yc, xend=pc, yend=0), colour="red") +
  geom_segment(mapping=aes(x=pc, y=0, xend=1, yend=0), colour="red")

print(p2)

