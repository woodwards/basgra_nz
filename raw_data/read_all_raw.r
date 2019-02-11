# read in data set
options(error=NULL)
options(error = function(){.rs.recordTraceback(TRUE)})

# load libraries
library(tidyverse)
library(readxl)
library(stringr)
library(lubridate)

# load useful functions
source("raw_data/utility_functions.r")

#### set up and loop ####

# sites to loop through
sites <- c("Northland", "Scott", "Lincoln")
# sites <- sites[2]
for (site in sites){
print(site)
  
# choose output data
acultivar <- c("Alto", "Commando", "Halo", "Nui") 
acultivar <- acultivar[1]
aseed_rate <- c("6kg", "12kg", "18kg", "24kg", "30kg")
aseed_rate <- aseed_rate[3] 
# calib_start <- ymd("20110401") # period for data weight = 1
calib_start <- ymd("20120401") # period for data weight = 1
calib_end <- ymd("20171231") # period for data weight = 1

# seed rates in order
seed_rate_levels <- c("6kg", "12kg", "18kg", "24kg", "30kg")

#### read and normalise data ####

if (site=="Scott"){
  
  suppressMessages({ # suppress .name_repair messages
    raw_file_name <- "raw_data/FD1004 Data For Modelling.xlsx"
    data_rpm <- read_xlsx(raw_file_name, sheet="Waikato RPM Height data") %>% autosnake() 
    data_cut1 <- read_xlsx(raw_file_name, sheet="Cut Yield Data Year1 2011to12") %>% autosnake() %>% 
      mutate(mean_dm=NA_real_)
    data_cut2 <- read_xlsx(raw_file_name, sheet="Cut Yield Data Year2 onwards") %>% autosnake() %>% 
      mutate(yield_kg_dm_ha=if_else(yield_kg_dm_ha>0, yield_kg_dm_ha, NA_real_))
    data_till <- read_xlsx(raw_file_name, sheet="Tiller density data Waikato") %>% autosnake() 
    data_bot <- read_xlsx(raw_file_name, sheet="Botanical Composition data ", skip=1) %>% autosnake() 
    data_li <- read_xlsx(raw_file_name, sheet="Waikato LightInterception") %>% autosnake() 
    data_sm <- read_xlsx(raw_file_name, sheet="Soil moisture data") %>% autosnake() %>% 
      group_by(date_measured_d, block, cultivar, seed_rate) %>% 
      summarise(soil_moisture = mean(soil_moisture, na.rm=TRUE))
  })
  
} else if (site=="Lincoln"){
  
  raw_file_name <- "raw_data/FD1004_C Data For Modelling.xlsx"
  data_rpm <- read_xlsx(raw_file_name, sheet="Canterbury RPM Heights") %>% autosnake() %>% 
    rename(
      date_pre_rpm_d  = date_pre_rpm,
      date_post_rpm_d = date_post_rpm,
      grazing         = graze_no
    )
  data_cut1 <- read_xlsx(raw_file_name, sheet='Cut Yield Data Year1&2') %>% autosnake() %>% 
    mutate(mean_dm=NA_real_)
  data_cut2 <- read_xlsx(raw_file_name, sheet='Cut Yield Data Year3on Canty') %>% autosnake() %>% 
    mutate(mean_dm=dm_1)
  data_till <- read_xlsx(raw_file_name, sheet='Tiller density Canterbury') %>% autosnake() 
  data_bot <- read_xlsx(raw_file_name, sheet='Botanical Composition Data', skip=1) %>% autosnake() 
  data_sm <- read_xlsx(raw_file_name, sheet='SoilMoist', skip=9, 
                       col_names = c("date_measured_d", "rn", 
                                     "sm5cm_av", "sm10cm_av", 
                                     "sm5cm_max", "sm10cm_max", 
                                     "sm5cm_min", "sm10cm_min")) %>% 
    mutate(
      dom = day(date_measured_d),
      soil_moisture = sm10cm_av * 100 
    ) %>% 
    select(date_measured_d, dom, soil_moisture) %>% 
    filter(dom==1) %>%  # only take one point per month
    mutate(
      cultivar = acultivar,
      seed_rate = aseed_rate,
      block = rep_len(c(2,3,5), n()) # split soil moisture data across blocks 2,3,5
    )

} else if (site=="Northland"){

  raw_file_name <- "raw_data/FD1004_N Data For Modelling.xlsx"
  data_rpm <- read_xlsx(raw_file_name, sheet="Northland RPM Height data") %>% autosnake() 
  data_cut1 <- read_xlsx(raw_file_name, sheet="Northland Yield Data") %>% autosnake() %>%
    mutate(mean_dm=dm) %>% 
    filter(graze_no<=5) # split into two bits to match Scott
  data_cut2 <- read_xlsx(raw_file_name, sheet="Northland Yield Data") %>% autosnake() %>% 
    mutate(mean_dm=dm) %>% 
    filter(graze_no>5)
  data_till <- read_xlsx(raw_file_name, sheet="Tiller density data Northland") %>% autosnake() 
  data_bot <- read_xlsx(raw_file_name, sheet="Botanical Composition data ", skip=1) %>% autosnake() 
  data_sm <- read_xlsx(raw_file_name, sheet="SoilMoist", skip=1, col_names = c("excel_date", "soil_moisture")) %>% 
    mutate(
      date_measured_d = as.POSIXct(floor(excel_date)*(60*60*24), origin="1899-12-30"),
      dom = day(date_measured_d)
    ) %>% 
    group_by(date_measured_d, dom) %>%
    summarise(soil_moisture = mean(soil_moisture, na.rm=TRUE)) %>% 
    ungroup() %>% 
    filter(dom==1) %>%  # %>%  only take one point per month
    mutate(
      cultivar = acultivar,
      seed_rate = aseed_rate,
      block = rep_len(c(2,3,4), n()) # split soil moisture data across blocks 2,3,4
    )  
  
}

#### Rising Plate Meter Data ####

# deal with a few missing values (there are others but they don't affect us at this stage)
data_rpm <- data_rpm %>% 
  arrange(grazing, plot)
i <- which(is.na(data_rpm$pregrazing_mass_kg_dm_ha)) # find missing 
data_rpm$pregrazing_mass_kg_dm_ha[i] <- data_rpm$pregrazing_mass_kg_dm_ha[i-100]*0.5+data_rpm$pregrazing_mass_kg_dm_ha[i+100]*0.5
i <- which(is.na(data_rpm$pregrazing_mass_kg_dm_ha)) # find missing 
data_rpm$pregrazing_mass_kg_dm_ha[i] <- mean(data_rpm$pregrazing_mass_kg_dm_ha[data_rpm$grazing==data_rpm$grazing[i]], na.rm=TRUE)
i <- which(is.na(data_rpm$pregrazing_mass_kg_dm_ha)) # find missing 
stopifnot(length(i)==0)

data_rpm <- data_rpm %>%
  mutate(seed_rate = factor(seed_rate, levels=seed_rate_levels)) %>%
  rename(
       date_pre = date_pre_rpm_d,
       mass_pre = pregrazing_mass_kg_dm_ha,
       date_post = date_post_rpm_d,
       mass_post = postgrazing_mass_kg_dm_ha,
       date_grazed = date_grazed_d
       ) %>%
  select(date_grazed, block, cultivar, seed_rate, grazing, date_pre, mass_pre, date_post, mass_post) %>%
  group_by(grazing) %>% 
  mutate( 
    date_pre = as_date(mean(date_pre)),
    date_post = as_date(mean(date_post)),
    date_grazed = as_date(mean(date_grazed))
  ) %>% 
  ungroup() %>% 
  group_by(block, cultivar, seed_rate, grazing) %>% 
  mutate( 
    samples = n(),
    mass_pre = mean(mass_pre, na.rm=TRUE),
    mass_post = mean(mass_post, na.rm=TRUE)
  ) %>% 
  ungroup()

#### Pasture Cut Data ####

data_cut1 <- data_cut1 %>%
  mutate(seed_rate = factor(seed_rate, levels=seed_rate_levels)) %>%
  rename(
    date_cut = date_d,
    yield_cut = yield_kg_dm_ha,
    dm_cut = mean_dm,
    grazing = graze_no
  ) %>%
  select(date_cut, block, seed_rate, cultivar, yield_cut, dm_cut, grazing) %>% 
  group_by(block, seed_rate, cultivar, grazing) %>%
  summarise(
    samples = n(),
    date_cut = mean(date_cut),
    yield_cut = mean(yield_cut),
    dm_cut = mean(dm_cut)
  ) %>% 
  ungroup()

data_cut2 <- data_cut2 %>%
  mutate(seed_rate = factor(seed_rate, levels=seed_rate_levels)) %>%
  rename(
    date_cut = date_d,
    yield_cut = yield_kg_dm_ha,
    dm_cut = mean_dm,
    grazing = graze_no
  ) %>%
  select(date_cut, block, seed_rate, cultivar, yield_cut, dm_cut, grazing) %>% 
  group_by(block, seed_rate, cultivar, grazing) %>% 
  summarise(
    samples = n(),
    date_cut = mean(date_cut),
    yield_cut = mean(yield_cut),
    dm_cut = mean(dm_cut)
  ) %>% 
  ungroup()

data_cut <- suppressMessages(
  full_join(data_cut1, data_cut2) %>% 
    group_by(grazing) %>% 
    mutate(date_cut = as_date(mean(date_cut, na.rm=TRUE))) %>% 
    ungroup()
  )

# add grazing date from data_rpm
data_cut <- data_cut %>% 
  left_join(
    select(data_rpm, block, seed_rate, cultivar, grazing, date_grazed),
    by = c("block", "cultivar", "seed_rate", "grazing")
  )

#### Tiller Density Data ####

# tiller sampling is not associated with a grazing event
data_till <- data_till %>%
  mutate(seed_rate = factor(seed_rate, levels=seed_rate_levels)) %>%
  rename(
    date_till = date_d,
    tillers = ryegrass_tiller_density_tillers_m2
  ) %>% 
  mutate(
    date_till = as_date(date_till),
    sampling = paste(season, year)
    )

# create separate summary df
data_till_sum <- data_till %>%
  group_by(block, cultivar, seed_rate, sampling) %>%
  summarise(
    samples = n(),
    mean_tillers = mean(tillers),
    sd_tillers = sd(tillers),
    date_till = mean(date_till, na.rm=TRUE),
    tillers = mean_tillers
    ) %>%
  ungroup() %>% 
  group_by(sampling) %>% 
  mutate(
    date_till = as_date(mean(date_till, na.rm=TRUE))
    ) %>% 
  ungroup()

#### Botanical Composition Data ####

# associated with pregrazing cuts
data_bot <- data_bot %>%
  mutate(seed_rate = factor(seed_rate, levels=seed_rate_levels)) %>%
  rename(
    date_bot = date_d,
    leaf = perennial_ryegrass_leaf, # (includes pseudostem) these are percentages, sum = 100
    stem = perennial_ryegrass_reproductive_stem,
    ann = annual_ryegrass,
    wc = white_clover,
    poa = poa_sp,
    ogx = other_grasses_excluding_poa,
    og = other_grasses_including_poa,
    weed = weeds,
    dead = dead
    ) %>%
  mutate(
    date_bot = as_date(date_bot),
    total = round(leaf + stem + ann + wc + og + weed + dead), # should always be 100
    rgfrac =  100 * (leaf + stem) / (leaf + stem + ann + wc + og + weed), # of green
    month = month(dmy(paste("01", month, "2011"))),
    # estimate botanicals composition below cutting height
    # Tozer data - from Ruakura only - tiller dissection only not sward
    # may be ok to apply these assumptions to all sites since residual very similar
    # be careful to ensure it gives reasonable estimates
    x = month %% 12 / 12, # could use yday? but would need to fix regression
    y = -51.56*x^5 + 119.17*x^4 - 93.30*x^3 + 26.71*x^2 - 1.53*x + 1.07,
    leaf_below = leaf / 1.09, # based on Tozer data (March/April difference ignored)
    stem_below = stem / 5.38, # based on Tozer data
    dead_below = (100 - leaf_below - stem_below) * rgfrac / 100 # FIXME Tozer data is for removed tiller not whole sward
  ) %>%
  select(leaf, stem, ann, wc, og, weed, dead, date_bot, total, block, seed_rate, cultivar, rgfrac, leaf_below, stem_below, dead_below) %>% 
  filter(!is.na(total)) 
stopifnot(all(data_bot$total==100))

# add cut info to data_bot 
closest_cut <- function(block_bot, cultivar_bot, seed_rate_bot, date_bot){
  temp <- data_cut %>% 
    filter(block==block_bot, cultivar==cultivar_bot, seed_rate==seed_rate_bot)
  ans <- as.POSIXct(rep(NA, length(date_bot)))
  if (nrow(temp)>0){ # there are some cut dates for this block-cultivar-seed_rate
    for (i in 1:length(date_bot)){
      delay <- as.double(difftime(date_bot[i], temp$date_cut, units="days")) # vector of delays
      min_delay <- min(abs(delay))
      if (min_delay<7){ # only cuts within 7 days of botanical
        ans[i] <- temp$date_cut[which(abs(delay)==min_delay)[[1]]]
      }
    }
  }
  return(ans)
}
  
data_bot <- data_bot %>% 
  group_by(block, cultivar, seed_rate, date_bot) %>%
  mutate(
    date_cut = as_date(closest_cut(block, cultivar, seed_rate, date_bot))
    ) %>% 
  filter(!is.na(date_cut)) %>% 
  ungroup %>% 
  left_join(
    select(data_cut, block, cultivar, seed_rate, grazing, date_grazed, date_cut, yield_cut),
    by = c("block", "cultivar", "seed_rate", "date_cut")
    ) %>% 
  filter(!is.na(date_grazed)) %>% 
  group_by(grazing) %>%
  mutate(
    date_cut_min = min(date_cut),
    date_cut_max = max(date_cut),
    date_cut_diff = interval(date_cut_max, date_cut_min) %/% days(1),
    date_bot_min = min(date_bot),
    date_bot_max = max(date_bot),
    date_bot_diff = interval(date_bot_max, date_bot_min) %/% days(1)
  ) %>% 
  ungroup
i <- data_bot$date_cut>data_bot$date_grazed
if (length(i)>0){
  print("Warning: Shifted Cut Date After Grazing Date")
  data_bot$date_cut[i] <- data_bot$date_grazed[i]
}

# check n samples
temp <- data_bot %>% 
  group_by(block, cultivar, seed_rate, grazing) %>% 
  summarise(
    samples = n()
  ) %>% 
  ungroup()
stopifnot(all(temp$samples==1))

#### harvest model ####

# add cut yield into data_rpm
data_rpm <- data_rpm %>% 
  left_join(
    select(data_cut, block, seed_rate, cultivar, grazing, date_cut, yield_cut),
    by = c("block", "seed_rate", "cultivar", "grazing")
  ) %>% 
  mutate(
    yield_below = mass_post, # simple assumption, quite stable, not actually calibrated though
    yield_total = yield_cut + yield_below, # yield_cut is not available at every grazing, so use this to create a model
    harv0 = 100 * yield_cut * 1.09 / (yield_cut * 1.09 + yield_below) # use this to create a model also
  )

# fill in harv values using lm
harv_lm <- lm(harv0 ~ mass_pre + mass_post, data_rpm)
harv_predict = predict(harv_lm, data_rpm)

# fill in yield values using lm
yield_lm <- lm(yield_total ~ mass_pre + mass_post, data_rpm)
yield_predict = predict(yield_lm, data_rpm)

# put yield and harv models into data_rpm
data_rpm <- data_rpm %>% 
  mutate(
    yield_model = yield_predict,
    yield_above = yield_model - yield_below,
    harv = harv_predict
  )

print(
ggplot(data_rpm) +
  labs(title=paste("Estimated Yield and Harvest for", site)) +
  geom_jitter(aes(x=grazing, y=mass_post), colour="blue") +
  geom_jitter(aes(x=grazing, y=mass_pre), colour="lightblue") +
  geom_jitter(aes(x=grazing, y=yield_total), colour="red", shape=1) +
  geom_jitter(aes(x=grazing, y=yield_model), colour="black", shape=1) +
  geom_jitter(aes(x=grazing, y=harv0*10), colour="lightgreen", shape=1) +
  geom_jitter(aes(x=grazing, y=harv*10), colour="darkgreen", shape=1) +
  geom_abline(slope=0, intercept=0, colour="orange") +
  geom_abline(slope=0, intercept=1000, colour="orange")
)

#### gather data_bot into data_bot2 ####

data_bot2 <- data_bot %>%
  select(leaf, stem, ann, wc, og, weed, dead, block, cultivar, seed_rate, date_bot) %>% 
  gather(leaf, stem, ann, wc, og, weed, dead, key="species", value="fraction") %>%
  mutate(species = factor(species, 
                          levels=c("leaf", "stem", "ann", "og", "wc", "weed", "dead")))

#### calculate botanical mass ####

# join cut and botanical data by grazing
data_bm <- data_bot %>%
  left_join(
    select(data_rpm, block, seed_rate, cultivar, grazing, yield_above, yield_below),
    by = c("block", "seed_rate", "cultivar", "grazing")
  )

#### gather data_bm into data_bm2 ####

data_bm2 <- data_bm %>%
  gather(leaf, stem, ann, wc, og, weed, dead, key="species", value="fraction") %>%
  mutate(
    species = factor(species, levels=c("leaf", "stem", "ann", "og", "wc", "weed", "dead")),
    species_mass = fraction / 100 * yield_above
    ) 

#### Soil Moisture Data ####

data_sm <- data_sm %>%
  mutate(seed_rate = factor(seed_rate, levels=seed_rate_levels)) %>%
  rename(
    date_sm = date_measured_d,
    sm = soil_moisture
  ) %>%
  ungroup() %>% 
  mutate(date_sm = as_date(date_sm)) %>% 
  group_by(block, cultivar, seed_rate, date_sm) %>%
  summarise(mean_sm = mean(sm)) %>% 
  ungroup()

#### make model input files ####

# move post harvest measurements if on grazing day since model does harvest at start of day
# data_rpm$date_post2 <- as.Date(with(data_rpm,
#                             ifelse(date_post > date_grazed, 
#                                    as.Date(date_post), 
#                                    as.Date(date_post) + 1)
#                             ), origin="1970-01-01")

# some calculations
# data_rpm <- data_rpm %>%
#   mutate(
#     year_pre = year(date_pre),
#     doy_pre = yday(date_pre),
#     year_post = year(date_post2),
#     doy_post = yday(date_post2),
#     year_grazed = year(date_grazed),
#     doy_grazed = yday(date_grazed), 
#     days_pre = difftime(date_grazed, date_pre),
#     days_post = difftime(date_post2, date_grazed)
#     )

print("Preparing model input files...")
ablock <- 3 # for testing
for (ablock in c(0, 1,2,3,4,5)) { # loop through blocks (0 = mean)

  # write harvest dates and harvest % for selected series
  harv_file_name <- paste("raw_data/harvest_", site, "_", ablock, ".txt", sep="")
  if (ablock==0){
    data_h <- data_rpm %>%
      filter(seed_rate %in% aseed_rate & cultivar %in% acultivar) %>% 
      group_by(grazing) %>% 
      summarise(
        date_range = max(date_grazed, na.rm=TRUE)-min(date_grazed, na.rm=TRUE),
        date_grazed = mean(date_grazed, na.rm=TRUE),
        harv = mean(harv, na.rm=TRUE)
      )
  } else {
    data_h <- data_rpm %>%
      filter(block==ablock & seed_rate %in% aseed_rate & cultivar %in% acultivar) %>% 
      group_by(grazing) %>% 
      summarise(
        date_range = max(date_grazed, na.rm=TRUE)-min(date_grazed, na.rm=TRUE),
        date_grazed = mean(date_grazed, na.rm=TRUE),
        harv = mean(harv, na.rm=TRUE)
      )
  }
  grazing_dates <- data_h$date_grazed
  stopifnot(!any(is.na(data_h$harv)))
  days_harvest <- matrix(as.integer(-1), nrow=100, ncol=3) # up to 100 harvests
  days_harvest[1:nrow(data_h),] <- c(year(data_h$date_grazed), yday(data_h$date_grazed), data_h$harv)
  write.table(days_harvest, file=harv_file_name, row.names=FALSE, col.names=FALSE, sep="\t")
  
  # write calibration data
  data_file_name <- paste("raw_data/data_calibration_", site, "_", ablock, ".txt", sep="")
  
  # collect the data in this list
  data_c <- vector("list", 10) 
  err_c <- c(TILTOT=2000, CLV=30, CST=10, CLVD=20, CRT=60, WCLM=10, BASAL=10)
  
  # pre and post mass (but this includes other species!)
  # temp <- data_rpm %>%
  #   select(block, seed_rate, cultivar, mass_pre, year_pre, doy_pre, 
  #          mass_post, year_post, doy_post) %>%
  #   filter(block==ablock & seed_rate %in% aseed_rate & cultivar %in% acultivar) 
  # data_c[[1]] <- with(temp, tibble(var="DM", year=year_pre, 
  #                                  doy=doy_pre, data=mass_pre/10) %>% drop_na())
  # data_c[[2]] <- with(temp, tibble(var="DM", year=year_post, 
  #                                  doy=doy_post, data=mass_post/10) %>% drop_na())
  
  # ryegrass tillers
  if (ablock==0){
    temp <- data_till_sum %>%
      select(block, seed_rate, cultivar, sampling, date_till, tillers) %>%
      filter(seed_rate %in% aseed_rate & cultivar %in% acultivar) %>% 
      mutate(
        date_till2=if_else(date_till %in% grazing_dates, date_till-days(1), date_till), # avoid grazing dates
        incalib=((date_till2>=calib_start)&(date_till2<=calib_end)),
        weight=as.numeric(incalib)
      ) %>% 
      group_by(sampling) %>% 
      summarise(
        date_range = max(date_till2, na.rm=TRUE)-min(date_till2, na.rm=TRUE),
        date_till2 = as_date(mean(date_till2, na.rm=TRUE)),
        tillers = mean(tillers, na.rm=TRUE),
        weight = mean(weight, na.rm=TRUE)
      )
  } else {
    temp <- data_till_sum %>%
      select(block, seed_rate, cultivar, sampling, date_till, tillers) %>%
      filter(block==ablock & seed_rate %in% aseed_rate & cultivar %in% acultivar) %>% 
      mutate(
        date_till2=if_else(date_till %in% grazing_dates, date_till-days(1), date_till), # avoid grazing dates
        incalib=((date_till2>=calib_start)&(date_till2<=calib_end)),
        weight=as.numeric(incalib)
      ) %>% 
      group_by(sampling) %>% 
      summarise(
        date_range = max(date_till2, na.rm=TRUE)-min(date_till2, na.rm=TRUE),
        date_till2 = as_date(mean(date_till2, na.rm=TRUE)),
        tillers = mean(tillers, na.rm=TRUE),
        weight = mean(weight, na.rm=TRUE)
      )
  }
  # stopifnot(all(temp$date_till2==temp$date_till))
  data_c[[3]] <- with(temp, tibble(var="TILTOT", year=year(date_till2), doy=yday(date_till2), 
                                   data=tillers, sd=err_c["TILTOT"], type="sd",
                                   weight=weight) %>% drop_na())
  
  # ryegrass mass (total or above cutting height? depending on definition of yield_bot)
  if (ablock==0){
    temp <- data_bm %>%
      select(block, seed_rate, cultivar, grazing, leaf, stem, dead, leaf_below, stem_below, dead_below, yield_above, yield_below, date_bot, date_grazed) %>%
      filter(seed_rate %in% aseed_rate & cultivar %in% acultivar) %>% 
      mutate(
        date_bot2=if_else(date_bot >= date_grazed, date_grazed-days(1), date_bot), # avoid grazing dates
        incalib=((date_bot2>=calib_start)&(date_bot2<=calib_end)),
        weight=as.numeric(incalib)
      ) %>% 
      group_by(grazing) %>% 
      summarise(
        date_range = max(date_bot2, na.rm=TRUE)-min(date_bot2, na.rm=TRUE),
        date_bot2 = mean(date_bot2, na.rm=TRUE),
        leaf = mean(leaf, na.rm=TRUE),
        leaf_below = mean(leaf_below, na.rm=TRUE),
        stem = mean(stem, na.rm=TRUE),
        stem_below = mean(stem_below, na.rm=TRUE),
        dead = mean(dead, na.rm=TRUE),
        dead_below = mean(dead_below, na.rm=TRUE),
        yield_above = mean(yield_above, na.rm=TRUE),
        yield_below = mean(yield_below, na.rm=TRUE),
        weight = mean(weight, na.rm=TRUE)
      )
  } else {
    temp <- data_bm %>%
      select(block, seed_rate, cultivar, grazing, leaf, stem, dead, leaf_below, stem_below, dead_below, yield_above, yield_below, date_bot, date_grazed) %>%
      filter(block==ablock & seed_rate %in% aseed_rate & cultivar %in% acultivar) %>% 
      mutate(
        date_bot2=if_else(date_bot >= date_grazed, date_grazed-days(1), date_bot), # avoid grazing dates
        incalib=((date_bot2>=calib_start)&(date_bot2<=calib_end)),
        weight=as.numeric(incalib)
      ) %>% 
      group_by(grazing) %>% 
      summarise(
        date_range = max(date_bot2, na.rm=TRUE)-min(date_bot2, na.rm=TRUE),
        date_bot2 = mean(date_bot2, na.rm=TRUE),
        leaf = mean(leaf, na.rm=TRUE),
        leaf_below = mean(leaf_below, na.rm=TRUE),
        stem = mean(stem, na.rm=TRUE),
        stem_below = mean(stem_below, na.rm=TRUE),
        dead = mean(dead, na.rm=TRUE),
        dead_below = mean(dead_below, na.rm=TRUE),
        yield_above = mean(yield_above, na.rm=TRUE),
        yield_below = mean(yield_below, na.rm=TRUE),
        weight = mean(weight, na.rm=TRUE)
      )
  }
  # stopifnot(all(temp$date_bot2==temp$date_bot))
  data_c[[4]] <- with(temp, tibble(var="CLV", year=year(date_bot2), doy=yday(date_bot2), 
                                   data=leaf/100*yield_above/10*0.45+leaf_below/100*yield_below/10*0.45, 
                                   sd=err_c["CLV"], type="sd",
                                   weight=weight) %>% drop_na())
  data_c[[5]] <- with(temp, tibble(var="CST", year=year(date_bot2), doy=yday(date_bot2), 
                                   data=stem/100*yield_above/10*0.45+stem_below/100*yield_below/10*0.45, 
                                   sd=err_c["CST"], type="sd",
                                   weight=weight) %>% drop_na())
  data_c[[6]] <- with(temp, tibble(var="CLVD", year=year(date_bot2), doy=yday(date_bot2), 
                                   data=dead/100*yield_above/10*0.45+dead_below/100*yield_below/10*0.45, 
                                   sd=err_c["CLVD"], type="sd",
                                   weight=weight) %>% drop_na())
  data_c[[7]] <- with(temp, tibble(var="CRT", year=year(date_bot2), doy=yday(date_bot2),
                                   data=leaf/100*yield_above/10*0.45+leaf_below/100*yield_below/10*0.45,
                                   sd=err_c["CRT"], type="sd",
                                   weight=weight) %>% drop_na())
  
  # light interception (but this includes all species!)
  # temp <- data_li %>%
  #   select(block, seed_rate, cultivar, li, date) %>%
  #   filter(block==ablock & seed_rate %in% aseed_rate & cultivar %in% acultivar) 
  # data_c[[5]] <- with(temp, tibble(var="LINT", year=year(date), 
  #                                  doy=yday(date), data=li) %>% drop_na())
  
  # soil moisture
  if (ablock==0){
    temp <- data_sm %>%
      select(block, seed_rate, cultivar, mean_sm, date_sm) %>%
      filter(seed_rate %in% c(aseed_rate, NA_character_) & cultivar %in% c(acultivar, NA_character_)) %>% 
      mutate(
        # soil moisture is not affected by grazing so no need to adjust date
        incalib=((date_sm>=calib_start)&(date_sm<=calib_end)),
        weight=as.numeric(incalib)
      ) 
  } else {
    temp <- data_sm %>%
      select(block, seed_rate, cultivar, mean_sm, date_sm) %>%
      filter(block==ablock & seed_rate %in% c(aseed_rate, NA_character_) & cultivar %in% c(acultivar, NA_character_)) %>% 
      mutate(
        # soil moisture is not affected by grazing so no need to adjust date
        incalib=((date_sm>=calib_start)&(date_sm<=calib_end)),
        weight=as.numeric(incalib)
      ) 
  }
  data_c[[8]] <- with(temp, tibble(var="WCLM", year=year(date_sm), doy=yday(date_sm), 
                                   data=mean_sm, sd=err_c["WCLM"], type="sd",
                                   weight=weight) %>% drop_na())
  
  # basal area
  # temp <- data_till_sum %>%
  #   select(block, seed_rate, cultivar, p2000, date_till) %>%
  #   filter(block==ablock & seed_rate %in% aseed_rate & cultivar %in% acultivar)
  # data_c[[7]] <- with(temp, tibble(var="BASAL", year=year(date_till),
  #                                  doy=yday(date_till), data=p2000*100, sd=err_c["BASAL"], type="sd",
  #                                  weight=ifelse(((date_till>=calib_start)&(date_till<=calib_end)), 1, 0)) %>% drop_na())

  # bind list and write file
  data_calib <- bind_rows(data_c)
  data_calib <- arrange(data_calib, var, year, doy)
  write.table(data_calib, file=data_file_name, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

} # end data set loop

# save data
rdata_file_name <- paste("raw_data/", site, ".RData", sep="")
save.image(file=rdata_file_name)

# produce report if desired
print("Preparing report... (please ensure output file is not open in Word)")
time_stamp <- today()
rmarkdown::render(input = "raw_data/report_any.rmd",
                  # output_format = "word_document",
                  output_file = paste("report_", str_to_lower(site), ".docx", sep=""),
                  output_dir = "raw_data")

} # next site
