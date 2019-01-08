# read in data set

# load libraries
library(tidyverse)
library(readxl)
library(stringr)
library(lubridate)

# load useful functions
source("raw_data/utility_functions.r")

#### set options ####

# choose site
sites <- c("Northland", "Scott", "Lincoln")
for (site in sites){
print(site)
  
# choose output data
acultivar <- "Alto" # only Alto and Halo have light interception data
aseed_rate <- "18kg"
calib_start <- ymd("20120401") # period for data weight = 1
calib_end <- ymd("20171231") # period for data weight = 1

# seed rates in order
seed_rate_levels <- c("6kg", "12kg", "18kg", "24kg", "30kg")

#### read and normalise data ####
if (site=="Scott"){
  
  raw_file_name <- "raw_data/FD1004 Data For Modelling.xlsx"
  data_rpm <- read_xlsx(raw_file_name, sheet="Waikato RPM Height data") %>% autosnake() 
  data_cut1 <- read_xlsx(raw_file_name, sheet="Cut Yield Data Year1 2011to12") %>% autosnake() %>% 
    mutate(mean_dm=NA_real_)
  data_cut2 <- read_xlsx(raw_file_name, sheet="Cut Yield Data Year2 onwards") %>% autosnake() %>% 
    mutate(yield_kg_dm_ha=if_else(yield_kg_dm_ha>0, yield_kg_dm_ha, NA_real_))
  data_till <- read_xlsx(raw_file_name, sheet="Tiller density data Waikato") %>% autosnake() 
  data_bot <- read_xlsx(raw_file_name, sheet="Botanical Composition data ", skip=1) %>% autosnake() 
  data_li <- read_xlsx(raw_file_name, sheet="Waikato LightInterception") %>% autosnake() 
  data_sm <- read_xlsx(raw_file_name, sheet="Soil moisture data") %>% autosnake() 

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
      seed_rate = factor(aseed_rate, levels=seed_rate_levels),
      cultivar = acultivar,
      dom = day(date_measured_d),
      soil_moisture = sm10cm_av * 100 
    ) %>% 
    filter(dom==1) %>%  # only take one point per month
    mutate(
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
      dom = day(date_measured_d),
      seed_rate = factor(aseed_rate, levels=seed_rate_levels),
      cultivar = acultivar
    ) %>% 
    group_by(cultivar, seed_rate, date_measured_d, dom) %>%
    summarise(soil_moisture = mean(soil_moisture)) %>% 
    ungroup() %>% 
    filter(dom==1) %>%  # %>%  only take one point per month
    mutate(
      block = rep_len(c(2,3,4), n()) # split soil moisture data across blocks 2,3,4
    )  
  
}

#### Rising Plate Meter Data ####
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
  group_by(block, cultivar, seed_rate, grazing) %>% 
  mutate( 
    samples = n(),
    date_pre = mean(date_pre),
    mass_pre = mean(mass_pre),
    date_post = mean(date_post),
    mass_post = mean(mass_post),
    date_grazed = mean(date_grazed)
    ) %>% 
  ungroup()

data_rpm <- data_rpm %>%
  group_by(block, cultivar, seed_rate) %>%
  arrange(grazing) %>%
  mutate(
    lag_grazing = grazing-lag(grazing,1),
    date_post_last = lag(date_post ,1), # FIXME could fail if missing grazings
    mass_post_last = lag(mass_post ,1), # FIXME could fail if missing grazings 
    growth_days_pre = as.numeric(difftime(date_pre, date_post_last), units="days"),
    growth_pre = mass_pre - mass_post_last,
    growth_rate_pre = growth_pre / as.double(growth_days_pre),
    growth_rate_post = lead(growth_rate_pre, 1)
    ) %>% 
  ungroup()

stopifnot(all(data_rpm$lag_grazing %in% c(NA,1)))

data_rpm <- data_rpm %>%
  group_by(block, cultivar, seed_rate, grazing) %>%
  mutate(
    delay_pre = as.numeric(difftime(date_grazed, date_pre), units="days"),
    delay_post = as.numeric(difftime(date_post, date_grazed), units="days"),
    mass_grazed = mass_pre, # + growth_rate_pre * delay_pre * 0,
    mass_resid = mass_post, # - growth_rate_post * delay_post * 0, 
    harv = (1 - mass_resid / mass_grazed) * 100 # proportion harvested
    ) %>% 
  ungroup()

## Pasture Cut Data ####
data_cut1 <- data_cut1 %>%
  mutate(seed_rate = factor(seed_rate, levels=seed_rate_levels)) %>%
  rename(
    date_cut = date_d,
    yield = yield_kg_dm_ha,
    dm_pc = mean_dm,
    grazing = graze_no
  ) %>%
  select(date_cut, block, seed_rate, cultivar, yield, dm_pc, grazing) %>% 
  group_by(block, seed_rate, cultivar, grazing) %>%
  summarise(
    samples = n(),
    date_cut = mean(date_cut),
    yield = mean(yield),
    dm_pc = mean(dm_pc)
  ) %>% 
  ungroup()

data_cut2 <- data_cut2 %>%
  mutate(seed_rate = factor(seed_rate, levels=seed_rate_levels)) %>%
  rename(
    date_cut = date_d,
    yield = yield_kg_dm_ha,
    dm_pc = mean_dm,
    grazing = graze_no
  ) %>%
  select(date_cut, block, seed_rate, cultivar, yield, dm_pc, grazing) %>% 
  group_by(block, seed_rate, cultivar, grazing) %>% 
  summarise(
    samples = n(),
    date_cut = mean(date_cut),
    yield = mean(yield),
    dm_pc = mean(dm_pc)
  ) %>% 
  ungroup()

# combine tables by row keeping common variables
data_cut <- suppressMessages(full_join(data_cut1, data_cut2))

#### Estimate Mass Below Cut ####
# join cut and rpm data by grazing no. and estimate mass on cutting date
# mass below cutting height is this minus cut yield
data_bc <- data_rpm %>%
  filter(!is.na(date_pre)) %>% 
  left_join(data_cut, by=c("block", "seed_rate", "cultivar", "grazing")) %>%
  select(date_pre, date_cut, block, seed_rate, cultivar, mass_pre, yield, grazing) %>%
  mutate(delay = as.integer(as.Date(date_cut) - as.Date(date_pre))) %>% 
  # filter(abs(delay)<4) %>%
  # drop_na() %>%
  mutate(
    mass_cut = mass_pre, # + growth_rate_pre * delay * 0,
    below = mass_cut - yield
    )
unique(data_bc$delay) # check delay between cut and rpm

# copy below back into data_cut
temp <- data_cut %>%
  left_join(data_bc, by=c("block", "seed_rate", "cultivar", "grazing")) %>%
  select(date_cut.x, block, seed_rate, cultivar, yield.x, dm_pc, grazing,
         date_pre, mass_pre, delay, below) %>%
  rename(date_cut = date_cut.x,
         yield = yield.x)
data_cut <- temp

## Tiller Density Data ####
data_till <- data_till %>%
  mutate(seed_rate = factor(seed_rate, levels=seed_rate_levels)) %>%
  rename(
    date_till = date_d,
    tillers = ryegrass_tiller_density_tillers_m2
  ) %>% 
  mutate(
    sampling_month = floor_date(date_till, "months")
    )

# create separate summary df
data_till_sum <- data_till %>%
  group_by(block, cultivar, seed_rate, sampling_month) %>%
  summarise(
    samples = n(),
    mean_tillers = mean(tillers),
    sd_tillers = sd(tillers),
    date_till = mean(date_till),
    tillers = mean_tillers
    ) %>%
  ungroup()

#### Botanical Composition Data ####
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
    total = leaf + stem + ann + wc + og + weed + dead, # should always be 100
    rgfrac = 100 * (leaf + stem) / (leaf + stem + ann + wc + og + weed), # of green
    month = month(dmy(paste("01", month, "2011"))),
    # estimate botanicals composition below cutting height
    # based on Tozer data from Ruakura only
    # seems dangerous to apply these assumptioins to all sites!
    # be careful to ensure it gives reasonable estimates
    # weighted average seems safe?
    leaf_below = 0*leaf / (0.79*leaf + 0.56*stem + 0.79*ann + 0.50*wc + 0.79*og + weed + 3.30*dead) * 100,
    stem_below = 0*stem / (0.79*leaf + 0.56*stem + 0.79*ann + 0.50*wc + 0.79*og + weed + 3.30*dead) * 100,
    dead_below = 0*dead / (0.79*leaf + 0.56*stem + 0.79*ann + 0.50*wc + 0.79*og + weed + 3.30*dead) * 100,
    total_below = leaf_below + stem_below + dead_below
  ) %>%
  select(leaf, stem, ann, wc, og, weed, dead, date_bot, total, block, seed_rate, cultivar, rgfrac, leaf_below, stem_below, dead_below, total_below) 
  
# stop()

# check samples
temp <- data_bot %>% 
  group_by(block, cultivar, seed_rate, date_bot) %>% 
  summarise(
    samples = n()
  ) %>% 
  ungroup()

#### add grazing number to data_bot ####
temp_cut <- data_cut %>% 
  # distinct(block, cultivar, seed_rate, date_cut) %>% 
  group_by(block, cultivar, seed_rate, grazing) %>% 
  summarise(
    number = n(),
    date_cut = mean(date_cut)
  ) %>% 
  ungroup()

closest_cut <- function(block_bot, cultivar_bot, seed_rate_bot, date_bot){
  temp <- temp_cut %>% 
    filter(block==block_bot, cultivar==cultivar_bot, seed_rate==seed_rate_bot)
  ans <- as.POSIXct(rep(NA, length(date_bot)))
  if (nrow(temp)>0){ # there are some cut dates for this block-cultivar-seed_rate
    for (i in 1:length(date_bot)){
      delay <- as.double(difftime(date_bot[i], temp$date_cut, units="days")) # vector of delays
      ans[i] <- temp$date_cut[which(abs(delay)==min(abs(delay)))[[1]]]
    }
  }
  return(ans)
}

temp_bot <- data_bot %>% 
  group_by(block, cultivar, seed_rate, date_bot) %>%
  mutate(date_cut = closest_cut(block, cultivar, seed_rate, date_bot)) %>% 
  ungroup %>% 
  left_join(select(data_cut, block, cultivar, seed_rate, date_cut, grazing),
            by = c("block", "cultivar", "seed_rate", "date_cut")) %>% 
  mutate(delay = as.integer(as.Date(date_bot) - as.Date(date_cut)),
         date_cut = if_else(abs(delay) < 7, date_cut, as.POSIXct(NA)),
         grazing = if_else(abs(delay) < 7, grazing, NA_real_))
  
data_bot <- temp_bot
  
#### gather data_bot into data_bot2 ####
data_bot2 <- data_bot %>%
  select(leaf, stem, ann, wc, og, weed, dead, block, cultivar, seed_rate, date_bot) %>% 
  gather(leaf, stem, ann, wc, og, weed, dead, key="species", value="fraction") %>%
  mutate(species = factor(species, 
                          levels=c("leaf", "stem", "ann", "og", "wc", "weed", "dead")))

#### calculate botanical mass ####
# join cut and botanical data by grazing
# data_cut$date_cut <- data_cut$date
data_bm <- data_bot %>%
  # rename(date_bot = date) %>%
  left_join(data_cut, by=c("grazing", "date_cut", "block", "seed_rate", "cultivar")) %>%
  select(date_bot, block, seed_rate, cultivar, 
         leaf, stem, ann, wc, og, weed, dead, leaf_below, stem_below, dead_below,
         date_cut, grazing, yield, below) %>%
  rename(yield_cut = yield) %>%
  mutate(yield_bot = yield_cut #+ below * 0 + growth_rate_pre * delay * 0 # assumed on date_bot, not including below mass
    )

#### gather data_bm into data_bm2 ####
data_bm2 <- data_bm %>%
  gather(leaf, stem, ann, wc, og, weed, dead, key="species", value="fraction") %>%
  mutate(
    species = factor(species, levels=c("leaf", "stem", "ann", "og", "wc", "weed", "dead")),
    species_mass = fraction / 100 * yield_bot # include mass below cutting?
    ) 

#### Soil Moisture Data ####
data_sm <- data_sm %>%
  mutate(seed_rate = factor(seed_rate, levels=seed_rate_levels)) %>%
  rename(
    date_sm = date_measured_d,
    sm = soil_moisture
  ) %>%
  group_by(block, cultivar, seed_rate, date_sm) %>%
  summarise(mean_sm = mean(sm)) %>% 
  ungroup()

#### make model input files ####

# move post harvest measurements if on grazing day since model does harvest at start of day
data_rpm$date_post2 <- as.Date(with(data_rpm,
                            ifelse(date_post > date_grazed, 
                                   as.Date(date_post), 
                                   as.Date(date_post) + 1)
                            ), origin="1970-01-01")

# some calculations
data_rpm <- data_rpm %>%
  mutate(
    year_pre = year(date_pre),
    doy_pre = yday(date_pre),
    year_post = year(date_post2),
    doy_post = yday(date_post2),
    year_grazed = year(date_grazed),
    doy_grazed = yday(date_grazed), 
    days_pre = difftime(date_grazed, date_pre),
    days_post = difftime(date_post2, date_grazed)
    )

ablock <- 3 # for testing
for (ablock in c(1,2,3,4,5)) { # loop through blocks

  # write harvest dates and harvest % for selected series
  harv_file_name <- paste("raw_data/harvest_", site, "_", ablock, ".txt", sep="")
  data_h <- data_rpm %>%
    select(block, seed_rate, cultivar, year_grazed, doy_grazed, harv) %>%
    filter(block == ablock & seed_rate==aseed_rate & cultivar==acultivar) %>%
    drop_na()
  days_harvest <- matrix(as.integer(-1), nrow=100, ncol=3) # up to 100 harvests
  days_harvest[1:nrow(data_h),] <- c(data_h$year_grazed, data_h$doy_grazed, data_h$harv)
  write.table(days_harvest, file=harv_file_name, row.names=FALSE, col.names=FALSE, sep="\t")
  
  # write calibration data
  data_file_name <- paste("raw_data/data_calibration_", site, "_", ablock, ".txt", sep="")
  
  # collect the data in this list
  data_c <- vector("list", 10) 
  err_c <- c(TILTOT=1000, CLV=25, CST=2.5, CLVD=40, WCLM=5, BASAL=10)
  
  # pre and post mass (but this includes other species!)
  # temp <- data_rpm %>%
  #   select(block, seed_rate, cultivar, mass_pre, year_pre, doy_pre, 
  #          mass_post, year_post, doy_post) %>%
  #   filter(block==ablock & seed_rate==aseed_rate & cultivar==acultivar) 
  # data_c[[1]] <- with(temp, tibble(var="DM", year=year_pre, 
  #                                  doy=doy_pre, data=mass_pre/10) %>% drop_na())
  # data_c[[2]] <- with(temp, tibble(var="DM", year=year_post, 
  #                                  doy=doy_post, data=mass_post/10) %>% drop_na())
  
  # ryegrass tillers
  temp <- data_till_sum %>%
    select(block, seed_rate, cultivar, tillers, date_till) %>%
    filter(block==ablock & seed_rate==aseed_rate & cultivar==acultivar) 
  data_c[[3]] <- with(temp, tibble(var="TILTOT", year=year(date_till), doy=yday(date_till), 
                                   data=tillers, sd=err_c["TILTOT"], type="sd",
                                   weight=ifelse(((date_till>=calib_start)&
                                                    (date_till<=calib_end)), 1, 0))
                      %>% drop_na())
  
  # ryegrass mass (total or above cutting height? depending on definition of yield_bot)
  temp <- data_bm %>%
    # rename(date = date_cut) %>%
    select(block, seed_rate, cultivar, leaf, stem, dead, leaf_below, stem_below, dead_below, yield_bot, below, date_bot) %>%
    filter(block==ablock & seed_rate==aseed_rate & cultivar==acultivar) 
  data_c[[4]] <- with(temp, tibble(var="CLV", year=year(date_bot), 
                                   doy=yday(date_bot), data=leaf/100*yield_bot/10*0.45+leaf_below/100*below/10*0.45, 
                                   sd=err_c["CLV"], type="sd",
                                   weight=ifelse(((date_bot>=calib_start)&(date_bot<=calib_end)), 1, 0)) %>% drop_na())
  data_c[[5]] <- with(temp, tibble(var="CST", year=year(date_bot), 
                                   doy=yday(date_bot), data=stem/100*yield_bot/10*0.45+stem_below/100*below/10*0.45, 
                                   sd=err_c["CST"], type="sd",
                                   weight=ifelse(((date_bot>=calib_start)&(date_bot<=calib_end)), 1, 0)) %>% drop_na())
  # data_c[[6]] <- with(temp, tibble(var="CLVD", year=year(date_bot), 
  #                                  doy=yday(date_bot), data=dead/100*yield_bot/10*0.45+dead_below/100*below/10*0.45, 
  #                                  sd=err_c["CLVD"], type="sd",
  #                                  weight=ifelse(((date_bot>=calib_start)&(date_bot<=calib_end)), 1, 0)) %>% drop_na())
  
  # light interception (but this includes all species!)
  # temp <- data_li %>%
  #   select(block, seed_rate, cultivar, li, date) %>%
  #   filter(block==ablock & seed_rate==aseed_rate & cultivar==acultivar) 
  # data_c[[5]] <- with(temp, tibble(var="LINT", year=year(date), 
  #                                  doy=yday(date), data=li) %>% drop_na())
  
  # soil moisture
  temp <- data_sm %>%
    select(block, seed_rate, cultivar, mean_sm, date_sm) %>%
    filter(block==ablock & seed_rate==aseed_rate & cultivar==acultivar) 
  data_c[[8]] <- with(temp, tibble(var="WCLM", year=year(date_sm), 
                                   doy=yday(date_sm), data=mean_sm, sd=err_c["WCLM"], type="sd",
                                   weight=ifelse(((date_sm>=calib_start)&(date_sm<=calib_end)), 1, 0)) %>% drop_na())
  
  # basal area
  # temp <- data_till_sum %>%
  #   select(block, seed_rate, cultivar, p2000, date_till) %>%
  #   filter(block==ablock & seed_rate==aseed_rate & cultivar==acultivar)
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

} # next site
