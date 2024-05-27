#### Script Info ####

# Authors: Dario San Segundo Molina, Ignacio Morales Castilla
# 
# GitHub repo: github.com/dario-ssm/PhenoBrassicaPests

# Aim: apply phyisiological models from phenology studies of Brassica pests 
#      to spatio-temporal contexts. Specifically to obtain long-term temporal trends
#      of phenological cues of a Brassica pest species ("Pieris rapae") in Spain with different methodologies

# Description: we use Spain02 climatic database and thermal traits (i.e. DDs, LDTs) of insect pest to predict
#              dates of emergence and voltinism across Spain and years. We analyse temporal trends for yearly variation
#              of thermal traits obtained by linear degree-days modeling.
#
# Species: Pieris rapae
# Climate dataset: Spain02 v5  (Herrera et al. 2016 and Kotlarsky et al. 2017; see http://www.meteo.unican.es/datasets/spain02 )
## Aknowledgements: The authors thank AEMET and UC by the data provided for this work (Spain02v5 gridded temperature data set). 

# 0. Load ----
library(raster)
library(tidyverse)
library(lubridate)
library(here)
library(nlme)
library(ggthemes)
library(sf)
library(viridis)
library(mapSpain)
source(here("Scripts/1-functions_phenobraspests.R"))
# source(here("Scripts/3-box_pieris_rapae.R"))

load(here("Data/daily_tmax_df.RData"))
# print(daily_tmax_df)
load(here("Data/daily_tmin_df.RData"))
# print(daily_tmin_df)
pieris_data <- read_delim(here("Data/pieris_devdata/gilbert_pupa.csv"),
                          delim = ";") %>% 
  rename(stage = ...3,
         parasite = ...4, 
         reference = interact)


# Spatio-Temporal Forecasts  ---------------------------------------
## after model selection, best candidates were: 
#   gaussian, oneill, briere1, weibull and ssi_low

## after comparing with observational data, 
## we use gaussian oneill with daily data.

##### a) daily temp at broccoli regions -------------------------------------------------------
## downloaded data from http://www.meteo.unican.es/datasets/spain02 on Jan-24th-2022
#  (now in 2023 it is not available anymore)

tmax_spain <- brick("~/Spain02_v5.0_010reg_aa3d/Spain02_v5.0_DD_010reg_aa3d_tasmax.nc")
tmin_spain <- brick("~/Spain02_v5.0_010reg_aa3d/Spain02_v5.0_DD_010reg_aa3d_tasmin.nc")

broccoli_spain_regions <- mapSpain::esp_get_prov(prov = c("Murcia", "Alicante",
                                                          "Almería", "Granada",
                                                          "Cádiz", "Badajoz", 
                                                          "Albacete", "Zaragoza",
                                                          "Navarra")) #regions with > 10,000 tons of production each year
broccoli_spain <- st_as_sf(broccoli_spain_regions$geometry) %>%  # transform to sf object to operate
  st_transform(crs = st_crs(tmax_spain)) #transform to cooordinate system of the raster

## extract climate into broccoli regions

#tmin
spain_broccoli_tmin_mask <- mask(tmin_spain, broccoli_spain)
spain_broccoli_tbl_tmin <- as(spain_broccoli_tmin_mask, "SpatialPixelsDataFrame") %>% 
  as_tibble() %>% 
  drop_na() 
spain_broccoli_tbl_tmin_long <- spain_broccoli_tbl_tmin %>% 
  pivot_longer(cols = -c(x,y),
               names_to = "date",
               values_to = "daily_tmin") %>% 
  mutate(date = (str_sub(date, 2, -1L)),
         date = ymd(date))
## tmax
spain_broccoli_tmax_mask <- mask(tmax_spain, broccoli_spain)
spain_broccoli_tbl_tmax <- as(spain_broccoli_tmax_mask, "SpatialPixelsDataFrame") %>% 
  as_tibble() %>% 
  drop_na() 
spain_broccoli_tbl_tmax_long <- spain_broccoli_tbl_tmax %>% 
  pivot_longer(cols = -c(x,y),
               names_to = "date",
               values_to = "daily_tmax") %>% 
  mutate(date = (str_sub(date, 2, -1L)),
         date = ymd(date))

## join and average
spain_broccoli_temps <- inner_join(spain_broccoli_tbl_tmin_long,
                                   spain_broccoli_tbl_tmax_long) %>% 
  mutate(temperature = map2_dbl(.x = daily_tmin,
                                .y = daily_tmax,
                                .f = ~mean(c(.x, .y)))
  )
head(spain_broccoli_temps) # <- check it out
save(spain_broccoli_temps, file = here("Data/spain_broccoli_temps.RData"))

# load(file = here("Data/spain_broccoli_temps.RData"))
head(spain_broccoli_temps)

##### b) Broccoli model -------------------------------------------------------
## model from Tan et al. (2000) Scientia Horticulturae for Marathon broccoli
values_Tan2000 <- function(){
  tbase= c(0, 0)
  tupper = c(25, 25)
  gdds = c(627, 678)
  crop_stage = c("head initiation","harvest")
  tan2000 <- tibble(crop_stage, tbase, tupper, gdds)
  return(print(tan2000))
}
traits_tan2000 <- values_Tan2000()
## Note: according to some papers, Marathon does not support high temperatures
## (~25 for t_avg)
## first set planting date on 15th August 

broccoli_planting <- function(curr_year){
  planting_date = ifelse(leap_year(curr_year),
                         yday("2023-08-16"),
                         yday("2023-08-15")
  )
  return(planting_date)
}

crop_calendar <- function(doy, planting_date){
  day_after_planting = case_when(doy >= planting_date ~ doy - planting_date + 1,
                                 doy < planting_date ~ max(doy) - planting_date + doy + 1)
  return(day_after_planting) # <- dap
}


broccoli_temps_cropcal <- spain_broccoli_temps %>% 
  group_by(x, y) %>% 
  mutate(id_cell = cur_group_id()) %>%  #add a cell identifier for (x, y) coordinates.
  mutate(doy = yday(date),
         year = year(date),
         dop = broccoli_planting(year)
  ) %>% # <- day of planting (dop)
  group_by(year) %>% 
  mutate(crop_day = case_when(leap_year(year) == TRUE & doy < dop ~ crop_calendar(doy, dop),
                              leap_year(year) == TRUE & doy >= dop ~ crop_calendar(doy, dop),
                              leap_year(year) == FALSE ~ crop_calendar(doy, dop)
  )
  ) %>% 
  ungroup() %>% 
  mutate(day_after_planting = ifelse(year == 2015 & doy >= dop,
                                     NA,
                                     ifelse(year == 1950 & doy < dop,
                                            NA,
                                            crop_day))
  ) %>%
  mutate(season = if_else(doy < dop,
                          year - 1,
                          year)) %>% 
  drop_na()

# view(broccoli_temps_cropcal %>% filter(id_cell == 1) %>% filter(year %in% c(1950:1955)))

## now we compute broccoli head initiation and harvest day
daily_avg_dds_broccoli <- function(tavg_day, t_base, t_high){
  daily_dds = ifelse(tavg_day >= t_base & tavg_day <= t_high,
                     tavg_day - t_base,
                     0)
  
}

broccoli_gdds <- broccoli_temps_cropcal %>% 
  mutate(daily_avg_dds = map_dbl(.x = temperature,
                                 .f = ~daily_avg_dds_broccoli(tavg_day = .x, 
                                                              t_base = mean(traits_tan2000$tbase),
                                                              t_high = mean(traits_tan2000$tupper)
                                 )
  )
  ) %>% 
  group_by(x, y, id_cell, season) %>% 
  summarise(gdds_season = cumsum(daily_avg_dds),
            date = date) %>% 
  ungroup() %>% 
  mutate(crop_remains = map_dbl(.x = gdds_season,
                                .f = ~logic_dd(heat_units = sum(traits_tan2000$gdds),
                                               rate_cumsum = .x)
  )
  )
broccoli_gdds_harvest <- broccoli_gdds %>% 
  group_by(x, y, season) %>% 
  summarise(harvest_day = sum(crop_remains))
save(broccoli_gdds, file = "~/broccoli_gdds.RData")
save(broccoli_gdds_harvest, file = here("Data/broccoli_gdds_harvest.RData")) 


### LOAD IF RE-RUNNING THIS SCRIPT ---
# load(file = here("Data/broccoli_gdds.RData"))
# load(file = "~/broccoli_gdds.RData")
# load(file = here("Data/broccoli_gdds_harvest.RData"))


# and plot it
broccoli_gdds_average <- broccoli_gdds_harvest %>% 
  group_by(x, y) %>% 
  summarise(harvest_day = round(mean(harvest_day))) %>% 
  filter(harvest_day < 150) # <- avoid high altitudes
spain_map_sf <- st_read("~/Dario Investigacion/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp")
xrange <- c(-10, 3)
yrange <- c(35, 44)
province_borders <- broccoli_spain_regions$geometry
ocean <- st_read("/Users/dario-ssm/Downloads/ne_10m_ocean/ne_10m_ocean.shp")


spain_broccoli_phenomap <- ggplot()+  
  geom_sf(data = ocean, 
          fill = "#98c1d9")+
  geom_sf(data = spain_map_sf)+
  geom_tile(data = broccoli_gdds_average,
            aes(x, y, fill = harvest_day))+
  scale_fill_viridis(name = "Harvest day (Day after common transplant on Aug 15th",
                     direction = -1)+
  geom_sf(data = province_borders, 
          fill = NA,
          color = "gray23")+
  coord_sf(xlim = xrange, ylim = yrange)+
  theme_map()+
  theme(legend.position = "bottom")
spain_broccoli_phenomap
ggsave(here("Data/spain_broccoli_phenomap.png"),
       width = 16, height = 16, units = "cm")

# broccoli phenology shift
#first pooling
broccoli_pheno_trends_pooled <- broccoli_gdds %>%
  group_by(x, y, season) %>% 
  mutate(cell = cur_group_id()) %>% 
  summarise(harvest_day = sum(crop_remains)) %>% 
  ungroup() %>% 
  group_by(x, y) %>% 
  mutate(id_cell = cur_group_id())
broccoli_pheno_trends_pooled
phenoshift_broccoli_glmer <- lme4::glmer(harvest_day ~ season + (1|id_cell),
                                         data = broccoli_pheno_trends_pooled,
                                         na.action = na.exclude,
                                         family = Gamma(link = "log")
)
phenoshift_broccoli_glmer
visreg::visreg(phenoshift_broccoli_glmer)
sum_phenoshift_broccoli <- summary(phenoshift_broccoli_glmer)
shift_phenoshift_broccoli <- tibble(intercept = sum_phenoshift_broccoli$coefficients[1,1],
                                    slope_season = sum_phenoshift_broccoli$coefficients[2,1]
) %>% 
  mutate(across(c(1:2), 
                ~exp(.x))
  ) %>% 
  rowwise() %>% 
  summarise(shift_decade = (intercept*slope_season -intercept)*10)
print(shift_phenoshift_broccoli) # <- -0.471 days per decade

phenoint_broccoli_glmer <- lme4::glmer(harvest_day ~ 1 + (1|id_cell),
                                       data = broccoli_pheno_trends_pooled,
                                       na.action = na.exclude,
                                       family = Gamma(link = "log")
)
sum_phenoint_broccoli <- summary(phenoint_broccoli_glmer)
intercept_broccoli <- exp(sum_phenoint_broccoli$coefficients[1,1])
se_broccoli <- exp(sum_phenoint_broccoli$coefficients[1,2])
intercept_broccoli # <- mean day of harvest = day 85, i.e. 8th-Nov
se_broccoli # <- +/- 1.

##### c) Pieris rapae full development all year -------------------------------------------------------
## Note that this is a REALLY SLOW script. Life-stages' emergences and durations are computed
## sequentially. This script should better be functionalized for automation, but it is 
## currently working in a reproducible way.


###### i. Linear models ----
# from Jones & Ives (1987)
values_pieris_all<- function(){
  ldt= c(9.8, 9.8, 6.7, 10)
  heat_units = c(51, 171, 140, 60)
  life_stage = c("eggs","larvae","pupae", "adults")
  reference = c(rep("jones1987", 3), "gossard1977")
  pieris_all_stages <- tibble(life_stage,ldt,heat_units, reference)
  return(print(pieris_all_stages))
}

traits_pieris_allstages <- values_pieris_all() 

## first we compute voltinism shifts across the whole year
load(file = here("~/spain_broccoli_temps.RData")) # <- reload all data temps
head(spain_broccoli_temps)
test_pieris_pheno_spatial <- spain_broccoli_temps %>% 
  group_by(x, y) %>% 
  mutate(id_cell = cur_group_id()) %>% 
  filter(year(date) %in% c(1950:1959)) 


###### ii.  50's ----
pieris_pheno_spatial  <- test_pieris_pheno_spatial %>% 
  ## pupa (1st generation)
  mutate(year = year(date),
         daily_rate  = map2_dbl(.x = daily_tmax,  
                                .y = daily_tmin,
                                ~daily_avg(tmax = .x,
                                           tmin = .y, 
                                           LDT = as_vector(traits_pieris_allstages[3,2])
                                )
         )
  ) %>% 
  group_by(id_cell, year) %>% 
  summarise(phys_time = cumsum(daily_rate),
            date,
            daily_tmin,
            daily_tmax) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]),
                                   1,
                                   0),
         life_stage = if_else(phys_time_logic == 1,
                              "pupa",
                              "unknown"),
         complete = if_else(life_stage == "pupa" & phys_time_logic == 1,
                            "yes",
                            "no")
  ) %>% 
  ungroup() %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  mutate(life_stage = if_else(phys_time_logic == 0 & complete  == "no",
                              "pupa",
                              life_stage))


save(pieris_pheno_spatial, file = "~/50s_pieris_pheno_spatial.RData") 

###### iii. 60's ----
pieris_60s <-  spain_broccoli_temps %>% 
  group_by(x, y) %>% 
  mutate(id_cell = cur_group_id()) %>% 
  filter(year(date) %in% c(1960:1969)) 

pieris_pheno_spatial_60s  <- pieris_60s %>% 
  ## pupa (1st generation)
  mutate(year = year(date),
         daily_rate  = map2_dbl(.x = daily_tmax,  
                                .y = daily_tmin,
                                ~daily_avg(tmax = .x,
                                           tmin = .y, 
                                           LDT = as_vector(traits_pieris_allstages[3,2])
                                )
         )
  ) %>% 
  group_by(id_cell, year) %>% 
  summarise(phys_time = cumsum(daily_rate),
            date,
            daily_tmin,
            daily_tmax) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]),
                                   1,
                                   0),
         life_stage = if_else(phys_time_logic == 1,
                              "pupa",
                              "unknown"),
         complete = if_else(life_stage == "pupa" & phys_time_logic == 1,
                            "yes",
                            "no")
  ) %>% 
  ungroup() %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  
  mutate(life_stage = if_else(phys_time_logic == 0 & complete  == "no",
                              "pupa",
                              life_stage))


save(pieris_pheno_spatial_60s, file = "~/60s_pieris_pheno_spatial.RData") 

###### iv. 70's ----
pieris_70s <-  spain_broccoli_temps %>% 
  group_by(x, y) %>% 
  mutate(id_cell = cur_group_id()) %>% 
  filter(year(date) %in% c(1970:1979)) 

pieris_pheno_spatial_70s  <- pieris_70s %>% 
  ## pupa (1st generation)
  mutate(year = year(date),
         daily_rate  = map2_dbl(.x = daily_tmax,  
                                .y = daily_tmin,
                                ~daily_avg(tmax = .x,
                                           tmin = .y, 
                                           LDT = as_vector(traits_pieris_allstages[3,2])
                                )
         )
  ) %>% 
  group_by(id_cell, year) %>% 
  summarise(phys_time = cumsum(daily_rate),
            date,
            daily_tmin,
            daily_tmax) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]),
                                   1,
                                   0),
         life_stage = if_else(phys_time_logic == 1,
                              "pupa",
                              "unknown"),
         complete = if_else(life_stage == "pupa" & phys_time_logic == 1,
                            "yes",
                            "no")
  ) %>% 
  ungroup() %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  )

mutate(life_stage = if_else(phys_time_logic == 0 & complete  == "no",
                            "pupa",
                            life_stage))


save(pieris_pheno_spatial_70s, file = "~/70s_pieris_pheno_spatial.RData") 


###### v. 80's ----
pieris_80s <-  spain_broccoli_temps %>% 
  group_by(x, y) %>% 
  mutate(id_cell = cur_group_id()) %>% 
  filter(year(date) %in% c(1980:1989)) 

pieris_pheno_spatial_80s  <- pieris_80s %>% 
  ## pupa (1st generation)
  mutate(year = year(date),
         daily_rate  = map2_dbl(.x = daily_tmax,  
                                .y = daily_tmin,
                                ~daily_avg(tmax = .x,
                                           tmin = .y, 
                                           LDT = as_vector(traits_pieris_allstages[3,2])
                                )
         )
  ) %>% 
  group_by(id_cell, year) %>% 
  summarise(phys_time = cumsum(daily_rate),
            date,
            daily_tmin,
            daily_tmax) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]),
                                   1,
                                   0),
         life_stage = if_else(phys_time_logic == 1,
                              "pupa",
                              "unknown"),
         complete = if_else(life_stage == "pupa" & phys_time_logic == 1,
                            "yes",
                            "no")
  ) %>% 
  ungroup() %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  mutate(life_stage = if_else(phys_time_logic == 0 & complete  == "no",
                              "pupa",
                              life_stage))


save(pieris_pheno_spatial_80s, file = "~/80s_pieris_pheno_spatial.RData") 

###### vi. 90's ----
pieris_90s <-  spain_broccoli_temps %>% 
  group_by(x, y) %>% 
  mutate(id_cell = cur_group_id()) %>% 
  filter(year(date) %in% c(1990:1999)) 

pieris_pheno_spatial_90s  <- pieris_90s %>% 
  ## pupa (1st generation)
  mutate(year = year(date),
         daily_rate  = map2_dbl(.x = daily_tmax,  
                                .y = daily_tmin,
                                ~daily_avg(tmax = .x,
                                           tmin = .y, 
                                           LDT = as_vector(traits_pieris_allstages[3,2])
                                )
         )
  ) %>% 
  group_by(id_cell, year) %>% 
  summarise(phys_time = cumsum(daily_rate),
            date,
            daily_tmin,
            daily_tmax) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]),
                                   1,
                                   0),
         life_stage = if_else(phys_time_logic == 1,
                              "pupa",
                              "unknown"),
         complete = if_else(life_stage == "pupa" & phys_time_logic == 1,
                            "yes",
                            "no")
  ) %>% 
  ungroup() %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  mutate(life_stage = if_else(phys_time_logic == 0 & complete  == "no",
                              "pupa",
                              life_stage))


save(pieris_pheno_spatial_90s, file = "~/90s_pieris_pheno_spatial.RData") 

###### vii. 2000's ----
pieris_2000s <-  spain_broccoli_temps %>% 
  group_by(x, y) %>% 
  mutate(id_cell = cur_group_id()) %>% 
  filter(year(date) %in% c(2000:2009)) 

pieris_pheno_spatial_2000s  <- pieris_2000s %>% 
  ## pupa (1st generation)
  mutate(year = year(date),
         daily_rate  = map2_dbl(.x = daily_tmax,  
                                .y = daily_tmin,
                                ~daily_avg(tmax = .x,
                                           tmin = .y, 
                                           LDT = as_vector(traits_pieris_allstages[3,2])
                                )
         )
  ) %>% 
  group_by(id_cell, year) %>% 
  summarise(phys_time = cumsum(daily_rate),
            date,
            daily_tmin,
            daily_tmax) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]),
                                   1,
                                   0),
         life_stage = if_else(phys_time_logic == 1,
                              "pupa",
                              "unknown"),
         complete = if_else(life_stage == "pupa" & phys_time_logic == 1,
                            "yes",
                            "no")
  ) %>% 
  ungroup() %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  mutate(life_stage = if_else(phys_time_logic == 0 & complete  == "no",
                              "pupa",
                              life_stage))


save(pieris_pheno_spatial_2000s, file = "~/2000s_pieris_pheno_spatial.RData") 

###### viii. 2010's ----
pieris_2010s <-  spain_broccoli_temps %>% 
  group_by(x, y) %>% 
  mutate(id_cell = cur_group_id()) %>% 
  filter(year(date) %in% c(2010:2015)) 

pieris_pheno_spatial_2010s  <- pieris_2010s %>% 
  ## pupa (1st generation)
  mutate(year = year(date),
         daily_rate  = map2_dbl(.x = daily_tmax,  
                                .y = daily_tmin,
                                ~daily_avg(tmax = .x,
                                           tmin = .y, 
                                           LDT = as_vector(traits_pieris_allstages[3,2])
                                )
         )
  ) %>% 
  group_by(id_cell, year) %>% 
  summarise(phys_time = cumsum(daily_rate),
            date,
            daily_tmin,
            daily_tmax) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]),
                                   1,
                                   0),
         life_stage = if_else(phys_time_logic == 1,
                              "pupa",
                              "unknown"),
         complete = if_else(life_stage == "pupa" & phys_time_logic == 1,
                            "yes",
                            "no")
  ) %>% 
  ungroup() %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## adult (1st flight)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[4,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[4,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "adult",
                              life_stage),
         complete = if_else(life_stage == "adult" & complete == "no",
                            "yes",
                            complete)
  ) %>%
  ## eggs (1st generation eggs)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[1,2])
                                       )
                              ),
                              0)
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[1,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "egg",
                              life_stage),
         complete = if_else(life_stage == "egg" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## larva (1st generation larva)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[2,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[2,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "larva",
                              life_stage),
         complete = if_else(life_stage == "larva" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  ## pupae (2nd generation pupae)
  mutate(daily_rate = if_else(complete == "no",
                              map2_dbl(.x = daily_tmax,
                                       .y = daily_tmin,
                                       ~daily_avg(tmax = .x,
                                                  tmin = .y, 
                                                  LDT = as_vector(traits_pieris_allstages[3,2])
                                       )
                              ),
                              0
  )
  ) %>%
  group_by(id_cell, year) %>% 
  summarise(phys_time = if_else(complete == "no",
                                cumsum(daily_rate),
                                phys_time),
            date, daily_tmin, daily_tmax,
            phys_time_logic, life_stage, complete) %>% 
  mutate(phys_time_logic = if_else(phys_time < as_vector(traits_pieris_allstages[3,3]) &
                                     complete == "no",
                                   1,
                                   phys_time_logic),
         life_stage = if_else(phys_time_logic == 1 & complete == "no",
                              "pupa",
                              life_stage),
         complete = if_else(life_stage == "pupa" & complete == "no",
                            "yes",
                            complete)
  ) %>% 
  mutate(life_stage = if_else(phys_time_logic == 0 & complete  == "no",
                              "pupa",
                              life_stage))


save(pieris_pheno_spatial_2010s, file = "~/2010s_pieris_pheno_spatial.RData") 



####### ix. pieris phenology... --------------------------------------------------------
#load subsets
for(decade in c(50, 60, 70, 80, 90, 2000, 2010)){
  load(file = paste0("~/", decade, "s_pieris_pheno_spatial.RData"))
  print(decade)
}
#and merge them
pieris_pheno_spatial_all <- pieris_pheno_spatial %>% 
  bind_rows(pieris_pheno_spatial_60s,
            pieris_pheno_spatial_70s,
            pieris_pheno_spatial_80s,
            pieris_pheno_spatial_90s,
            pieris_pheno_spatial_2000s,
            pieris_pheno_spatial_2010s) %>% 
  mutate(doy = yday(date),
         photoperiod = map_dbl(.x = doy,
                               ~unlist(daylength(latitude = 40, JDay = .x))[3]) # add photoperiod
  )  
save(pieris_pheno_spatial_all, file = "~/pieris_pheno_spatial_all.RData")
#load(file = "~/broccoli_gdds.RData")
head(pieris_pheno_spatial_all)

#obtain id_cell locations
broccoli_gdds_locations <- broccoli_gdds %>% 
  rename(year = season) %>% 
  group_by(id_cell) %>% 
  summarise(x = mean(x),
            y = mean(y))
voltinism_lag <- pieris_pheno_spatial_all %>% 
  ungroup() %>% 
  mutate(life_stage = if_else(photoperiod < 9.5 & # avoid development under pupal overwintering conditions (see Shapiro 1984)
                                life_stage != "pupa",
                              "pupa",
                              life_stage)
  ) %>% 
  mutate(genlag =  if_else(life_stage == "pupa" & 
                             lag(life_stage, default = "larva") == "larva",
                           1,
                           0) #to index each generation transition starting from pupae
  )

## past voltinism (70s)

pieris_pheno_all_past <- voltinism_lag %>% 
  filter(year %in% 1966:1975) %>% 
  inner_join(broccoli_gdds_locations, by = "id_cell") %>% 
  group_by(x, y, year) %>% 
  mutate(generation = cumsum(genlag)) %>%  # compute generations
  summarise(generation = max(generation)) %>% 
  group_by(x,y) %>% 
  summarise(generation = round(mean(generation)))# calculate voltinism (i.e. number of generations)

## "present" voltinism (2010s)
pieris_pheno_all_present <- voltinism_lag %>% 
  filter(year %in% 2005:2014) %>% 
  inner_join(broccoli_gdds_locations, by = "id_cell") %>% 
  group_by(x, y, year) %>% 
  mutate(generation = cumsum(genlag)) %>%  # compute generations
  summarise(generation = max(generation)) %>% 
  group_by(x,y) %>% 
  summarise(generation = round(mean(generation)))# calculate voltinism (i.e. number of generations)

## shift_voltinism (present-past)
pieris_pheno_all_shift <- pieris_pheno_all_present %>% 
  inner_join(pieris_pheno_all_past, by = c("x", "y")) %>% 
  rename(voltinism_present = generation.x,
         voltinism_past = generation.y) %>% 
  mutate(voltinism_shift = voltinism_present - voltinism_past) %>% 
  select(x, y, voltinism_shift) 

## compute voltinism average all period
voltinism_all_period <- voltinism_lag %>% 
  inner_join(broccoli_gdds_locations, by = "id_cell") %>% 
  group_by(x, y, year) %>% 
  mutate(generation = cumsum(genlag)) %>%  # compute generations
  summarise(generation = max(generation)) %>%  # calculate voltinism (i.e. number of generations)
  group_by(x, y) %>% 
  summarise(generation = mean(generation)) %>% 
  mutate(generation = as_factor(str_extract(generation,pattern = ".")
  )
  )


voltinism_pieris_year <- voltinism_lag %>% 
  left_join(broccoli_gdds_locations, by = "id_cell") %>% 
  group_by(id_cell, year) %>% 
  mutate(generation = cumsum(genlag)) %>%  # compute generations
  summarise(generation = max(generation)) %>% 
  ungroup() %>% 
  group_by(id_cell, year) %>% 
  summarise(generation = mean(generation))

###### ... and export rasters ----
pieris_voltinism_past_rast <- terra::rast(pieris_pheno_all_past)
terra::writeRaster(pieris_voltinism_past_rast,
                   here::here(paste0("Data/pieris_voltinism_past_rast.tif")), 
                   overwrite = TRUE)

pieris_voltinism_present_rast <- terra::rast(pieris_pheno_all_present)
terra::writeRaster(pieris_voltinism_present_rast,
                   here::here(paste0("Data/pieris_voltinism_present_rast.tif")), 
                   overwrite = TRUE)

pieris_voltinism_shift_rast <- terra::rast(pieris_pheno_all_shift)
terra::writeRaster(pieris_voltinism_shift_rast,
                   here::here(paste0("Data/pieris_voltinism_shift_rast.tif")), 
                   overwrite = TRUE)


###### x. ggmap ----
####### a. present ----
spain_map_sf <- st_read("~/Dario Investigacion/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp")
xrange <- c(-10, 3)
yrange <- c(35, 44)
province_borders <- broccoli_spain_regions$geometry
ocean <- st_read("/Users/dario-ssm/Downloads/ne_10m_ocean/ne_10m_ocean.shp")
palette_voltinism_present <- colorRampPalette(c("#001219","#005f73", "#0a9396", "#94d2bd",
                                                "#e9d8a6", "#ee9b00", "#ca6702", "#bb3e03",
                                                "#ae2012", "#9b2226","#6a040f"))(max(pieris_pheno_all_present$generation))

spain_pieris_phenomap_present <- ggplot()+  
  geom_sf(data = ocean, 
          fill = "#98c1d9")+
  geom_sf(data = spain_map_sf)+
  geom_tile(data = pieris_pheno_all_present,
            aes(x, y, fill = as_factor(generation)))+
  scale_fill_manual(values = palette_voltinism_present)+
  geom_sf(data = province_borders, 
          fill = NA,
          color = "gray23")+
  coord_sf(xlim = xrange, ylim = yrange)+
  theme_map()+
  labs(fill ="Voltinism (# generations per year)",
       title = "Present (2005-2014)")+
  guides(fill = guide_legend(ncol = 14,nrow = 1, byrow = TRUE))
spain_pieris_phenomap_present
ggsave(filename = here("Data/spain_pieris_phenomap_present.png"),
       height = 14, width = 14, units = "cm", dpi = 300)
####### b. past ----
palette_voltinism_past <- colorRampPalette(c("#001219","#005f73", "#0a9396", "#94d2bd",
                                             "#e9d8a6", "#ee9b00", "#ca6702", "#bb3e03",
                                             "#ae2012", "#9b2226","#6a040f"))(max(pieris_pheno_all_present$generation))

spain_pieris_phenomap_past <- ggplot()+  
  geom_sf(data = ocean, 
          fill = "#98c1d9")+
  geom_sf(data = spain_map_sf)+
  geom_tile(data = pieris_pheno_all_past,
            aes(x, y, fill = as_factor(generation)))+
  scale_fill_manual(values = palette_voltinism_past)+
  geom_sf(data = province_borders, 
          fill = NA,
          color = "gray23")+
  coord_sf(xlim = xrange, ylim = yrange)+
  theme_map()+
  labs(fill ="Voltinism (# generations per year)",
       title = "Past (1966-1975)")+
  guides(fill = guide_legend(ncol = 14,nrow = 1, byrow = TRUE))
spain_pieris_phenomap_past
ggsave(filename = here("Data/spain_pieris_phenomap_past.png"),
       height = 14, width = 14, units = "cm", dpi = 300)
####### c. shift ----
palette_voltinism_shift <- c("#01665E","#80CDC1", "#F6E8C3","#DFC27D", "#BF812D", "#8C510A") 

spain_pieris_phenomap_shift <- ggplot()+  
  geom_sf(data = ocean, 
          fill = "#98c1d9")+
  geom_sf(data = spain_map_sf)+
  geom_tile(data = pieris_pheno_all_shift,
            aes(x, y, fill = as_factor(voltinism_shift)))+
  scale_fill_manual(values = palette_voltinism_shift)+
  geom_sf(data = province_borders, 
          fill = NA,
          color = "gray23")+
  coord_sf(xlim = xrange, ylim = yrange)+
  theme_map()+
  labs(fill ="Voltinism shift (Change in number of generations per year from past to present)",
       title = "Present Generations (2005-2014) - Past Generations (1966-1975)")+
  guides(fill = guide_legend(ncol = 14,nrow = 1, byrow = TRUE))
spain_pieris_phenomap_shift
ggsave(filename = here("Data/spain_pieris_phenomap_shift.png"),
       height = 16, width = 16, units = "cm", dpi = 300)

###### xi. glmer ----
voltinism_pieris_year
lm_voltinism <- lm(generation ~ year,
                   data = voltinism_pieris_year) 
plot(lm_voltinism)
performance::check_model(lm_voltinism) # <- okay
hist(resid(lm_voltinism), 50)  # <- okay
phenoshift_pieris_lmer <-lme4::lmer(generation ~ year + (1|id_cell),
                                    data = voltinism_pieris_year,
                                    na.action = na.exclude)

sum_phenoshift_pieris_lmer <- summary(phenoshift_pieris_lmer)
set.seed(2023)
sim_voltinism_slope <- rnorm(1000, 
                             mean = sum_phenoshift_pieris_lmer$coefficients[2,1],
                             sd = sum_phenoshift_pieris_lmer$coefficients[2,2])
sim_voltinism_intercept <- rnorm(1000, 
                                 mean = sum_phenoshift_pieris_lmer$coefficients[1,1],
                                 sd = sum_phenoshift_pieris_lmer$coefficients[1,2])
sim_voltinism <- tibble(sim_slope = sim_voltinism_slope,
                        sim_intercept = sim_voltinism_intercept)

n_draws <- 1000
alpha_level <- 0.15
col_draw <- "grey72"
col_median <- "darkcyan"
species_name <- expression(~italic("Pieris rapae"))
label_shift <- paste0(round(sum_phenoshift_pieris_lmer$coefficients[2,1]*66, 2)," generations increase")
voltinism_all_model_plot <- ggplot(voltinism_pieris_year, aes(x = year, y = generation))+
  geom_point(alpha = 0.01, color = "darkcyan")+
  geom_abline(aes(slope = sim_slope,
                  intercept = sim_intercept),
              data = slice_sample(sim_voltinism, n = n_draws),
              color = col_draw,
              alpha = alpha_level)+
  geom_abline(aes(slope = sum_phenoshift_pieris_lmer$coefficients[2,1],
                  intercept = sum_phenoshift_pieris_lmer$coefficients[1,1]),
              color = col_median,
              size = 1.4)+
  theme_few()+
  labs(title = "Voltinism ~ Year",
       subtitle = species_name,
       x = "Year",
       y = "Voltinism (# generations per year)")+
  annotate(geom = "text", x = 1992, y = 14, label = label_shift)
voltinism_all_model_plot
ggsave(filename = here("Data/voltinism_all_model_plot.png"),
       width = 12, height = 12, units = "cm")

##### d) Pieris rapae voltinism at broccoli growing season -------------------------------------------------------
load(file = here("Data/broccoli_gdds.RData"))
load(file = "~/pieris_pheno_spatial_all.RData")

hist(broccoli_gdds_harvest$harvest_day) # <- subset with harvest lower than 188 (doy 365) given the data distribution
broccoli_grow_season <- broccoli_gdds %>% 
  filter(crop_remains == 1 &
           yday(date) <= 365) %>% 
  rename(year = season)

pieris_pheno_broccoli <- broccoli_grow_season %>% 
  inner_join(pieris_pheno_spatial_all, by = c("id_cell", "year", "date")) %>% 
  mutate(life_stage = if_else(photoperiod < 9.5 & # avoid development under pupal overwintering conditions (see Shapiro 1984)
                                life_stage != "pupa",
                              "pupa",
                              life_stage)
  ) %>% 
  mutate(genlag =  if_else(life_stage == "larva" & 
                             lag(life_stage, default = "egg") == "egg",
                           1,
                           0) #to index each generation transition starting from pupae
  )

## past voltinism (70s)
pieris_pheno_broccoli_past <- pieris_pheno_broccoli %>% 
  filter(year %in% 1966:1975) %>% 
  group_by(x, y, year) %>% 
  mutate(generation = cumsum(genlag)) %>%  # compute generations
  summarise(generation = max(generation)) # calculate voltinism (i.e. number of generations)

## "present" voltinism (2010s)
pieris_pheno_broccoli_present <- pieris_pheno_broccoli %>% 
  filter(year %in% 2005:2014) %>% 
  group_by(x, y, year) %>% 
  mutate(generation = cumsum(genlag)) %>%  # compute generations
  summarise(generation = max(generation)) # calculate voltinism (i.e. number of generations)

## shift_voltinism (present-past)
pieris_pheno_broccoli_shift <- pieris_pheno_broccoli_present %>% 
  inner_join(pieris_pheno_broccoli_past, by = c("x", "y")) %>% 
  rename(voltinism_present = generation.x,
         voltinism_past = generation.y) %>% 
  mutate(voltinism_shift = voltinism_present - voltinism_past) %>% 
  select(x, y, voltinism_shift) 
###### ... and export rasters ----
pieris_pheno_broccoli_past_rast <- terra::rast(pieris_pheno_broccoli_past)
terra::writeRaster(pieris_pheno_broccoli_past_rast,
                   here::here(paste0("Data/pieris_pheno_broccoli_past_rast.tif")), 
                   overwrite = TRUE)

pieris_pheno_broccoli_present_rast <- terra::rast(pieris_pheno_broccoli_present)
terra::writeRaster(pieris_pheno_broccoli_present_rast,
                   here::here(paste0("Data/pieris_pheno_broccoli_present_rast.tif")), 
                   overwrite = TRUE)

pieris_pheno_broccoli_shift_rast <- terra::rast(pieris_pheno_broccoli_shift)
terra::writeRaster(pieris_pheno_broccoli_shift_rast,
                   here::here(paste0("Data/pieris_pheno_broccoli_shift_rast.tif")), 
                   overwrite = TRUE) 
###### ii. ggmap ----
####### a. present ----
spain_map_sf <- st_read("~/Dario Investigacion/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp")
xrange <- c(-10, 3)
yrange <- c(35, 44)
province_borders <- broccoli_spain_regions$geometry
ocean <- st_read("/Users/dario-ssm/Downloads/ne_10m_ocean/ne_10m_ocean.shp")
palette_voltinism_pres_broc <- c("#94d2bd", colorRampPalette(brewer.pal(8,"YlOrRd"))(max(pieris_pheno_broccoli_present$generation)))

spain_pieris_phenomap_broccoli <- ggplot()+  
  geom_sf(data = ocean, 
          fill = "#98c1d9")+
  geom_sf(data = spain_map_sf)+
  geom_tile(data = pieris_pheno_broccoli_present,
            aes(x, y, fill = as_factor(generation)))+
  scale_fill_manual(values = palette_voltinism_pres_broc)+
  geom_sf(data = province_borders, 
          fill = NA,
          color = "gray23")+
  coord_sf(xlim = xrange, ylim = yrange)+
  theme_map()+
  labs(fill ="Present Voltinism (# generations per year during broccoli growth)",
       title = "Present (2005-2014)")+
  guides(fill = guide_legend(ncol = 14,nrow = 1, byrow = TRUE))
spain_pieris_phenomap_broccoli
ggsave(filename = here("Data/pieris_phenomap_broccoli_present.png"),
       height = 14, width = 14, units = "cm", dpi = 300)
####### b. past ----
palette_voltinism_past_broc <- c("#94d2bd", colorRampPalette(brewer.pal(8,"YlOrRd"))(max(pieris_pheno_broccoli_past$generation)))

spain_pieris_phenomap_broccoli_past <- ggplot()+  
  geom_sf(data = ocean, 
          fill = "#98c1d9")+
  geom_sf(data = spain_map_sf)+
  geom_tile(data = pieris_pheno_broccoli_past,
            aes(x, y, fill = as_factor(generation)))+
  scale_fill_manual(values = palette_voltinism_past_broc)+
  geom_sf(data = province_borders, 
          fill = NA,
          color = "gray23")+
  coord_sf(xlim = xrange, ylim = yrange)+
  theme_map()+
  labs(fill ="Past Voltinism (# generations per year during broccoli growth)",
       title = "Past (1966-1975)")+
  guides(fill = guide_legend(ncol = 14,nrow = 1, byrow = TRUE))
spain_pieris_phenomap_broccoli_past
ggsave(filename = here("Data/pieris_phenomap_broccoli_past.png"),
       height = 14, width = 14, units = "cm", dpi = 300)

######## c. shift ----
palette_voltinism_shift_broc <- rev(colorRampPalette(brewer.pal(9,"BrBG"))(length(unique(pieris_pheno_broccoli_shift$voltinism_shift))
)
)
spain_pieris_phenomap_broccoli_shift <- ggplot()+  
  geom_sf(data = ocean, 
          fill = "#98c1d9")+
  geom_sf(data = spain_map_sf)+
  geom_tile(data = pieris_pheno_broccoli_shift,
            aes(x, y, fill = as_factor(voltinism_shift)))+
  scale_fill_manual(values = palette_voltinism_shift_broc)+
  geom_sf(data = province_borders, 
          fill = NA,
          color = "gray23")+
  coord_sf(xlim = xrange, ylim = yrange)+
  theme_map()+
  labs(fill ="Voltinism shift (Change in number of generations during broccoli growing season from past to present",
       title = "Present Generations (2005-2014) - Past Generations (1966-1975)")+
  guides(fill = guide_legend(ncol = 14,nrow = 1, byrow = TRUE))
spain_pieris_phenomap_broccoli_shift
ggsave(filename = here("Data/pieris_phenomap_broccoli_shift.png"),
       height = 16, width = 16, units = "cm", dpi = 300)

###### d. glmer ----
voltinism_pieris_broccoli_year <- pieris_pheno_broccoli %>% 
  group_by(id_cell, year) %>% 
  mutate(generation = cumsum(genlag)) %>%  # compute generations
  summarise(generation = max(generation)) %>% 
  ungroup() %>% 
  group_by(id_cell, year) %>% 
  summarise(generation = mean(generation))

lm_voltinism_broccoli <- lm(generation ~ year,
                            data = voltinism_pieris_broccoli_year) 
hist(resid(lm_voltinism_broccoli), 20)  # <- okay
# performance::check_model(lm_voltinism_broccoli)

phenoshift_pieris_broccoli_lmer <-lme4:: lmer(generation ~ year + (1|id_cell),
                                              data = voltinism_pieris_broccoli_year,
                                              na.action = na.exclude)
sum_phenoshift_pieris_broccoli_lmer <- summary(phenoshift_pieris_broccoli_lmer)
set.seed(2023)
sim_voltinism_broc_slope <- rnorm(1000, 
                                  mean = sum_phenoshift_pieris_broccoli_lmer$coefficients[2,1],
                                  sd = sum_phenoshift_pieris_broccoli_lmer$coefficients[2,2])
sim_voltinism_broc_intercept <- rnorm(1000, 
                                      mean = sum_phenoshift_pieris_broccoli_lmer$coefficients[1,1],
                                      sd = sum_phenoshift_pieris_broccoli_lmer$coefficients[1,2])
sim_voltinism_broc <- tibble(sim_slope = sim_voltinism_broc_slope,
                             sim_intercept = sim_voltinism_broc_intercept)

n_draws <- 1000
alpha_level <- 0.15
col_draw <- "grey72"
col_median <- "firebrick4"
species_name <- expression(~italic("Pieris rapae"))
label_shift <- paste0(round(sum_phenoshift_pieris_broccoli_lmer$coefficients[2,1]*66, 2)," generations increase")

voltinism_broc_model_plot <- ggplot(voltinism_pieris_broccoli_year, aes(x = year, y = generation))+
  geom_point(alpha = 0.01, color = "firebrick2")+
  geom_abline(aes(slope = sim_slope,
                  intercept = sim_intercept),
              data = slice_sample(sim_voltinism_broc, n = n_draws),
              color = col_draw,
              alpha = alpha_level)+
  geom_abline(aes(slope = sum_phenoshift_pieris_broccoli_lmer$coefficients[2,1],
                  intercept = sum_phenoshift_pieris_broccoli_lmer$coefficients[1,1]),
              color = col_median,
              size = 1.4)+
  theme_few()+
  labs(title = "Voltinism ~ Year",
       subtitle = species_name,
       x = "Year",
       y = "Voltinism (# generations per broccoli growing season each year)")+
  annotate(geom = "text", x = 1986, y = 5, label = label_shift)
voltinism_broc_model_plot
ggsave(filename = here("Data/voltinism_broc_model_plot.png"),
       width = 12, height = 12, units = "cm")



