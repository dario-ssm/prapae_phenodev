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
library(rTPC)
library(nls.multstart)
library(nlme)
library(ggthemes)
library(chillR)
library(sf)
library(viridis)
source(here("Scripts/1-functions_phenobraspests.R"))
load(here("Data/daily_tmax_df.RData"))
print(daily_tmax_df)
load(here("Data/daily_tmin_df.RData"))
print(daily_tmin_df)
pieris_data <- read_delim(here("Data/pieris_devdata/gilbert_pupa.csv"),
                          delim = ";") %>% 
  rename(stage = ...3,
         parasite = ...4, 
         reference = interact)

# 1. Models and parameters ----
##### a) Given degree-days ------------------------------------------------------
values_Samson1984_linear <- function(){
  ldt= 9
  udt=30
  heat_units = 192
  life_stage = "larvae"
  samson1984 <- tibble(life_stage, ldt, udt, heat_units)
  return(print(samson1984))
}
values_Samson1984_linear()

values_Davies1985a <- function(){
  ldt= c(10,10,9.3)
  heat_units = c(54,157,116)
  life_stage = c("eggs","larvae","pupae")
  davies1985 <- tibble(life_stage,ldt,heat_units)
  return(print(davies1985))
}
values_Davies1985a()

#### b) Fitting curves to devdata ----
###### i. pupa data ----
pieris_data_pupa <- pieris_data %>% 
  filter(stage == "pupa")
plot(pieris_data_pupa$temperature,
     pieris_data_pupa$devrate)
## for start values, we will use the estimates of the previous model fitting
#briere-2
start_vals_br2 <- rTPC::get_start_vals(x = pieris_data_pupa$temperature,
                                       y = pieris_data_pupa$devrate,
                                       model_name = "briere2_1999")
briere2_pieris <- gnls(devrate ~ briere2(temp = temperature,
                                         tmin, tmax, a, b),
                       data = pieris_data_pupa,
                       start = start_vals_br2,
                       na.action = na.exclude,
                       control = gnlsControl(nlsTol = 1e-06)
)

summary(briere2_pieris)

#briere-1
start_vals_br1 <- c(tmin = 8, tmax = 33, a = 2e-04)
briere1_pieris <- gnls(devrate ~ briere1(temp = temperature,
                                         tmin, tmax, a),
                       data = pieris_data_pupa,
                       start = start_vals_br1,
                       na.action = na.exclude,
                       control = gnlsControl(nlsTol = 1e-06)
)
summary(briere1_pieris)

#lactin-2
start_vals_vague_lactin2 <- devRate::devRateEqStartVal[["lactin2_95"]] # take literature start values from devRate
names(start_vals_vague_lactin2) <- c("a", "tmax", "delta_t", "b")
start_upper_vals <- purrr::map(.x = start_vals_vague_lactin2,
                               .f = ~.x + abs(.x))
start_lower_vals <- purrr::map(.x = start_vals_vague_lactin2,
                               .f = ~.x - abs(.x))

multstart_vals_fit_lactin2 <- nls.multstart::nls_multstart(formula = devrate~lactin2(temperature, a, tmax, delta_t, b),
                                                        data = pieris_data_pupa,
                                                        iter = 500,
                                                        start_lower = start_lower_vals,
                                                        start_upper = start_upper_vals,
                                                        supp_errors = "Y")
sum_start_vals_fit_lactin2 <- summary(multstart_vals_fit_lactin2)
start_lactin2 <- sum_start_vals_fit_lactin2$parameters[,1]
lactin2_pieris <- gnls(devrate ~ lactin2(temperature, a, tmax, delta_t, b),
                    data = pieris_data_pupa,
                    start = start_lactin2,
                    na.action = na.exclude,
                    control = gnlsControl(nlsTol = 1e-03))
summary(lactin2_pieris)

#lactin-1
start_vals_lactin1 <- c(a = 0.16, tmax = 36.49, delta_t = 6.22) # <-  devRate::lactin1_95$startVal

lactin1_pieris<- gnls(devrate ~ lactin1(temp = temperature,
                                         a, tmax, delta_t),
                       data = pieris_data_pupa,
                       start = start_vals_lactin1,
                       control = gnlsControl(nlsTol = 1e-05))
summary(lactin1_pieris)


#rezende

start_vals_rezende <- c(q_10 = 3.789, cte = 0.008, thr = 20.96, decay = 0.0051)
rezende_pieris <- gnls(devrate ~ rezende(temp = temperature,
                                         q_10, cte, thr, decay),
                       data = pieris_data_pupa,
                       start = start_vals_rezende,
                       na.action = na.exclude,
                       weights = varExp(form = ~temperature),
                       control = gnlsControl(nlsTol = 1e-06))
summary(rezende_pieris)
#ssilow
start_vals_vague_ssilow <- devRate::devRateEqStartVal[["schoolfieldLow_81"]] # take literature start values from devRate
names(start_vals_vague_ssilow) <- c("p25", "a", "b", "c")
start_upper_vals <- purrr::map(.x = start_vals_vague_ssilow,
                               .f = ~.x + abs(.x/0.5))
start_lower_vals <- purrr::map(.x = start_vals_vague_ssilow,
                               .f = ~.x - abs(.x/0.5))
multstart_vals_fit_ssilow <- nls.multstart::nls_multstart(formula = devrate~ssilow(temperature,  p25, a, b, c),
                                                           data = pieris_data_pupa,
                                                           iter = 500,
                                                           start_lower = start_lower_vals,
                                                           start_upper = start_upper_vals,
                                                           supp_errors = "Y")
sum_start_vals_fit_ssilow <- summary(multstart_vals_fit_ssilow)
start_ssilow <- sum_start_vals_fit_ssilow$parameters[,1]
ssilow_pieris <- gnls(devrate ~ssilow(temperature,  p25, a, b, c),
                       data = pieris_data_pupa,
                       start = start_ssilow,
                       na.action = na.exclude,
                       weights = varExp(form = ~temperature),
                       control = gnlsControl(nlsTol = 1e-06))
summary(ssilow_pieris)

#weibull
start_vals_weibull <- rTPC::get_start_vals(x = pieris_data_pupa$temperature,
                                           y = pieris_data_pupa$devrate,
                                           model_name = "weibull_1995")
weibull_pieris <- gnls(devrate ~ weibull_1995(temp = temperature,
                                              a, topt, b, c),
                       data = pieris_data_pupa,
                       start = start_vals_weibull,
                       na.action = na.exclude,
                       control = gnlsControl(nlsTol = 1e-02)
)
summary(weibull_pieris)

#beta
start_vals_beta <- rTPC::get_start_vals(x = pieris_data_pupa$temperature,
                                        y = pieris_data_pupa$devrate,
                                        model_name = "beta_2012")

beta_pieris <- gnls(devrate ~ beta_2012(temp = temperature,
                                        a, b, c, d, e),
                    data = pieris_data_pupa,
                    start = start_vals_beta,
                    na.action = na.exclude,
                    control = gnlsControl(nlsTol = 1e-03))
summary(beta_pieris)

#gaussian
start_vals_gaussian <- rTPC::get_start_vals(x = pieris_data_pupa$temperature,
                                            y = pieris_data_pupa$devrate,
                                            model_name = "gaussian_1987")
gaussian_pieris <- gnls(devrate ~ gaussian_1987(temp = temperature,
                                                rmax, topt, a),
                        data = pieris_data_pupa,
                        start = start_vals_gaussian,
                        na.action = na.exclude,
                        control = gnlsControl(nlsTol = 1e-06))
summary(gaussian_pieris)

#oneill
start_vals_oneill <- rTPC::get_start_vals(x = pieris_data$temperature,
                                          y = pieris_data$devrate,
                                          model_name = "oneill_1972")

oneill_pieris <- gnls(devrate ~ oneill_1972(temp = temperature,
                                            rmax, ctmax, topt, q10),
                      data = pieris_data,
                      start = start_vals_oneill,
                      na.action = na.exclude,
                      control = gnlsControl(nlsTol = 1e-01))
summary(oneill_pieris)

#wang
start_vals_vague_wang <- devRate::devRateEqStartVal[["wang_82"]] # take literature start values from devRate
names(start_vals_vague_wang) <- c("k", "r", "topt", "tmin", "tmax", "a")
start_upper_vals <- purrr::map(.x = start_vals_vague_wang,
                               .f = ~.x + abs(.x/.1))
start_lower_vals <- purrr::map(.x = start_vals_vague_wang,
                               .f = ~.x - abs(.x/.1))

multstart_vals_fit_wang <- nls.multstart::nls_multstart(formula = devrate~wang(temperature, k, r, topt, tmin, tmax, a),
                                                        data = pieris_data_pupa,
                                                        iter = 500,
                                                        start_lower = start_lower_vals,
                                                        start_upper = start_upper_vals,
                                                        supp_errors = "Y")
sum_start_vals_fit_wang <- summary(multstart_vals_fit_wang)
start_wang <- sum_start_vals_fit_wang$parameters[,1]
# wang_pieris <- gnls(devrate ~ wang(temperature, k, r, topt, tmin, tmax, a),
#                     data = pieris_data_pupa,
#                     start = start_wang,
#                     na.action = na.exclude,
#                     control = gnlsControl(nlsTol = 1e-03))
# summary(wang_pieris) #not converging

#mod_polynomial
start_vals_vague_poly <- devRate::devRateEqStartVal[["poly4"]] # take literature start values from devRate
names(start_vals_vague_poly) <- c("a_0", "a_1", "a_2", "a_3", "a_4")
start_upper_vals <- purrr::map(.x = start_vals_vague_poly,
                               .f = ~.x + abs(.x/.1))
start_lower_vals <- purrr::map(.x = start_vals_vague_poly,
                               .f = ~.x - abs(.x/.1))

multstart_vals_fit_poly <- nls.multstart::nls_multstart(formula = devrate~mod_polynomial(temperature, a_0, a_1, a_2, a_3, a_4),
                                                        data = pieris_data_pupa,
                                                        iter = 500,
                                                        start_lower = start_lower_vals,
                                                        start_upper = start_upper_vals,
                                                        supp_errors = "Y")
sum_start_vals_fit_poly <- summary(multstart_vals_fit_poly)
start_poly <- sum_start_vals_fit_poly$parameters[,1]
mod_poly_pieris <- gnls(devrate ~ mod_polynomial(temp = temperature,
                                                 a_0, a_1, a_2, a_3, a_4),
                        data = pieris_data_pupa,
                        start = start_poly,
                        na.action = na.exclude,
                        control = gnlsControl(nlsTol = 1e-06))
summary(mod_poly_pieris)

#stinner
start_vals_vague_stinner <- unlist(devRate::devRateEqStartVal[["stinner_74"]]) # take literature start values from devRate
names(start_vals_vague_stinner) <- c("c", "k1", "k2", "topt")
start_upper_vals <- purrr::map(.x = start_vals_vague_stinner,
                               .f = ~.x + abs(.x))
start_lower_vals <- purrr::map(.x = start_vals_vague_stinner,
                               .f = ~.x - abs(.x))

multstart_vals_fit_stinner <- nls.multstart::nls_multstart(formula = devrate~stinner1974(temperature, c, k1, k2, topt),
                                                          data = pieris_data_pupa,
                                                          iter = 500,
                                                          start_lower = start_lower_vals,
                                                          start_upper = start_upper_vals,
                                                          supp_errors = "Y")
sum_start_vals_fit_stinner <- summary(multstart_vals_fit_stinner)
start_stinner <- sum_start_vals_fit_stinner$parameters[,1]
stinner_pieris <- gnls(devrate ~ stinner1974(temp = temperature,
                                             c, k1, k2, topt),
                      data = pieris_data_pupa,
                      start = start_stinner,
                      na.action = na.exclude,
                      weights = varExp(form = ~temperature),
                      control = gnlsControl(nlsTol = 1e-02))
summary(stinner_pieris) #large errors

#shi2011
start_vals_vague_shi2011 <- devRate::devRateEqStartVal[["shi_11"]] # take literature start values from devRate
names(start_vals_vague_shi2011) <- c("c", "k1",  "t1", "k2", "t2")
start_upper_vals <- purrr::map(.x = start_vals_vague_shi2011,
                               .f = ~.x + abs(.x/.1))
start_lower_vals <- purrr::map(.x = start_vals_vague_shi2011,
                               .f = ~.x - abs(.x/.1))

multstart_vals_fit_shi2011 <- nls.multstart::nls_multstart(formula = devrate~shi2011(temperature, c, k1, t1, k2, t2),
                                                          data = pieris_data_pupa,
                                                          iter = 500,
                                                          start_lower = start_lower_vals,
                                                          start_upper = start_upper_vals,
                                                          supp_errors = "Y")
sum_start_vals_fit_shi2011 <- summary(multstart_vals_fit_shi2011)
start_shi2011 <- sum_start_vals_fit_shi2011$parameters[,1]
shi2011_pieris <- gnls(devrate ~ shi2011(temp = temperature,
                                         c, k1, t1, k2, t2),
                      data = pieris_data_pupa,
                      start = start_shi2011,
                      na.action = na.exclude,
                      weights = varExp(form = ~temperature),
                      control = gnlsControl(nlsTol = 1e-03))
summary(shi2011_pieris) # large errors



#linear
linear_pieris <- lm(devrate~temperature, 
                    data = pieris_data_pupa)

##comparisons

anova(briere1_pieris, briere2_pieris, lactin2_pieris, lactin1_pieris, ssilow_pieris,
      weibull_pieris, gaussian_pieris, mod_poly_pieris, rezende_pieris, linear_pieris)

## combine them
ssilow_fit <- tibble(model_name = "ssilow",
                  nls_model = list(ssilow_pieris))
br1_fit <- tibble(model_name = "briere1",
                  nls_model = list(briere1_pieris))
br2_fit <- tibble(model_name = "briere2",
                  nls_model = list(briere2_pieris))
lactin2_fit <- tibble(model_name = "lactin2",
                      nls_model = list(lactin2_pieris))
lactin1_fit <- tibble(model_name = "lactin1",
                      nls_model = list(lactin1_pieris))
weibull_fit <- tibble(model_name = "weibull",
                        nls_model = list(weibull_pieris))
mod_poly_fit <- tibble(model_name = "mod_poly",
                      nls_model = list(mod_poly_pieris))
rezende_fit <- tibble(model_name = "rezende",
                     nls_model = list(rezende_pieris))
gaussian_fit <- tibble(model_name = "gaussian",
                       nls_model = list(gaussian_pieris))
linear_fit <- tibble(model_name = "linear",
                     nls_model = list(linear_pieris))

extract_se_model <- function(model_fit){
  model_summary <- summary(model_fit)
  model_se <- model_summary$tTable[,2]
  return(model_se)
}

model_outputs <- bind_rows(ssilow_fit,
                           br1_fit,
                           br2_fit,
                           lactin2_fit,
                           lactin1_fit,
                           weibull_fit,
                           mod_poly_fit,
                           rezende_fit,
                           gaussian_fit,
                           linear_fit) |>  
  mutate(est = map(.x = nls_model,
                   .f = coef),
         param = map(.x = nls_model, 
                     .f = param_names_extractor),
         se = map(.x = nls_model,
                  .f = extract_se_model)) |> 
  unnest(c(est, se, param))

seq_temps <- seq(0, 48, 0.01)

model_preds <- bind_rows(ssilow_fit,
                         br1_fit,
                         br2_fit,
                         lactin2_fit,
                         lactin1_fit,
                         weibull_fit,
                         mod_poly_fit,
                         rezende_fit,
                         gaussian_fit,
                         linear_fit) %>% 
  mutate(preds = case_when(model_name == "ssilow" ~ list(ssilow(seq_temps, 
                                                                p25 = ssilow_pieris$coefficients[1],
                                                                a = ssilow_pieris$coefficients[2],
                                                                b = ssilow_pieris$coefficients[3],
                                                                c = ssilow_pieris$coefficients[4]
                                                               )
                                                        ),
                           model_name == "briere1"  ~ list(briere1(seq_temps,
                                                                   tmin = briere1_pieris$coefficients[1],
                                                                   tmax = briere1_pieris$coefficients[2],
                                                                   a = briere1_pieris$coefficients[3]
                                                                   )
                                                           ),
                           model_name == "briere2"  ~ list(briere2(seq_temps,
                                                                   tmin = briere2_pieris$coefficients[1],
                                                                   tmax = briere2_pieris$coefficients[2],
                                                                   a = briere2_pieris$coefficients[3],
                                                                   b = briere2_pieris$coefficients[4]
                                                                   )
                                                           ),
                           model_name == "lactin2" ~ list(lactin2(seq_temps,
                                                                  a = lactin2_pieris$coefficients[1],
                                                                  tmax = lactin2_pieris$coefficients[2],
                                                                  delta_t = lactin2_pieris$coefficients[3],
                                                                  b = lactin2_pieris$coefficients[4]
                                                                  )
                                                          ),
                           model_name == "lactin1" ~ list(lactin1(seq_temps,
                                                                  a = lactin1_pieris$coefficients[1],
                                                                  tmax = lactin1_pieris$coefficients[2],
                                                                  delta_t = lactin1_pieris$coefficients[3]
                                                                  )
                                                          ),
                           model_name == "mod_poly" ~ list(mod_polynomial(seq_temps,
                                                                          a_0 = mod_poly_pieris$coefficients[1],
                                                                          a_1 = mod_poly_pieris$coefficients[2],
                                                                          a_2 = mod_poly_pieris$coefficients[3],
                                                                          a_3 = mod_poly_pieris$coefficients[4],
                                                                          a_4 = mod_poly_pieris$coefficients[5]
                                                                          )
                                                           ),
                           model_name == "rezende" ~ list(rezende(seq_temps,
                                                                  q_10 = rezende_pieris$coefficients[1],
                                                                  cte = rezende_pieris$coefficients[2],
                                                                  thr = rezende_pieris$coefficients[3],
                                                                  decay = rezende_pieris$coefficients[4]
                                                                  )
                                                          ),
                           model_name == "gaussian" ~ list(gaussian_1987(seq_temps,
                                                                         rmax = gaussian_pieris$coefficients[1],
                                                                         topt = gaussian_pieris$coefficients[2],
                                                                         a = gaussian_pieris$coefficients[3])),
                           model_name == "weibull" ~ list(weibull_1995(seq_temps,
                                                                       a = weibull_pieris$coefficients[1],
                                                                       topt = weibull_pieris$coefficients[2],
                                                                       b = weibull_pieris$coefficients[3],
                                                                       c = weibull_pieris$coefficients[4])),
                          model_name == "linear" ~ list(linear(seq_temps,
                                                                intercept = linear_pieris$coefficients[1],
                                                                slope = linear_pieris$coefficients[2])
                           )
  )
  ) %>% 
  unnest(preds) %>%
  group_by(model_name) %>% 
  mutate(temperature = seq_temps) %>%
  filter(preds >= 0 & preds <1)


aic_scores <- model_preds |> 
  pull(nls_model) |> 
  map_dbl(AIC)

model_preds_aic <- model_preds |> 
  ungroup() |> 
  mutate(aic_score = aic_scores)

aic_text <-  model_preds_aic  |>
  group_by(model_name)  |>
  summarise(aic = unique(aic_score)) |> 
  arrange(aic)
aic_order <- aic_text  |>
  pull(model_name)
aic_values <- aic_text |>
  mutate(aic =   paste("AIC =",
                       round(aic, 2)),
         temp = min(pieris_data_pupa$temperature),
         preds = 1.5*max(pieris_data_pupa$devrate))


## and ggplot-it
label_facets_num <- function(string){
  len <- length(string)
  string = paste('(', 1:len, ') ', string, sep = '')
  return(string)
}
species_name <- expression(~italic("Pieris rapae"))
my_pallette <- colorRampPalette(colors = mordechai <- c("#E9C86B", "#808D5C", "#9C4E64", "#3A848B", "#E85038", "#F9B69C")
)
pieris_all_curves_pupa  <-   ggplot()+
  geom_point(data = pieris_data_pupa, 
             aes(x = temperature,
                 y = devrate),
             color = "gray62", 
             alpha = .5)+
  geom_line(data = model_preds_aic |>
              filter(preds < 0.3),
            aes(x = temperature, y = preds, color = model_name),
            linewidth = 1.3)+
  facet_wrap(~factor(model_name, levels = aic_order))+
  theme(legend.position = "none")+
  labs(x = "Temperature (ºC)", y = "Development rate (1/days)")+
  ggthemes::theme_few()+
  scale_color_manual(values = my_pallette(length(unique(model_preds_aic$model_name))))+
  scale_fill_manual(values = my_pallette(length(unique(model_preds_aic$model_name))))+
  theme(legend.position = "none")+
  labs(x = "Temperature (ºC)", 
       y = "Development rate (1/days)",
       title = species_name,
       subtitle = "Pupa")+
  geom_label(data = aic_values,
             aes(label = aic,
                 x = temp,
                 y = preds,
                 fill = model_name),
             color = "white",
             size = 3)
pieris_all_curves_pupa
ggsave(filename = here("Data/pieris_all_curves_pupa.png"),
       width = 20,
       height = 15, 
       units = "cm",
       dpi = 300)
ggsave(filename = here("Data/pieris_all_curves_pupa.svg"),
       width = 20,
       height = 15, 
       units = "cm",
       dpi = 300)

pieris_all_curves_nogrid <-  ggplot()+
  geom_point(data = pieris_data_pupa, 
             aes(x = temperature,
                 y = devrate),
             color = "gray62", 
             alpha = .4)+
  geom_line(data = model_preds_aic |>
              filter(preds < (1.5*max(model_preds_aic$preds))),
            aes(x = temperature, y = preds, color = model_name),
            linewidth = 1.3)+
  theme(legend.position = "none")+
  labs(x = "Temperature (ºC)", y = "Development rate (1/days)")+
  theme_light()+
  theme(legend.position = "bottom")+
  labs(x = "Temperature (ºC)", 
       y = "Development rate (1/days)",
       title = species_name,
       subtitle = "Pupa")+
  scale_color_viridis_d()
pieris_all_curves_nogrid
ggsave(here("Data/pieris_all_curves_overlap.png"),
       width = 16, height = 16, units= "cm")
# 2. Linear rate summation trends ----
dds_pupa_fitted <-  model_outputs %>% 
  filter(model_name == "linear" &
           param == "temperature") %>% 
  select(est) %>% 
  mutate(est = 1/est) %>% 
  as_vector()

ldt_pupa_fitted <-  model_outputs %>% 
  filter(model_name == "linear") %>%
  pivot_wider(values_from = est,
              names_from = param) %>%
  rename(intercept = `(Intercept)`, 
         slope = temperature) %>% 
  mutate(ldt = -intercept/slope) %>% 
  pull(ldt) 


dds_pupa_davies1985 <- values_Davies1985a() %>% 
  filter(life_stage == "pupae") %>% 
  select(heat_units) %>% 
  as_vector()
ldt_pupa_davies1985 <- values_Davies1985a() %>% 
  filter(life_stage == "pupae") %>% 
  select(ldt) %>% 
  as_vector()

##### a) heat_units ------------------------------------------------------
daily_dds_pupa <- inner_join(daily_tmin_df, daily_tmax_df) %>% 
  mutate(year = year(date),
         temperature = map2_dbl(.x = daily_tmin,
                               .y = daily_tmax,
                               .f = ~mean(c(.x, .y))
         )
  ) %>% 
  group_by(year) %>% 
  mutate(ratesum_dailyavg_fitted = map2_dbl(.x = daily_tmin,
                                            .y = daily_tmax,
                                            .f = ~ daily_avg(.x, .y,
                                                             LDT = ldt_pupa_fitted)
  ),
  ratesum_sinewave_fitted = map2_dbl(.x = daily_tmin,
                                     .y = daily_tmax,
                                     .f = ~ sinewave_ldt(.x, .y,
                                                         LDT = ldt_pupa_fitted)
  ),
  ratesum_dailyavg_davies = map2_dbl(.x = daily_tmin,
                                     .y = daily_tmax,
                                     .f = ~ daily_avg(.x, .y,
                                                      LDT = ldt_pupa_davies1985)
  ),
  
  ratesum_sinewave_davies = map2_dbl(.x = daily_tmin,
                                     .y = daily_tmax,
                                     .f = ~ sinewave_ldt(.x, .y,
                                                         LDT = ldt_pupa_davies1985)
  )
  ) %>% 
  mutate(across(c(5:8), ~cumsum(.x)))
##### b) day of emergence linear ------------------------------------------------------
doe_linear_pupa <- daily_dds_pupa %>%
  mutate(doy = yday(date)) %>% 
  mutate(across(c(5,6), 
                ~logic_dd(heat_units = dds_pupa_fitted,
                          rate_cumsum = .x)
  ),
  across(c(7,8),
         ~logic_dd(heat_units = dds_pupa_davies1985,
                   rate_cumsum = .x))
  ) %>% 
  group_by(year) %>% 
  summarise(across(c(5:8),
                   ~sum(.x))
  ) %>% 
  rename(doe_dailyavg_fitted = ratesum_dailyavg_fitted,
         doe_sinewave_fitted = ratesum_sinewave_fitted,
         doe_dailyavg_davies = ratesum_dailyavg_davies,
         doe_sinewave_davies = ratesum_sinewave_davies) %>% 
  pivot_longer(cols = -1,
               names_to = "method",
               values_to = "doe") %>% 
  filter(str_sub(method,-6) == "fitted") %>%  
  mutate(method = str_sub(method, 5, -8))

# and gg-plot it
doe_linear_trends <- ggplot(doe_linear_pupa,
                            aes(x = year,
                                y = doe,
                                color = method,
                                fill = method))+
  geom_point()+
  geom_line(linetype = "dashed")+
  geom_smooth()+
  theme_clean()+
  labs(x = "Year",
       y = "Day of first adult emergence")
doe_linear_trends
ggsave(filename = here("Data/doe_linear_trends.png"),
       width = 16,
       height = 16,
       units = "cm")


# 3. Nonlinear rate summation trends ---------------------------------------

## using gaussian, ssilow, oneill, weibull, briere1, ratkowsky, lactin2, briere2 and linear, 
params_briere1 <- model_outputs %>% 
  filter(model_name == "briere1") %>% 
  pull(est)
params_briere2 <- model_outputs %>% 
  filter(model_name == "briere2") %>% 
  pull(est)
params_weibull <- model_outputs %>% 
  filter(model_name == "weibull") %>% 
  pull(est)
params_gaussian <- model_outputs %>% 
  filter(model_name == "gaussian") %>% 
  pull(est)
params_ssilow <- model_outputs %>% 
  filter(model_name == "ssilow") %>% 
  pull(est)
params_rezende<- model_outputs %>% 
  filter(model_name == "rezende") %>% 
  pull(est)
params_mod_poly<- model_outputs %>% 
  filter(model_name == "mod_poly") %>% 
  pull(est)
params_lactin2 <- model_outputs %>% 
  filter(model_name == "lactin2") %>% 
  pull(est)
params_lactin1 <- model_outputs %>% 
  filter(model_name == "lactin1") %>% 
  pull(est)

##### a) daily temps ----------------------------------------------------------

daily_ratesum_pupa <- inner_join(daily_tmin_df, daily_tmax_df) %>% 
  mutate(year = year(date),
         temperature = map2_dbl(.x = daily_tmin,
                               .y = daily_tmax,
                               .f = ~mean(c(.x, .y))
         )
  ) %>% 
  mutate(rate_summation_briere1 = map_dbl(.x = temperature,
                                          .f = ~ briere1(temp = .x,
                                                         tmin = params_briere1[1],
                                                         tmax = params_briere1[2],
                                                         a = params_briere1[3])
  ),
  rate_summation_briere2 = map_dbl(.x = temperature,
                                   .f = ~ briere2(temp = .x,
                                                  tmin = params_briere2[1],
                                                  tmax = params_briere2[2],
                                                  a = params_briere2[3],
                                                  b = params_briere2[4])
  ),
  rate_summation_lactin2 = map_dbl(.x = temperature,
                                   .f = ~ lactin2(temp = .x,
                                                  a = params_lactin2[1],
                                                  tmax = params_lactin2[2],
                                                  delta_t = params_lactin2[3],
                                                  b = params_lactin2[4])
  ),
  rate_summation_ssilow = map_dbl(.x = temperature,
                                  .f = ~ ssilow(temp = .x,
                                                p25 = params_ssilow[1],
                                                a = params_ssilow[2],
                                                b = params_ssilow[3],
                                                c = params_ssilow[4]
                                                )
                                  ),
  rate_summation_weibull = map_dbl(.x = temperature,
                                   .f = ~ weibull_1995(temp = .x, 
                                                       a = params_weibull[1],
                                                       topt = params_weibull[2],
                                                       b = params_weibull[3],
                                                       c = params_weibull[4]
                                   )
  ),
  rate_summation_gaussian = map_dbl(.x = temperature,
                                    .f = ~ gaussian_1987(temp = .x, 
                                                         rmax = params_gaussian[1],
                                                         topt = params_gaussian[2],
                                                         a = params_gaussian[3]
                                    )
  ),
  rate_summation_lactin1 = map_dbl(.x = temperature,
                               .f = ~ lactin1(temp = .x,
                                              a = params_lactin1[1],
                                              tmax = params_lactin1[2],
                                              delta_t = params_lactin1[3]
                                              )
                               ),
  rate_summation_rezende = map_dbl(.x = temperature,
                                   .f = ~ rezende(temp = .x,
                                                  q_10 = params_rezende[1],
                                                  cte = params_rezende[2],
                                                  thr = params_rezende[3],
                                                  decay = params_rezende[4]
                                                  )
                                   )
  )  |>  
  mutate(across(c(6:13), 
                ~if_else(condition = .x <0, true = 0, false = .x)),
         across(c(6:13), 
                ~ 100*.x))

rate_summation_emergence <- daily_ratesum_pupa %>% 
  group_by(year) %>% 
  mutate(across(c(5:12),
                ~cumsum(.x)),
         across(c(5:12), 
                ~logic_ratesum(.x)),
         doy = yday(date)
  ) %>% 
  group_by(year) %>% 
  summarise(across(c(5:12),
                   ~sum(.x))
  ) %>% 
  rename(briere1 = rate_summation_briere1,
         briere2 = rate_summation_briere2,
         lactin1 = rate_summation_lactin1,
         lactin2 = rate_summation_lactin2,
         ssilow = rate_summation_ssilow,
         rezende = rate_summation_rezende,
         weibull = rate_summation_weibull,
         gaussian = rate_summation_gaussian) %>% 
  pivot_longer(cols = -1,
               names_to = "model",
               values_to = "doe") 
rate_summation_emergence

# gg-plot it
doe_nonlinear_trends <- ggplot(rate_summation_emergence,
                               aes(x = year, 
                                   y = doe,
                                   color = model,
                                   fill = model))+
  geom_point()+
  geom_line(linetype = "dashed")+
  geom_smooth()+
  facet_wrap(~model)+
  theme_clean()+
  labs(x = "Year",
       y = "Day of first adult emergence")
doe_nonlinear_trends

##### b) hourly temps ----------------------------------------------------------
##  we use chillR package to simulate temperature variation across hours within days.
daily_temps_doy <- daily_tmax_df %>% 
  inner_join(daily_tmin_df) %>% 
  mutate(JDay = yday(date)) %>% 
  rename(Tmin = daily_tmin,
         Tmax = daily_tmax) 

hourly_temp_doy <- chillR::make_hourly_temps(latitude = 40.5,
                                             year_file = daily_temps_doy,
                                             keep_sunrise_sunset = TRUE ) %>%
  pivot_longer(cols = 8:31,
               names_to = "hour_of_day",
               values_to = "temperature") %>% 
  rename(tmax = Tmax,
         tmin = Tmin,
         doy = JDay,
         sunrise = Sunrise,
         sunset = Sunset,
         daylength = Daylength) 

hourly_ratesum_pupa <- hourly_temp_doy %>% 
  mutate(year = year(date),
        rate_summation_briere1 = map_dbl(.x = temperature,
                                         .f = ~ briere1(temp = .x,
                                                        tmin = params_briere1[1],
                                                        tmax = params_briere1[2],
                                                        a = params_briere1[3])
                                                 ),
         rate_summation_briere2 = map_dbl(.x = temperature,
                                          .f = ~ briere2(temp = .x,
                                                         tmin = params_briere2[1],
                                                         tmax = params_briere2[2],
                                                         a = params_briere2[3],
                                                         b = params_briere2[4])
         ),
         rate_summation_lactin2 = map_dbl(.x = temperature,
                                          .f = ~ lactin2(temp = .x,
                                                         a = params_lactin2[1],
                                                         tmax = params_lactin2[2],
                                                         delta_t = params_lactin2[3],
                                                         b = params_lactin2[4])
         ),
         rate_summation_ssilow = map_dbl(.x = temperature,
                                         .f = ~ ssilow(temp = .x,
                                                       p25 = params_ssilow[1],
                                                       a = params_ssilow[2],
                                                       b = params_ssilow[3],
                                                       c = params_ssilow[4]
                                         )
         ),
         rate_summation_weibull = map_dbl(.x = temperature,
                                          .f = ~ weibull_1995(temp = .x, 
                                                              a = params_weibull[1],
                                                              topt = params_weibull[2],
                                                              b = params_weibull[3],
                                                              c = params_weibull[4]
                                          )
         ),
         rate_summation_gaussian = map_dbl(.x = temperature,
                                           .f = ~ gaussian_1987(temp = .x, 
                                                                rmax = params_gaussian[1],
                                                                topt = params_gaussian[2],
                                                                a = params_gaussian[3]
                                           )
         ),
         rate_summation_lactin1 = map_dbl(.x = temperature,
                                          .f = ~ lactin1(temp = .x,
                                                         a = params_lactin1[1],
                                                         tmax = params_lactin1[2],
                                                         delta_t = params_lactin1[3]
                                          )
         ),
         rate_summation_rezende = map_dbl(.x = temperature,
                                          .f = ~ rezende(temp = .x,
                                                         q_10 = params_rezende[1],
                                                         cte = params_rezende[2],
                                                         thr = params_rezende[3],
                                                         decay = params_rezende[4]
                                                         )
                                          )
         ) %>% 
  mutate(across(c(11:18), 
                ~if_else(condition = .x <0, true = 0, false = .x)),
         across(c(11:18), 
                ~ .x*100/24))

rate_summation_emergence_hourly <- hourly_ratesum_pupa %>% 
  group_by(year) %>% 
  mutate(across(c(10:17),
                ~cumsum(.x)
  ),
  across(c(10:17), 
         ~logic_ratesum(.x)
  )
  ) %>% 
  summarise(across(c(10:17),
                   ~round(sum(.x, na.rm = TRUE)/24)
  )
  ) %>% 
  rename(briere1 = rate_summation_briere1,
         briere2 = rate_summation_briere2,
         lactin1 = rate_summation_lactin1,
         lactin2 = rate_summation_lactin2,
         ssilow = rate_summation_ssilow,
         rezende = rate_summation_rezende,
         weibull = rate_summation_weibull,
         gaussian = rate_summation_gaussian) %>% 
  pivot_longer(cols = -1,
               names_to = "model",
               values_to = "doe") 
rate_summation_emergence_hourly

# gg-plot it
doe_nonlinear_hourly_trends <- ggplot(rate_summation_emergence_hourly,
                                      aes(x = year,
                                          y = doe,
                                          color = model,
                                          fill = model))+
  geom_point()+
  geom_line(linetype = "dashed")+
  geom_smooth()+
  facet_wrap(~model)+
  theme_clean()+
  labs(x = "Year",
       y = "Day of first adult emergence")
doe_nonlinear_hourly_trends

# 4. Comparisons (long-term trends) ---------------------------------------
##### a) Visual comparison ----
doe_linear_daily <- doe_linear_pupa %>% 
  filter(method == "dailyavg") %>% 
  mutate(method = "linear") %>% 
  rename(model = method,
         daily = doe)

doe_linear_hourly <- doe_linear_pupa %>% 
  filter(method == "sinewave") %>% 
  mutate(method = "linear") %>% 
  rename(model = method,
         hourly = doe)
doe_linear_scale_compare <- inner_join(doe_linear_daily,
                                       doe_linear_hourly)

rate_summation_scale_compare <- rate_summation_emergence %>% 
  inner_join(rate_summation_emergence_hourly, 
             by = c("year", "model")) %>%
  rename(daily = doe.x,
         hourly = doe.y) %>% 
  bind_rows(doe_linear_scale_compare) %>% 
  pivot_longer(cols = 3:4,
               names_to = "time_res",
               values_to = "doe")
# incorporate observational data from Gordo & Sanz (2006)
pieris_obs_trends <- read_csv(here("Data/pieris_devdata/gordo_sanz_2006_doe.csv")) 
rate_summation_scale_obs <- rate_summation_scale_compare %>% 
  bind_rows(pieris_obs_trends) %>% 
  mutate()
scale_compare_trends_obs <- ggplot()+
  geom_point(data = rate_summation_scale_compare,
             aes(x = year, y = doe, color = time_res))+
  geom_line(data = rate_summation_scale_compare,
            aes(x = year, y = doe, color = time_res),
            linetype = "dashed")+
  geom_smooth(data = rate_summation_scale_compare,
              aes(x = year, y = doe, 
                  color = time_res, fill = time_res))+
  geom_point(data = pieris_obs_trends, 
             aes(x = year, y = doe),
             color = "#e9c46a")+
  geom_line(data = pieris_obs_trends, 
            aes(x = year, y = doe),
            color = "#e9c46a",
            linetype = "dashed")+
  geom_smooth(data = pieris_obs_trends, 
              aes(x = year, y = doe),
              color = "#e9c46a", fill = "#e9c46a")+
  facet_wrap(~model)+
  theme_clean()+
  labs(x = "Year",
       y = "Day of first adult emergence",
       color = "Resolution (model)",
       fill = "Resolution (model)")
scale_compare_trends_obs
ggsave(filename = here("Data/scale_compare_trends.png"),
       width = 20,
       height = 20,
       units = "cm")

##### b) RMSE ----
preds_vs_obs_trends <- rate_summation_scale_compare %>% 
  rename(doe_pred = doe) %>% 
  full_join(pieris_obs_trends) %>% #combine with observations from Gordo & Sanz (2006)
  drop_na() %>% # remove years of predictions with no observational data
  rename(doe_obs = doe) %>% 
  group_by(model, time_res) %>% 
  summarise(rmsep = chillR::RMSEP(predicted = doe_pred,
                                  observed = doe_obs)) %>% 
  arrange(rmsep) %>% 
  print() #gaussian, oneill and weibull are better models, followed by ssilow, hourly resolution does not improve accuracy
#for worse models, such as linear, briere, ssi and lactin, hourly resolution clearly improve accuracy.

transfer_rank_resolution <- preds_vs_obs_trends %>% 
  group_by(time_res) %>% 
  summarise(rmsep = mean(rmsep)) %>% 
  arrange(rmsep) %>% 
  print()# overall hourly rmsep is lower

transfer_rank_model <- preds_vs_obs_trends %>% 
  group_by(model) %>% 
  summarise(rmsep = mean(rmsep)) %>% 
  arrange(rmsep) %>% 
  print() # overall gaussian, oneill and weibull are better

#re-plot the graph (manually) annotating rmsep

preds_obs_compare_rmsep <- rate_summation_scale_compare %>% 
  full_join(preds_vs_obs_trends) %>% 
  group_by(year, model) %>% 
  summarise(rmsep = min(rmsep)) %>% 
  arrange(rmsep) %>%
  print()
rmsep_text <-  preds_obs_compare_rmsep %>% 
  group_by(model) %>% 
  summarise(rmsep = mean(rmsep)) %>%
  arrange(rmsep)
rmsep_order <- rmsep_text %>%
  select(model) %>% 
  as_vector()
rmsep_values <- rmsep_text %>% 
  mutate(rmsep =   paste("RMSEP =",
                         round(rmsep,2)),
         year = 1995,
         doe = 148)
preds_obs_compare_plot <- ggplot()+
  geom_point(data = rate_summation_scale_compare %>% filter(year %in% pieris_obs_trends$year),
             aes(x = year, y = doe, color = time_res),
             alpha = .66)+
  geom_line(data = rate_summation_scale_compare %>% filter(year %in% pieris_obs_trends$year),
            aes(x = year, y = doe, color = time_res),
            linetype = "longdash",
            alpha = .66)+
  geom_smooth(data = rate_summation_scale_compare %>% filter(year %in% pieris_obs_trends$year),
              aes(x = year, y = doe, 
                  color = time_res, fill = time_res))+
  geom_point(data = pieris_obs_trends, 
             aes(x = year, y = doe),
             color = "#e9c46a",
             fill = "#e9c46a",
             alpha = .66)+
  geom_line(data = pieris_obs_trends, 
            aes(x = year, y = doe),
            color = "#e9c46a",
            linetype = "longdash",
            alpha = .66)+
  geom_smooth(data = pieris_obs_trends, 
              aes(x = year, y = doe),
              color = "#e9c46a",
              fill = "#e9c46a",)+
  labs(x = "Year",
       y = "Day of first adult emergence",
       color = "Resolution (model)",
       fill = "Resolution (model)")+
  facet_wrap(~factor(model, levels = rmsep_order))+
  theme_few()+
  theme(strip.text.x = element_text(size = 10, 
                                    color = "grey42", 
                                    face = "bold")
  )+
  geom_text(data = rmsep_values, 
            aes(label = rmsep,
                x = year,
                y = doe), 
            size = 3)
preds_obs_compare_plot  

## save files into .svg and .png
ggsave(filename = here("Data/preds_obs_compare_plot.png"),
       width = 25,
       height = 25,
       units = "cm")
ggsave(filename = here("Data/preds_obs_compare_plot.svg"),
       width = 25,
       height = 25,
       units = "cm")


