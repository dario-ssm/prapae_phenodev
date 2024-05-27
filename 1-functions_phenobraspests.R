### Script Info ####

# Authors: Darío San Segundo Molina, Ignacio Morales Castilla
# 
# GitHub repo: github.com/dario-ssm/PhenoBrassicaPests

# Aim: functions script
#
# Species: Pieris rapae
# Source study: Davies, 1985
# Climate dataset: Spain02 v5  (Herrera et al. 2016 and Kotlarsky et al. 2017; see http://www.meteo.unican.es/datasets/spain02 )
## Aknowledgements: The authors thank AEMET and UC by the data provided for this work (Spain02v5 gridded temperature data set). 


# 1. Sinewave Functions ----
# We define within each emergence function a sinewave function to compute Degree Days accumulated
# based on minimum and maximum temperatures. An imput of climatic
# data with daily Minimum Temperatures (Tmin) and Maximum temperatures (Tmax) is
# required
# define sinewave function
sinewave_lower <- function(Tmin,Tmax,LDT){
  a= (Tmax-Tmin)/2
  b= (Tmax+Tmin)/2
  theta=asin((LDT-b)/a)
  theta[theta =="NaN"] <- 0
  DDs <- (1/pi)*((b-LDT)*((pi/2)-theta)+a*cos(theta))
  DDs[DDs<0] <- 0
  return(DDs)
}
sinewave_upper <- function(Tmin,Tmax,LDT,UDT){
  a = (Tmax-Tmin)/2
  b = (Tmax+Tmin)/2
  theta2=asin((UDT-b)/a)
  theta2[theta2 =="NaN"] <- 0
  DDs=(1/pi)*((b-LDT)*(theta2+(pi/2))+(UDT-LDT)*((pi/2)-theta2)-(a*cos(theta2)))
  DDs[DDs<0] <- 0
  return(DDs)
}
sinewave_both <- function(Tmin,Tmax,LDT,UDT){
  a= (Tmax-Tmin)/2
  b= (Tmax+Tmin)/2
  theta=asin((LDT-b)/a)
  theta2=asin((UDT-b)/a)
  theta[theta=="NaN"] <- 0
  theta2[theta2 =="NaN"] <- 0
  DDs=(1/pi)*((b-LDT)*(theta2-theta)+a*(cos(theta)-cos(theta2))+(UDT-LDT)*((pi/2)-theta2))
  DDs[DDs<0] <- 0
  return(DDs)
}
sinewave_bet <- function(Tmin,Tmax,LDT){
  a= (Tmax-Tmin)/2
  b= (Tmax+Tmin)/2
  DDs=b-LDT
  DDs[DDs<0] <- 0
  return(DDs)
}
sinewave_ldt <- function(Tmin,Tmax,LDT){
  ifelse((Tmin<LDT & Tmax<=LDT),0,
         ifelse((Tmin<LDT & Tmax>LDT),sinewave_lower(Tmin,Tmax,LDT),
                ifelse((Tmin>=LDT),sinewave_bet(Tmin,Tmax,LDT),NA)))
  
}

sinewave <- function(Tmin,Tmax,LDT,UDT){
  a= (Tmax-Tmin)/2
  b= (Tmax+Tmin)/2
  theta=asin((LDT-b)/a)
  theta2=asin((UDT-b)/a)
  ifelse((Tmin<LDT & Tmax<=LDT),0,
         ifelse((Tmin<LDT & UDT>Tmax & Tmax>LDT),sinewave_lower(Tmin,Tmax,LDT),
                ifelse((UDT>Tmin & Tmin>LDT & Tmax>UDT), sinewave_upper(Tmin,Tmax,LDT,UDT),
                       ifelse((Tmin<LDT & Tmax > UDT),
                              sinewave_both(Tmin,Tmax,LDT,UDT),
                              ifelse((Tmin>=UDT & Tmax > UDT),UDT-LDT,
                                     ifelse((UDT>Tmin & Tmin>=LDT & LDT<Tmax & Tmax<=UDT),
                                            sinewave_bet(Tmin,Tmax,LDT),NA))))))
}

# 2. Triangle Method ----
# We define within each emergence function a triangle function to compute Degree Days accumulated
# based on minimum and maximum temperatures. An imput of climatic
# data with daily Minimum Temperatures (Tmin) and Maximum temperatures (Tmax) is
# required
# define triangle function (Lindsay 1956?)

triangle <- function(tmin, tmax, LDT){
  DDs = (12*(tmax-LDT)^2)/(tmax-tmin)
  DDs[DDs<0] <- 0
  return(DDs)
}

# 3. daily avg ----
# We define within each emergence function a daily avg function to compute Degree Days accumulated
# based on minimum and maximum temperatures. An imput of climatic
# data with daily Minimum Temperatures (Tmin) and Maximum temperatures (Tmax) is
# required
# define dailyavg function (Lindsay 1956?)

daily_avg <- function(tmax,tmin,LDT){
  DDs = mean(c(tmax,tmin)) - LDT
  DDs[DDs<0] <- 0
  return(DDs)
}


# 4. Development rate models ----------------------------------------------
library(rTPC)
# get parameters
param_names_extractor <- function(nls_object){
  parameter_est <- coef(nls_object)
  param_names <- names(parameter_est)
  return(param_names)
}

briere1 <- function (temp, tmin, tmax, a) {
  est <- a * temp * (temp - tmin) * (tmax - temp)^(1/2)
  return(est)
}

briere2 <- function (temp, tmin, tmax, a, b) {
  est <- a * temp * (temp - tmin) * (tmax - temp)^(1/b)
  return(est)
}

## rTPC::gaussian_1987()

lactin2 <- function (temp, a, tmax, delta_t, b) {
  est <- exp(a * temp) - exp(a * tmax - ((tmax - temp)/delta_t)) +  b
  return(est)
}

lactin1 <- function (temp, a, tmax, delta_t) {
  est <- exp(a * temp) - exp(a * tmax - (tmax - temp)/delta_t)
  return(est)
}

mod_polynomial <- function(temp, a_0, a_1, a_2, a_3, a_4){
  est <- a_0 + a_1 * temp + a_2 * temp^2 + a_3 * temp^3 + a_4 * temp^4
  return(est)
}

ratkowsky1983 <- function (temp, a, tmin, tmax, b) 
{
  est <- ((a * (temp - tmin)) * (1 - exp(b * (temp - tmax))))^2
  return(est)
}

regniere <- function(temp, phi, b, btmin, tmax, delta_b, delta_m) {
  est <- phi* (exp(b * (temp - tmin)) - ((tmax - temp)/(tmax - tmin)) * exp(-b *
                                                                              (temp - tmin)/delta_b) - ((temp - tmin)/(tmax - tmin)) * exp(b * (tmax - tmin) - (tmax - temp)/delta_m))
  return(est)
}

rezende <- function (temp, q_10, cte, thr, decay) {
  est <- {
    ifelse(temp < thr, 
           (cte * 10^(log10(q_10)/(10/temp))), 
           (cte * 10^(log10(q_10)/(10/temp))) * (1 - decay * (thr - temp)^2))
  }
  return(est)
}

shi2011 <- function (temp, c, k1, t1, k2, t2) {
  est <- c * (1 - exp(-k1 * (temp - t1))) * (1 - exp(k2 * (temp - t2)))
  return(est)
}

ssi1981 <- function (temp, r_tref, e, el, tl, eh, th, tref) {
  tref <- 273.15 + tref
  k <- 8.62e-05
  boltzmann.term <- r_tref * exp(e/k * (1/tref - 1/(temp +273.15)))
  inactivation.term <- 1/(1 + exp(-el/k * (1/(tl + 273.15) - 1/(temp + 273.15))) 
                          + exp(eh/k * (1/(th + 273.15) -   1/(temp + 273.15))))
  return(boltzmann.term * inactivation.term)
}

ssilow <- function(temp, p25, a, b, c) {
  est <- (p25 * (temp + 273.16)/298 * exp(a/1.987 * (1/298 - 1/(temp + 273.16))))/(1 + exp(b/1.987 * (1/c - 1/(temp + 273.16))))
  return(est)
}

stinner1974 <- function (temp, c, k1, k2, topt){
  est <- c/(1 + exp(k1 + k2 * (2 * topt - temp)))
  return(est)
}

wang <- function(temp, k, r, topt, tmin, tmax, a){
  est <- (k/(1 + exp(-r * (temp - topt)))) * (1 - exp(-(temp - tmin)/a)) *
    (1 - exp(-(tmax - temp)/a))
  return(est)
}

linear <- function (temp, intercept, slope) 
{
  est <- slope*temp + intercept
  return(est)
}


devRate::poly4$eq


# 5. doe_dds ----------------------------------------------
logic_dd <- function(heat_units, rate_cumsum){
  logic_dd = if_else(rate_cumsum < heat_units,
                     1,
                     0)
  return(logic_dd)
}

# 6. logic rate summation nonlinear ----------------------------------------------
logic_ratesum <- function(rate_cumsum){
  logic_rate_summation = if_else(rate_cumsum < 100,
                     1,
                     0)
  return(logic_rate_summation)
}


