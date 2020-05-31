# MSc Project 2020: Comparison of Spatial Gaussian Copula and Random Forests in Spatial Prediction
# Sun, N
#
# Parallelized simulations of gaussian copula and rfsp for spatial prediction
# Nonparallelized simulations are found in "simulated_data_comparison.R"

setwd("~/MS_Project/simulation_study/")
simdata <- readRDS("../simulated_Y.rds")
dat <- read.csv("../ForestDataFuzzed.csv")

source("foreach_func.R")

library(gstat)
library(tidyverse)
library(rgdal)
library(sp)
library(automap)
library(doParallel)

alb.xy <- project(with(dat,cbind(lon_fuzzed,lat_fuzzed)),
                  "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
colnames(alb.xy) <- c("x","y")
dat <- cbind(alb.xy,dat)

cl <- makeCluster(2)
registerDoParallel(cl)
sim.sets.n <- 100

start <- Sys.time()
foreach(k = 1:2) %dopar% {
    run.sim(index = k, dat = dat, train.n = 300)
}
end <- Sys.time()
end - start
