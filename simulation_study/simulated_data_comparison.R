setwd("~/MS_Project/simulation_study/")
# simdata <- readRDS("../simulated_Y.rds")
simdata <- readRDS("../simulated_tshevol.rds")
dat <- read.csv("../ForestDataFuzzed.csv")

source("../ZIGFunctions.R")
source("../metrics.R")
source("preprocess.R")
source("copula_func.R")
source("kriging_func.R")
source("rfsp_func.R")

library(gstat)
library(tidyverse)
library(rgdal)
library(sp)
library(automap)

alb.xy <- project(with(dat,cbind(lon_fuzzed,lat_fuzzed)),
                  "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
colnames(alb.xy) <- c("x","y")
dat <- cbind(alb.xy,dat)

train.n <- 300
sim.sets.n <- 1000

obs_responses <- vector(mode = "list",length = sim.sets.n)
# gauscop_predictions_alg <- vector(mode = "list",
#                               length =sim.sets.n)
gauscop_predictions_krige <- vector(mode = "list", length =sim.sets.n)
rfsp_predictions150 <- vector(mode = "list", length = sim.sets.n)
rfsp_predictions100 <- vector(mode = "list", length = sim.sets.n)
rfsp_predictions50 <- vector(mode = "list", length = sim.sets.n)
ok_predictions <- vector(mode = "list", length = sim.sets.n)
rfsp_predictions_tr0 <- vector(mode = "list", length = sim.sets.n)
ok_predictions_tr0 <- vector(mode = "list", length = sim.sets.n)

gauscop_times <- vector(mode = "numeric", length = sim.sets.n)
rfsp_times150 <- vector(mode = "numeric", length = sim.sets.n)
rfsp_times100 <- vector(mode = "numeric", length = sim.sets.n)
rfsp_times50 <- vector(mode = "numeric", length = sim.sets.n)
ok_times <- vector(mode = "numeric", length = sim.sets.n)

pic90s_gauscop <- vector(mode = "numeric", length = sim.sets.n)
pic90s_rfsp <- vector(mode = "numeric", length = sim.sets.n)
pic90s_ok <- vector(mode = "numeric", length = sim.sets.n)

start <- Sys.time()

for (i in 1:sim.sets.n) {
    
    set.seed(182)
    full_simdata <- data.frame(
        x = dat$x,
        y = dat$y,
        resp = simdata[i,]
    )
   
    tmp <- split_data(full_simdata, n = train.n)
    if (sum(tmp$train$resp) == 0) {
        obs_responses[[i]] <- tmp$test$resp
        # gauscop_predictions_alg[[i]] <- yhat_gauscop_alg
        gauscop_predictions_krige[[i]] <- tmp$test$resp
        rfsp_predictions150[[i]] <- tmp$test$resp
        rfsp_predictions100[[i]] <- tmp$test$resp
        rfsp_predictions50[[i]] <- tmp$test$resp
        ok_predictions[[i]] <- tmp$test$resp
        rfsp_predictions_tr0[[i]] <- tmp$test$resp
        ok_predictions_tr0[[i]] <- tmp$test$resp
        next
    }
    
    # Create spatial random forest models
    rfsp_start150 <- Sys.time()
    rfsp_model150 <- calculate_rfsp_predictions_rasterpd(training_data = tmp$train, test_points = tmp$test)
    rfsp_time150 <- Sys.time() - rfsp_start150
    
    rfsp_start100 <- Sys.time()
    rfsp_model100 <- calculate_rfsp_predictions_rasterpd(training_data = tmp$train, test_points = tmp$test, num.trees = 100)
    rfsp_time100 <- Sys.time() - rfsp_start100
    
    rfsp_start50 <- Sys.time()
    rfsp_model50 <- calculate_rfsp_predictions_rasterpd(training_data = tmp$train, test_points = tmp$test, num.trees = 50)
    rfsp_time50 <- Sys.time() - rfsp_start50
    
    # Create spatial Gaussian copula model
    gauscop_start <- Sys.time()
    transformed.training.set <- cube_root_fix_zeros(tmp$train)
    
    theta <- calculate_spatial_params(simdata = transformed.training.set$transformed_data,
                                      testdata = tmp$test,
                                      zeros = transformed.training.set$zeros)
    krige_zhat_gauscop(resp = transformed.training.set$transformed_data$resp,
                       mu = theta$mu,
                       beta = theta$beta,
                       epsilon = transformed.training.set$epsilon,
                       Pi = transformed.training.set$Pi,
                       training_data = tmp$train,
                       test_data = tmp$test) -> zhat.zigs.df.krige
    gauscop_time <- Sys.time() - gauscop_start
    
    # calculate_zhat_gauscop(resp = transformed.training.set$transformed_data$resp,
    #                        mu = theta$mu,
    #                        beta = theta$beta,
    #                        epsilon = transformed.training.set$epsilon,
    #                        Pi = transformed.training.set$Pi,
    #                        Sigma_obs = theta$Sigma_obs,
    #                        Sigma_obs_noobs = theta$Sigma_obs_noobs) -> zhat.zigs.df.alg
    
    yhat_gauscop_krige <- backtransform_gauscop(zhat.zigs.df.krige$zigs, transformed.training.set$epsilon)
    # yhat_gauscop_alg <- backtransform_gauscop(zhat.zigs.df.alg$zigs, transformed.training.set$epsilon)
    
    # Create ordinary kriging model
    ok_start <- Sys.time()
    training_data.cpy <- tmp$train
    test_data.cpy <- tmp$test
    coordinates(training_data.cpy) <- ~x+y
    coordinates(test_data.cpy) <- ~x+y
    
    ak_model <- autoKrige(resp ~ 1, training_data.cpy, test_data.cpy)
    ok_time <- Sys.time() - ok_start
    
    obs_responses[[i]] <- tmp$test$resp
    # gauscop_predictions_alg[[i]] <- yhat_gauscop_alg
    gauscop_predictions_krige[[i]] <- yhat_gauscop_krige
    rfsp_predictions150[[i]] <- rfsp_model150$predictions
    rfsp_predictions100[[i]] <- rfsp_model100$predictions
    rfsp_predictions50[[i]] <- rfsp_model50$predictions
    ok_predictions[[i]] <- ak_model$krige_output$var1.pred
    rfsp_predictions_tr0[[i]] <- small_preds_to_zeros(rfsp_model150$predictions, transformed.training.set$epsilon)
    ok_predictions_tr0[[i]] <- small_preds_to_zeros(ak_model$krige_output$var1.pred, transformed.training.set$epsilon)
    
    gauscop_times[i] <- gauscop_time
    rfsp_times150[i] <- rfsp_time150
    rfsp_times100[i] <- rfsp_time100
    rfsp_times50[i] <- rfsp_time50
    ok_times[i] <- ok_time
    
    pic90s_gauscop[i] <- pic90(yhat_gauscop_krige, tmp$test$resp)
    pic90s_rfsp[i] <- pic90(rfsp_model150$predictions, tmp$test$resp)
    pic90s_ok[i] <- pic90(ak_model$krige_output$var1.pred, tmp$test$resp)
    
    print(paste0("Completed set: ", toString(i)))
}

runtime <- Sys.time() - start

sim.final.df <- data.frame(
    obs = unlist(obs_responses),
    gauscop_krige = unlist(gauscop_predictions_krige),
    # gauscop_alg = unlist(gauscop_predictions_alg),
    rfsp150 = unlist(rfsp_predictions150),
    rfsp100 = unlist(rfsp_predictions100),
    rfsp50 = unlist(rfsp_predictions50),
    ok = unlist(ok_predictions),
    rfsp_tr0 = unlist(rfsp_predictions_tr0),
    ok_tr0 = unlist(ok_predictions_tr0)
)

sim.final.df$gauscop_krige[is.na(sim.final.df$gauscop_krige)] <- 0
sim.final.df$ok[is.na(sim.final.df$ok)] <- 0
sim.final.df$ok_tr0[is.na(sim.final.df$ok_tr0)] <- 0

# RMSPE of the three methods
rmspe(sim.final.df$gauscop_krige, sim.final.df$obs)
# rmspe(sim.final.df$gauscop_alg, sim.final.df$obs)
rmspe(sim.final.df$rfsp150, sim.final.df$obs)
rmspe(sim.final.df$rfsp100, sim.final.df$obs)
rmspe(sim.final.df$rfsp50, sim.final.df$obs)
rmspe(sim.final.df$ok, sim.final.df$obs)
rmspe(sim.final.df$rfsp_tr0, sim.final.df$obs)
rmspe(sim.final.df$ok_tr0, sim.final.df$obs)

# signed relative bias
srb(sim.final.df$gauscop_krige, sim.final.df$obs)
srb(sim.final.df$rfsp150, sim.final.df$obs)
srb(sim.final.df$rfsp100, sim.final.df$obs)
srb(sim.final.df$rfsp50, sim.final.df$obs)
srb(sim.final.df$ok, sim.final.df$obs)
srb(sim.final.df$rfsp_tr0, sim.final.df$obs)
srb(sim.final.df$ok_tr0, sim.final.df$obs)

# PIC90
mean(pic90s_gauscop[!is.na(pic90s_gauscop)])
mean(pic90s_rfsp[!is.na(pic90s_rfsp)])
mean(pic90s_ok[!is.na(pic90s_ok)])

# saveRDS(sim.final.df,
#         file = "sim_final_df.rds")

# summary(gauscop_times)
# summary(ok_times)
# summary(rfsp_times150)
# summary(rfsp_times100)
# summary(rfsp_times50)

# Histograms of predictions vs observed values
# par(mfrow=c(2,2))
# hist(sim.final.df$gauscop_krige,
#      main = "Gaussian Copula (kriging) predictions",
#      breaks = seq(0, 3000, 100),
#      xaxt = 'n')
# hist(sim.final.df$gauscop_alg[sim.final.df$gauscop_alg <= 3000],
#      main = "Gaussian Copula (matrix algebra) predictions",
#      breaks = seq(0, 3000, 100),
#      xaxt = 'n')
# hist(sim.final.df$ok,
#      main = "Ordinary Kriging predictions",
#      breaks = seq(min(sim.final.df$ok), 3000, 100),
#      xaxt = 'n')
# hist(sim.final.df$rfsp150,
#      main = "Random forest predictions",
#      breaks = seq(0, 3000, 100),
#      xaxt = 'n')
# hist(sim.final.df$obs,
#      main = "Observed values",
#      breaks = seq(0, max(sim.final.df$obs)+100, 100))
# par(mfrow=c(1,1))

# How are zeros showing up in the predictions
# sum(sim.final.df$gauscop == 0) / nrow(sim.final.df)
# sum(sim.final.df$rfsp_tr0 == 0) / nrow(sim.final.df)
# sum(sim.final.df$obs == 0) / nrow(sim.final.df)

# Let's do some residual plots
# par(mfrow=c(1,3))
# plot(x = sim.final.df$obs,
#      y = (sim.final.df$gauscop - sim.final.df$obs),
#      main = "Residual plot of Gaussian Copula",
#      xlab = "Observed Values",
#      ylab = "Residuals")
# plot(x = sim.final.df$obs,
#      y = (sim.final.df$ok - sim.final.df$obs),
#      main = "Residual plot of Ordinary Kriging",
#      xlab = "Observed Values",
#      ylab = "Residuals")
# plot(x = sim.final.df$obs,
#      y = (sim.final.df$rfsp150 - sim.final.df$obs),
#      main = "Residual plot of Random Forest",
#      xlab = "Observed Values",
#      ylab = "Residuals")
# par(mfrow=c(1,1))

resids <- data.frame(obs = sim.final.df$obs,
                     cop_res = (sim.final.df$gauscop_krige - sim.final.df$obs),
                     rfsp_res = (sim.final.df$rfsp150 - sim.final.df$obs),
                     ok_res = (sim.final.df$ok - sim.final.df$obs))

ggplot(resids) +
    geom_point(mapping = aes(x = obs, y= cop_res),
               color ="#013d09", alpha = .2) +
    geom_line(mapping = aes(x = obs, y = -obs),
              linetype = 2, color = "black") +
    labs(
        title = "Copula Residuals",
        x = element_blank(),
        y = "Residuals"
    ) +
    theme_minimal() -> g1

ggplot(resids) +
    geom_point(mapping = aes(x = obs, y= rfsp_res),
               color ="#013d09", alpha = .2) +
    geom_line(mapping = aes(x = obs, y = -obs),
              linetype = 2, color = "black") +
    labs(
        title = "RFsp Residuals",
        x = element_blank(),
        y = element_blank()
    ) +
    theme_minimal() -> g2

ggplot(resids) +
    geom_point(mapping = aes(x = obs, y= ok_res),
               color ="#013d09", alpha = .2) +
    geom_line(mapping = aes(x = obs, y = -obs),
              linetype = 2, color = "black") +
    labs(
        title = "Kriging Residuals",
        x = element_blank(),
        y = element_blank()
    ) +
    theme_minimal() -> g3

library(gridExtra)
grid.arrange(g1, g3, g2, nrow = 1,
