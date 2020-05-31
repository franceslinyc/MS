setwd("~/MS_Project/resampling_study/")
dat <- read.csv("../ForestDataFuzzed.csv")

source("../ZIGFunctions.R")
source("../metrics.R")
source("preprocess.R")
source("copula_func.R")
source("kriging_func.R")
source("rfsp_func.R")

suppressMessages({
    library(gstat)
    library(tidyverse)
    library(rgdal)
    library(sp)
    library(automap)
})

alb.xy <- project(with(dat,cbind(lon_fuzzed,lat_fuzzed)),
                  "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
colnames(alb.xy) <- c("x","y")
dat <- cbind(alb.xy,dat)

train.n <- 500
sim.sets.n <- 1000

obs_responses <- vector(mode = "list", length = sim.sets.n)
gauscop_predictions_krige <- vector(mode = "list", length =sim.sets.n)
rfsp_predictions150 <- vector(mode = "list", length = sim.sets.n)
rfsp_predictions100 <- vector(mode = "list", length = sim.sets.n)
rfsp_predictions50 <- vector(mode = "list", length = sim.sets.n)
uk_predictions <- vector(mode = "list", length = sim.sets.n)
rfsp_predictions_tr0 <- vector(mode = "list", length = sim.sets.n)
uk_predictions_tr0 <- vector(mode = "list", length = sim.sets.n)

test.x <- vector(mode = "list", length = sim.sets.n)
test.y <- vector(mode = "list", length = sim.sets.n)

gauscop_times <- vector(mode = "numeric", length = sim.sets.n)
rfsp_times150 <- vector(mode = "numeric", length = sim.sets.n)
rfsp_times100 <- vector(mode = "numeric", length = sim.sets.n)
rfsp_times50 <- vector(mode = "numeric", length = sim.sets.n)
uk_times <- vector(mode = "numeric", length = sim.sets.n)

pic90s_gauscop <- vector(mode = "numeric", length = sim.sets.n)
pic90s_rfsp <- vector(mode = "numeric", length = sim.sets.n)
pic90s_uk <- vector(mode = "numeric", length = sim.sets.n)

start <- Sys.time()

full_simdata <- data.frame(
    x = dat$x,
    y = dat$y,
    annpre = dat$annpre,
    resp = dat$tshevol
)

set.seed(182)

for (i in 1:sim.sets.n) {
    
    tmp <- split_data(full_simdata, n = train.n)
    
    # Create spatial random forest models
    rfsp_start150 <- Sys.time()
    rfsp_model150 <- calculate_rfsp_predictions_rasterpd(training_data = tmp$train, test_data = tmp$test)
    rfsp_time150 <- Sys.time() - rfsp_start150
    
    rfsp_start100 <- Sys.time()
    rfsp_model100 <- calculate_rfsp_predictions_rasterpd(training_data = tmp$train, test_data = tmp$test, num.trees = 100)
    rfsp_time100 <- Sys.time() - rfsp_start100
    
    rfsp_start50 <- Sys.time()
    rfsp_model50 <- calculate_rfsp_predictions_rasterpd(training_data = tmp$train, test_data = tmp$test, num.trees = 50)
    rfsp_time50 <- Sys.time() - rfsp_start50
    
    # Create spatial Gaussian copula model
    gauscop_start <- Sys.time()
    theta <- calculate_glm_params(data = tmp$train)
    krige_zhat_gauscop(training_data = tmp$train,
                       test_data = tmp$test,
                       resp = theta$transformed_data$resp,
                       mu = theta$mu,
                       beta = theta$beta,
                       Pi = theta$Pi,
                       epsilon = theta$epsilon,
                       nonzero.glm.coe = theta$nonzero.glm.coe,
                       zero.glm.coe = theta$zero.glm.coe,
                       s2 = theta$s2) -> zhat.zigs.df.krige
    yhat_gauscop_krige <- backtransform_gauscop(zhat.zigs.df.krige$zigs, theta$epsilon)
    gauscop_time <- Sys.time() - gauscop_start
    
    # Create ordinary kriging model
    uk_start <- Sys.time()
    training_data.cpy <- tmp$train
    test_data.cpy <- tmp$test
    coordinates(training_data.cpy) <- ~x+y
    coordinates(test_data.cpy) <- ~x+y
    
    uk_preds <- tryCatch({
        suppressWarnings({
            uk_preds <- autoKrige(resp ~ annpre, training_data.cpy, test_data.cpy)$krige_output$var1.pred
        })
    }, error = function(e){
        print(e)
        uk_preds <- rep(0, nrow(tmp$test))
        return(uk_preds)
    })
    uk_time <- Sys.time() - uk_start
   
    # Add the predictions to our list objects 
    obs_responses[[i]] <- tmp$test$resp
    gauscop_predictions_krige[[i]] <- yhat_gauscop_krige
    rfsp_predictions150[[i]] <- rfsp_model150$predictions
    rfsp_predictions100[[i]] <- rfsp_model100$predictions
    rfsp_predictions50[[i]] <- rfsp_model50$predictions
    uk_predictions[[i]] <- uk_preds
    rfsp_predictions_tr0[[i]] <- small_preds_to_zeros(rfsp_model150$predictions, theta$epsilon)
    uk_predictions_tr0[[i]] <- small_preds_to_zeros(uk_preds, theta$epsilon)
    
    test.x[[i]] <- tmp$test$x
    test.y[[i]] <- tmp$test$y
   
    # Add the runtimes to the appropriate vectors 
    gauscop_times[i] <- gauscop_time
    rfsp_times150[i] <- rfsp_time150
    rfsp_times100[i] <- rfsp_time100
    rfsp_times50[i] <- rfsp_time50
    uk_times[i] <- uk_time
    
    # PIC90 of the methods
    pic90s_gauscop[i] <- pic90(yhat_gauscop_krige, tmp$test$resp)
    pic90s_rfsp[i] <- pic90(rfsp_model150$predictions, tmp$test$resp)
    pic90s_uk[i] <- pic90(uk_preds, tmp$test$resp)
    
    print(paste0("Completed set: ", toString(i)))
}

runtime <- Sys.time() - start
runtime

sim.final.df <- data.frame(
    x = unlist(test.x),
    y = unlist(test.y),
    obs = unlist(obs_responses),
    gauscop_krige = unlist(gauscop_predictions_krige),
    rfsp150 = unlist(rfsp_predictions150),
    rfsp100 = unlist(rfsp_predictions100),
    rfsp50 = unlist(rfsp_predictions50),
    uk = unlist(uk_predictions),
    rfsp_tr0 = unlist(rfsp_predictions_tr0),
    uk_tr0 = unlist(uk_predictions_tr0)
)

# saveRDS(sim.final.df,
#         file = "sim_final_df.rds")

# list(
#     obs = obs_responses,
#     gauscop = gauscop_predictions,
#     rfsp = rfsp_predictions,
#     ok = ok_predictions
# ) -> all.output.data
# 
# saveRDS(all.output.data,
#         file = "sim_final_lists.rds")

summary(gauscop_times)
summary(uk_times)
summary(rfsp_times150)
summary(rfsp_times100)
summary(rfsp_times50)

# RMSPE of the three methods
rmspe(sim.final.df$gauscop_krige, sim.final.df$obs)
rmspe(sim.final.df$rfsp150, sim.final.df$obs)
rmspe(sim.final.df$rfsp100, sim.final.df$obs)
rmspe(sim.final.df$rfsp50, sim.final.df$obs)
rmspe(sim.final.df$uk, sim.final.df$obs)
rmspe(sim.final.df$rfsp_tr0, sim.final.df$obs)
rmspe(sim.final.df$uk_tr0, sim.final.df$obs)

# Signed rank bias of the methods
srb(sim.final.df$gauscop_krige, sim.final.df$obs)
srb(sim.final.df$rfsp150, sim.final.df$obs)
srb(sim.final.df$rfsp100, sim.final.df$obs)
srb(sim.final.df$rfsp50, sim.final.df$obs)
srb(sim.final.df$uk, sim.final.df$obs)
srb(sim.final.df$rfsp_tr0, sim.final.df$obs)
srb(sim.final.df$uk_tr0, sim.final.df$obs)

# Prediction inclusion interval metrics
mean(pic90s_gauscop)
mean(pic90s_rfsp)
mean(pic90s_uk)

# colnames(sim.final.df)
# sim.final.df %>% 
#     group_by(obs) %>% 
#     summarise(
#         gauscop.rmspe = rmspe(gauscop_krige, obs),
#         rfsp.rmspe = rmspe(rfsp150, obs),
#         uk.rmspe = rmspe(uk, obs),
#     ) %>% 
#     ggplot() +
#     geom_line(mapping = aes(x = obs, y = gauscop.rmspe),
#               color = "#013d09") +
#     theme_minimal()

# Histograms of predictions vs observed values
# par(mfrow=c(2,2))
# hist(sim.final.df$gauscop_krige,
#      main = "Gaussian Copula (kriging) predictions",
#      breaks = seq(0, 3000, 100))
# hist(sim.final.df$uk,
#      main = "Universal Kriging predictions",
#      breaks = seq(min(sim.final.df$uk), 3000, 100))
# hist(sim.final.df$rfsp150,
#      main = "Random forest predictions",
#      breaks = seq(0, 3000, 100))
# hist(sim.final.df$obs,
#      main = "Observed values",
#      breaks = seq(0, max(sim.final.df$obs)+100, 100))
# par(mfrow=c(1,1))

# How are zeros showing up in the predictions
sum(sim.final.df$gauscop == 0) / nrow(sim.final.df)
sum(sim.final.df$rfsp_tr0 == 0) / nrow(sim.final.df)
sum(sim.final.df$obs == 0) / nrow(sim.final.df)

colnames(sim.final.df)

sim.final.df %>% 
    group_by(obs) %>% 
    summarise(
        count = n(),
        mean.gauscop = median(gauscop_krige),
        mean.rfsp = median(rfsp150),
        mean.uk = median(uk),
        mean.rfsp_tr0 = median(rfsp_tr0),
        mean.uk_tr0 = median(uk_tr0),
        rmspe.gauscop = rmspe(gauscop_krige, obs),
        rmspe.rfsp = rmspe(rfsp150, obs),
        rmspe.uk = rmspe(uk, obs),
        rmspe.rfsp_tr0 = rmspe(rfsp_tr0, obs),
        rmspe.uk_tr0 = rmspe(uk_tr0, obs)
    ) %>% 
    filter(obs == 0)

sim.final.df %>% 
    group_by(obs) %>% 
    summarise(
        count = n(),
        gauscop.pzero = sum(gauscop_krige == 0)/count,
        rfsp.pzero = sum(rfsp150 == 0)/count,
        uk.pzero = sum(uk == 0)/count
    )

sim.final.df %>% 
    filter(obs == 0) %>% 
    dplyr::select(rfsp_tr0) %>% 
    ggplot() +
    geom_histogram(mapping = aes(x = rfsp_tr0),
                   fill = "#013d09") +
    labs(
        title = "Random forest(150)",
        x = element_blank()
    ) +
    theme_minimal() -> g1

sim.final.df %>% 
    filter(obs == 0) %>% 
    dplyr::select(gauscop_krige) %>% 
    ggplot() +
    geom_histogram(mapping = aes(x = gauscop_krige),
                   fill = "#013d09") +
    labs(
        title = "Gaussian copula",
        x = element_blank()
    ) +
    theme_minimal() -> g2

sim.final.df %>% 
    filter(obs == 0) %>% 
    dplyr::select(uk) %>% 
    ggplot() +
    geom_histogram(mapping = aes(x = uk),
                   fill = "#013d09") +
    labs(
        title = "Universal kriging",
        x = element_blank()
    ) +
    theme_minimal() -> g3

grid.arrange(g2, g1, g3, nrow = 1)

sim.final.df %>%
    group_by(x, y) %>%
    summarise(
        mse.gauscop = sd(gauscop_krige),
        mse.rfsp = sd(rfsp150),
        mse.uk = sd(uk),
    ) -> plot.mse

ggplot(dat, aes(x=x/1000,y=y/1000)) +
    geom_point(mapping = aes(color = totvol),
               size=2) +
    geom_point(stroke=1,
               shape = 21,
               color = "black"
               ) +
    scale_color_gradient(low = "white", high = "#013d09") +
    scale_x_continuous(breaks = round(seq(-2150, -2000, by = 50),2)) +
    labs(
        title = "Total Timber Volume in Sampled Plots",
        x = "X (kilometers)",
        y = "Y (kilometers)",
        color = "Total Volume"
    ) +
    theme(
        axis.text.x=element_text(angle=90,hjust=1),
        panel.background = element_rect(fill = 'white'),
        panel.grid = element_line(color = "grey", size = .1),
        plot.caption = element_text(hjust = 0, face = "italic")
    ) -> g1

ggplot(plot.mse, aes(x=x/1000,y=y/1000)) +
    geom_point(mapping = aes(color = mse.gauscop),
               size=2) +
    geom_point(stroke=1,
               shape = 21,
               color = "black"
               ) +
    scale_color_gradient(low = "white", high = "#013d09") +
    scale_x_continuous(breaks = round(seq(-2150, -2000, by = 50),2)) +
    labs(
        title = "RMSPE of Gaussian Copula Model",
        x = "X (kilometers)",
        y = "Y (kilometers)",
        color = "RMSPE"
    ) +
    theme(
        axis.text.x=element_text(angle=90,hjust=1),
        panel.background = element_rect(fill = 'white'),
        panel.grid = element_line(color = "grey", size = .1),
        plot.caption = element_text(hjust = 0, face = "italic")
    ) -> g2

ggplot(plot.mse, aes(x=x/1000,y=y/1000)) +
    geom_point(mapping = aes(color = mse.rfsp),
               size=2) +
    geom_point(stroke=1,
               shape = 21,
               color = "black"
               ) +
    scale_color_gradient(low = "white", high = "#013d09") +
    scale_x_continuous(breaks = round(seq(-2150, -2000, by = 50),2)) +
    labs(
        title = "RMSPE of Random Forest(150)",
        x = "X (kilometers)",
        y = "Y (kilometers)",
        color = "RMSPE"
    ) +
    theme(
        axis.text.x=element_text(angle=90,hjust=1),
        panel.background = element_rect(fill = 'white'),
        panel.grid = element_line(color = "grey", size = .1),
        plot.caption = element_text(hjust = 0, face = "italic")
    ) -> g3

ggplot(plot.mse, aes(x=x/1000,y=y/1000)) +
    geom_point(mapping = aes(color = mse.uk),
               size=2) +
    geom_point(stroke=1,
               shape = 21,
               color = "black"
               ) +
    scale_color_gradient(low = "white", high = "#013d09") +
    scale_x_continuous(breaks = round(seq(-2150, -2000, by = 50),2)) +
    labs(
        title = "RMSPE of Universal Kriging",
        x = "X (kilometers)",
        y = "Y (kilometers)",
        color = "RMSPE"
    ) +
    theme(
        axis.text.x=element_text(angle=90,hjust=1),
        panel.background = element_rect(fill = 'white'),
        panel.grid = element_line(color = "grey", size = .1),
        plot.caption = element_text(hjust = 0, face = "italic")
    ) -> g4

require(gridExtra)
grid.arrange(g1, g2, g3, g4, nrow =2)

# Let's do some residual plots
# par(mfrow=c(1,3))
# plot(x = sim.final.df$obs,
#      y = (sim.final.df$gauscop - sim.final.df$obs),
#      ylim = c(min((sim.final.df$gauscop - sim.final.df$obs)), max((sim.final.df$gauscop - sim.final.df$obs))),
#      main = "Residual plot of Gaussian Copula",
#      xlab = "Observed Values",
#      ylab = "Residuals")
# plot(x = sim.final.df$obs,
#      y = (sim.final.df$uk - sim.final.df$obs),
#      ylim = c(min((sim.final.df$gauscop - sim.final.df$obs)), max((sim.final.df$gauscop - sim.final.df$obs))),
#      main = "Residual plot of Universal Kriging",
#      xlab = "Observed Values",
#      ylab = "Residuals")
# plot(x = sim.final.df$obs,
#      y = (sim.final.df$rfsp150 - sim.final.df$obs),
#      ylim = c(min((sim.final.df$gauscop - sim.final.df$obs)), max((sim.final.df$gauscop - sim.final.df$obs))),
#      main = "Residual plot of Random Forest",
#      xlab = "Observed Values",
#      ylab = "Residuals")
# par(mfrow=c(1,1))

# resids <- data.frame(obs = sim.final.df$obs,
#                      cop_res = (sim.final.df$gauscop_krige - sim.final.df$obs),
#                      rfsp_res = (sim.final.df$rfsp150 - sim.final.df$obs),
#                      uk_res = (sim.final.df$uk - sim.final.df$obs))
# 
# ggplot(resids) +
#     geom_point(mapping = aes(x = obs, y= cop_res),
#                color ="#013d09", alpha = .2) +
#     geom_line(mapping = aes(x = obs, y = -obs),
#               linetype = 2, color = "black") +
#     labs(
#         title = "Copula Residuals",
#         x = element_blank(),
#         y = "Residuals"
#     ) +
#     theme_minimal() -> g1
# 
# ggplot(resids) +
#     geom_point(mapping = aes(x = obs, y= rfsp_res),
#                color ="#013d09", alpha = .2) +
#     geom_line(mapping = aes(x = obs, y = -obs),
#               linetype = 2, color = "black") +
#     labs(
#         title = "RFsp Residuals",
#         x = element_blank(),
#         y = element_blank()
#     ) +
#     theme_minimal() -> g2
# 
# ggplot(resids) +
#     geom_point(mapping = aes(x = obs, y= uk_res),
#                color ="#013d09", alpha = .2) +
#     geom_line(mapping = aes(x = obs, y = -obs),
#               linetype = 2, color = "black") +
#     labs(
#         title = "Kriging Residuals",
#         x = element_blank(),
#         y = element_blank()
#     ) +
#     theme_minimal() -> g3
# 
# grid.arrange(g1, g3, g2, nrow = 1,
#              bottom = textGrob("Observed Values"))
