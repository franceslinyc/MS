# Playing around with splitting the area into three regions
# We assume that each region can be modeled using a ZIG distribution
# We then use the Gaussian copula with these marginal distributions

setwd("/mnt/c/Users/njnic/Documents/MS_Project_copy")
# simdata <- readRDS("../simulated_tshevol.rds")
dat <- read.csv("ForestDataFuzzed.csv")

source("ZIGFunctions.R")
source("metrics.R")
source("simulation_study/preprocess.R")
setwd("/mnt/c/Users/njnic/Documents/MS_Project_copy/simulation_study")
source("copula_func.R")
setwd("/mnt/c/Users/njnic/Documents/MS_Project_copy")

library(gstat)
library(tidyverse)
library(rgdal)
library(sp)
library(automap)
library(ggplot2)

alb.xy <- project(with(dat,cbind(lon_fuzzed,lat_fuzzed)),
                  "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
colnames(alb.xy) <- c("x","y")
dat <- cbind(alb.xy,dat)

ggplot(dat, aes(x=x/1000,y=y/1000,color=totvol)) +
    geom_point(size=2) +
    geom_point(stroke=1,
               shape = 21,
               color = "black"
               ) +
	geom_vline(xintercept = -2125,
			   color = "red",
			   lty = 2) +
	geom_vline(xintercept = -2050,
			   color = "red",
			   lty = 2) +
    scale_color_gradient(low = "white", high = "#013d09") +
    scale_x_continuous(breaks = round(seq(-2150, -2000, by = 50),2)) +
    labs(
        title = "Total Timber Volume in Sampled Plots",
        subtitle = "Northwest Oregon",
        x = "X (kilometers)",
        y = "Y (kilometers)",
        color = "Total Volume"
    ) +
    theme(
        axis.text.x=element_text(angle=90,hjust=1),
        panel.background = element_rect(fill = 'white'),
        panel.grid = element_line(color = "grey", size = .1),
        plot.caption = element_text(hjust = 0, face = "italic")
    )

x11()

dim(dat[dat$x/1000 <= -2125,])
dim(dat[dat$x/1000 >= -2050,])

hist(dat[dat$x/1000 <= -2125,"totvol"], main = "left region")
hist(dat[dat$x/1000 >= -2050,"totvol"], main = "right region")
hist(dat[dat$x/1000 > -2125 & dat$x/1000 < -2050,"totvol"], main = "center region")

summary(dat[dat$x/1000 <= -2125,"totvol"])
summary(dat[dat$x/1000 >= -2050,"totvol"])
summary(dat[dat$x/1000 > -2125 & dat$x/1000 < -2050,"totvol"])

mean(dat[dat$x/1000 <= -2125,"totvol"] == 0)
mean(dat[dat$x/1000 >= -2050,"totvol"] == 0)
mean(dat[dat$x/1000 > -2125 & dat$x/1000 < -2050,"totvol"] == 0)

full_simdata <- data.frame(
	x = dat$x,
	y = dat$y,
	resp = dat$totvol
)
full_simdata %>% 
	mutate(region = case_when(x/1000 <= -2125 ~ 1,
							  x/1000 >= -2050 ~ 3,
							  x/1000 > -2125 & x/1000 < -2050 ~ 2)) -> full_simdata

set.seed(182)
tmp <- split_data(full_simdata, n = 300)
colnames(tmp$train) <- c("x", "y", "resp", "region")
colnames(tmp$test) <- c("x", "y", "resp", "region")

full_simdata %>% 
	count(region)

tmp$train %>% 
	count(region)
tmp$test %>% 
	count(region)

n.regions <- 3
mu <- vector(mode = "list", length = n.regions)
beta <- vector(mode = "list", length = n.regions)

# section to be put into a loop
current.training.set <- subset(tmp$train, region == 1)
current.testing.set <- subset(tmp$test, region == 1)

transformed.training.set <- cube_root_fix_zeros(current.training.set)
theta <- calculate_spatial_params(simdata = transformed.training.set$transformed_data,
								  testdata = current.testing.set,
								  zeros = transformed.training.set$zeros)
z <- qnorm(pzig(y = transformed.training.set$transformed_data$resp,
			    mu = theta$mu,
			    beta = theta$beta,
			    epsilon = transformed.training.set$epsilon,
			    Pi = transformed.training.set$Pi))
transformed.training.set$transformed_data$z <- z
transformed.training.set$transformed_data$region <- 1

zhat <- ord_kriging_z(training_data = training_data, test_points = test_data, test_z = z)$krige_output$var1.pred

pz <- pnorm(zhat)
pz[which(pz==1)] <- 1-.Machine$double.eps
pz[which(pz==0)] <- .Machine$double.eps

zigs <- qzig(u = pz,
		   mu = rep(mu[1], length(pz)),
		   beta = beta,
		   epsilon = epsilon,
		   Pi = rep(Pi[1], length(pz)))

krige_zhat_gauscop(resp = transformed.training.set$transformed_data$resp,
				   mu = theta$mu,
				   beta = theta$beta,
				   epsilon = transformed.training.set$epsilon,
				   Pi = transformed.training.set$Pi,
				   training_data = transformed.training.set$transformed_data,
				   test_data = current.testing.set) -> zhat.zigs.df.krige

names(theta)
names(zhat.zigs.df.krige)
zhat.zigs.df.krige$zigs

yhat_gauscop_krige <- backtransform_gauscop(zhat.zigs.df.krige$zigs, transformed.training.set$epsilon)
