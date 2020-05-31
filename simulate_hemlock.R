library(rgdal)
library(tidyverse)
# Actual forest inventory data but with fuzzed locations.
dat <- read.csv("ForestDataFuzzed.csv")
# Perform Alber's equal area map projection. 
# Resulting xy coordinates are in meters.
alb.xy <- project(with(dat,cbind(lon_fuzzed,lat_fuzzed)),
                  "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
colnames(alb.xy) <- c("x","y")
dat <- cbind(alb.xy,dat)

colnames(dat)
sum(dat$psmevol == 0)
sum(dat$tshevol == 0)
sum(dat$totvol == 0)

# Plot total volume, a potential response variable of interest.
library(ggplot2)
ggplot(dat, aes(x=x/1000,y=y/1000,color=tshevol)) +
  geom_point(size=2) +
  scale_color_gradient(low = "white", high = "red") +
  scale_x_continuous(breaks = round(seq(-2150, -2000, by = 50),2)) +
  theme(axis.text.x=element_text(angle=90,hjust=1))

# Simulate data that resembles totvol (total volume). For now, don't consider 
# covariates,but tc3, forind, annpre, anntmp, smrtp, ndvi, and cover are all
# potential covariates. totvol, psmevol, tshevol, totbiom, and numtree are 
# all potential response variables.

# The model for totvol will be a zero-inflated gamma model.
# Explore the data to find reasonable parameters to use in a simulation.

zeros <- which(dat$tshevol==0) # Which observations are 0?
(epsilon <- min(dat$tshevol[-zeros])) # Smallest non-0
Pi <- length(zeros)/length(dat$tshevol) # Proportion of 0's
# Anticipating modeling Pi as a function of parameters, so each obs. will have its own Pi.
Pi <- rep(Pi, times=length(dat$tshevol)) 


# Estimate the gamma parameters via Method of Moments on the non-0 data.
mu <- mean(dat$totvol[-zeros])
beta <- var(dat$totvol[-zeros])/mu
# Each obs. may have its own mu also, but beta remains scalar.
mu <- rep(mu, times=length(dat$totvol))

# Create the smoothed-out zero-inflated gamma response. The zeros are replaced by 
# random unif(0,epsilon)'s. This makes everything continuous, which is nicer
# in copula-land.

resp <- dat$tshevol
resp[zeros] <- runif(length(zeros),0,epsilon)

# Create distance matrix
H <- as.matrix(dist(alb.xy))
# write.table(H,
#             file = "distance_matrix.txt",
#             row.names = FALSE,
#             col.names = FALSE)

source("ZIGFunctions.R") # Functions for zero-inflated gamma.

# Estimate the spatial parameters. Create Nresp_df, a data frame of smoothed response
# transformed to standard normal via the Gaussian copula double-transformation.
# Fit a variogram to that. Without covariates, the Gaussian copula model fits,
# but the exponential and spherical don't.
Nresp_df <- data.frame(x=dat$x,y=dat$y,
                       z=qnorm(pzig(resp,mu,beta,epsilon,Pi)))
coordinates(Nresp_df) <- ~x+y
library(automap)
v_fit<-autofitVariogram(z~1, input_data = Nresp_df)
plot(v_fit)

names(v_fit)
v_fit$var_model$psill

# Translate to spatial correlogram parameters, since we will be simulating 
# marginally standard normal data.
(alphaN<-v_fit$var_model$psill[2]/sum(v_fit$var_model$psill)) #psill/(nugget+psill) from fitted variogram to ranked counts
(alphaR<-v_fit$var_model$range[2])

# create copula correlation matrix
Sigma<-alphaN*exp(-(H/alphaR)^2)
diag(Sigma)<-1

# Simulate data
set.seed(183)
nsim <- 1000 # How many data sets.
Usim <- matrix(0,nrow=nsim, ncol=ncol(Sigma))
Ysim <- Usim

require(mvtnorm)
Usim <- pnorm(rmvnorm(nsim,sigma=Sigma)) # Simulate normals, then transform to uniform.
Ysim <- t(apply(Usim,1,function(x) qzig(x,mu = rep(mu, ncol(Sigma)),beta,epsilon,Pi))) # transform uniforms to zigs
Ysim[which(Ysim<epsilon)] <- 0 # Replace the data smaller than epsilon with 0.

saveRDS(Ysim,
        file = "simulated_tshevol.rds")

# Ys <- readRDS("simulated_Y.rds")
# write.table(Ysim,
#             file="simulated_data_sets.txt",
#             row.names = FALSE,
#             col.names = FALSE)

# Map one of the simulated data sets.
Y <- Ysim[2,]
ggplot(dat, aes(x=x/1000,y=y/1000,color=Y)) +
  geom_point(size=2) +
  scale_color_gradient(low = "white", high = "red") +
  scale_x_continuous(breaks = round(seq(-2150, -2000, by = 50),2)) +
  theme(axis.text.x=element_text(angle=90,hjust=1))

hist(apply(Ysim, 1, mean),
     main = "Distribution of Simulated Means vs. Actual Mean")
abline(v = mean(resp),
       col = "red",
       lty = 2)

# Now we have a simulated dataset -- great!
# Let's start by trying transformations to make the response "more normal"
# Lisa recommended using a cubic root transformation

train_idx <- sample(1:1224, 300, replace = FALSE)
train_idx
Y <- Ysim[1, train_idx]
mean(Y); mean(resp)

data.frame(x=dat$x[train_idx],
           y=dat$y[train_idx],
           z=Y) -> sim_data

ggplot(sim_data, aes(x=x/1000,y=y/1000,color=z)) +
  geom_point(size=2) +
  scale_color_gradient(low = "white", high = "red") +
  scale_x_continuous(breaks = round(seq(-2150, -2000, by = 50),2)) +
  theme(axis.text.x=element_text(angle=90,hjust=1))

hist(Y, breaks = 20)
hist(resp, breaks = 20)
Y <- Y^(1/3)
hist(Y^(1/3))

zeros <- which(Y==0) 
(epsilon <- min(Y[-zeros])) 
(Pi <- length(zeros)/length(Y))
Pi <- rep(Pi, length(Y))
Y[zeros] <- runif(length(zeros),0,epsilon)

# method of moment estimators
mu <- mean(Y[-zeros])
beta <- var(Y[-zeros])/mu
mu <- rep(mu, times=length(Y))

tr_std_norms <- qnorm(pzig(Y, mu, beta, epsilon, Pi))
hist(pzig(Y, mu, beta, epsilon, Pi))
hist(tr_std_norms)

sim_data$z <- Y
emp_var<-variogram(z~1,
                   loc=~x+y,
                   sim_data)
plot(emp_var)
v_fit<-fit.variogram(emp_var,vgm("Gau"))
plot(emp_var,v_fit)
v_fit

(alphaN<-v_fit$psill[2]/sum(v_fit$psill))
(alphaR<-v_fit$range[2])

# Let's try following Lisa's work
# This method uses per-observation estimates of the GLM parameters instead of MOM
#

H <- as.matrix(dist(sim_data[,c("x", "y")]))

Sigma<-alphaN*exp(-(H/alphaR)^2)
diag(Sigma)<-1

X <- matrix(cbind(rep(1, 300),
                  dat$annpre[train_idx]),
            nrow = 300)
zero_resp_train <- ifelse(Y < epsilon, 1, 0)
zero_glm_train <- glm(zero_resp_train ~ dat$annpre[train_idx],
                      family=binomial)
xi_Z <- zero_glm_train$coefficients

nonzero_resp_idx <- which(Y>=epsilon)
nonzero_glm_train <- glm(Y[nonzero_resp_idx] ~ X[nonzero_resp_idx, 2],
                         family=Gamma("log"))
xi_NZ <- nonzero_glm_train$coefficients
summary(nonzero_glm_train)

Pi <- alogit(X %*% xi_Z)
mu <- exp(X %*% xi_NZ)
s2 <- var(Y[Y > 0])
beta <- mean(s2/mu)

ZY <- data.frame(x = sim_data$x,
                 y = sim_data$y,
                 z = qnorm(pzig(Y, mu, beta, epsilon, Pi)))

ZY_var <- variogram(z~1,
                    loc=~x+y,
                    ZY)
plot(ZY_var)
ZY_fit <- fit.variogram(ZY_var,
                       vgm(psill = .6, "Gau", range = 80000, nugget = .2))
plot(ZY_var, ZY_fit)

(alphaN <- ZY_fit$psill[2]/sum(ZY_fit$psill))
(alphaR <- ZY_fit$range[2])

theta <- c(logit(alphaN),
           log(alphaR),
           log(1/beta),
           xi_Z,
           xi_NZ)
names(theta) <- NULL

start <- Sys.time()
result <- optim(theta, NL, 
                method = "Nelder-Mead",
                Y = Y,
                epsilon = epsilon,
                H = H,
                X = X)
Sys.time() - start

data.frame(
  c(alphaN, alphaR, beta, xi_Z, xi_NZ),
  result$par
)

# Let's use these parameters to predict the values for the other test points

test_idx <- c(1:1224)[-train_idx]
Y2 <- Ysim[1, test_idx]
hist(Y2)

X <- matrix(cbind(rep(1, length(dat$annpre)),
                  dat$annpre),
            nrow = length(dat$annpre))

xi_Z <- matrix(result$par[4:(3+ncol(X))], ncol=1)
xi_NZ <- matrix(result$par[(4+ncol(X)):(length(result$par))], ncol=1)
Pi <- alogit(X %*% xi_Z)
mu <- exp(X %*% xi_NZ)
beta <- 1/(exp(result$par[3]))
alphaN <- alogit(result$par[1])
alphaR <- exp(result$par[2])

H <- as.matrix(dist(dat[,c("x", "y")]))
test_dist <- data.frame(x = dat$x[test_idx],
                        y = dat$y[test_idx])
H <- as.matrix(dist(test_dist))
Sigma<-alphaN*exp(-(H/alphaR)^2)
diag(Sigma)<-1



ZY_test <- data.frame(x = sim_data$x,
                      y = sim_data$y,
                      z = qnorm(pzig(Y, mu, beta, epsilon, Pi)))

