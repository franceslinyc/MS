dat <- read.csv("ForestDataFuzzed.csv")
str(dat)

library(rgdal)
library(gridExtra)
alb.xy <- project(with(dat,cbind(lon_fuzzed,lat_fuzzed)),
                  "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
colnames(alb.xy) <- c("x","y")
dat <- cbind(alb.xy,dat)


# Plot total volume, a potential response variable of interest.
library(ggplot2)
ggplot(dat, aes(x=x/1000,y=y/1000,color=totvol)) +
    geom_point(size=2) +
    geom_point(stroke=1,
               shape = 21,
               color = "black"
               ) +
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
    ) -> g1

ggplot(dat, aes(x=x/1000,y=y/1000)) +
    geom_point(mapping = aes(color = annpre),
               size=2) +
    geom_point(stroke=1,
               shape = 21,
               color = "black"
               ) +
    scale_color_gradient(low = "white", high = "#013d09") +
    scale_x_continuous(breaks = round(seq(-2150, -2000, by = 50),2)) +
    labs(
        title = "Annual Precipitation in Sampled Plots",
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
    ) -> g2

grid.arrange(g1, g2, nrow =1)

ggplot(dat, aes(x = totvol)) +
    geom_histogram(fill = "#49ad57", color = "black",
                   bins = 20) +
    labs(
        title = "Total Volume",
        x = element_blank(),
        y = "Number of Plots"
    ) +
    theme(
        panel.background = element_rect(fill = 'white'),
        panel.grid = element_line(color = "grey", size = .1)
    ) -> g1

ggplot(dat, aes(x = totbiom)) +
    geom_histogram(fill = "#49ad57", color = "black",
                   bins = 20) +
    labs(
        title = "Total Biomass",
        x = element_blank(),
        y = "Number of Plots"
    ) +
    theme(
        panel.background = element_rect(fill = 'white'),
        panel.grid = element_line(color = "grey", size = .1)
    ) -> g2

ggplot(dat, aes(x = numtree)) +
    geom_histogram(fill = "#49ad57", color = "black",
                   bins = 20) +
    labs(
        title = "Number of Trees",
        x = element_blank(),
        y = "Number of Plots"
    ) +
    theme(
        panel.background = element_rect(fill = 'white'),
        panel.grid = element_line(color = "grey", size = .1)
    ) -> g3

ggplot(dat, aes(x = psmevol)) +
    geom_histogram(fill = "#49ad57", color = "black",
                   bins = 20) +
    labs(
        title = "PSME Volume",
        x = element_blank(),
        y = "Number of Plots"
    ) +
    theme(
        panel.background = element_rect(fill = 'white'),
        panel.grid = element_line(color = "grey", size = .1)
    ) -> g4

grid.arrange(g1, g2, g3, g4, ncol = 2)

# Explore some covariates

str(dat)
hist(dat$annpre)
hist(dat$anntmp)
hist(dat$tc3)
hist(dat$tc3for)

plot(dat$tc3, dat$annpre)

pairs(totvol ~ tc3 + annpre + anntmp + ndvi,
      data = dat,
      upper.panel = NULL)

pairs(numtree ~ tc3 + annpre + anntmp + ndvi,
      data = dat,
      upper.panel = NULL)
pairs(totvol ~ tc3 + annpre + anntmp + ndvi,
      data = dat,
      upper.panel = NULL)

ggplot(dat, aes(x=x/1000,y=y/1000,color=totvol)) +
    geom_point(size=2) +
    scale_color_gradient(low = "#c5e3ca", high = "#013d09") +
    scale_x_continuous(breaks = round(seq(-2150, -2000, by = 50),2)) +
    labs(
        title = "Total Timber Volume in Sampled Plots",
        subtitle = "Northwest Oregon",
        caption = "Alber's Equal Area Projection used here",
        x = "X (kilometers)",
        y = "Y (kilometers)",
        color = "Total Volume"
    ) +
    theme(
        axis.text.x=element_text(angle=90,hjust=1),
        panel.background = element_rect(fill = 'white'),
        panel.grid = element_line(color = "grey", size = .1)
    )
