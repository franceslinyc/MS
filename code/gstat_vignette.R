library(sp)
data(meuse)

class(meuse)
str(meuse)

colnames(meuse)

coordinates(meuse) <- ~x+y
head(coordinates(meuse))

summary(meuse)

bubble(meuse,
       "zinc",
       main = "zinc concentrations (ppm)")

# ?bubble
# ?meuse.grid

data("meuse.grid")
summary(meuse.grid)
coordinates(meuse.grid) <- ~x+y
class(meuse.grid)

gridded(meuse.grid) <- TRUE
class(meuse.grid["dist"])
str(meuse.grid)

# ?image
# meuse.grid["dist"]

image(meuse.grid["dist"])
title("distance to river (red = 0)",
      xlab = "x",
      ylab = "y")
# no legend?

library(gstat)
zinc.idw <- idw(zinc~1, meuse, meuse.grid)
class(zinc.idw)

spplot(zinc.idw["var1.pred"],
       main = "zinc inverse distance weighted interpolations")
?spplot
head(meuse)
plot(log(zinc) ~ sqrt(dist),
     meuse)
abline(lm(log(zinc) ~ sqrt(dist), meuse),
       col = "red",
       lty = 2)

lzn.vgm <- variogram(log(zinc) ~ 1, meuse)
lzn.vgm

lzn.fit <- fit.variogram(lzn.vgm,
                         model = vgm(1, "Sph", 900, 1))
lzn.fit
plot(lzn.vgm, lzn.fit,
     main = "variogram")

lznr.vgm <- variogram(log(zinc) ~ sqrt(dist), meuse)
lznr.fit <- fit.variogram(lznr.vgm,
                          model = vgm(1, "Exp", 300, 1))
lznr.fit
plot(lznr.vgm, lznr.fit)

lzn.kriged <- krige(log(zinc) ~ 1,
                    meuse,
                    meuse.grid,
                    model = lzn.fit)
spplot(lzn.kriged["var1.pred"])

lzn.condsim <- krige(log(zinc) ~ 1,
                     meuse,
                     meuse.grid,
                     model = lzn.fit,
                     nmax =30,
                     nsim = 4)
spplot(lzn.condsim,
       main = "four conditional simulations")

lzn.condsim2 <- krige(log(zinc) ~ sqrt(dist),
                      meuse,
                      meuse.grid,
                      model = lznr.fit,
                      nmax =30,
                      nsim = 4)
spplot(lzn.condsim,
       main = "four UK conditional simulations")

lzn.dir <- variogram(log(zinc) ~ 1,
                     meuse,
                     alpha = c(0, 45, 90, 135))
lzndir.fit <- vgm(.59,
                  "Sph",
                  1200,
                  .05,
                  anis =c(45, .4))
plot(lzn.dir,
     lzndir.fit,
     as.table = TRUE)

# North is 0 degrees
# East is 90 degrees
# Point pairs are assigned to the 
# directional variogram with the closest direction.
# After being assigned a direction, point pairs are then
# binned by distance.

# anisotropic, anise top pic. ah nice, so tropic!

lznr.dir <- variogram(log(zinc) ~ sqrt(dist),
                      meuse,
                      alpha = c(0, 45, 90, 135))
plot(lznr.dir, lznr.fit, as.table = TRUE)

# Instead of classifying point pairs by direction and
# distance separately, we can classify them jointly.
# This map is a two dimensional heat map where x and y
# are the direction and distance of the point pairs.

vgm.map <- variogram(log(zinc) ~ sqrt(dist),
                     meuse,
                     cutoff = 1500,
                     width = 100,
                     map = TRUE)
plot(vgm.map, threshold = 5)

# linear model of coregionalization

g <- gstat(NULL, "log(zn)", log(zinc) ~ sqrt(dist), meuse)
g <- gstat(g, "log(cd)", log(cadmium) ~ sqrt(dist), meuse)
g <- gstat(g, "log(pb)", log(lead) ~ sqrt(dist), meuse)
g <- gstat(g, "log(cu)", log(copper) ~ sqrt(dist), meuse)
v <- variogram(g)

g <- gstat(g, model = vgm(1, "Exp", 300, 1), fill.all = TRUE)
g.fit <- fit.lmc(v, g)
g.fit

plot(v, g.fit)
vgm.map <- variogram(g, cutoff = 1500, width = 100, map = TRUE)
plot(vgm.map,
     threshold = 5,
     col.regions = bpy.colors(),
     xlab = "",
     ylab = "")

















