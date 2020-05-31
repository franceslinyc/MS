# Should be running R version 3.6.x

install.packages(c("sp", 
                   "gstat",
                   "ranger",
                   "rgdal",
                   "geoR",
                   "GSIF",
                   "raster",
                   "fields",
                   "automap"))

# For plotting and data cleaning:

install.packges("tidyverse")

# For parallel programming:

install.packages(c("doParallel",
                   "foreach"))