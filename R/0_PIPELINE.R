#Load the necessary packages
library(filesstrings)
library(data.table)
library(lidR)
source('R/angle_FUN.R')

setwd("input")

#Calculate normals for TLS point cloud and convert to leaf angle
source('../R/00_calculateNormals.R')
source('../R/01_calculateAngle.R')

#Correct topography and calculate the LAD and vertical LAD
correct.topography = TRUE
source('../R/02_calculateLAD.R')
source('../R/03_calculateVerticalLAD.R')

#Remove unecessary files created during processing and free memory
file.remove(output_file)
setwd("..")
rm(dat,ground,las,leaf_angle,r,slope,topo,topo.df,topo.las,topo.las.r)
gc()
