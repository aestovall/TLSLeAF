
#### EDIT THE SETUP FILE FOR YOUR SYSTEM ####

#Make sure the setup file is correct and save
# The most important step is to have Cloud Compare installed
# and add the executable into the setup file

# file.edit('R/000_setup.R')

source('R/000_setup.R')

###### RUN THE PIPELINE #########

#Calculate normals for TLS point cloud and convert to leaf angle
source('../R/00_calculateNormals.R')
source('../R/01_calculateAngle.R')
source('../R/02_classify_wood_leaves.R')

#Correct topography and calculate the LAD and vertical LAD
correct.topography = TRUE
source('../R/02_calculateLAD.R')
source('../R/03_calculateVerticalLAD.R')

#Remove unecessary files created during processing and free memory
file.remove(output_file)
setwd("..")
rm(dat,ground,las,leaf_angle,r,slope,topo,topo.df,topo.las,topo.las.r)
gc()
