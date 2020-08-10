#### EDIT THE SETUP FILE FOR YOUR SYSTEM ####

#Make sure the setup file is correct and save
# The most important step is to have Cloud Compare installed
# and add the executable into the setup file
file.edit('R/000_setup.R')

# Run the setup file
source('R/000_setup.R')

###### RUN THE PIPELINE #########

###input file ####
files<-list.files(pattern = "ptx", 
                  recursive = TRUE, 
                  full.names = TRUE)
i=1

input_file = files[i]
output_file = gsub(".ptx",".asc", input_file)

#Find center coordinates. Modify according to file format
center<-fread(input_file, nrows=4, header=FALSE)[1,]

#Calculate normals for gridded TLS point cloud
source('../R/00_calculateNormals.R')

#Calculate scattering angle and leaf angle
source('../R/01_calculateAngle.R')

#Classify wood and leaf from random forest classifier
source('../R/02_classify_wood_leaves.R')

#Correct topography and calculate the LAD and vertical LAD
source('../R/03_normalize_topography.R')

# leaf angle voxelation and density normalization
source('../R/04_voxelize.R')

#simulate LAD from voxel statistics
source('../R/05_simulate_LAD.R')

#fit beta function and get beta parameters from LAD
source('../R/06_fit_beta.R')

#Remove unecessary files created during processing and free memory
file.remove(output_file, gsub(".asc","_angles.asc",output_file),
            gsub(".asc","_class.asc",c2c.file),
            gsub(".asc","_0_10_NORM.asc",c2c.file),
            gsub(".asc","_0_50_NORM.asc",c2c.file),
            gsub(".asc","_0_75_NORM.asc",c2c.file))
gc()

setwd("..")
