#### EDIT THE SETUP FILE FOR YOUR SYSTEM ####

#Make sure the setup file is correct and save
# The most important step is to have Cloud Compare installed
# and add the executable into the setup file
file.edit('R/000_setup.R')

# Run the setup file
source('R/000_setup.R')

###input file ####
files<-list.files("input",pattern = "ptx", 
                  recursive = TRUE, 
                  full.names = TRUE)


input_file = files[i]

#Find center coordinates. Modify according to file format
center<-fread(input_file, nrows=4, header=FALSE)[1,]
colnames(center)<-c("x","y","z")

df<-TLSLeAF(input_file, 
        overwrite=FALSE,
        center, 
        scatterLim=80,
        SS=0.02, 
        scales=c(0.1,0.5,0.75),
        rf_model,
        correct.topography = TRUE,
        voxRes=0.1,
        minVoxDensity=5,
        superDF=TRUE)



