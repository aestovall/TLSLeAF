#### EDIT THE SETUP FILE FOR YOUR SYSTEM ####

#Make sure the setup file is correct and save
# The most important step is to have Cloud Compare installed
# and add the executable into the setup file
file.edit('R/000_setup_windows.R')

# Run the setup file
source('R/000_setup_windows.R')

#read the randomForest leaf/wood classification model
rf_model<-readRDS("leaf_wood_class_RF.rds")

###input file ####
files<-list.files("input", full.names = TRUE)
i<-1
input_file = files[i]
print(input_file)

#Find center coordinates. Modify according to file format (this is for *.ptx)
center<-fread(input_file, nrows=4, header=FALSE)[1,]
colnames(center)<-c("x","y","z")

#Run TLSLeAF
TLSLeAF.df<-TLSLeAF(input_file,
                    overwrite=FALSE,
                    center,
                    scatterLim=80,
                    SS=0.02,
                    scales=c(0.1,0.5,0.75),
                    rf_model=rf_model,
                    correct.topography = TRUE,
                    voxRes=0.1,
                    minVoxDensity=5,
                    superDF=TRUE,
                    clean=TRUE,
                    OS='windows')


library(ggplot2)
ggplot(TLSLeAF.df@LAD, 
       aes(a))+
  geom_density(fill="lightgreen")

ggplot(TLSLeAF.df@G, 
       aes(inc_bin, G_p, group=assumption, color=assumption))+
  facet_wrap(~density)+
  geom_line()

