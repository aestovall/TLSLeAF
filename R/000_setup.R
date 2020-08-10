#Setup
#Load the necessary packages
library(filesstrings)
library(data.table)
library(lidR)
library(geometry)
library(ggplot2)
library(randomForest)

rf_model<-readRDS("leaf_wood_class_RF.rds")
source('R/angle_FUN.R')

# create input directory
if(!dir.exists("input")) dir.create("input")
setwd("input")

# create input directory
if(!dir.exists("output")) dir.create("output")
if(!dir.exists("figures")) dir.create("figures")

#Set optional parameters
SCATTER_LIM<-85 #threshold of scattering angle to remove'
correct.topography = TRUE #topographically normalize the voxels?

OS<-"mac" #What is your operating system?

# Add the path to your cloud compare executable
cloudcompare<-"/Applications/CloudCompare.app/Contents/MacOS/CloudCompare"
# On a Windows machine your path to CC may look like:
# cloudcompare<-"C:/Program Files/CloudCompare/cloudcompare.exe"

if(OS=="mac") run<-function(x) rstudioapi::terminalExecute(x) else shell(x)
