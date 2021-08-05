#####Setup#####

#Install the necessary packages
# install.packages(c('filesstrings',
#                    'data.table',
#                    'lidR',
#                    'geometry',
#                    'ggplot2',
#                    'randomForest',
#                    'fitdistrplus'))

#### Set optional parameters ####
# scatterLim = 85          #Threshold of scattering angle to remove'
# correct.topography = TRUE #topographically normalize the voxels?
# SS = 0.02                 #Spatial subsampling resolution
# scales = c(0.1,0.5,0.75)  #3 scales of normal computation. Set for RF model
# voxRes = 0.1             #voxel resolution for LAD normalization
# minVoxDensity = 5         #minimum number of measurments per voxel
# superDF = TRUE            #Merges output into a single TLSLeAF.class
# clean = TRUE              #removes temporary files created during processing
# overwrite = FALSE
# SCATTER_LIM = 80
# 
# overwrite=FALSE
# scatterLim=80
# SS=0.02
# scales=c(0.1,0.5,0.75)
# # rf_model,
# # rf_model_path=rf_model_path
# correct.topography = TRUE
# voxRes=0.1
# minVoxDensity=5
# superDF=TRUE
# clean=FALSE

#### SETUP CloudCompare ####
#What is your operating system?
#OS<-"mac" 
OS<-"windows" 

#Is CloudCompare already in the PATH?
PATH = FALSE

#If not in the PATH, what directory? 
#Edit path to CloudCompare in quotes.

#Mac
if(!PATH & OS=='mac') cc_dir = 
  
  #"/Applications/CloudCompare.app/Contents/MacOS/CloudCompare"
  '/Users/aestoval/Desktop/CloudCompare.app/Contents/MacOS/CloudCompare'
#Windows
if(!PATH & !OS=='mac') cc_dir = 
  
  shQuote('C:\\Program Files\\CloudCompare\\CloudCompare.exe')



#### NO EDITS REQUIRED BELOW THIS LINE ####
#
#
#

library(filesstrings)
library(data.table)
library(lidR)
# install.packages('raster')
library(raster)
library(geometry)
library(ggplot2)
library(randomForest)
library(fitdistrplus)
library(vroom)

#load the randomForest model
rf_model<-readRDS("leaf_wood_class_RF.rds")
rf_model_path<-"leaf_wood_class_RF.rds"

#load TLSLeAF functions
source('R/TLSLeAF_FUN_update_v2.R')

# create directories
if(!dir.exists("input")) dir.create("input")
if(!dir.exists("output")) dir.create("output")
if(!dir.exists("figures")) dir.create("figures")

# Add the path to your cloud compare executable
if(PATH & OS=='mac') cloudcompare<-'cloudcompare' else cloudcompare <- cc_dir
if(PATH & OS!='mac') cloudcompare<-'cloudcompare' else cloudcompare <- cc_dir
if(OS=="mac") run<-function(x) rstudioapi::terminalExecute(x) else run<-function(x) shell(x, intern=FALSE, wait=TRUE)

