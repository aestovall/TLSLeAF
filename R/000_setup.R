#Setup
#Load the necessary packages

# install.packages(c('filesstrings',
#                    'data.table',
#                    'lidR',
#                    'geometry',
#                    'ggplot2',
#                    'randomForest',
#                    'fitdistrplus'))

library(filesstrings)
library(data.table)
library(lidR)
library(geometry)
library(ggplot2)
library(randomForest)
library(fitdistrplus)

rf_model<-readRDS("leaf_wood_class_RF.rds")
source('R/angle_FUN.R')

# create directories
if(!dir.exists("input")) dir.create("input")
if(!dir.exists("output")) dir.create("output")
if(!dir.exists("figures")) dir.create("figures")

#Set optional parameters
SCATTER_LIM = 85 #threshold of scattering angle to remove'
correct.topography = TRUE #topographically normalize the voxels?
SS=0.02
scales=c(0.1,0.5,0.75)
vox.res=0.1
superDF=TRUE
clean=TRUE

#What is your operating system?
#OS<-"mac" 
OS<-"windows"

#Is CloudCompare already in the PATH?
PATH = FALSE

#If not in the PATH, what directory?
#Mac
#cc_dir = "/Applications/CloudCompare.app/Contents/MacOS/CloudCompare"
#Windows
cc_dir = shQuote('C:\\Program Files\\CloudCompare\\CloudCompare.exe')

# Add the path to your cloud compare executable
if(PATH & OS=='mac') cloudcompare<-'cloudcompare' else cloudcompare <- cc_dir

# On a Windows machine your path to CC may look like:
if(PATH & OS!='mac') cloudcompare<-'cloudcompare' else cloudcompare <- cc_dir


if(OS=="mac") run<-function(x) rstudioapi::terminalExecute(x) else run<-function(x) shell(x)
