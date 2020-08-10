#Setup
#Load the necessary packages
library(filesstrings)
library(data.table)
library(lidR)
library(geometry)
library(ggplot2)
source('R/angle_FUN.R')

if(!dir.exists("input")) dir.create("input")
setwd("input")

SCATTER_LIM<-85 #threshold of scattering angle to remove
rf_model<-readRDS("./leaf_wood_class_RF.rds")

OS<-"mac"
if(OS=="mac") run<-function(x) rstudioapi::terminalExecute(x) else shell(x)

# Add the path to your cloud compare executable
cloudcompare<-"/Applications/CloudCompare.app/Contents/MacOS/CloudCompare"

# On a Windows machine your path to CC may look like:
# cloudcompare<-"C:/Program Files/CloudCompare/cloudcompare.exe"

