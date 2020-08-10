c2c.file<-gsub(".asc","_angles.asc",output_file)
rm(dat)
gc()

print(paste("Processing", input_file))
termId<-run(paste(cloudcompare, # call Cloud Compare. The .exe file folder must be in the system PATH
                  "-SILENT",
                  "-C_EXPORT_FMT", "ASC", "-PREC", 6, #Set asc as export format
                  "-NO_TIMESTAMP",
                  "-AUTO_SAVE OFF",
                  "-O", c2c.file, #open the subsampled file
                  "-SS SPATIAL 0.02",
                  "-OCTREE_NORMALS", 0.1,
                  "-SAVE_CLOUDS","FILE", gsub(".asc","_0_10_NORM.asc",c2c.file),
                  "-OCTREE_NORMALS", 0.5,
                  "-SAVE_CLOUDS","FILE", gsub(".asc","_0_50_NORM.asc",c2c.file),
                  "-OCTREE_NORMALS", 0.75,
                  "-SAVE_CLOUDS", "FILE", gsub(".asc","_0_75_NORM.asc",c2c.file),
                  sep = " "))

while (is.null(rstudioapi::terminalExitCode(termId))) {
  Sys.sleep(0.1)
}

result <- rstudioapi::terminalBuffer(termId)

# Delete the buffer and close the session in the IDE
rstudioapi::terminalKill(termId)

#######IMPORT NORMALS#########
dat<-data.table::fread(gsub(".asc","_0_10_NORM.asc",c2c.file),header = FALSE)
colnames(dat)[1:7]<-c("X","Y","Z","angle","nX10","nY10","nZ10")

dat1<-data.table::fread(gsub(".asc","_0_50_NORM.asc",c2c.file), select = c(5,6,7),header = FALSE)
colnames(dat1)[1:3]<-c("nX50","nY50","nZ50")

dat2<-data.table::fread(gsub(".asc","_0_75_NORM.asc",c2c.file), select = c(5,6,7),header = FALSE)
colnames(dat2)[1:3]<-c("nX75","nY75","nZ75")

dat<-cbind(dat,dat1,dat2)
remove(dat1,dat2)
gc()

#ensures all imports are numeric, avoiding errors
dat<-as.data.frame(dat)
for(ii in 1:length(dat)) if (class(dat[[1,ii]])=="character") dat[,ii]<- as.numeric( dat[,ii] )

#######Steps used to train the RF model#######
# library(randomForest)
# require(caret)
# 
# train_rows<-sample(nrow(dat),nrow(dat)*0.25)
# train<-dat[train_rows,6:15]
# validate<-dat[!train_rows, 6:15]
# 
# train$class<-as.factor(train$class)
# 
# rf_model<-randomForest(train[,2:10], train$class, do.trace = TRUE, ntree=100)
# print(rf_model)
# 
# validate$predict<-predict(rf_model, validate[,1:10])
# validate$class<-as.factor(validate$class)

# save the model to disk
# saveRDS(rf_model, "./leaf_wood_class_RF.rds")

#RF model import is now at the beginning of the script
# rf_model<-readRDS("./leaf_wood_class_RF.rds")

#### Use RF Model to classify Leaves #####
dat<-na.omit(dat) #RF cannot accept NA values

time<-Sys.time()
dat$predict<-predict(rf_model, dat)
print(Sys.time()-time)#

dat<-cbind(dat[,1:4], dat$predict)
colnames(dat)[4:5]<-c("dip_deg", "class")

fwrite(dat, file = gsub(".asc","_class.asc",c2c.file), sep = " ")
