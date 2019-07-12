#Calculate leaf orientation and leaf angle

dat<-data.table::fread(output_file, header = FALSE)
colnames(dat)[1:10]<-c("X","Y","Z", "R","G","B","I","nX","nY","nZ")

#Convert normals to leaf orientation and leaf angle
dat<-cbind(dat,ConvertNormalToDipAndDipDir(dat[,8:10]))
