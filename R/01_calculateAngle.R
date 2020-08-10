#Calculate leaf orientation and leaf angle
dat<-data.table::fread(output_file, header = FALSE)
colnames(dat)[1:10]<-c("X","Y","Z", "R","G","B","I","nX","nY","nZ")

#Adjust coordinates to scanner center
dat[,1:3]<-dat[,1:3]-t(c(center))

#ensures all imports are numeric, avoiding errors
dat<-as.data.frame(dat)
for(ii in 1:length(dat)) if (class(dat[[1,ii]])=="character") dat[,ii]<- as.numeric( dat[,ii] )

#calculate the radius, zenith, and azimuth angles from XYZ coordinates
dat<-RZA(dat)

# estimate scattering angle and filter, removing steep angles
dat$scatter<-scatter(dat)
dat<-dat[dat$scatter<=SCATTER_LIM,]
dat<-na.omit(dat)

#Convert normals to leaf orientation and leaf angle
dat<-cbind(dat[,1:3],
           ConvertNormalToDipAndDipDir(dat[,8:10])[,2])

fwrite(dat, file = gsub(".asc","_angles.asc",output_file), sep = " ")
