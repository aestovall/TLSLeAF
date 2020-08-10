# Normalize topography
dat<-data.table::fread(gsub(".asc","_class.asc",c2c.file), header = TRUE)
colnames(dat)[1:4]<-c("X","Y","Z", 'dip_deg')

#Create a surface model and topographically correct the TLS data
time<-Sys.time()
if(correct.topography == TRUE) dat$z_cor<-normalize_topography(dat) else dat$z_cor<-dat$Z
print(Sys.time()-time)#

fwrite(na.exclude(dat), file = gsub(".asc","_angles_topocorrect.asc",output_file), sep = " ")
