dat<-data.table::fread(gsub(".asc","_angles_topocorrect.asc",output_file), header = TRUE)
dat<-dat[dat$class==0,]

#normalize measurment density with voxels
time<-Sys.time()
voxels<-LAvoxel(dat)
print(Sys.time()-time)#

fwrite(na.exclude(voxels), file = gsub(".asc","_angles_topocorrect_vox.asc",output_file), sep = " ")
rm(las, voxels)