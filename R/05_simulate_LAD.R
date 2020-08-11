#Import LA Voxels
dat<-data.table::fread(gsub("input","output",gsub(".asc","_angles_topocorrect_vox.asc",output_file)), 
                       header = TRUE)

#Simulate density-equalized LAD -- each voxel has 10 measurments simulated 
#that represent the mean and standard deviation of that voxel.
LAD<-sim_LAD(dat$dip_dir,dat$dip_dir_sd)

LAD<-na.exclude(LAD[LAD$a>=0 & LAD$a<=90,])
LAD<-data.frame(a=LAD)

fwrite(LAD, 
       file = gsub("input","output",gsub(".asc","_LAD.txt",output_file)),
       sep = " ")

rm(LAD)
