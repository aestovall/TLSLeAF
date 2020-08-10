#Import LA Voxels
dat<-data.table::fread(gsub(".asc","_angles_topocorrect_vox.asc",output_file), header = TRUE)

#Simulate density-equalized LAD -- each voxel has 10 measurments simulated 
#that represent the mean and standard deviation of that voxel.
LAD<-sim_LAD(dat$dip_dir,dat$dip_dir_sd)
