#Estimate the average vertical leaf angle distribution 

leaf_angle<-data.frame(z=dat$z_cor,
                       z_bin = floor(dat$z_cor),
                       leaf = dat$dip_deg,
                       count = 1)

la_dist <- aggregate(leaf~z_bin, FUN = "mean", data = leaf_angle)
la_dist$sd <- aggregate(leaf~z_bin, FUN = "sd", data = leaf_angle)$leaf
n_dist <- aggregate(count~z_bin, FUN = "sum", data = leaf_angle)

la_dist$sderr<-la_dist$sd/sqrt(n_dist$count)
la_dist$ci_max<-la_dist$leaf + 1.96*la_dist$sderr
la_dist$ci_min<-la_dist$leaf - 1.96*la_dist$sderr

la_dist$ID<-gsub(".asc","",input_file)

plot(la_dist$z_bin,la_dist$leaf, col = "white")
lines(la_dist$z_bin,la_dist$leaf, col = "blue")
