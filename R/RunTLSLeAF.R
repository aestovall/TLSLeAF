#### EDIT THE SETUP FILE FOR YOUR SYSTEM ####

#Make sure the setup file is correct and save
# The most important step is to have Cloud Compare installed
# and add the executable into the setup file
file.edit('R/000_setup.R')

# Run the setup file
source('R/000_setup.R')

###input file ####
files<-list.files(pattern = "ptx", 
                  recursive = TRUE, 
                  full.names = TRUE)
input_file = files[1]

#Find center coordinates. Modify according to file format
center<-fread(input_file, nrows=4, header=FALSE)[1,]
colnames(center)<-c("x","y","z")

df<-TLSLeAF(input_file, 
            overwrite=TRUE,
            center, 
            SCATTER_LIM=85,
            SS=0.02, 
            scales=c(0.1,0.5,0.75),
            rf_model,
            vox.res=5,
            minVoxDensity=5,
            superDF=TRUE)


param<-df@Beta_parameters
beta<-data.frame(a= seq(0.01,0.98, 0.01),y=dbeta(seq(0.01,0.98, 0.01), 
                                                 param$alpha0,param$beta0))

param$label<-paste0('alpha ==',round(param$alpha0, 1))
param$label1<-paste0('beta ==',round(param$beta0, 1))

ggplot(df@LAD,aes(x=a))+
  geom_density(fill="green", color=NA, alpha=0.5)+
  geom_path(data=beta, aes(x=a*90,y=y/90), linetype=2, size=1)+
  theme_bw()+
  geom_text(data=param,aes(x=5,y=max(beta$y/90)-0.002,label= label ), parse=TRUE, hjust=0)+
  geom_text(data=param,aes(x=5,y=max(beta$y/90)-0.003,label= label1 ), parse=TRUE, hjust=0)+
  xlim(0,90)+
  xlab("TLS-Measured Leaf Angle (degrees)")+
  ylab("Density")

clean.temp(output_file, c2c.file)



