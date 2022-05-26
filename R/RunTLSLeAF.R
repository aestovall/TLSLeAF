#### EDIT THE SETUP FILE FOR YOUR SYSTEM ####

#Make sure the setup file is correct and save
# The most important step is to have Cloud Compare installed
# and add the executable into the setup file
file.edit('R/000_setup.R')

# Run the setup file
source('R/000_setup.R')

#read the randomForest leaf/wood classification model
rf_model<-readRDS("leaf_wood_class_RF.rds")

###input file ####
files<-list.files("input", full.names = TRUE, pattern = ".ptx")
i<-2
input_file = files[i]
print(input_file)

#Find center coordinates. Modify according to file format (this is for *.ptx)
center<-fread(input_file, nrows=4, header=FALSE)[1,]
colnames(center)<-c("x","y","z")

#Run TLSLeAF
TLSLeAF.df<-TLSLeAF(input_file,
                    overwrite=TRUE,
                    center,
                    scatterLim=80,
                    SS=0.02,
                    scales=c(0.1,0.5,0.75),
                    rf_model=rf_model,
                    correct.topography = FALSE,
                    voxRes=0.1,
                    minVoxDensity=5,
                    superDF=TRUE)


library(ggplot2)
ggplot(TLSLeAF.df@LAD, 
       aes(a))+
  geom_density(fill="lightgreen")

ggplot(TLSLeAF.df@G, 
       aes(inc_bin, G_p, group=assumption, color=assumption))+
  facet_wrap(~density)+
  geom_line()

library(data.table)

tls.scan<-readTLS(input_file)
header<-ptx.header(input_file)

pgap.all<-pgap.angle(tls.scan,
                     z_res=1,
                     adj_z=1.6,
                     a_range=c(10,60),
                     a_bin = 5,
                     plane = FALSE,
                     header)
  
plot(pgap.all$list,pgap.all$pgap)

pai<-pgap2PAI(pgap_all = pgap.all, 
              method="A", 
              G = TLSLeAF.df@G[TLSLeAF.df@G$assumption=="non-random",])

plot(pai$list, pai$f.prime)
lines(pai$list, pai$f.prime, col="red")

