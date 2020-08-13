#FIT BETA DISTRIBUTION
LAD_lim<-data.table::fread(gsub("input","output",gsub(".asc","_LAD.txt",output_file)), 
                                header = TRUE)
m<-fitdist(LAD_lim$a/90,
           'beta',
           method='mme')

# Get alpha and beta parametrs
alpha0 <- m$estimate[1] # parameter 1
beta0 <- m$estimate[2] # parameter 2
beta<-data.frame(a= seq(0.01,0.98, 0.01),y=dbeta(seq(0.01,0.98, 0.01), 
                                 alpha0,beta0))
param<-data.frame(alpha0,beta0)
param$label<-paste0('alpha ==',round(param$alpha0, 1))
param$label1<-paste0('beta ==',round(param$beta0, 1))

ggplot(LAD_lim,aes(x=a))+
  geom_density(fill="green", color=NA, alpha=0.5)+
  geom_path(data=beta, aes(x=a*90,y=y/90), linetype=2, size=1)+
  theme_bw()+
  geom_text(data=param,aes(x=5,y=max(beta$y/90)-0.002,label= label ), parse=TRUE, hjust=0)+
  geom_text(data=param,aes(x=5,y=max(beta$y/90)-0.003,label= label1 ), parse=TRUE, hjust=0)+
  xlim(0,90)+
  xlab("TLS-Measured Leaf Angle (degrees)")+
  ylab("Density")

ggsave(gsub("input","figures",gsub(".asc","_LAD.png",output_file)), 
       width = 3.75, height = 3.75, units = "in", dpi=300)

print(rbind(paste( "The alpha and beta parameters are:"),
paste("alpha=", param[1,1]),
paste("beta=", param[1,2])))
