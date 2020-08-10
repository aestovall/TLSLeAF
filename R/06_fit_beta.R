#FIT BETA DISTRIBUTION
library(fitdistrplus)
m<-fitdist(LAD$a/90,
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

ggplot(beta,aes(x=a))+
  geom_density(color="black", alpha=0.5)+
  geom_path(data=site_beta, aes(x=a*90,y=y/85), linetype=2)+
  scale_fill_viridis(discrete=TRUE)+
  theme_bw()+
  geom_text(data=param,aes(x=5,y=0.025, group=site, label= label ), parse=TRUE, hjust=0)+
  geom_text(data=param,aes(x=5,y=0.022, group=site, label= label1 ), parse=TRUE, hjust=0)+
  xlim(0,90)+
  xlab("TLS-Measured Leaf Angle (degrees)")+
  ylab("Density")

ggsave(paste("LAD_Site_examples", ".pdf", sep=""),path = file.path(getwd(),"output/figures"), 
       width = 3.75*0.8, height = 3.75*0.8, units = "in")

