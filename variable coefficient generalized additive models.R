# This script fits variable coefficient GAMs to test the hypothesis
# of era-dependent inverse production regimes / responses to SST variability

library(mgcv)
library(nlme)
library(maps)
library(plotfunctions) # for alpha transparency
source('/Users/lciannel/Documents/MyDocuments/R/distance.function.R')

setwd("/Users/lciannel/Documents/MyDocuments/Funding/NSF/goa/CalCOFI/coastwide-salmon-master")
data.dump<-read.csv('coastwide salmon data with different sst groupings.csv', header=T,stringsAsFactors=F);
data.dump$our.region<-data.dump$region;data.dump$loc.sst<-data.dump$sst
head(data.dump);dim(data.dump)
table(data.dump$species)

#Chum data
chum<-data.dump[data.dump$species=='Chum',];
chum$ln.rs<-log(chum$recruits/chum$spawners);
chum$era<-as.factor(chum$era);chum$stock<-as.factor(chum$stock)
#Chum least square line and distance from southernmost point
lon.chum<--tapply(chum$long,chum$stock,mean)
lat.chum<-tapply(chum$lat,chum$stock,mean)
lm.chum<-lm(lat.chum~lon.chum)
summary(lm.chum)
#Project each stock onto the least square line
lat.chum.prj<-coef(lm.chum)[1]+coef(lm.chum)[2]*lon.chum
#Calculate shortest distance (km) from northwest origin of least square line to projection point of each stock on the least square line
stock.dist<-lon.chum*NA
for(i in 1:length(lon.chum)){
	stock.dist[i]<-distance.function(lat.chum.prj[i],lon.chum[i],max(lat.chum.prj),min(lon.chum))/1000
}
range(stock.dist)
#Make new column on the file with the stock distance from the origin of the least square line
index<-match(chum$stock,names(stock.dist))
chum$stock.dist<-stock.dist[index]
head(chum);dim(chum)

#Pink data
pink<-data.dump[data.dump$species=='Pink',];
pink$ln.rs<-log(pink$recruits/pink$spawners);
pink$era<-as.factor(pink$era);pink$stock<-as.factor(pink$stock)
#Least square line and distance from southernmost point
lon.pink<--tapply(pink$long,pink$stock,mean)
lat.pink<-tapply(pink$lat,pink$stock,mean)
lm.pink<-lm(lat.pink~lon.pink)
summary(lm.pink)
#Project each stock onto the least square line
lat.pink.prj<-coef(lm.pink)[1]+coef(lm.pink)[2]*lon.pink
#Calculate shortest distance (km) from northwest origin of least square line to projection point of each stock on the least square line
stock.dist<-lon.pink*NA
for(i in 1:length(lon.pink)){
	stock.dist[i]<-distance.function(lat.pink.prj[i],lon.pink[i],max(lat.pink.prj),min(lon.pink))/1000
}
range(stock.dist)
#Make new column on the file with the stock distance from the origin of the least square line
index<-match(pink$stock,names(stock.dist))
pink$stock.dist<-stock.dist[index]
head(pink);dim(pink)

#Sockeye data
sock<-data.dump[data.dump$species=='Sockeye',];
sock$ln.rs<-log(sock$recruits/sock$spawners);
sock$era<-as.factor(sock$era);sock$stock<-as.factor(sock$stock)
#Least square line and distance from southernmost point
lon.sock<--tapply(sock$long,sock$stock,mean)
lat.sock<-tapply(sock$lat,sock$stock,mean)
lm.sock<-lm(lat.sock~lon.sock)
summary(lm.sock)
#Project each stock onto the least square line
lat.sock.prj<-coef(lm.sock)[1]+coef(lm.sock)[2]*lon.sock
#Calculate shortest distance (km) from northwest origin of least square line to projection point of each stock on the least square line
stock.dist<-lon.sock*NA
for(i in 1:length(lon.sock)){
	stock.dist[i]<-distance.function(lat.sock.prj[i],lon.sock[i],max(lat.sock.prj),min(lon.sock))/1000
}
range(stock.dist)
#Make new column on the file with the stock distance from the origin of the least square line
index<-match(sock$stock,names(stock.dist))
sock$stock.dist<-stock.dist[index]
head(sock);dim(sock)


#Check time series of r/s (note that for pink must use 'full.stock')
xyplot(ln.rs ~ entry.yr | factor(stock), type = "l", col = 1,
  xlab = "Entry year",
  ylab = 'ln(r/s)',
  strip = function(bg = 'white', ...)
  strip.default(bg = c('white'), ...),
  data = pink)

apply(1*(tapply(pink$long,list(pink$stock,pink$region),mean,na.rm=T)>0),2,sum,na.rm=T)

range(tapply(pink$entry.yr,pink$stock,min,na.rm=T))
range(tapply(pink$entry.yr,pink$stock,max,na.rm=T))

#Check location
dev.new(width=13,height=5)
par(mfrow=c(1,3))
plot(-sock$long,sock$lat,ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),ylim=range(pink$lat),xlim=range(-pink$long),pch=16,main='Sockeye',col='red')
map("world",fill=T,col="lightblue4",add=T)
points(-sock$long,sock$lat,col="red",pch=16)
N<-length(unique(sock$stock))
text(-132,64,paste('N = ',as.character(N)),col='white',cex=2)
abline(lm.sock)
points(lon.sock,lat.sock.prj,pch='+',col='blue')
plot(-chum$long,chum$lat,ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),ylim=range(pink$lat),xlim=range(-pink$long),pch=16,main='Chum',col='red')
map("world",fill=T,col="lightblue4",add=T)
points(-chum$long,chum$lat,col='red',pch=16)
N<-length(unique(chum$stock))
text(-132,64,paste('N = ',as.character(N)),col='white',cex=2)
abline(lm.chum)
points(lon.chum,lat.chum.prj,pch='+',col='blue')
plot(-pink$long,pink$lat,ylab=expression(paste("Latitude ("^0,'N)')),xlab=expression(paste("Longitude ("^0,'E)')),ylim=range(pink$lat),xlim=range(-pink$long),pch=16,main='Pink',col='red')
map("world",fill=T,col="lightblue4",add=T)
points(-pink$long,pink$lat,col='red',pch=16)
N<-length(unique(pink$stock))
text(-132,64,paste('N = ',as.character(N)),col='white',cex=2)
abline(lm.pink)
points(lon.pink,lat.pink.prj,pch='+',col='blue')

#Sockeye: VC GAM, spatial 

#Early period
sock$loc.sst2<-sock$sst1+1+abs(min(sock$sst1))
sock$long2<--sock$long
#Early
tvc.early.sock<-gamm(ln.rs~s(spawners,k=4)+ s(stock.dist,k=5)+s(stock.dist,by=loc.sst2,k=5)+ s(stock,bs='re'),data=sock[sock$era==1,],correlation=corAR1(form= ~ 1 | stock))  
summary(tvc.early.sock$gam)
#R-sq.(adj) =  0.477   
# Scale est. = 0.6357    n = 729
AIC(tvc.early.sock$lme)
#1719.157
plot(tvc.early.sock $gam,pages=1,scale=0)

sock$loc.sst2<-sock$sst1+1+abs(min(sock$sst1))
gam.early.sock<-gamm(ln.rs~s(spawners,k=4)+s(stock.dist,k=5) +loc.sst2+ s(stock,bs='re'),data=sock[sock$era==1,],correlation=corAR1(form= ~ 1 | stock))  
summary(gam.early.sock$gam)
#R-sq.(adj) =  0.444   
# Scale est. = 0.67625   n = 729
AIC(gam.early.sock$lme)
#1741.922
plot(gam.early.sock$gam,pages=1,scale=0)

anova(tvc.early.sock$lme,gam.early.sock$lme)
#tvc.early.sock$lme     1 11 1719.157 1769.665 -848.5785                       
#gam.early.sock$lme     2  9 1741.922 1783.247 -861.9611 1 vs 2 26.7652  <.0001

#Late
tvc.late.sock<-gamm(ln.rs~s(spawners,k=4)+s(stock.dist,k=5) +s(stock.dist,by=loc.sst2,k=5)+ s(stock,bs='re'),data=sock[sock$era==2,],correlation=corAR1(form= ~ 1 | stock))  
summary(tvc.late.sock$gam)
#R-sq.(adj) =  0.281   
# Scale est. = 0.63239   n = 640
AIC(tvc.late.sock$lme)
#1531.496
plot(tvc.late.sock $gam,pages=1,scale=0)

sock$loc.sst2<-sock$sst1+1+abs(min(sock$sst1))
gam.late.sock<-gamm(ln.rs~s(spawners,k=4)+ s(stock.dist,k=5) +loc.sst2+ s(stock,bs='re'),data=sock[sock$era=='2',],correlation=corAR1(form= ~ 1 | stock))  
summary(gam.late.sock$gam)
#R-sq.(adj) =  0.279   
#  Scale est. = 0.6358    n = 640
AIC(gam.late.sock$lme)
#1531.55
plot(gam.late.sock$gam,pages=1,scale=0)

anova(tvc.late.sock$lme,gam.late.sock$lme)
#tvc.late.sock$lme     1 11 1531.496 1580.572 -754.7480                        
#gam.late.sock$lme     2  9 1531.550 1571.703 -756.7748 1 vs 2 4.053682  0.1318


#Plotting results on a map
covariate<-sock$loc.sst2[sock$era==1]
pred<-predict(tvc.early.sock$gam,type='terms')  [,3]
pred.se<-predict(tvc.early.sock$gam,type='terms',se.fit=T)[[2]][,3]
slope.early<-pred/covariate
covariate<-sock$loc.sst2[sock$era==2]
pred<-predict(tvc.late.sock$gam,type='terms')  [,3]
pred.se<-predict(tvc.late.sock$gam,type='terms',se.fit=T)[[2]][,3]
slope.late<-pred/covariate
range(slope.early)
range(slope.late)
max.slope<-max(abs(c(slope.early,slope.late)))

dev.new(width=9,height=7)
par(mfrow=c(2,2),mai=c(0.7,0.7,0.4,0.2))
plot(1,1,,ylim=range(pink$lat),xlim=range(-pink$long),main='Sockeye - early',xlab="",ylab="")
map('world',fill=T,col='grey',add=T)
symbols(sock$long2[sock$era==1][slope.early>0],sock$lat[sock$era==1][slope.early>0],circles=slope.early[slope.early>0],inches=0.1*max(abs(slope.early[slope.early>0]))/max.slope,bg='red',add=T)
symbols(sock$long2[sock$era==1][slope.early<0],sock$lat[sock$era==1][slope.early<0],circles=abs(slope.early)[slope.early<0],inches=0.1*max(abs(slope.early[slope.early<0]))/max.slope,bg='blue',add=T)
legend.bubble<-c(seq(0.4,0.1,by=-0.1))
symbols(rep(-160,4),c(52,50.8,49.8,48.9),circle=legend.bubble,add=T,inches=0.1*max(legend.bubble)/max.slope,bg='grey')
text(rep(-157,4),c(52,50.8,49.8,48.9),labels=as.character(legend.bubble),cex=1)
lines(sort(lon.sock),lat.sock.prj[order(lon.sock)],lwd=2,lty=2)
points(lon.sock,lat.sock.prj,pch='+',col='blue')


plot(1,1,,ylim=range(pink$lat),xlim=range(-pink$long),main='Sockeye - late',xlab="",ylab="")
map('world',fill=T,col='grey',add=T)
symbols(sock$long2[sock$era==2][slope.late>0],sock$lat[sock$era==2][slope.late>0],circles=slope.late[slope.late>0],inches=0.1*max(abs(slope.late[slope.late>0]))/max.slope,bg='red',add=T)
symbols(sock$long2[sock$era==2][slope.late<0],sock$lat[sock$era==2][slope.late<0],circles=abs(slope.late)[slope.late<0],inches=0.1*max(abs(slope.late[slope.late<0]))/max.slope,bg='blue',add=T)

#Slope of sst as a function of distance
plot(tvc.early.sock$gam,select=3,xlab='Distance (km)',ylab='Slope of SST',ylim=c(-0.5,0.5),lwd=2)
abline(h=0,lty=2)
plot(tvc.late.sock$gam,select=3,xlab='Distance (km)',ylab='Slope of SST',scale=0,ylim=c(-0.5,0.5),lwd=2)
abline(h=0,lty=2)
dev.copy(jpeg,'Sockeye.jpg',width=9,height=7,units='in',res=100)
dev.off()


#Pink VC GAM, spatial 
#Early period
pink$loc.sst2<-pink$sst1+1+abs(min(pink$sst1))
pink$long2<--pink$long
#Early
tvc.early.pink<-gamm(ln.rs~s(spawners,k=4)+ s(stock.dist,k=5)+s(stock.dist,by=loc.sst2,k=5)+ s(stock,bs='re'),data=pink[pink$era==1,],correlation=corAR1(form= ~ 1 | stock))  
summary(tvc.early.pink$gam)
#R-sq.(adj) =  0.187   
#Scale est. = 1.1089    n = 455
AIC(tvc.early.pink$lme)
#1406.374
plot(tvc.early.pink$gam,pages=1,scale=0)

pink$loc.sst2<-pink$sst1+1+abs(min(pink$sst1))
gam.early.pink<-gamm(ln.rs~s(spawners,k=4)+s(stock.dist,k=5) +loc.sst2+ s(stock,bs='re'),data=pink[pink$era==1,],correlation=corAR1(form= ~ 1 | stock))  
summary(gam.early.pink$gam)
#R-sq.(adj) =   0.17   
#Scale est. = 1.1348    n = 455
AIC(gam.early.pink$lme)
#1408.801
plot(gam.early.pink$gam,pages=1,scale=0)

anova(tvc.early.pink$lme, gam.early.pink$lme)
#                   Model df      AIC      BIC    logLik   Test  L.Ratio p-value
#tvc.early.pink$lme     1 11 1406.374 1451.697 -692.1870                       
#gam.early.pink$lme     2  9 1408.801 1445.884 -695.4006 1 vs 2 6.42728  0.0402

#Late
tvc.late.pink<-gamm(ln.rs~s(spawners,k=4)+ s(stock.dist,k=4)+s(stock.dist,by=loc.sst2,k=5)+ s(stock,bs='re'),data=pink[pink$era==2,],correlation=corAR1(form= ~ 1 | stock),)  
summary(tvc.late.pink$gam)
#R-sq.(adj) =  0.0499   
#Scale est. = 1.136     n = 373
AIC(tvc.late.pink $lme)
#1142.632
plot(tvc.late.pink $gam,pages=1,scale=0)

pink$loc.sst2<-pink$loc.sst+1+abs(min(pink$loc.sst))
gam.late.pink<-gamm(ln.rs~s(spawners,k=4)+s(stock.dist,k=5) +loc.sst2+ s(stock,bs='re'),data=pink[pink$era==2,],correlation=corAR1(form= ~ 1 | stock))  
summary(gam.late.pink$gam)
#R-sq.(adj) =  0.0501   
 # Scale est. = 1.1399    n = 373
AIC(gam.late.pink$lme)
#1140.028
plot(gam.late.pink$gam,pages=1,scale=0)

anova(tvc.late.pink$lme, gam.late.pink$lme)
#                  Model df      AIC      BIC    logLik   Test  L.Ratio p-value
#tvc.late.pink$lme     1 11 1142.631 1185.769 -560.3158                        
#gam.late.pink$lme     2  9 1140.028 1175.322 -561.0142 1 vs 2 1.396791  0.4974

#Plotting results on a map
covariate<-pink$loc.sst2[pink$era==1]
pred<-predict(tvc.early.pink$gam,type='terms')  [,3]
pred.se<-predict(tvc.early.pink$gam,type='terms',se.fit=T)[[2]][,3]
slope.early<-pred/covariate
covariate<-pink$loc.sst2[pink$era==2]
pred<-predict(tvc.late.pink$gam,type='terms')  [,3]
pred.se<-predict(tvc.late.pink$gam,type='terms',se.fit=T)[[2]][,3]
slope.late<-pred/covariate
range(slope.early)
range(slope.late)
max.slope<-max(abs(c(slope.early,slope.late)))

dev.new(width=9,height=7)
par(mfrow=c(2,2),mai=c(0.7,0.7,0.4,0.2))
plot(1,1,,ylim=range(pink$lat),xlim=range(-pink$long),main='Pink - early',xlab="",ylab="")
map('world',fill=T,col='grey',add=T)
symbols(pink$long2[pink$era==1][slope.early>0],pink$lat[pink$era==1][slope.early>0],circles=slope.early[slope.early>0],inches=0.1*max(abs(slope.early[slope.early>0]))/max.slope,bg='red',add=T)
symbols(pink$long2[pink$era==1][slope.early<0],pink$lat[pink$era==1][slope.early<0],circles=abs(slope.early)[slope.early<0],inches=0.1*max(abs(slope.early[slope.early<0]))/max.slope,bg='blue',add=T)
legend.bubble<-c(seq(0.5,0.2,by=-0.1))
symbols(rep(-160,4),c(52,50.8,49.8,48.9),circle=legend.bubble,add=T,inches=0.1*max(legend.bubble)/max.slope,bg='grey')
text(rep(-157,4),c(52,50.8,49.8,48.9),labels=as.character(legend.bubble),cex=1)
lines(sort(lon.pink),lat.pink.prj[order(lon.pink)],lwd=2,lty=2)
points(lon.pink,lat.pink.prj,pch='+',col='blue')

plot(1,1,,ylim=range(pink$lat),xlim=range(-pink$long),main='Pink - late',xlab="",ylab="")
map('world',fill=T,col='grey',add=T)
symbols(pink$long2[pink$era==2][slope.late>0],pink$lat[pink$era==2][slope.late>0],circles=slope.late[slope.late>0],inches=0.1*max(abs(slope.late[slope.late>0]))/max.slope,bg='red',add=T)
symbols(pink$long2[pink$era==2][slope.late<0],pink$lat[pink$era==2][slope.late<0],circles=abs(slope.late)[slope.late<0],inches=0.1*max(abs(slope.late[slope.late<0]))/max.slope,bg='blue',add=T)

#Slope of sst as a function of distance
plot(tvc.early.pink$gam,select=3,xlab='Distance (km)',ylab='Slope of SST',ylim=c(-0.5,1),lwd=2)
abline(h=0,lty=2)
plot(tvc.late.pink$gam,select=3,xlab='Distance (km)',ylab='Slope of SST',scale=0,ylim=c(-0.5,1),lwd=2)
abline(h=0,lty=2)
dev.copy(jpeg,'Pink.jpg',width=9,height=7,units='in',res=100)
dev.off()


##############################CHUM########################################
#VC GAM, spatial 

#Early period
chum$loc.sst2<-chum$sst1+1+abs(min(chum$sst1))
chum$long2<--chum$long
#Early
tvc.early.chum<-gamm(ln.rs~s(spawners,k=4)+ s(stock.dist,k=5)+s(stock.dist,by=loc.sst2,k=5)+ s(stock,bs='re'),data=chum[chum$era=='1',],correlation=corAR1(form= ~ 1 | stock))  
summary(tvc.early.chum$gam)
#R-sq.(adj) =  0.271   
#Scale est. = 0.61769   n = 457
AIC(tvc.early.chum$lme)
#1102.241
plot(tvc.early.chum$gam,pages=1,scale=0)

chum$loc.sst2<-chum$sst1+1+abs(min(chum$sst1))
gam.early.chum<-gamm(ln.rs~s(spawners,k=4)+s(stock.dist,k=5) +loc.sst2+ s(stock,bs='re'),data=chum[chum$era=='1',],correlation=corAR1(form= ~ 1 | stock))  
summary(gam.early.chum$gam)
#R-sq.(adj) =  0.274   
#Scale est. = 0.616     n = 457
AIC(gam.early.chum$lme)
#1100.351
plot(gam.early.chum$gam,pages=1,scale=0)

anova(tvc.early.chum$lme,gam.early.chum$lme)
#tvc.early.chum$lme     1 11 1102.241 1147.613 -540.1205                        
#gam.early.chum$lme     2  9 1100.351 1137.473 -541.1756 1 vs 2 2.110059  0.3482

#Late
tvc.late.chum<-gamm(ln.rs~s(spawners,k=4)+ s(stock.dist,k=5)+s(stock.dist,by=loc.sst2,k=5)+ s(stock,bs='re'),data=chum[chum$era=='2',],correlation=corAR1(form= ~ 1 |stock))  
summary(tvc.late.chum $gam)
#R-sq.(adj) =  0.423   
#Scale est. = 0.46621   n = 398
AIC(tvc.late.chum $lme)
#890.2313
plot(tvc.late.chum $gam,pages=1,scale=0)

chum$loc.sst2<-chum$sst1+1+abs(min(chum$sst1))
gam.late.chum<-gamm(ln.rs~s(spawners,k=4)+s(stock.dist,k=5) +loc.sst2+ s(stock,bs='re'),data=chum[chum$era=='2',],correlation=corAR1(form= ~ 1 |stock))  
summary(gam.late.chum$gam)
#R-sq.(adj) =  0.406   
 #Scale est. = 0.48095   n = 398
AIC(gam.late.chum$lme)
#895.2719
plot(gam.late.chum$gam,pages=1,scale=0)

anova(tvc.late.chum$lme,gam.late.chum$lme)
#tvc.late.chum$lme     1 11 890.2311 934.0821 -434.1156                        
#gam.late.chum$lme     2  9 895.2719 931.1500 -438.6359 1 vs 2 9.040771  0.0109

#Plotting results on a map
covariate<-chum$loc.sst2[chum$era=='1']
pred<-predict(tvc.early.chum$gam,type='terms')  [,3]
pred.se<-predict(tvc.early.chum$gam,type='terms',se.fit=T)[[2]][,3]
slope.early<-pred/covariate
covariate<-chum$loc.sst2[chum$era=='2']
pred<-predict(tvc.late.chum$gam,type='terms')  [,3]
pred.se<-predict(tvc.late.chum$gam,type='terms',se.fit=T)[[2]][,3]
slope.late<-pred/covariate
range(slope.early)
range(slope.late)
max.slope<-max(abs(c(slope.early,slope.late)))

dev.new(width=9,height=7)
par(mfrow=c(2,2),mai=c(0.7,0.7,0.4,0.2))
plot(1,1,,ylim=range(chum$lat),xlim=range(-chum$long),main='Chum - early',xlab="",ylab="")
map('world',fill=T,col='grey',add=T)
symbols(chum$long2[chum$era=='1'][slope.early>0],chum$lat[chum$era=='1'][slope.early>0],circles=slope.early[slope.early>0],inches=0.1*max(abs(slope.early[slope.early>0]))/max.slope,bg='red',add=T)
symbols(chum$long2[chum$era=='1'][slope.early<0],chum$lat[chum$era=='1'][slope.early<0],circles=abs(slope.early)[slope.early<0],inches=0.1*max(abs(slope.early[slope.early<0]))/max.slope,bg='blue',add=T)
legend.bubble<-c(seq(0.4,0.1,by=-0.1))
symbols(rep(-160,4),c(49.8,48.9,48.0,47.2),circle=legend.bubble,add=T,inches=0.1*max(legend.bubble)/max.slope,bg='grey')
text(rep(-157,4),c(49.8,48.9,48.0,47.2),labels=as.character(legend.bubble),cex=1)
lines(sort(lon.chum),lat.chum.prj[order(lon.chum)],lwd=2,lty=2)
points(lon.chum,lat.chum.prj,pch='+',col='blue')

plot(1,1,,ylim=range(chum$lat),xlim=range(-chum$long),main='Chum - late',xlab="",ylab="")
map('world',fill=T,col='grey',add=T)
symbols(chum$long2[chum$era=='2'][slope.late>0],chum$lat[chum$era=='2'][slope.late>0],circles=slope.late[slope.late>0],inches=0.1*max(abs(slope.late[slope.late>0]))/max.slope,bg='red',add=T)
symbols(chum$long2[chum$era=='2'][slope.late<0],chum$lat[chum$era=='2'][slope.late<0],circles=abs(slope.late)[slope.late<0],inches=0.1*max(abs(slope.late[slope.late<0]))/max.slope,bg='blue',add=T)

#Slope of sst as a function of distance
plot(tvc.early.chum$gam,select=3,xlab='Distance (km)',ylab='Slope of SST',ylim=c(-0.6,0.5),lwd=2)
abline(h=0,lty=2)
plot(tvc.late.chum$gam,select=3,xlab='Distance (km)',ylab='Slope of SST',scale=0,ylim=c(-0.6,0.5),lwd=2)
abline(h=0,lty=2)
dev.copy(jpeg,'Chum.jpg',width=9,height=7,units='in',res=100)
dev.off()

###Figure for CJFAS paper
dev.new(width=4,height=8)
par(mfrow=c(3,1))
#pink
plot(tvc.early.pink$gam,select=3,xlab='Distance (km)',ylab='Slope of SST',ylim=c(-0.5,1),lwd=2,shade=T,shade.col=alpha('orange',0.5),scale=0,main='Pink')
par(new=T)
plot(tvc.late.pink$gam,select=3,scale=0,ylim=c(-0.5,1),lwd=2,shade=T,shade.col=alpha('blue',0.5),xlab='',ylab='')
abline(h=0,lty=2)
#sockeye
plot(tvc.early.sock$gam,select=3,xlab='Distance (km)',ylab='Slope of SST',ylim=c(-0.5,0.5),lwd=2,main='Sockeye',scale=0,shade.col=alpha('orange',0.5),shade=T)
par(new=T)
plot(tvc.late.sock$gam,select=3,xlab='',ylab='',scale=0,ylim=c(-0.5,0.5),lwd=2,shade.col=alpha('blue',0.5),shade=T)
abline(h=0,lty=2)
#chum
plot(tvc.early.chum$gam,select=3,xlab='Distance (km)',ylab='Slope of SST',ylim=c(-0.5,0.5),lwd=2,shade=T,shade.col=alpha('orange',0.5),scale=0,main='Chum')
par(new=T)
plot(tvc.late.chum$gam,select=3,xlab='',ylab='',scale=0,ylim=c(-0.5,0.5),lwd=2,shade.col=alpha('blue',0.5),shade=T)
abline(h=0,lty=2)
dev.copy(jpeg,'Fig4.jpg',width=4,height=8,units='in',res=100)
dev.off()

