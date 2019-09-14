
# This script calculates era-specific (pre/post 1988/89) correlations between salmon productivity 
# and coastal SST at various lags, and also produces the study site figure for the CJFAS manuscript.

# load required packages
library(maps)
library(mapdata)
library(fields)
library(tidyr)
library(ncdf4)
library(chron)
library(ggplot2)
library(dplyr)
library(mapproj)
library(sinkr)

# to load sinkr, use package devtools and uncomment the following:
# devtools::install_github("marchtaylor/sinkr")
options(max.print = 99999)

# First, download ERSST v5 data and make study site figure. 
# Data are being downloaded on 1/11/19. 44º-68ºN, 186º-236ºE, 1/1950-11/2018.

# uncomment following line to download the SST data

# download.file("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1950-01-01):1:(2018-11-01T00:00:00Z)][(0.0):1:(0.0)][(44):1:(68)][(186):1:(236)]", "~all.sst")

nc <- nc_open("~all.sst")

# get lat/long
x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# assign dates
raw <- ncvar_get(nc, "time") # seconds since January 1, 1970
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

# save month and year for processing later
m <- months(d)
yr <- years(d)

# get required sst data
SST <- ncvar_get(nc, "sst")

# need to change to matrix for easier use
SST <- aperm(SST, 3:1) # transpose array

SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  # Change to matrix

# set names
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# define coastal cells - those touching land, excluding the Aleutians
# note that 56N200E is excluded b/c it straddles the GOA and EBS and so can't be assigned

coast <- c("N68E196", "N68E194", "N66E194", "N66E192", "N64E200", "N64E198", "N64E196",
           "N64E194", "N62E194", "N60E194", "N60E196", "N60E198", "N58E202", "N58E200",
           "N58E198", "N56E198", "N56E196",
           "N52E228", "N52E230", "N52E232", "N54E196",
           "N54E198", "N54E200","N54E226", "N54E228", "N54E230", "N56E202", "N56E204",
           "N56E206", "N56E224", "N56E226", "N56E228", "N58E204", "N58E206", "N58E208",
           "N58E222", "N58E224", "N58E226", "N60E208", "N60E210", "N60E212", "N60E214",
           "N60E216", "N60E218", "N60E220", "N50E232", "N50E234", "N50E236" ,"N48E236",
           "N48E234", "N46E236", "N44E236")

coast.sst <- SST
keep <- colnames(SST) %in% coast
coast.sst[,!keep] <- NA # blank out all non-coastal cells in the data

# now make the study site fig.
# load salmon data!
run.dat <- read.csv("coastwide salmon data.csv", row.names = 1)

pink.runs <- run.dat %>%
  filter(species=="Pink") %>%
  group_by(stock) %>%
  summarise(lat=mean(lat), long=mean(360-long))

sock.runs <- run.dat %>%
  filter(species=="Sockeye") %>%
  group_by(stock) %>%
  summarise(lat=mean(lat), long=mean(360-long))

chum.runs <- run.dat %>%
  filter(species=="Chum") %>%
  group_by(stock) %>%
  summarise(lat=mean(lat), long=mean(360-long))

# produce study site figure
png("study site fig.png", 8,5, units="in", res=300)
par(mfrow=c(2,2), mar=c(0.5,0.5,1.5,1), las=1)

# define projection details
my.proj <- "orthographic"
my.orien <- c(55,220,340)

map("world2Hires", 
    proj=my.proj, parameters = NULL, orient=my.orien,
    xlim=c(195,241), ylim=c(48,64.5),
    fill=FALSE, lforce="e")

map('world2Hires', 'usa', fill=T, add=T, 
    lwd=0.5, col="lightyellow3", proj=my.proj, parameters = NULL, orient=my.orien)
map('world2Hires', 'Canada',fill=T, add=T, 
    lwd=0.5, col="lightyellow3", proj=my.proj, parameters = NULL, orient=my.orien)

plot.points <- mapproject(pink.runs$long, pink.runs$lat,projection=my.proj, orientation=my.orien) 
points(plot.points$x, plot.points$y, pch=21, bg ="red", cex=1)
box()  

mtext("a) Pink", adj=0)

map("world2Hires", 
    proj=my.proj, parameters = NULL, orient=my.orien,
    xlim=c(195,240), ylim=c(49,60),
    fill=FALSE, lforce="e")

map('world2Hires', 'usa', fill=T, add=T, 
    lwd=0.5, col="lightyellow3", proj=my.proj, parameters = NULL, orient=my.orien)
map('world2Hires', 'Canada',fill=T, add=T, 
    lwd=0.5, col="lightyellow3", proj=my.proj, parameters = NULL, orient=my.orien)

plot.points <- mapproject(sock.runs$long, sock.runs$lat,projection=my.proj, orientation=my.orien) 
points(plot.points$x, plot.points$y, pch=21, bg ="red", cex=1)
box()  
mtext("b) Sockeye", adj=0)

map("world2Hires", 
    proj=my.proj, parameters = NULL, orient=my.orien,
    xlim=c(195,240), ylim=c(49,60),
    fill=FALSE, lforce="e")

map('world2Hires', 'usa', fill=T, add=T, 
    lwd=0.5, col="lightyellow3", proj=my.proj, parameters = NULL, orient=my.orien)
map('world2Hires', 'Canada',fill=T, add=T, 
    lwd=0.5, col="lightyellow3", proj=my.proj, parameters = NULL, orient=my.orien)

plot.points <- mapproject(chum.runs$long, chum.runs$lat,projection=my.proj, orientation=my.orien) 
points(plot.points$x, plot.points$y, pch=21, bg ="red", cex=1)
box()  
mtext("c) Chum", adj=0)

# Re-shape mean SST data  to a matrix with latitudes in columns, longitudes in rows
# get mean value for each cell
coast.mean <- colMeans(coast.sst)
# turn into matrix for plotting
z <- t(matrix(coast.mean,length(y)))  
# and change to a set of polygons for projection
polys <- matrixPoly(x, y, z)

map("world2Hires", 
    proj=my.proj, parameters = NULL, orient=my.orien,
    xlim=c(193,235.5), ylim=c(35,65),
    fill=FALSE, lforce="e")

COLS <- val2col(z, col = tim.colors(64))
for(i in seq(polys)){
  tmp <- mapproject(polys[[i]],
                    proj=my.proj, parameters = NULL, orient=my.orien)
  polygon(tmp$x, tmp$y, col=COLS[i], border=COLS[i], lwd=0.1)
}

# make grid for cells used
plot.lat <- lat[keep]
plot.long <- lon[keep]

for(i in 1:length(plot.lat)){
  xgrid <- c(plot.long[i]-1, plot.long[i]+1, plot.long[i]+1, plot.long[i]-1, plot.long[i]-1)
  ygrid <- c(plot.lat[i]+1, plot.lat[i]+1, plot.lat[i]-1, plot.lat[i]-1, plot.lat[i]+1)
  proj.lines <- mapproject(xgrid, ygrid, projection=my.proj, orientation=my.orien) 
  lines(proj.lines$x, proj.lines$y, lwd=0.5)
}

map('world2Hires', 'usa', fill=T, add=T, 
    lwd=0.5, col="lightyellow3", proj=my.proj, parameters = NULL, orient=my.orien)
map('world2Hires', 'Canada',fill=T, add=T, 
    lwd=0.5, col="lightyellow3", proj=my.proj, parameters = NULL, orient=my.orien)
box()

# add legend strip
mt.cex <- 1.1
l.mar <- 3
l.cex <- 0.8
l.l <- 1.2
tc.l <- -0.2
image.plot(z, legend.only=TRUE, horizontal =TRUE,  legend.lab = "Mean annual SST (ºC)", 
           smallplot = c(0.25,0.78,0.2,0.23), 
           legend.cex=0.8,
           legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0))) 
mtext("d) SST data", adj=0)

dev.off()

# Now process the SST data for analysis.

# drop all non-coastal cells from the data!
coast.sst <- coast.sst[,keep]

# and set a couple functions for standardizing
f1 <- function(x) tapply(x, m, mean)
f2 <- function(x) tapply(x, m, sd)

# get anomalies for each cell (difference from overall monthly mean, in units of standard deviation) 
mu <- apply(coast.sst, 2, f1)	# Compute monthly means for each cell
mu <- mu[rep(1:12, floor(length(d)/12)),]  # Replicate means matrix for each year at each location 

#now need to account for trailing months (i.e., fractions of a year that are available with updated data)
add <- length(d)-12*floor(length(d)/12)

# add in additional required months
mu <- rbind(mu, mu[1:add,])

std <- apply(coast.sst, 2, f2)	# Compute monthly sd for each cell

# and stack as for mu
std <- std[rep(1:12, floor(length(d)/12)),] 

# add in additional required months
std <- rbind(std, std[1:add,])

# now calculate anomalies...
coast.anom <- (coast.sst - mu)/std # difference from seasonal mean as SD!

# now assign each run to a set of coastal cells
# begin by defining an object with the lat and long for each cell!
coast.cells <- data.frame(cell=colnames(coast.sst), lat=lat[keep], long=lon[keep])

# add EBS/GOA factor!
# this is a convenience for identifying local cells without including cells on the other
# side of the Alaska Peninsula from the ocean entry point
# so in this case, "GOA" will include both GOA and South runs...
coast.cells$region <- "GOA"
change <- coast.cells$lat >= 57 & coast.cells$long <=202
coast.cells$region[change] <- "EBS"

change <- coast.cells$lat == 56 & coast.cells$long <=200
coast.cells$region[change] <- "EBS"

# now we need an object with the lat longs of each run!

# simplify spawner-recruit time series for each species

pink <- run.dat %>%
  filter(species=="Pink") %>%
  select(species, stock, brood.yr, spawners, recruits, lat, long, region)

sock <- run.dat %>%
  filter(species=="Sockeye") %>%
  select(species, stock, brood.yr, spawners, recruits, lat, long, region)

chum <- run.dat %>%
  filter(species=="Chum") %>%
  select(species, stock, brood.yr, spawners, recruits, lat, long, region)

pink.cells <- data.frame(species="Pink", stock=names(tapply(pink$lat, pink$stock, mean)), 
                         lat=tapply(pink$lat, pink$stock, mean), 
                         long=tapply(pink$long, pink$stock, mean))

sock.cells <- data.frame(species="Sockeye", stock=names(tapply(sock$lat, sock$stock, mean)), 
                         lat=tapply(sock$lat, sock$stock, mean), 
                         long=tapply(sock$long, sock$stock, mean))

chum.cells <- data.frame(species="Chum", stock=names(tapply(chum$lat, chum$stock, mean)), 
                         lat=tapply(chum$lat, chum$stock, mean), 
                         long=tapply(chum$long, chum$stock, mean))

run.cells <- rbind(pink.cells, sock.cells, chum.cells)
run.cells <- na.omit(run.cells)
run.cells$cells <- NA
rownames(run.cells) <- 1:nrow(run.cells)

# quickly add region
run.cells$region <- NA

for(i in 1:nrow(run.cells)){
  
  run.cells$region[i] <- unique(as.character(run.dat$region[run.dat$species==run.cells$species[i] & run.dat$stock==run.cells$stock[i]]))
  
}

# change coast.cells long back to degrees W!
coast.cells$long <- 360-coast.cells$long

# get cells <= 400 km away for ebs runs...
ebs.x1 <- as.matrix(cbind(run.cells$long[run.cells$region=="EBS"], 
                          run.cells$lat[run.cells$region=="EBS"]))
ebs.x2 <- as.matrix(cbind(coast.cells$long[coast.cells$region=="EBS"], 
                          coast.cells$lat[coast.cells$region=="EBS"]))

dist.ebs <- rdist.earth(ebs.x1, ebs.x2, miles=F)

colnames(dist.ebs) <- coast.cells$cell[coast.cells$region=="EBS"]
rownames(dist.ebs) <- run.cells$stock[run.cells$region=="EBS"]

for(i in 1:nrow(dist.ebs)){
  # i <- 1
  run <- rownames(dist.ebs)[i]
  use <- dist.ebs[i,] <= 400 
  cells <- as.vector(colnames(dist.ebs)[use])
  cc <- cells[1]
  
  for(j in 2:length(cells)){
    cc <- c(cc,cells[j])
    
  }
  
  run.cells$cells[run.cells$stock==run] <- paste(cells, collapse=",")
}

###
# get cells <= 400 km away for goa runs...
goa.x1 <- as.matrix(cbind(run.cells$long[run.cells$region %in% c("South", "GOA")], 
                          run.cells$lat[run.cells$region %in% c("South", "GOA")]))
goa.x2 <- as.matrix(cbind(coast.cells$long[coast.cells$region %in% c("South", "GOA")], 
                          coast.cells$lat[coast.cells$region %in% c("South", "GOA")]))

dist.goa <- rdist.earth(goa.x1, goa.x2, miles=F)

colnames(dist.goa) <- coast.cells$cell[coast.cells$region %in% c("South", "GOA")]
rownames(dist.goa) <- run.cells$stock[run.cells$region %in% c("South", "GOA")]

for(i in 1:nrow(dist.goa)){
  # i <- 1
  run <- rownames(dist.goa)[i]
  use <- dist.goa[i,] <= 400
  cells <- as.vector(colnames(dist.goa)[use])
  cc <- cells[1]
  
  for(j in 2:length(cells)){
    cc <- c(cc,cells[j])
  }
  run.cells$cells[run.cells$stock==run] <- paste(cells, collapse=",")
}

# here we produce a set of maps to check our selection of cells for each run

pdf("runs and local sst cells.pdf", 10,8)
par(mfrow=c(4,3), las=1, mar=c(2,2,2,0.5))

for(i in 1:nrow(run.cells)){
  
  sub <- run.cells[i,]
  sub.sst <- SST
  temp <- strsplit(sub$cells,",")
  use <- matrix(unlist(temp), ncol=1, byrow=TRUE)
  
  keep <- colnames(sub.sst) %in% use
  sub.sst[,!keep] <- NA
  
  lim <- range(colMeans(SST), na.rm = T)
  
  sub.mean <- colMeans(sub.sst)
  z <- t(matrix(sub.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
  image(x,y,z, col=tim.colors(64), zlim=lim, xlim=c(160,240), ylim=c(44,70), ylab="", xlab="", yaxt="n", xaxt="n")
  #contour(x, y, z, add=T)  
  map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="lightyellow3")
  map('world2Hires', 'Canada',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="lightyellow3")
  map('world2Hires', 'ussr',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="lightyellow3")
  points(360-sub$long, sub$lat, pch=21, bg="red", cex=1.5)
  mtext(paste(sub$species, sub$stock))
}
dev.off()

# Now the processing is done and we can begin the analysis. 
# First step is calculating era-specific corrrelations between productivity and SST at various lags for each run.

# create matrix of correlations for early and late period
run.corrs.early <- run.cells %>%
  select(species, stock, region)

run.corrs.early[,4:45] <- NA
names(run.corrs.early)[4:45] <- 1:42

run.corrs.late <- run.corrs.early

# set up an entry year column to use in partitioning by era
run.dat$entry.yr <- run.dat$brood.yr + 1
for(i in 1:nrow(run.dat)){
  if(run.dat$species[i] == "Sockeye") {run.dat$entry.yr[i]<-run.dat$brood.yr[i]+2} 
}

# fit a single SR model to the entire time series and get residuals
run.dat$res <- NA

for(i in 1:nrow(run.cells)){
  
  temp <- run.dat %>%
    filter(stock==run.cells$stock[i])
  
  # get residuals from Ricker model
  run.dat$res[run.dat$stock==run.cells$stock[i]] <- resid(lm(log(recruits/spawners) ~ spawners, data=temp)) 
}

# now calculate correlations with SST for each era
# first, the early period
for(i in 1:nrow(run.cells)){ # loop through each run
  
  temp <- run.dat %>%
    filter(stock==run.cells$stock[i], entry.yr <=1988) # select run of interest
  
  # and get the relevant SST data
  cells <- unlist(strsplit(run.cells$cell[i], ","))
  temp.sst <- coast.anom[,colnames(coast.anom) %in% cells]
  time0 <- d[m=="Jun" & yr %in% temp$brood.yr]
  
  count <- 1:nrow(temp.sst)
  
  count <- count[rownames(temp.sst) %in% as.character(time0)]
  
  for(j in 1:42){ # loop through each of the lags
    
    run.corrs.early[(i),(j+3)] <- cor(temp$res, rowMeans(temp.sst[count+j,]), use="p")
  }
}

# and now the late era
for(i in 1:nrow(run.cells)){ # loop through each run
  
  temp <- run.dat %>%
    filter(stock==run.cells$stock[i], entry.yr >1988) # select run of interest
  
  # and get the relevant SST data
  cells <- unlist(strsplit(run.cells$cell[i], ","))
  temp.sst <- coast.anom[,colnames(coast.anom) %in% cells]
  time0 <- d[m=="Jun" & yr %in% temp$brood.yr]
  
  count <- 1:nrow(temp.sst)
  
  count <- count[rownames(temp.sst) %in% as.character(time0)]
  
  for(j in 1:42){ # loop through each of the lags
    
    run.corrs.late[(i),(j+3)] <- cor(temp$res, rowMeans(temp.sst[count+j,]), use="p")
    
  }
}

# now get averages for each species in each era
pink.early <- pink.late <- sock.early <- sock.late <- chum.early <- chum.late <- matrix(nrow=3, ncol=42)

for(j in 1:42){
  
  pink.early[,j] <- tapply(run.corrs.early[run.corrs.early$species=="Pink",(j+3)],
                           run.corrs.early$region[run.corrs.early$species=="Pink"], mean, na.rm=T)
  pink.late[,j] <- tapply(run.corrs.late[run.corrs.late$species=="Pink",(j+3)],
                          run.corrs.late$region[run.corrs.late$species=="Pink"], mean, na.rm=T)
  
  chum.early[,j] <- tapply(run.corrs.early[run.corrs.early$species=="Chum",(j+3)],
                           run.corrs.early$region[run.corrs.early$species=="Chum"], mean, na.rm=T)
  chum.late[,j] <- tapply(run.corrs.late[run.corrs.late$species=="Chum",(j+3)],
                          run.corrs.late$region[run.corrs.late$species=="Chum"], mean, na.rm=T)
  
  sock.early[,j] <- tapply(run.corrs.early[run.corrs.early$species=="Sockeye",(j+3)],
                           run.corrs.early$region[run.corrs.early$species=="Sockeye"], mean, na.rm=T)
  sock.late[,j] <- tapply(run.corrs.late[run.corrs.late$species=="Sockeye",(j+3)],
                          run.corrs.late$region[run.corrs.late$species=="Sockeye"], mean, na.rm=T)
}

pink.early <- as.data.frame(pink.early)
pink.late <- as.data.frame(pink.late)

chum.early <- as.data.frame(chum.early)
chum.late <- as.data.frame(chum.late)

sock.early <- as.data.frame(sock.early)
sock.late <- as.data.frame(sock.late)

pink.early$region <- pink.late$region <- chum.early$region <- chum.late$region <- sock.early$region <- sock.late$region <- c("Bering Sea", "Gulf of Alaska", "South")

pink.early$era <- chum.early$era <- sock.early$era <- "Before 1988/89"
pink.late$era <- chum.late$era <- sock.late$era <- "After 1988/89"

pink.early$species <- pink.late$species <- "Pink"
sock.early$species <- sock.late$species <- "Sockeye"
chum.early$species <- chum.late$species <- "Chum"

plot.corrs <- rbind(pink.early, pink.late, chum.early, chum.late, sock.early, sock.late)

plot.corrs <- plot.corrs %>%
  gather(lag, correlation, -region, -era, -species) %>%
  mutate(lag.plot=rep(1:42, each=18))

# add a month column for plot labels
month.lab <- c("Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# change Bering Sea to Eastern Bering Sea fpor plotting
change <- grep("Bering", plot.corrs$region)
plot.corrs$region[change] <- "Eastern Bering Sea"

plot.corrs$plot.era <- as.factor(ifelse(plot.corrs$era=="Before 1988/89", 1,2))

plot.corrs$species <- as.factor(plot.corrs$species)

plot.corrs$sp.order <- ifelse(plot.corrs$species=="Pink", 1,
                              ifelse(plot.corrs$species=="Sockeye", 2, 3))

plot.corrs$species <- reorder(plot.corrs$species, plot.corrs$sp.order)

# make column for indicating general times of ocean entry
plot.corrs$entry <- NA
plot.val1 <- -0.5
plot.val2 <- -0.3

plot.corrs$entry[plot.corrs$species=="Pink" 
                 & plot.corrs$region %in% c("Eastern Bering Sea", "Gulf of Alaska") 
                 & plot.corrs$lag.plot %in% 10:15] <- plot.val2

plot.corrs$entry[plot.corrs$species=="Pink" 
                 & plot.corrs$region == "South" 
                 & plot.corrs$lag.plot %in% 9:13] <- plot.val1

plot.corrs$entry[plot.corrs$species=="Sockeye" 
                 & plot.corrs$region %in% c("Eastern Bering Sea", "Gulf of Alaska") 
                 & plot.corrs$lag.plot %in% c(23:27, 35:39)] <- plot.val2

plot.corrs$entry[plot.corrs$species=="Sockeye" 
                 & plot.corrs$region == "South" 
                 & plot.corrs$lag.plot %in% 22:26] <- plot.val1

plot.corrs$entry[plot.corrs$species=="Chum" 
                 & plot.corrs$region %in% c("Eastern Bering Sea", "Gulf of Alaska") 
                 & plot.corrs$lag.plot %in% 10:15] <- plot.val2

plot.corrs$entry[plot.corrs$species=="Chum" 
                 & plot.corrs$region == "South" 
                 & plot.corrs$lag.plot %in% 9:13] <- plot.val1

# reduce the range of plotted months for pink and chum

pink.reduced <- plot.corrs %>%
  filter(species=="Pink", lag.plot <= 23)

chum.reduced <- plot.corrs %>%
  filter(species=="Chum", lag.plot <= 30)

plot.corrs <- rbind(pink.reduced, chum.reduced,
                    filter(plot.corrs, species == "Sockeye"))

# and add vertical lines for delineating years
lines.plot <- data.frame(species=c(rep("Pink",6), rep("Sockeye", 9), rep("Chum",6)),
                         region=c(rep(c("Eastern Bering Sea", "Gulf of Alaska", "South"), each=2),
                                  rep(c("Eastern Bering Sea", "Gulf of Alaska", "South"), each=3),
                                  rep(c("Eastern Bering Sea", "Gulf of Alaska", "South"), each=2)),
                         int=c(rep(c(6.5,18.5), 3),
                               rep(c(6.5, 18.5, 30.5), 3),
                               rep(c(6.5,18.5), 3)))

# define colorblind palette 
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

correlation.plot <- 
  ggplot(plot.corrs, aes(lag.plot, correlation, color=plot.era)) + theme_linedraw()+ geom_line(lwd=0.5) + 
  geom_point() + facet_grid(region~species, scale="free", space="free_x") +
  geom_hline(yintercept = 0, color="dark grey") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5), legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size=10), axis.title.x = element_blank()) + ylab("Correlation coefficient") + 
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks=seq(1,42, by=3), labels = month.lab[seq(1,42, by=3)]) + 
  scale_color_manual(labels = c("Before 1988/89", "After 1988/89"), values = cb[c(2,6)]) +
  geom_line(aes(lag.plot,entry), size=2, color=cb[1]) +
  geom_vline(data=lines.plot, aes(xintercept = int), lty=2, color="dark grey") 


# also make a table of runs used in analysis
run.table <- run.dat %>%
  select(species, region, stock, brood.yr) %>%
  group_by(region, species, stock) %>%
  summarize(min(brood.yr), max(brood.yr))

write.csv(run.table, "run table.csv")
