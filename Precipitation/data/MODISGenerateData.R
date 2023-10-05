rm(list=ls())

# Load and Process Data
dat<-read.csv("data.csv", header = F)
dat<-dat/1000 # Scale values 
z<-as.numeric(as.matrix(dat))

# Load and process locations
lon<-read.csv("longitude.csv", header = F)
lat<-read.csv("latitude.csv", header = F)
newLon<-seq(min(lon), max(lon), length.out=2030)
newLat<-seq(min(lat), max(lat), length.out=1354)
locs<-expand.grid(newLon, newLat)
# Remove Missing Values
rmInd<-which(z==-9.999)
z<-z[-rmInd]
locs<-locs[-rmInd,]

save(z,locs, file="../samples/water_vapor_20230622.RData")

# Function to create Plots
plotColorRF<-function(dat,rangeDat=dat,label="Plot",location,length.out=10,pch=16,cex=1, cex.main=1.5){
  breaks <- seq(range(rangeDat,na.rm = TRUE)[1],
                range(rangeDat,na.rm = TRUE)[2],
                length.out=length.out)
  pal <- tim.colors(length(breaks)-1,alpha = 1)
  fb <- classIntervals(dat, n = length(pal),
                       style = "fixed", fixedBreaks = breaks)
  col <- findColours(fb, pal)
  plot(x=location[,1],y=location[,2],col=col, pch=pch,cex=cex,cex.main=cex.main,
       main=label)
  return(pal)
}



# Sample 150k locations for plotting purposes
sampInd<-sample(1:length(z),150000)
useObs<-z[sampInd]
useGridLocation<-locs[sampInd,]

library(fields); library(classInt)
png(filename = "TotalPrecipitableWater.png", width = 1000, height = 600)
# Plot Observations and Probability Surface for Validation Locations (Figure 5)
layout(matrix(1:2,ncol=2,nrow=1), width = c(2,0.6),height = c(1,1))
par(mar=c(2.2,2,2,2))
# Observations
colA<-plotColorRF(dat=useObs , rangeDat = useObs,
                  location = useGridLocation, 
                  label = "Observations" , length.out = 1000,cex=0.5 , cex.main=1.5)

# Predicted Probability Surface
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Scale', cex.main=1.5)
legend_image <- as.raster(matrix(rev(colA), ncol=1))
text(x=1.6, y = seq(0,1,l=5), 
     labels = round(seq(min(useObs),max(useObs),l=5),2), cex = 1.5)
rasterImage(legend_image, -0.5, 0, 0.5,1)
dev.off()


