
# Load in Data and Libraries ----------------------------------------------

library(raster)
library(pracma)
library(RColorBrewer)
library(maptools)
library(rasterVis)
library(reshape)
library(ggplot2)
library(gridExtra)
library(graticule)
library(lineup)
# library(paran) # used to determine number of pc's to use

# save default plot margins to reset later
old_mar <- par()$mar

# import, mask and trim melt onset data:
melt_stack <- readRDS('data/Melt.RDS')
masks <- raster('data/SSMI_Regions2.tif')
masks[masks !=10] <- NA
laptev_outline <- rasterToPolygons(masks,dissolve=T)
laptev_masked <- mask(melt_stack,masks)
laptev <- trim(laptev_masked)

# get values from raster stack, find location with no missing values,
# detrend and perform PCA
values <- getValues(laptev)
values <- values[,1:37]  # for full melt years
val <- t(values)
indexicemelt <- complete.cases(t(val))
icemeltval <- val[,complete.cases(t(val))]
rm(values)

y <- detrend(icemeltval)
y_trend <- icemeltval-y

pcay <- prcomp(y, scale = TRUE)

# import, mask and trim melt onset data:
snow_stack <- readRDS('data/snowretreat.RDS')
masks <- raster('data/EASE2_N0_SnowRegions2.tif')
masks[masks !=24] <- NA
ob_outline <- rasterToPolygons(masks,dissolve=T)
ob_masked <- mask(snow_stack,masks)
ob <- trim(ob_masked)

ob[ob<=0] <- NA

# get values from raster stack, find location with no missing values,
# detrend and perform PCA
values <- getValues(ob)
val <- t(values)
indexsnowmelt <- complete.cases(t(val))
snowmeltval <- val[,complete.cases(t(val))]
rm(values)

x <- detrend(snowmeltval)
x_trend <- snowmeltval-x

pcax <- prcomp(x, scale = TRUE)

namex <- 'Snow Melt'
namey <- 'Ice Melt'

# Snow PoV
PoVx <- pcax$sdev^2/sum(pcax$sdev^2)
plot(1:37,PoVx[1:37],'b',
     main='Proportion ov Variance Captures by PCs\nSnow Melt',
     xlab="Modes", ylab="Frac. Var. explained")
grid()

# Ice PoV
PoVy <- pcay$sdev^2/sum(pcay$sdev^2)
plot(1:37,PoVy[1:37],'b',
     main='Proportion ov Variance Captures by PCs\nIce Melt',
     xlab="Modes", ylab="Frac. Var. explained")
grid()

plot(1979:2015, pcax$x[,1],type="b",main='PC1 Trend with Time, Ice Melt',xlab="Year",ylab="PC1")
plot(1979:2015, pcax$x[,2],type="b",main='PC2 Trend with Time, Ice Melt',xlab="Year",ylab="PC2")
plot(1979:2015, pcax$x[,3],type="b",main='PC3 Trend with Time, Ice Melt',xlab="Year",ylab="PC3")
plot(1979:2015, pcax$x[,4],type="b",main='PC4 Trend with Time, Ice Melt',xlab="Year",ylab="PC4")


# HPA
# paran(y, iterations=100, all=FALSE, graph=TRUE)

# paran(x, iterations=100, all=FALSE, graph=TRUE)

# Number of PCs to use determined from above HPA
xnpc <- 4
ynpc <- 4

# load in other atmospheric variables for later examination
sat <- readRDS('data/SAT_global.RDS') #surface air temperature
slp <- readRDS('data/SLP_global.RDS') #sea level pressure
gph <- readRDS('data/GPH500_global.RDS') #geopotential height at 500hpa
qv <- readRDS('data/Qv10_global.RDS') #specific humidity



# CCA Maps ---------------------------------------------------------------------

# Perform canonical correlation analysis on leading PCs
cca <- cancor(pcax$x[,1:xnpc],
              pcay$x[,1:ynpc])

# get CCA time series for x and y
x_canvar <- pcax$x[,1:xnpc] %*% cca$xcoef
y_canvar <- pcay$x[,1:ynpc] %*% cca$ycoef


i <- 1
# Calculate correlation between CCA x time-series and x values
xcca1corr <- cor(x,x_canvar[,i])

xraster <- raster(ob,layer=1)
new_vals <- rep(NA,ncell(xraster))
new_vals[indexsnowmelt] = xcca1corr
xcca_raster <- setValues(xraster, new_vals)
xcca_raster <- projectRaster(xcca_raster,crs='+proj=laea +lat_0=90 +lon_0=80')
extent(xcca_raster) <- extent(c(-1735277,700963.4,-4243553,-1478613 ))

data(wrld_simpl)
out <- crop(wrld_simpl, extent(-180, 180, 40, 83.57027))
wrld_snow <- spTransform(out,CRS('+proj=laea +lat_0=90 +lon_0=80'))

mapTheme <- rasterTheme(region=(brewer.pal(11,'RdBu')))

p1 <- levelplot(xcca_raster,margin=F,par.settings=mapTheme,
                main=paste0(namex,' CCA ',i),
                at=seq(-1,1, length.out=102),
                colorkey=T,scales=list(draw=F),
                xlim=c(-2400000,1400000),ylim=c(-4900000,-1000000)) +
  layer(sp.lines(wrld_snow, col="black",lwd=0.5,fill='grey90'),under=T) +
  layer(sp.lines(wrld_snow, col="black",lwd=0.5),under=F)

# Calculate correlation between CCA y time-series and y values
ycca1corr <- cor(y,y_canvar[,i])

yraster <- raster(laptev,layer=1)
new_vals <- rep(NA,ncell(yraster))
new_vals[indexicemelt] = ycca1corr
ycca_raster <- setValues(yraster, new_vals)
ycca_raster <- projectRaster(ycca_raster,crs='+proj=stere +lat_0=90 +lon_0=125')


out <- crop(wrld_simpl, extent(-180, 180, 50, 83.57027))
wrld_ice <- spTransform(out,CRS('+proj=stere +lat_0=90 +lon_0=125'))

p2 <- levelplot(ycca_raster,margin=F,par.settings=mapTheme,
                main=paste0(namey,' CCA ',i),
                at=seq(-1,1, length.out=102),
                colorkey=T,scales=list(draw=F)) +
  layer(sp.lines(wrld_ice, col="black", fill='grey90',lwd=0.5),under=F)


# Plot CCA mode 1 as time series
dates <- (2015-nrow(y)+1):2015

ts_cca1 <- as.data.frame(cbind(dates,x_canvar[,i],y_canvar[,i]))
colnames(ts_cca1) <- c('Year',namex,namey)
ts_cca1 <- melt(ts_cca1,id='Year')
library(ggplot2)
p3 <- ggplot(ts_cca1,aes(x=Year,y=value,color=variable,group=variable)) +
  geom_line() +
  labs(title=paste("CCA", i, "Variates; R =", round(cor(x_canvar[,i],y_canvar[,i]),2)),
       y='') +
  theme(legend.title=element_blank(),legend.justification=c(0.9,0.9),
        legend.position=c(0.98,0.98),legend.text=element_text(size = 8, face = "bold"),
        legend.margin = margin(0,0,0,0))
# ggplot2 masks a function in rasterVis, so need to detach it after use
detach('package:ggplot2',unload=T)

# Create figure of two maps and time series
lay <- rbind(c(1,2),c(3,3))
grid.arrange(p1, p2, p3, layout_matrix = lay,heights=c(3,2))



# Eigen Maps --------------------------------------------------------------

# Look at eigen maps to identify leading patterns of variability in variables
x_eig_list <- list()
x_lim <- abs(max(pcax$rotation[,1:4]))
mapTheme <- rasterTheme(region=(brewer.pal(11,'Spectral'))) 
for(i in 1:4){
  
  xraster <- raster(ob,layer=1)
  new_vals <- rep(NA,ncell(xraster))
  new_vals[indexsnowmelt] = pcax$rotation[,i]
  xcca_raster <- setValues(xraster, new_vals)
  xcca_raster <- projectRaster(xcca_raster,crs='+proj=laea +lat_0=90 +lon_0=80')
  extent(xcca_raster) <- extent(c(-1735277,700963.4,-4243553,-1478613 ))
  
  data(wrld_simpl)
  out <- crop(wrld_simpl, extent(-180, 180, 40, 83.57027))
  wrld_snow <- spTransform(out,CRS('+proj=laea +lat_0=90 +lon_0=80'))
  
  x_eig_list[[i]] <- levelplot(xcca_raster,margin=F,par.settings=mapTheme,
                  main=list(label=paste0('Eigenvector ',i),cex=0.8,font=1),
                  at=seq(-x_lim,x_lim, length.out=100),
                  colorkey=T,scales=list(draw=F),
                  xlim=c(-2400000,1400000),ylim=c(-4900000,-1000000)) +
    layer(sp.lines(wrld_snow, col="black",lwd=0.5,fill='grey90'),under=T) +
    layer(sp.lines(wrld_snow, col="black",lwd=0.5),under=F) 
  
}

lay <- rbind(c(1,1),c(2,3),c(4,5))
grid.arrange(textGrob('Snow Melt Eigen Decomposition Spatial Maps',gp=gpar(fontface='bold')),
             x_eig_list[[1]],x_eig_list[[2]],
             x_eig_list[[3]],x_eig_list[[4]], 
             layout_matrix = lay,heights=c(1,4,4))

y_eig_list <- list()
y_lim <- abs(max(pcay$rotation[,1:4]))
mapTheme <- rasterTheme(region=(brewer.pal(11,'Spectral'))) 
for(i in 1:4){
  
  yraster <- raster(laptev,layer=1)
  new_vals <- rep(NA,ncell(yraster))
  new_vals[indexicemelt] = pcay$rotation[,i]
  ycca_raster <- setValues(yraster, new_vals)
  ycca_raster <- projectRaster(ycca_raster,crs='+proj=stere +lat_0=90 +lon_0=125')
  # extent(xcca_raster) <- extent(c(-1735277,700963.4,-4243553,-1478613 ))

  y_eig_list[[i]] <- levelplot(ycca_raster,margin=F,par.settings=mapTheme,
                               main=list(label=paste0('Eigenvector ',i),cex=0.8,font=1),
                               at=seq(-y_lim,y_lim, length.out=100),
                               colorkey=T,scales=list(draw=F)) +
    layer(sp.lines(wrld_ice, col="black", fill='grey90',lwd=0.5),under=F) 
  
}

lay <- rbind(c(1,1),c(2,3),c(4,5))
grid.arrange(textGrob('Ice Melt Eigen Decomposition Spatial Maps',gp=gpar(fontface='bold')),
             y_eig_list[[1]],y_eig_list[[2]],
             y_eig_list[[3]],y_eig_list[[4]], 
             layout_matrix = lay,heights=c(1,4,4))


# Climate Maps ------------------------------------------------------------

# graticules
crs.longlat <- CRS("+init=epsg:4326")
prj <- CRS("+proj=stere +lat_0=90 +lon_0=100 +ellps=WGS84")
lons <- c(30,40,50,60,70,80,90,100,110,120,130,140,150,160,170)
lats <- c(40,50,60,70,80)

# set lat and lon limits
xl <-  range(lons) + c(-0.4, 0.4)
yl <- range(lats) + c(-0.4, 0.4)

grat <- graticule(lons, lats, proj = prj,
                  xlim = xl, ylim = yl)

grat <- graticule(lons, lats, proj = prj,
                  xlim = c(25, 175), ylim = c(35,85))

# Labels
labs <- graticule_labels(lons, lats,
                         xline = lons[2],
                         yline = lats[2],
                         proj = prj)
labsLon <- labs[labs$islon,]
labsLat <- labs[!labs$islon,]


# SAT:
satlatlon <- t(sat[1:2,])
satdata <- sat[c(-1,-2),]

# Extract values north of 40N
north_index <- which(satlatlon[,1]>40 & satlatlon[,2]>35 & satlatlon[,2]<170)
north_sat_normal <- satdata[,north_index]
north_sat <- detrend(north_sat_normal)
north_satlatlon <- satlatlon[north_index,]

# look at extreme years determined from CCA results
index <- c(17,18,19,20)
years <- c(1995,1996,1997,1998)

ylims <- max(abs(north_sat[17:20,]),na.rm=T)
y_ulim <- ylims
y_llim <- -ylims
rasmaps <- list()

for(i in 1:4){
  xraster <- as.data.frame(cbind(north_sat[index[i],],north_satlatlon))
  coordinates(xraster)=~lon+lat
  gridded(xraster) <- TRUE
  xraster <- raster(xraster)
  crs(xraster) <- crs(wrld_simpl)
  xraster <- projectRaster(xraster,
                           crs='+proj=stere +lat_0=90 +lon_0=100')
  out <- crop(wrld_simpl, extent(-180, 180, 30, 83.57027))
  wrld <- spTransform(out,crs(xraster))
  ob_l <- spTransform(ob_outline,crs(xraster))
  laptev_l <- spTransform(laptev_outline,crs(xraster))
  
  mapTheme <- rasterTheme(region=rev(brewer.pal(11,'Spectral'))) 
  rasmaps[[i]] <- levelplot(xraster,margin=F,par.settings=mapTheme,
            main=paste0('Detrended Surface Air Temperature, ',
                        years[i],'\nDegrees C'),
            at=seq(y_llim,y_ulim, length.out=100),
            colorkey=T,scales=list(draw=F)) +
    layer(sp.lines(wrld, col="black",lwd=0.5,fill='grey90'),under=T) +
    layer(sp.lines(wrld, col="black",lwd=0.5),under=F) +
    layer(sp.lines(ob_l, col="black",lwd=3),under=F) +
    layer(sp.lines(laptev_l, col="black",lwd=3),under=F) +
    layer(sp.lines(grat))  +
    layer(sp.text(coordinates(labsLon),
                  txt = parse(text = labsLon$lab),
                  adj = c(1.1, -0.25),
                  cex = 1)) +
    layer(sp.text(coordinates(labsLat),
                  txt = parse(text = labsLat$lab),
                  adj = c(-0.25, -0.25),
                  cex = 1))
}
rasmaps[[1]]
rasmaps[[2]]
rasmaps[[3]]
rasmaps[[4]]

# composite years
xraster <- as.data.frame(cbind(colMeans(rbind(north_sat[17,],
                                    north_sat[19,])),
                               north_satlatlon))
coordinates(xraster)=~lon+lat
gridded(xraster) <- TRUE
xraster <- raster(xraster)
crs(xraster) <- crs(wrld_simpl)
xraster <- projectRaster(xraster,
                         crs='+proj=stere +lat_0=90 +lon_0=100')

ylims <- max(abs(getValues(xraster)),na.rm=T)
y_ulim <- ylims
y_llim <- -ylims

mapTheme <- rasterTheme(region=rev(brewer.pal(11,'Spectral'))) 
levelplot(xraster,margin=F,par.settings=mapTheme,
          main='Detrended Surface Air Temperature, 1995-1997 Mean\nDegrees C',
          at=seq(y_llim,y_ulim, length.out=100),
          colorkey=T,scales=list(draw=F)) +
  layer(sp.lines(wrld, col="black",lwd=0.5,fill='grey90'),under=T) +
  layer(sp.lines(wrld, col="black",lwd=0.5),under=F) +
  layer(sp.lines(ob_l, col="black",lwd=3),under=F) +
  layer(sp.lines(laptev_l, col="black",lwd=3),under=F) +
  layer(sp.lines(grat))  +
  layer(sp.text(coordinates(labsLon),
                txt = parse(text = labsLon$lab),
                adj = c(1.1, -0.25),
                cex = 1)) +
  layer(sp.text(coordinates(labsLat),
                txt = parse(text = labsLat$lab),
                adj = c(-0.25, -0.25),
                cex = 1))


# SLP:
slplatlon <- t(slp[1:2,])
slpdata <- slp[c(-1,-2),]

# Extract values north of 40N
north_index <- which(slplatlon[,1]>40 & slplatlon[,2]>35 & slplatlon[,2]<170)
north_slp_normal <- slpdata[,north_index]
north_slp <- detrend(as.matrix(north_slp_normal))
north_slplatlon <- slplatlon[north_index,]

ylims <- max(abs(north_slp[17:20,]),na.rm=T)
y_ulim <- ylims
y_llim <- -ylims
rasmaps <- list()

for(i in 1:4){
  xraster <- as.data.frame(cbind(north_slp[index[i],],north_slplatlon))
  coordinates(xraster)=~lon+lat
  gridded(xraster) <- TRUE
  xraster <- raster(xraster)
  crs(xraster) <- crs(wrld_simpl)
  xraster <- projectRaster(xraster,
                           crs='+proj=stere +lat_0=90 +lon_0=100')
  out <- crop(wrld_simpl, extent(-180, 180, 30, 83.57027))
  wrld <- spTransform(out,crs(xraster))
  
  mapTheme <- rasterTheme(region=rev(brewer.pal(11,'Spectral'))) 
  rasmaps[[i]] <- levelplot(xraster,margin=F,par.settings=mapTheme,
            main=paste0('Detrended Sea Level Pressure, ',
                        years[i],'\nPascals'),
            at=seq(y_llim,y_ulim, length.out=100),
            colorkey=T,scales=list(draw=F)) +
    layer(sp.lines(wrld, col="black",lwd=0.5,fill='grey90'),under=T) +
    layer(sp.lines(wrld, col="black",lwd=0.5),under=F) +
    layer(sp.lines(ob_l, col="black",lwd=3),under=F) +
    layer(sp.lines(laptev_l, col="black",lwd=3),under=F)  +
    layer(sp.lines(grat))  +
    layer(sp.text(coordinates(labsLon),
                  txt = parse(text = labsLon$lab),
                  adj = c(1.1, -0.25),
                  cex = 1)) +
    layer(sp.text(coordinates(labsLat),
                  txt = parse(text = labsLat$lab),
                  adj = c(-0.25, -0.25),
                  cex = 1))
}
rasmaps[[1]]
rasmaps[[2]]
rasmaps[[3]]
rasmaps[[4]]

# composites
xraster <- as.data.frame(cbind(colMeans(rbind(north_slp[17,],
                                              north_slp[19,])),
                               north_satlatlon))
coordinates(xraster)=~lon+lat
gridded(xraster) <- TRUE
xraster <- raster(xraster)
crs(xraster) <- crs(wrld_simpl)
xraster <- projectRaster(xraster,
                         crs='+proj=stere +lat_0=90 +lon_0=100')

mapTheme <- rasterTheme(region=rev(brewer.pal(11,'Spectral'))) 
levelplot(xraster,margin=F,par.settings=mapTheme,
          main='Detrended Sea Level Pressure, 1995-1997 Mean\nPascals',
          at=seq(y_llim,y_ulim, length.out=100),
          colorkey=T,scales=list(draw=F)) +
  layer(sp.lines(wrld, col="black",lwd=0.5,fill='grey90'),under=T) +
  layer(sp.lines(wrld, col="black",lwd=0.5),under=F) +
  layer(sp.lines(ob_l, col="black",lwd=3),under=F) +
  layer(sp.lines(laptev_l, col="black",lwd=3),under=F) +
  layer(sp.lines(grat))  +
  layer(sp.text(coordinates(labsLon),
                txt = parse(text = labsLon$lab),
                adj = c(1.1, -0.25),
                cex = 1)) +
  layer(sp.text(coordinates(labsLat),
                txt = parse(text = labsLat$lab),
                adj = c(-0.25, -0.25),
                cex = 1))


# GPH500:
gphlatlon <- t(gph[1:2,])
gphdata <- gph[c(-1,-2),]

# Extract values north of 40N
north_index <- which(gphlatlon[,1]>40 & gphlatlon[,2]>35 & gphlatlon[,2]<170)
north_gph_normal <- gphdata[,north_index]
north_gph <- detrend(as.matrix(north_gph_normal))
north_gphlatlon <- gphlatlon[north_index,]

ylims <- max(abs(north_gph[17:20,]),na.rm=T)
y_ulim <- ylims
y_llim <- -ylims
rasmaps <- list()

for(i in 1:4){
  xraster <- as.data.frame(cbind(north_gph[index[i],],north_gphlatlon))
  coordinates(xraster)=~lon+lat
  gridded(xraster) <- TRUE
  xraster <- raster(xraster)
  crs(xraster) <- crs(wrld_simpl)
  xraster <- projectRaster(xraster,
                           crs='+proj=stere +lat_0=90 +lon_0=100')
  out <- crop(wrld_simpl, extent(-180, 180, 30, 83.57027))
  wrld <- spTransform(out,crs(xraster))

  mapTheme <- rasterTheme(region=rev(brewer.pal(11,'Spectral'))) 
  rasmaps[[i]] <- levelplot(xraster,margin=F,par.settings=mapTheme,
            main=paste0('Detrended Geopotential Height at 500mbar, ',
                        years[i],'\nMeters'),
            at=seq(y_llim,y_ulim, length.out=100),
            colorkey=T,scales=list(draw=F)) +
    layer(sp.lines(wrld, col="black",lwd=0.5,fill='grey90'),under=T) +
    layer(sp.lines(wrld, col="black",lwd=0.5),under=F) +
    layer(sp.lines(ob_l, col="black",lwd=3),under=F) +
    layer(sp.lines(laptev_l, col="black",lwd=3),under=F) +
    layer(sp.lines(grat))  +
    layer(sp.text(coordinates(labsLon),
                  txt = parse(text = labsLon$lab),
                  adj = c(1.1, -0.25),
                  cex = 1)) +
    layer(sp.text(coordinates(labsLat),
                  txt = parse(text = labsLat$lab),
                  adj = c(-0.25, -0.25),
                  cex = 1))
}
rasmaps[[1]]
rasmaps[[2]]
rasmaps[[3]]
rasmaps[[4]]

# composites
xraster <- as.data.frame(cbind(colMeans(rbind(north_gph[17,],
                                              north_gph[19,])),
                               north_satlatlon))
coordinates(xraster)=~lon+lat
gridded(xraster) <- TRUE
xraster <- raster(xraster)
crs(xraster) <- crs(wrld_simpl)
xraster <- projectRaster(xraster,
                         crs='+proj=stere +lat_0=90 +lon_0=100')

mapTheme <- rasterTheme(region=rev(brewer.pal(11,'Spectral'))) 
levelplot(xraster,margin=F,par.settings=mapTheme,
          main='Detrended Geopotential Height 500mbar, 1995-1997 Mean\nMeters',
          at=seq(y_llim,y_ulim, length.out=100),
          colorkey=T,scales=list(draw=F)) +
  layer(sp.lines(wrld, col="black",lwd=0.5,fill='grey90'),under=T) +
  layer(sp.lines(wrld, col="black",lwd=0.5),under=F) +
  layer(sp.lines(ob_l, col="black",lwd=3),under=F) +
  layer(sp.lines(laptev_l, col="black",lwd=3),under=F) +
  layer(sp.lines(grat))  +
  layer(sp.text(coordinates(labsLon),
                txt = parse(text = labsLon$lab),
                adj = c(1.1, -0.25),
                cex = 1)) +
  layer(sp.text(coordinates(labsLat),
                txt = parse(text = labsLat$lab),
                adj = c(-0.25, -0.25),
                cex = 1))

# Qv:
qvlatlon <- t(qv[1:2,])
qvdata <- qv[c(-1,-2),]

# Extract values north of 40N
north_index <- which(qvlatlon[,1]>40 & qvlatlon[,2]>35 & qvlatlon[,2]<170)
north_qv_normal <- qvdata[,north_index]
north_qv <- detrend(as.matrix(north_qv_normal))
north_qvlatlon <- qvlatlon[north_index,]

ylims <- max(abs(north_qv[17:20,]),na.rm=T)
y_ulim <- ylims
y_llim <- -ylims
rasmaps <- list()

for(i in 1:4){
  xraster <- as.data.frame(cbind(north_qv[index[i],],north_qvlatlon))
  coordinates(xraster)=~lon+lat
  gridded(xraster) <- TRUE
  xraster <- raster(xraster)
  crs(xraster) <- crs(wrld_simpl)
  xraster <- projectRaster(xraster,
                           crs='+proj=stere +lat_0=90 +lon_0=100')
  out <- crop(wrld_simpl, extent(-180, 180, 30, 83.57027))
  wrld <- spTransform(out,crs(xraster))

  mapTheme <- rasterTheme(region=(brewer.pal(11,'Spectral'))) 
  rasmaps[[i]] <- levelplot(xraster,margin=F,par.settings=mapTheme,
            main=paste0('Detrended Specific Humidity, ',
                        years[i],'\nkg H2O/kg air'),
            at=seq(y_llim,y_ulim, length.out=100),
            colorkey=T,scales=list(draw=F)) +
    layer(sp.lines(wrld, col="black",lwd=0.5,fill='grey90'),under=T) +
    layer(sp.lines(wrld, col="black",lwd=0.5),under=F) +
    layer(sp.lines(ob_l, col="black",lwd=3),under=F) +
    layer(sp.lines(laptev_l, col="black",lwd=3),under=F) +
    layer(sp.lines(grat))  +
    layer(sp.text(coordinates(labsLon),
                  txt = parse(text = labsLon$lab),
                  adj = c(1.1, -0.25),
                  cex = 1)) +
    layer(sp.text(coordinates(labsLat),
                  txt = parse(text = labsLat$lab),
                  adj = c(-0.25, -0.25),
                  cex = 1))
}
rasmaps[[1]]
rasmaps[[2]]
rasmaps[[3]]
rasmaps[[4]]

# composites
xraster <- as.data.frame(cbind(colMeans(rbind(north_qv[17,],
                                              north_qv[19,])),
                               north_satlatlon))
coordinates(xraster)=~lon+lat
gridded(xraster) <- TRUE
xraster <- raster(xraster)
crs(xraster) <- crs(wrld_simpl)
xraster <- projectRaster(xraster,
                         crs='+proj=stere +lat_0=90 +lon_0=100')

mapTheme <- rasterTheme(region=(brewer.pal(11,'Spectral'))) 
levelplot(xraster,margin=F,par.settings=mapTheme,
          main='Detrended Specific Humidity, 1995-1997 Mean\nkg H2O/kg air',
          at=seq(y_llim,y_ulim, length.out=100),
          colorkey=T,scales=list(draw=F)) +
  layer(sp.lines(wrld, col="black",lwd=0.5,fill='grey90'),under=T) +
  layer(sp.lines(wrld, col="black",lwd=0.5),under=F) +
  layer(sp.lines(ob_l, col="black",lwd=3),under=F) +
  layer(sp.lines(laptev_l, col="black",lwd=3),under=F) +
  layer(sp.lines(grat))  +
  layer(sp.text(coordinates(labsLon),
                txt = parse(text = labsLon$lab),
                adj = c(1.1, -0.25),
                cex = 1)) +
  layer(sp.text(coordinates(labsLat),
                txt = parse(text = labsLat$lab),
                adj = c(-0.25, -0.25),
                cex = 1))


# CCA Predictions ---------------------------------------------------------

xnpc <- 4
ynpc <- 4

pcaxall <- prcomp(x,scale.=T)$x
pcayall <- prcomp(y,scale.=T)$x

cca_crossval <- CCAForecast(x,y,xnpc,ynpc,pcaxall,pcayall)

cross_corr <- corbetw2mat(y,cca_crossval,what='paired')
cross_skill <- mean(cross_corr)
cca_pred <- CCAPred(y,x,xnpc,ynpc)
pred_corr <- corbetw2mat(y,cca_pred,what='paired')
pred_skill <- mean(pred_corr)

# cross map
yraster <- raster(laptev,layer=1)
new_vals <- rep(NA,ncell(yraster))
new_vals[indexicemelt] = cross_corr
cross_raster <- setValues(yraster, new_vals)
cross_raster <- projectRaster(cross_raster,crs='+proj=stere +lat_0=90 +lon_0=125')

mapTheme <- rasterTheme(region=(brewer.pal(11,'RdBu')))
p1 <- levelplot(cross_raster,margin=F,par.settings=mapTheme,
          main='Cross-Validation, Pred. Correlation w/ Obs.',
          at=seq(-1,1, length.out=100),
          colorkey=T,scales=list(draw=F))# +
  layer(sp.lines(wrld_ice, col="black", fill='grey90',lwd=0.5),under=F) 

# pred map
yraster <- raster(laptev,layer=1)
new_vals <- rep(NA,ncell(yraster))
new_vals[indexicemelt] = pred_corr
pred_raster <- setValues(yraster, new_vals)
pred_raster <- projectRaster(pred_raster,crs='+proj=stere +lat_0=90 +lon_0=125')

mapTheme <- rasterTheme(region=(brewer.pal(11,'RdBu')))
p2 <- levelplot(pred_raster,margin=F,par.settings=mapTheme,
                main='Fitted Model, Pred. Correlation w/ Obs.',
                at=seq(-1,1, length.out=100),
                colorkey=T,scales=list(draw=F)) +
  layer(sp.lines(wrld_ice, col="black", fill='grey90',lwd=0.5),under=F) 


sig_cross <- cross_corr
sig_cross[sig_cross <= 0.287] <- NA

percent_sig <- sum(!is.na(sig_cross))/length(sig_cross)*100

eig_raster <- raster(laptev,layer=1)
new_vals <- rep(NA,ncell(eig_raster))
new_vals[indexicemelt] = sig_cross
sig_raster <- setValues(eig_raster, new_vals)
sig_raster <- projectRaster(sig_raster,crs='+proj=stere +lat_0=90 +lon_0=125')


mapTheme <- rasterTheme(region=(brewer.pal(9,'Blues')))
p3 <- levelplot(sig_raster,margin=F,par.settings=mapTheme,
                main='Significant Correlations, 95% Confidence',
                at=seq(-1,1, length.out=100),
                colorkey=T,scales=list(draw=F)) +
  layer(sp.lines(wrld_ice, col="black", fill='grey90',lwd=0.5),under=F) 

lay <- c(1,2)
grid.arrange(p1,p3, layout_matrix = lay)
grid.arrange(p1,p3)

# 
# map showing the regions
xraster <- as.data.frame(cbind(north_gph[8,],north_gphlatlon))
coordinates(xraster)=~lon+lat
gridded(xraster) <- TRUE
xraster <- raster(xraster)
crs(xraster) <- crs(wrld_simpl)
xraster <- projectRaster(xraster,
                         crs='+proj=stere +lat_0=90 +lon_0=100')
out <- crop(wrld_simpl, extent(-180, 180, 30, 83.57027))
wrld <- spTransform(out,crs(xraster))


levelplot(xraster,margin=F,par.settings=mapTheme,
          main='Snow Melt in Western Siberia\nSea Ice in the Laptev Sea',
          at=seq(0,0.00001, length.out=100),
          colorkey=F,scales=list(draw=F)) +
  layer(sp.lines(wrld, col="black",lwd=0.5,fill='grey90'),under=T) +
  layer(sp.lines(wrld, col="black",lwd=0.5),under=F) +
  layer(sp.lines(ob_l, col="black",lwd=3),under=F) +
  layer(sp.lines(laptev_l, col="black",lwd=3),under=F) +
  layer(sp.lines(grat))  +
  layer(sp.text(coordinates(labsLon),
                txt = parse(text = labsLon$lab),
                adj = c(1.1, -0.25),
                cex = 1)) +
  layer(sp.text(coordinates(labsLat),
                txt = parse(text = labsLat$lab),
                adj = c(-0.25, -0.25),
                cex = 1))
