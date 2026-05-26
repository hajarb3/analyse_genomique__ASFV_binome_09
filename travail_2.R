setwd("C:/Users/hajar/ULB/MA1/épidémio/Binome_09")

##-----------Data preparation and formatting for TempEST analysis----##
install.packages("readr")
library(readr) 

file <- read.csv("ASFV_genomic_analyses_simulated_dataset_1-3.csv")

file$sample <- paste0("H-", file$ID)
str(file)

sample_date<- file[,c("sample","collection_date")]
write_tsv(sample_date,"ASFV_genomic_date.tsv") #file in tsv format used in TempEST

##The following code is highly inspired from codes of the Professor Simon Dellicour.

##---------Spatial data preparation: reprojecting motorways shapefile-------####
library(raster)
library(sf)

forest_areas = raster("Raster_Wallonian_forest_areas.asc")
motorways = shapefile("Motorways_lines_in_study_area.shp")
motorways_reprojected = sp::spTransform(motorways, "EPSG:3035")
plot(forest_areas); lines(motorways_reprojected, col="red")
st_write(st_as_sf(motorways_reprojected), "Motorways_lines_re-projected.shp")

##-----------Estimating the evolution of the effective reproduction number----------##

if (!require(EpiEstim)) install.packages("EpiEstim")
if (!require(lubridate)) install.packages("lubridate")
library(EpiEstim)
library(lubridate)


tab = read.csv("ASFV_genomic_analyses_simulated_dataset_1-3.csv", header = TRUE)
tab = as.data.frame(tab)


days = interval(min(ymd(tab[,"collection_date"])),ymd(tab[,"collection_date"]))%/%days(1)+1
total_number_of_days = interval(min(ymd(tab[,"collection_date"])),max(ymd(tab[,"collection_date"])))%/%days(1)-1
daily_cases = rep(NA, total_number_of_days)
for (i in 1:length(daily_cases)) {
  daily_cases[i] = sum(days==i)
}

n = 1000 # number of iterations for which the SI mean will be drawn from a uniform distribution
mean_range = c(13, 18) # uniform distribution for the SI mean value, ranging from 9 to 23 days a mean SI value 
#ranging from 13 to 18 days and a standard deviation ranging from
#5 to 7 for the simulated ASFV datasets
sd_range = c(5, 7) # uniform distribution for the SI standard deviation, ranging from 4 to 8 days
t_start = seq(2, length(daily_cases)-6) # to set a sliding window of 7 days (1 week)
t_end = seq(8, length(daily_cases)) # to set a sliding window of 7 days (1 week)
all_Rt = matrix(NA, nrow=length(t_start), ncol=n) # matrix in which all estimates will be stored
for (i in 1:n) {
  mean_si_i = runif(1, mean_range[1], mean_range[2]) # mean SI value drawn for the iteration i
  sd_si_i = runif(1, sd_range[1], sd_range[2]) # standard deviation drawn for the iteration i
  res_i = estimate_R(incid=daily_cases, method="parametric_si", config=make_config(list(
    mean_si=mean_si_i, std_si=sd_si_i, t_start=t_start, t_end=t_end)))
  all_Rt[,i] = res_i$R$`Mean(R)` # all Rt estimates obtained for the iteration i
}
R_median = apply(all_Rt, 1, median, na.rm=T); Rt_days = (t_start+t_end)/2
R_dates = as.Date(min(ymd(tab[,"collection_date"]))) + Rt_days
R_lower = apply(all_Rt, 1, quantile, probs=0.025, na.rm=T)
R_upper = apply(all_Rt, 1, quantile, probs=0.975, na.rm=T)


par(oma=c(0,0,0,0), mar=c(2.0,3.5,0.1,0.5), lwd=0.3, col="gray30", col.axis="gray30", fg="gray30") 
#to set general graphical parameters
plot(R_dates, R_median, lwd=0.7, type="l", cex.axis=0.8, cex.lab=0.8, col="gray30", axes=F,
     xlab=NA, ylab=NA, ylim=c(0, max(R_upper, na.rm=T)), xlim=range(R_dates, na.rm=T)) 
#to plot a "blank" plot

xx_l = c(R_dates,rev(R_dates)) 
yy_l = c(R_lower,rev(R_upper)) #lower and upper bonds
polygon(xx_l,yy_l,col=rgb(187/255,187/255,187/255,0.5),border=0)  #plotted as a grey ribbon
lines(R_dates, R_median, lwd=0.3, type="l", cex.axis=0.8, cex.lab=0.8, col="gray30") # median

abline(h=1, lty=3, lwd=1, col= "red") 
#horizontal dashed line displaying the threshold of 1

axis.Date(side=1, x=R_dates, format="%Y-%m-%d", lwd=0.5, cex.axis=0.7, mgp=c(0,0.17,0), lwd.tick=0.5, col="gray30", col.axis="gray30", tck=-0.03, las=1) 
# to specify the position of the tick marks on the x-axis
axis(side=2, lwd=0.5, cex.axis=0.7, mgp=c(0,0.6,0), lwd.tick=0.5, col="gray30", col.axis="gray30", tck=-0.03, las=1, padj=0.4) 
mtext("Effective reproduction number (Rt)", side=2, col="gray30", cex=0.9, line=1.7, las=3)
box(bty="l", col="gray30", lwd=0.5)

##------Visualising and analysing the wavefront progression of an outbreak--------##

if (!require(fields)) install.packages("fields")
if (!require(raster)) install.packages("raster")
if (!require(sf)) install.packages("magrittr")
if (!require(sf)) install.packages("sf")
if (!require(sp)) install.packages("sp")
if (!require(ks)) install.packages("ks")
if (!require(RColorBrewer)) install.packages("RColorBrewer")
if (!require(lubridate)) install.packages("lubridate")
library(fields); library(raster); library(sf); library(sp); library(ks); library(RColorBrewer); library(lubridate)

borders = shapefile("Vectorial_files_study_area/Country_borders_in_study_area.shp")
motorways = shapefile("Motorways_lines_re-projected.shp")

land_cover_variables = raster("Raster_Wallonian_forest_areas.asc")


data1 = read.csv("ASFV_genomic_analyses_simulated_dataset_1-3.csv",head=T); head(data1)
data1$days_death <- interval(data1$collection_date[1], data1$collection_date)%/%days(1)
head(data1)
data1 = data1[order(data1$collection_date),]


rast = land_cover_variables # just to get a template raster (we will not consider its cell values)
gridSize = c(rast@ncols, rast@nrows) # dimensions of the grid for kernel density computations
xyMin = c(rast@extent@xmin, rast@extent@ymin); xyMax = c(rast@extent@xmax, rast@extent@ymax)
H = ks::Hpi(data1[,c("longitude","latitude")]) # to select the bandwidth
data2 = data1[which(data1[,"days_death"]==0),] # will ultimately gather all filtered records
buffer = data1[which(data1[,"days_death"]==0),] # for a specific day, this table will gather all
# the occurrence records occurring until that day (included)
for (i in 1: max(data1[,"days_death"])) {
  cat("\tDay",i,"\n",sep=" ") # to display the progression of the loop
  indices1 = which(data1[,"days_death"]==i) # to only select occurrence data recorded on day i
  if (length(indices1) > 0) { # to only get a new kernel density computation if new cases that day
    if (dim(data2)[1] >= 3) { # ...and if there are at least three filtered cases already stored
      kde = ks::kde(buffer[,c("longitude","latitude")], compute.cont=T,
                    H=H, gridsize=gridSize, xmin=xyMin, xmax=xyMax) # kernel density computation
      r = raster(kde) # to convert the resulting kernel density polygon(s) in a raster format
      contour = rasterToContour(r, levels=kde$cont["5%"]) # to get the 95% contour polygon(s)
      # i.e. the polygon(s) gathering 95% of the occurrence records
      threshold = contourLevels(kde, 0.05) # to get the raster cell threshold value above
      # which we get 95% of the occurrence records
      r[r[]<threshold] = NA # to empty the raster cells with a value below that threshold
      r[!is.na(r[])] = i # to assign the number of days (equal to i) to the other cells
      crs(r) = crs(rast) # to assign to this raster the projection of the original raster
      file_name = paste0("KDE95_contour_polygons/KDE_95_contour_day_",i-1,".tif")
      writeRaster(r, file_name, overwrite=T) # to save 95% contour polygon(s) in a raster
      file_name = paste0("KDE95_contour_polygons/KDE_95_contour_day_",i-1,".shp")
      st_write(st_as_sf(contour), file_name, update=T) # ...and in a shapefile
      indices2 = c() # a vector to store all the filtered cases until day i
      for (j in 1:length(indices1)) { # a loop to check if each new occurrence records of day
        point_in_polygon = FALSE
        # i is located or not within the 95% contour polygon(s)
        pt.x = data1[indices1[j],"longitude"]
        pt.y = data1[indices1[j],"latitude"]
        for (k in 1:length(contour@lines[[1]]@Lines)) {
          pol.x = contour@lines[[1]]@Lines[[k]]@coords[,1]
          pol.y = contour@lines[[1]]@Lines[[k]]@coords[,2]
          if (point.in.polygon(pt.x, pt.y, pol.x, pol.y) != 0) {
            point_in_polygon = TRUE
          }
        }
        if (point_in_polygon == FALSE) {
          indices2 = c(indices2, indices1[j]) # the occurrence record is only conserved if
        }
        # not included within the 95% contour polygon(s)
      }
      # (i.e. if extending the wavefront extent)
      if (length(indices2) > 0) {
        data2 = rbind(data2, data1[indices2,])
      }
    }
    else {
      data2 = rbind(data2, data1[indices1,])
    }
    buffer = rbind(buffer, data1[indices1,])
  }
}
data2 = unique(data2) # to remove potential duplicates entries
write.csv(data2, "ASFV_filtered_cases_travail_2.csv", quote=F, row.names=F)


buffer = land_cover_variables; buffer[] = NA # just to get an empty raster
for (i in 1:max(data1[,"days_death"])) {
  if (file.exists(paste0("KDE95_contour_polygons/KDE_95_contour_day_",i,".tif"))) {
    kde = raster(paste0("KDE95_contour_polygons/KDE_95_contour_day_",i,".tif"))
    buffer = merge(buffer, kde)
  }
}
writeRaster(buffer, "KDE95_contour_raster_travail_2.tif", overwrite=T)

data2 = read.csv("ASFV_filtered_cases_travail_2.csv", head=T)
template = land_cover_variables; template[!is.na(template[])] = 0
H = Hpi(data2[,c("longitude","latitude")])
kde = kde(data2[,c("longitude","latitude")],
          H=H, compute.cont=T, gridsize=c(1000,1000))
rast1 = raster(kde); contour = rasterToContour(rast1, levels=kde$cont["5%"])
coords = data2[,c("longitude","latitude")]
threshold = min(raster::extract(rast1, coords))
p = Polygon(contour@lines[[1]]@Lines[[1]]@coords)
ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
rast2 = mask(rast1, sps); mask = raster::resample(rast2, template)

coords = data2[,c("longitude","latitude")] # filtered occurrence records
tps_model = Tps(x=coords, Y=data2[,"days_death"]) # to train the thin plate spline surface model
tps = interpolate(template, tps_model) # to conduct the actual interpolation using this model
tps_mask = mask(tps, mask) # using the mask to only keep the interpolation within the invaded area

raster_resolution = mean(c(res(tps_mask)[1],res(tps_mask)[2]))
f = matrix(1/raster_resolution, nrow=3, ncol=3)
f[c(1,3,7,9)] = 1/(sqrt(2)*raster_resolution); f[5] = 0

fun = function(x, ...) {
  sum(abs(x-x[5])*f)/8
}
friction = focal(tps_mask, w=matrix(1,nrow=3,ncol=3), fun=fun, pad=T, padValue=NA, na.rm=F)

myAvSize = 11
friction_sd10 = focal(friction,
                      w=matrix(1/(myAvSize^2), nrow=myAvSize, ncol=myAvSize),
                      pad=T, padValue=NA, na.rm=F)

wavefrontVelocity_sd10 = ((1/(friction_sd10))/1000)*7
print(mean(wavefrontVelocity_sd10[], na.rm=T)) # approximately 2 (2.079942) km/semaines

tps_mask[tps_mask[]<0] = 0 # to discard negative interpolated values for the first time of invasion
r = wavefrontVelocity_sd10; r[(!is.na(r[]))&(r[]>1)] = 1; wavefrontVelocity_sd10_truncated = r
# to get a truncated version of "wavefrontVelocity_sd10", in which all
# raster cell values >1 are set to 1 (just for a visualisation purpose)
cols_land_cover = cols_land_cover = c("white",                              # valeur 0 → blanc (non-forêt)
                                      rgb(68/255,165/255,68/255,0.2))       # valeur 1 → vert (forêt)
# colour scale for the land cover raster (to get the forest areas
# coloured in light green and the non-forest areas in white)
cols_invasion_times = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(161)[21:121])
# colour scale for the 3 first maps dedicated to first invasion times
cols_wavefront_velocity = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(161)[21:121])
# colour scale for the 4° map dedicated to the wavefront velocities

x_mid = mean(c(xmin(land_cover_variables), xmax(land_cover_variables)))
par(mfrow=c(2,2), mar=c(0.5,0.0,0.5,2.0), oma=c(0,0,3.0,0),
    mgp=c(0,0.4,0), lwd=0.2, bty="o") # general graphical parameters

# (1) To plot the spatio-temporal distribution of all occurrence records:
days = data1[,"days_death"] # to consider all occurrence records
cols_points = cols_invasion_times[(((days-min(days))/(max(days)-min(days)))*100)+1]
plot(land_cover_variables, col=cols_land_cover, box=F, axes=F, legend=F)
plot(borders, add=T, lwd=3, col="white", lty=1); plot(borders, add=T, lwd=0.5, col="black", lty=2)
plot(motorways, add=T, lwd=0.75, col=rgb(222,67,39,255,maxColorValue=255), lty=1)
for (i in dim(data1)[1]:1) {
  all_occurrence_records = data1[i,c("longitude","latitude")]
  points(all_occurrence_records, pch=16, cex=0.7, col=cols_points[i])
  points(all_occurrence_records, pch=1, cex=0.7, col="gray30", lwd=0.35)
}
mtext("1. First invasion times", side=3, line=-1.4, at=x_mid, cex=0.8, font=1, col="black")
mtext("(all, in days)", side=3, line=-2.2, at=x_mid, cex=0.8, font=1, col="black")
rect(xmin(tps_mask), ymin(tps_mask), xmax(tps_mask), ymax(tps_mask),
     xpd=T, lwd=0.2, border="gray30") # to add a rectangle around the map
legendRast = raster(as.matrix(seq(0,max(data1[,"days_death"]),1)))
plot(legendRast, legend.only=T, col=cols_invasion_times, legend.width=0.5,
     legend.shrink=0.3, smallplot=c(0.895,0.910,0.025,0.975), alpha=1,
     legend.args=list(text="", cex=0.5, line=0.5, col="gray30"),
     axis.args=list(cex.axis=0.75, lwd=0, lwd.tick=0.2, tck=-0.8, col.axis="gray30",
                    line=0, mgp=c(0,0.45,0), at=seq(0,360,30)))

# (2) To plot the spatio-temporal distribution of the filtered occurrence records:
days = data2[,"days_death"] # only filtered occurrence records
cols_points = cols_invasion_times[(((days-min(days))/(max(days)-min(days)))*100)+1]
plot(land_cover_variables, col=cols_land_cover, box=F, axes=F, legend=F)
plot(borders, add=T, lwd=3, col="white", lty=1); plot(borders, add=T, lwd=0.5, col="black", lty=2)
plot(motorways, add=T, lwd=0.75, col=rgb(222,67,39,255,maxColorValue=255), lty=1)
for (i in dim(data2)[1]:1) {
  filetered_occurrence_records = data2[i,c("longitude","latitude")]
  points(filetered_occurrence_records, pch=16, cex=0.7, col=cols_points[i])
  points(filetered_occurrence_records, pch=1, cex=0.7, col="gray30", lwd=0.35)
}
mtext("2. First invasion times", side=3, line=-1.4, at=x_mid, cex=0.8, font=1, col="black")
mtext("(filtered, in days)", side=3, line=-2.2, at=x_mid, cex=0.8, font=1, col="black")
rect(xmin(wavefrontVelocity_sd10), ymin(wavefrontVelocity_sd10), xmax(wavefrontVelocity_sd10),
     ymax(wavefrontVelocity_sd10), xpd=T, lwd=0.2, border="gray30")
legendRast = raster(as.matrix(seq(0,max(data2[,"days_death"]),1)))
plot(legendRast, legend.only=T, col=cols_invasion_times, legend.width=0.5,
     legend.shrink=0.3, smallplot=c(0.895,0.910,0.025,0.975), alpha=1,
     legend.args=list(text="", cex=0.5, line=0.5, col="gray30"),
     axis.args=list(cex.axis=0.75, lwd=0, lwd.tick=0.2, tck=-0.8, col.axis="gray30",
                    line=0, mgp=c(0,0.45,0), at=seq(0,360,30))) # colour legend

# (3) To plot the first invasion times interpolated within the invaded area:
plot(land_cover_variables,col=cols_land_cover, box=F, axes=F, legend=F)
plot(tps_mask,add = T, main="", cex.main=1, cex.axis=0.7, bty="n", box=F, axes=F, legend=F,
     axis.args=list(cex.axis=0.7), col=cols_invasion_times, colNA=NA)
points(data2[,c("longitude","latitude")], pch=3, cex=0.3, lwd=0.3, col="gray30")
plot(borders, add=T, lwd=3, col="white", lty=1); plot(borders, add=T, lwd=0.5, col="black", lty=2)
plot(motorways, add=T, lwd=0.75, col=rgb(222,67,39,255,maxColorValue=255), lty=1)
mtext("3. First invasion times", side=3, line=-1.4, at=x_mid, cex=0.8, font=1, col="black")
mtext("(interpolated, in days)", side=3, line=-2.2, at=x_mid, cex=0.8, font=1, col="black")
rect(xmin(tps_mask), ymin(tps_mask), xmax(tps_mask), ymax(tps_mask),
     xpd=T, lwd=0.2, border="gray30") # to add a rectangle around the map
plot(tps_mask, legend.only=T, add=T, col=cols_invasion_times, legend.width=0.5,
     legend.shrink=0.3, smallplot=c(0.895,0.910,0.025,0.975), alpha=1,
     legend.args=list(text="", cex=0.5, line=0.5, col="gray30"),
     axis.args=list(cex.axis=0.75, lwd=0, lwd.tick=0.2, tck=-0.8, col.axis="gray30",
                    line=0, mgp=c(0,0.45,0), at=seq(0,360,30))) # colour legend

# (4) To plot the local wavefront velocity values estimated within the invaded area:
plot(land_cover_variables, col=cols_land_cover, box=F, axes=F, legend=F)
plot(wavefrontVelocity_sd10_truncated,add = T, main="", cex.main=1, cex.axis=0.7, bty="n", box=F, axes=F,
     legend=F, axis.args=list(cex.axis=0.7), col=cols_wavefront_velocity, colNA=NA)
points(data2[,c("longitude","latitude")], pch=3, cex=0.3, lwd=0.3, col="gray30")
plot(borders, add=T, lwd=3, col="white", lty=1); plot(borders, add=T, lwd=0.5, col="black", lty=2)
plot(motorways, add=T, lwd=0.75, col=rgb(222,67,39,255,maxColorValue=255), lty=1)
mtext("4. Wavefront velocity", side=3, line=-1.4, at=x_mid, cex=0.8, font=1, col="black")
mtext("(km/week)", side=3, line=-2.2, at=x_mid, cex=0.8, font=1, col="black")
rect(xmin(wavefrontVelocity_sd10), ymin(wavefrontVelocity_sd10), xmax(wavefrontVelocity_sd10),
     ymax(wavefrontVelocity_sd10), xpd=T, lwd=0.2, border="gray30")
plot(wavefrontVelocity_sd10_truncated, legend.only=T, add=T, col=cols_wavefront_velocity,
     legend.width=0.5, legend.shrink=0.3, smallplot=c(0.895,0.910,0.025,0.975), alpha=1,
     legend.args=list(text="", cex=0.5, line=0.5, col="gray30"),
     axis.args=list(cex.axis=0.75, lwd=0, lwd.tick=0.2, tck=-0.8, col.axis="gray30", line=0,
                    mgp=c(0,0.45,0), at=seq(0,1,0.1), labels=c(seq(0,0.9,0.1),">=1")))

