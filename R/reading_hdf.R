
library(rgdal)
library(gdalUtils)
library(raster)
dev.off()
par(mfrow=c(1,1), mar = c(2.2,2.2,2.2,1))
file_loc = "../../../Downloads/MISR_AM1_CGAS_APR_01_2010_F15_0031.hdf"
file_loc = "../../../Downloads/MISR_AM1_CMV_T20180501002721_P100_O097695_F01_0001.hdf"
file_loc = "../../../Downloads/MISR_AM1_CMV_T20180501052359_P148_O097698_F01_0001.hdf"

sds = get_subdatasets(file_loc)
gdalinfo(file_loc)


?raster
gdal_translate(sds[3], dst_dataset = "cloud9.tif")
rast <- raster("cloud9.tif.aux.xml")
plot(rast)
rast
#calculate min and max raster values
rast = setMinMax(rast)
# check range of raster values
cellStats(rast, range)
#view coord system
rast@crs
# exact bounds of extent
rast@extent
hist(rast)
rast_df = as.data.frame(rast, xy=TRUE)
head(rast_df)
rast_df[!is.na(rast_df)]







######################################################################################
###################################   CMV   ##########################################
#####################################################################################

library(ncdf4)
ncin <- nc_open(filename = "../MISR_data/CMV/MISR_AM1_CMV_AUG_2015_F02_0002.nc")
print(ncin)
lon <- ncvar_get(ncin, "Longitude")
lat <- ncvar_get(ncin, "Latitude")
doy <- ncvar_get(ncin, "DayOfYear")

cm_east <- ncvar_get(ncin, "CloudMotionEast")
cm_north <- ncvar_get(ncin, "CloudMotionNorth")
t <- ncvar_get(ncin, "Time")
dim(t)
tunits <- ncatt_get(ncin, "Time", "units")

nt <- dim(t)
ncvar_get(ncin)
t[1:100]
lat[1:100]

t-t[1]


idx_plot = which(doy==213)#c(which(doy==214), which(doy==213))
idx_plot=idx_plot[1600:2000]
quilt.plot(lon[idx_plot], lat[idx_plot], (cm_east[idx_plot]), nx = 100, ny=100,main ="Cloud Motion Vector, East")

#pdf("East_winds.pdf", height = 8.5, width = 11)
quilt.plot(lon,lat,cm_east, nx = 1000, ny=1000, main ="Cloud Motion Vector, East")
#dev.off()


# split the time units string into fields
tustr <- strsplit(tunits$value, " ")
tdstr <- strsplit(unlist(tustr)[3], "-")
tmonth = as.integer(unlist(tdstr)[2])
tday = as.integer(unlist(tdstr)[3])
tyear = as.integer(unlist(tdstr)[1])
chron(t, origin = c(tmonth, tday, tyear))



######################################################################################
################################   Aerosol   ########################################
#####################################################################################


library(ncdf4)
get_dats = function(filename){
  ncin <- nc_open(filename)

  lon <- ncvar_get(ncin, "4.4_KM_PRODUCTS/Longitude")
  ydim <- ncvar_get(ncin, "4.4_KM_PRODUCTS/Y_Dim")
  xdim <- ncvar_get(ncin, "4.4_KM_PRODUCTS/X_Dim")
  lat <- ncvar_get(ncin, "4.4_KM_PRODUCTS/Latitude")
  elevation <- ncvar_get(ncin, "4.4_KM_PRODUCTS/Elevation")
  year <- ncvar_get(ncin, "4.4_KM_PRODUCTS/Year")
  doy <- ncvar_get(ncin, "4.4_KM_PRODUCTS/Day_Of_Year")

  AOD <- ncvar_get(ncin, "4.4_KM_PRODUCTS/Aerosol_Optical_Depth")
  SAOD <- ncvar_get(ncin, "4.4_KM_PRODUCTS/Small_Mode_Aerosol_Optical_Depth")
  file_data = list("AOD" = AOD[!is.na(AOD)],
                   "lon"=lon[!is.na(AOD)],
                   "lat"= lat[!is.na(AOD)]
                   )

  return(file_data)

}

f1 = "../MISR_data/Aerosol_7_1_2015/MISR_AM1_AS_AEROSOL_P014_O082632_F13_0023.nc"




loc = "../MISR_data/Aerosol_7_1_2015_week/"
loc = "../MISR_data/Aerosol_7_1_2015/"
aero_files = grep(pattern = ".*\\.nc(?!\\.xml)",
                  list.files(loc),
                  value = T, perl = T)

aero_data = c()
for( fi in aero_files){
  fp = file.path(loc,fi)
  fi_dats = get_dats(fp)
  file_matrix = cbind(fi_dats$lon, fi_dats$lat,fi_dats$AOD)
  aero_data = rbind(aero_data, file_matrix)
}
negi = which(aero_data[,1] < 0)

pdf("AOD_1week.pdf", height = 6, width = 8)
quilt.plot(aero_data[,1], aero_data[,2],aero_data[,3], nx = 1000, ny = 1000, main = "Aerosol Optical Depth, July week 2015")
dev.off()



dim(aero_data)
head(aero_data)



######################################################################################
################################   VL with Aerosol   #################################
#####################################################################################



data.distr = 'gamma' # options: "gaussian","logistic", "poisson", "gamma"
spatial.dim = 2 # number of spatial dimensions
n=length(aero_data[,1]) # number of observed locs
default_lh_params = list("alpha"=2, "sigma"=.3)
locs = aero_data[,1:2]
z = aero_data[,3]

m=3
vecchia.approx=vecchia_specify(z, locs, m)
covparms=c(1, .2, 2)
# Perform inference on latent mean with Vecchia Laplace approximation
posterior = calculate_posterior_VL(vecchia.approx, likelihood_model=data.distr, covparms,likparms = default_lh_params, max.iter = 50)


