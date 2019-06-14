library(ncdf4)

library(gdalUtils)
library(raster)

#setwd("../MODIS/")

sds = get_subdatasets("MODIS/MYD05_L2.A2019087.1345.061.NRT.hdf" )  # water vapor set
gdalinfo("MODIS/MYD05_L2.A2019087.1345.061.NRT.hdf")



gdal_translate(sds[7], "test.tiff")
rast <- raster("test.tiff")
plot(rast)
rast_df_na = as.data.frame(rast, xy=TRUE)
sum(!is.na(rast_df_na$test))
rast_df = rast_df_na[!is.na(rast_df_na$test),]
#pdf(file="MODIS/water_vapor.pdf", width = 8, height = 6)
#quilt.plot(rast_df$x, rast_df$y, rast_df$test, nx = 1354, ny=2030)
#dev.off()
length(rast_df$x)
colnames(rast_df)<- c("x","y","vapor")
dim(rast_df)
2030*1354

write.csv( rast_df, "MODIS/water_vapor_data.csv")



# remove data with some pattern, for prediction later (not used)
sub_idx = sample(length(rast_df$x), 250000)

sub_rast_df = rast_df[sub_idx,]
striped_rast_df = rast_df
striped_rast_df[which(striped_rast_df$x > 800 & striped_rast_df$x < 850),3] <- NA
striped_rast_df = striped_rast_df[which(striped_rast_df$y<1500),]
striped_rast_df = striped_rast_df[which(striped_rast_df$x>454),]
striped_rast_df[which(striped_rast_df$y > 200 & striped_rast_df$y < 250),3] <- NA
removed_idx=c()
for(y in 50*2:29){
  for( x in 450+50*2:17){
    #striped_rast_df[which(striped_rast_df$y > (y-10) & striped_rast_df$y < (y+10)
                   #       & striped_rast_df$x > (x-10) & striped_rast_df$x < (x+10)),3] <- NA
    removed_idx = c(removed_idx, which(striped_rast_df$y > (y-10) & striped_rast_df$y < (y+10)
                                       & striped_rast_df$x > (x-10) & striped_rast_df$x < (x+10)))
  }
}

removed_idx = c(removed_idx,
                which(striped_rast_df$x > 800 & striped_rast_df$x < 850),
                which(striped_rast_df$y>1500),
                which(striped_rast_df$x<454),
                which(striped_rast_df$y > 200 & striped_rast_df$y < 250))

removed_idx = unique(removed_idx)
water_vapor_pred= rast_df[removed_idx,]

write.csv( striped_rast_df, "water_vapor_data_missing.csv")


pdf(file="MODIS/water_vapor_set_strp.pdf", width = 8, height = 6)
quilt.plot(striped_rast_df$x, striped_rast_df$y, striped_rast_df$vapor, nx = 900, ny=1500)#2030)
dev.off()

save(water_vapor_pred, file="MODIS/pred_water_vapor_20190328.RData")

# 12 (latitude) & 13 (longitude)
gdal_translate(sds[13], "lat.tiff")
rast <- raster("lat.tiff")
plot(rast)
lat_df = as.data.frame(rast, xy=TRUE)
min(lat_df$lat)
max(lat_df$lat)

"  SOUTHBOUNDINGCOORDINATE=45.1262553445153"
"  NORTHBOUNDINGCOORDINATE=67.4762902910701"
"  EASTBOUNDINGCOORDINATE=4.44383248127601"
"  WESTBOUNDINGCOORDINATE=-42.7070983072549"
