library(ncdf4)
nc <- nc_open("autoregressive_predictions.nc")

lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")
time <- ncvar_get(nc, "time")
init_time <- ncvar_get(nc, "init_time")

## Restrict to generous bounding box around CONUS.
lats <- which(lat >= 20 & lat <= 50)
lons <- which(lon >= (360-130) & lon <= (360-60))

start <- c(lons[1], lats[1], 1, 1)
count <- c(length(lons), length(lats), -1, -1)
system.time(temp <- ncvar_get(nc, "TMP2m", start, count))  
system.time(prec <- ncvar_get(nc, "PRATEsfc", start, count))

lon <- lon[lons]
lat <- lat[lats]

system.time(save(temp,prec,lat,lon,time,init_time, file = paste0("preds_", month, ".Rda")))
