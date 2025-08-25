library(climextRemes)

var <- commandArgs(trailingOnly=TRUE)[1]

cat("Processing var: ", var, "\n")

load('era5.Rda')

source("init.R")
nBlocks <- length(1940:2022)

leaps <- 365*4+31+29
for(i in 1:length(seq(1948,2022,by=4)))
    leaps <- c(leaps, leaps[length(leaps)]+365*3+366)


if(var == "temp") daily <- temp_era5 else daily <- prec_era5

shapes_era5 <- array(as.numeric(NA), c(nlon, nlat, 2))
se_shapes_era5 <- array(as.numeric(NA), c(nlon, nlat, 2))
failures_era5 <- array(as.numeric(NA), c(nlon, nlat, 2))
rvs_era5 <- array(as.numeric(NA), c(nlon, nlat, nrp, 2))
se_rvs_era5 <- array(as.numeric(NA), c(nlon, nlat, nrp, 2))
ns_era5 <- array(as.numeric(NA), c(nlon, nlat, 2))

shapes_gev_era5 <- array(as.numeric(NA), c(nlon, nlat))
se_shapes_gev_era5 <- array(as.numeric(NA), c(nlon, nlat))
failures_gev_era5 <- array(as.numeric(NA), c(nlon, nlat))
rvs_gev_era5 <- array(as.numeric(NA), c(nlon, nlat, nrp))
se_rvs_gev_era5 <- array(as.numeric(NA), c(nlon, nlat, nrp))
ns_gev_era5 <- array(as.numeric(NA), c(nlon, nlat))


for(focal_lon in seq_len(nlon)) {
    cat("Working on longitude", focal_lon, "\n")
    for(focal_lat in seq_len(nlat)) {
        cat("Working on latitude", focal_lat, "\n")
        focal <- daily[focal_lon, focal_lat, ]

        thresh <- quantile(focal, c(1-.0025,.999))

        for(k in 1:2) {
            ns_era5[focal_lon, focal_lat,k] <- sum(focal > thresh[k])
            fit <- fit_pot(focal[focal > thresh[k]], threshold = thresh[k], nBlocks = nBlocks, returnPeriod = returnPeriods, getParams=TRUE)
            
            failures_era5[focal_lon, focal_lat,k] <- fit$info$failure
            if(!fit$info$failure) {
                shapes_era5[focal_lon, focal_lat,k] <- fit$mle['shape']
                se_shapes_era5[focal_lon, focal_lat,k] <- fit$se_mle['shape']
                rvs_era5[focal_lon, focal_lat, seq_len(nrp),k] <- fit$returnValue
                se_rvs_era5[focal_lon, focal_lat, seq_len(nrp),k] <- fit$se_returnValue
            }
        }

        ## GEV
        yearly <- matrix(focal[-leaps], nrow = 365)

        mx <- apply(yearly, 2, max)
        fit_gev <- fit_gev(mx, returnPeriod = returnPeriods, getParams=TRUE)
        failures_gev_era5[focal_lon, focal_lat] <- fit_gev$info$failure
        if(!fit_gev$info$failure) {
            shapes_gev_era5[focal_lon, focal_lat] <- fit_gev$mle['shape']
            se_shapes_gev_era5[focal_lon, focal_lat] <- fit_gev$se_mle['shape']
            rvs_gev_era5[focal_lon, focal_lat, seq_len(nrp)] <- fit_gev$returnValue
            se_rvs_gev_era5[focal_lon, focal_lat, seq_len(nrp)] <- fit_gev$se_returnValue
        }

        

    }
}

save(ns_era5, failures_era5, shapes_era5, se_shapes_era5, rvs_era5, se_rvs_era5,
     ns_gev_era5, failures_gev_era5, shapes_gev_era5, se_shapes_gev_era5, rvs_gev_era5, se_rvs_gev_era5,
     file = paste0(var, '_era5_fits.Rda'))
