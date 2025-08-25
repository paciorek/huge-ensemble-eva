## This takes about 4 hours for precip.

library(climextRemes)

var <- commandArgs(trailingOnly=TRUE)[1]

cat("Processing var: ", var, "\n")

load(paste0(var, '.Rda')) # 520 sec 

source("init.R")

nt <- dim(daily)[3]
last <- (nt-365+1):nt
## Indices of leap year days (Feb. 29).
leaps <- c(365*3+31+29, 365*6+366+31+29, 365*9+366*2+31+29, 365*12+366*3+31+29,
           365*15+366*4+31+29)

## POT

shapes <- array(as.numeric(NA), c(nlon, nlat, nthr))
failures <- array(as.numeric(NA), c(nlon, nlat, nthr))
rvs <- array(as.numeric(NA), c(nlon, nlat, nrp, nthr))
se_rvs <- array(as.numeric(NA), c(nlon, nlat, nrp, nthr))
ns <- array(as.numeric(NA), c(nlon, nlat, nthr))

fits <- list(); length(fits) <- nthr

## GEV

shapes_gev <- array(as.numeric(NA), c(nlon, nlat))
failures_gev <- array(as.numeric(NA), c(nlon, nlat))
rvs_gev <- array(as.numeric(NA), c(nlon, nlat, nrp))
se_rvs_gev <- array(as.numeric(NA), c(nlon, nlat, nrp))


## seasonal POT

seas <- list(); length(seas) <- 4
nthr_seas <- nthr-2  # Not enough stratified data for highest thresholds when using full-year thresholds.
fits_seas1 <- list();  length(fits_seas1) <- nthr_seas

shapes_seas1 <- array(as.numeric(NA), c(nlon, nlat, 4, nthr_seas))
failures_seas1 <- array(as.numeric(NA), c(nlon, nlat, 4, nthr_seas))
rvs_seas1 <- array(as.numeric(NA), c(nlon, nlat, 4, nrp, nthr_seas))
se_rvs_seas1 <- array(as.numeric(NA), c(nlon, nlat, 4, nrp, nthr_seas))
ns_seas1 <- array(as.numeric(NA), c(nlon, nlat, 4, nthr_seas))

shapes_seas2 <- array(as.numeric(NA), c(nlon, nlat, 4, nthr))
failures_seas2 <- array(as.numeric(NA), c(nlon, nlat, 4, nthr))
rvs_seas2 <- array(as.numeric(NA), c(nlon, nlat, 4, nrp, nthr))
se_rvs_seas2 <- array(as.numeric(NA), c(nlon, nlat, 4, nrp, nthr))
ns_seas2 <- array(as.numeric(NA), c(nlon, nlat, 4, nthr))

fits_seas2 <- list();  length(fits_seas2) <- nthr

emp_quants <- array(as.numeric(NA), c(nlon, nlat, length(emp_qs)))


for(focal_lon in seq_len(nlon)) {
    cat("Working on longitude", focal_lon, "\n")
    for(focal_lat in seq_len(nlat)) {
        cat("Working on latitude", focal_lat, "\n")
        focal <- daily[focal_lon, focal_lat, ]

        quants <- quantile(focal, c(.99, emp_qs, qs))

        ## empirical quantiles
        emp_quants[focal_lon, focal_lat,] <- quants[2:3]

        ## POT
        
        thresh <- quants[4:(length(qs)+3)]  # Thresholds for POT.

        sub <- focal[focal > quants[1]]  # Initial coarse thresholding to reduce computation.
        for(i in seq_along(thresh)) {
            ns[focal_lon, focal_lat,i] <- sum(sub > thresh[i]) 
            fits[[i]] <- fit_pot(sub[sub > thresh[i]], threshold = thresh[i], nBlocks = nBlocks, returnPeriod = returnPeriods, getParams=TRUE)
        }
        
        failures[focal_lon, focal_lat,] <- sapply(fits, function(x) x$info$failure)
        wh <- !failures[focal_lon, focal_lat,]
        shapes[focal_lon, focal_lat, wh] <- sapply(fits[wh], function(x) x$mle['shape'])
        rvs[focal_lon, focal_lat, seq_len(nrp), wh] <- sapply(fits[wh], function(x) x$returnValue)
        se_rvs[focal_lon, focal_lat, seq_len(nrp), wh] <- sapply(fits[wh], function(x) x$se_returnValue)

        ## GEV

        ## Shift full last replicate of year 2001 to replace partial initial year (previously discarded).
        blocks <- matrix(c(focal[last], focal[1:(last[1]-1)]), nrow = ndays)
        yearly <- matrix(blocks[-leaps,], nrow = 365) # Omit leap days for simplicity.

        mx <- apply(yearly, 2, max)
        fit_gev <- fit_gev(mx, returnPeriod = returnPeriods, getParams=TRUE)

        failures_gev[focal_lon, focal_lat] <- fit_gev$info$failure
        if(!fit_gev$info$failure) {
            shapes_gev[focal_lon, focal_lat] <- fit_gev$mle['shape']
            rvs_gev[focal_lon, focal_lat, ] <- fit_gev$returnValue
            se_rvs_gev[focal_lon, focal_lat, ] <- fit_gev$se_returnValue
        }

        ## seasonal POT
        seas[[1]] <- yearly[c(1:(31+28),(365-31+1):365),]
        seas[[2]] <- yearly[(31+28+1):(31+28+31+30+31),]
        seas[[3]] <- yearly[(31+28+31+30+31+1):(31+28+31+30+31+30+31+31),]
        seas[[4]] <- yearly[(31+28+31+30+31+30+31+31+1):(365-31),]

        ## Approach 1: use same threshold as for full-year analysis,
        ## so that seasonal analysis uses as much data.
        ## There can be many failures due to few observations (or no observations) in certain seasons in locations with seasonality.
        for(j in 1:4) {
            for(i in seq_len(nthr_seas)) {
                fits_seas1[[i]] <- fit_pot(seas[[j]][seas[[j]] > thresh[i]], threshold = thresh[i], nBlocks = nBlocks,
                                         returnPeriod = returnPeriods, getParams=TRUE)
                ns_seas1[focal_lon, focal_lat,j,i] <- sum(seas[[j]] > thresh[i])
            }
            failures_seas1[focal_lon, focal_lat,j,] <- sapply(fits_seas1, function(x) x$info$failure)
            wh <- !failures_seas1[focal_lon, focal_lat,j,]
            if(sum(wh)) {
                shapes_seas1[focal_lon, focal_lat, j, wh] <- sapply(fits_seas1[wh], function(x) x$mle['shape'])
                rvs_seas1[focal_lon, focal_lat, j, seq_len(nrp), wh] <- sapply(fits_seas1[wh], function(x) x$returnValue)
                se_rvs_seas1[focal_lon, focal_lat, j, seq_len(nrp), wh] <- sapply(fits_seas1[wh], function(x) x$se_returnValue)
            }
        }

        ## Approach 2: use season-specific thresholds based on full-year qunatiles
        ## so that seasonal analysis uses 4x as much data.
        for(j in 1:4) {
            thresh_seas <- quantile(seas[[j]], qs)
            
            for(i in seq_len(nthr)) {
                fits_seas2[[i]] <- fit_pot(seas[[j]][seas[[j]] > thresh_seas[i]], threshold = thresh_seas[i], nBlocks = nBlocks,
                                         returnPeriod = returnPeriods, getParams=TRUE)
                ns_seas2[focal_lon, focal_lat,j,i] <- sum(seas[[j]] > thresh_seas[i])
            }
            failures_seas2[focal_lon, focal_lat,j,] <- sapply(fits_seas2, function(x) x$info$failure)
            wh <- !failures_seas2[focal_lon, focal_lat,j,]
            if(sum(wh)) {
                shapes_seas2[focal_lon, focal_lat, j, wh] <- sapply(fits_seas2[wh], function(x) x$mle['shape'])
                rvs_seas2[focal_lon, focal_lat, j, seq_len(nrp), wh] <- sapply(fits_seas2[wh], function(x) x$returnValue)
                se_rvs_seas2[focal_lon, focal_lat, j, seq_len(nrp), wh] <- sapply(fits_seas2[wh], function(x) x$se_returnValue)
            }
        }
    }
}

rvs_max_seas1 <- apply(rvs_seas1, c(1,2,4,5), max, na.rm = TRUE)
## Season the max occurs in (for diagnostics).
max_seas1 <- apply(rvs_seas1, c(1,2,4,5), function(x) if(sum(!is.na(x))) which.max(x) else as.numeric(NA))

rvs_max_seas2 <- apply(rvs_seas2, c(1,2,4,5), max, na.rm = TRUE)
## Season the max occurs in (for diagnostics).
max_seas2 <- apply(rvs_seas2, c(1,2,4,5), function(x) if(sum(!is.na(x))) which.max(x) else as.numeric(NA))

## Empirical extremes
max100 <- apply(daily, c(1,2), function(x) max(x[1:(100*365)]))
max1000 <- apply(daily, c(1,2), function(x) max(x[1:(1000*365)]))
max10000 <- apply(daily, c(1,2), max)

save(qs, returnPeriods, emp_qs,
    ns, failures, shapes, rvs, se_rvs,
    emp_quants,
    failures_gev, shapes_gev, rvs_gev, se_rvs_gev,
    ns_seas1, failures_seas1, shapes_seas1, rvs_seas1, se_rvs_seas1, rvs_max_seas1, max_seas1,
    ns_seas2, failures_seas2, shapes_seas2, rvs_seas2, se_rvs_seas2, rvs_max_seas2, max_seas2,
    max100, max1000, max10000,
    file = paste0(var, '_fits.Rda'))


