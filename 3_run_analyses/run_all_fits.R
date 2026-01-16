## ~12 hours to run. Was about 4 hours before finding aggregated seasonal return values
## by optimization.

var <- commandArgs(trailingOnly=TRUE)[1]

library(climextRemes)

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
se_shapes <- array(as.numeric(NA), c(nlon, nlat, nthr))
failures <- array(as.numeric(NA), c(nlon, nlat, nthr))
rvs <- array(as.numeric(NA), c(nlon, nlat, nrp, nthr))
se_rvs <- array(as.numeric(NA), c(nlon, nlat, nrp, nthr))
ns <- array(as.numeric(NA), c(nlon, nlat, nthr))

fits <- list(); length(fits) <- nthr

## GEV

shapes_gev <- array(as.numeric(NA), c(nlon, nlat, nlens))
se_shapes_gev <- array(as.numeric(NA), c(nlon, nlat, nlens))
failures_gev <- array(as.numeric(NA), c(nlon, nlat, nlens))
rvs_gev <- array(as.numeric(NA), c(nlon, nlat, nrp, nlens))
se_rvs_gev <- array(as.numeric(NA), c(nlon, nlat, nrp, nlens))

fits_gev <- list(); length(fits_gev) <- nlens

## seasonal POT

seas <- list(); length(seas) <- 4
nthr_seas <- nthr-2  # Not enough stratified data for highest thresholds when using full-year thresholds.
fits_seas1 <- list();  length(fits_seas1) <- 4

shapes_seas1 <- array(as.numeric(NA), c(nlon, nlat, 4, nthr_seas))
se_shapes_seas1 <- array(as.numeric(NA), c(nlon, nlat, 4, nthr_seas))
failures_seas1 <- array(as.numeric(NA), c(nlon, nlat, 4, nthr_seas))
rvs_seas1 <- array(as.numeric(NA), c(nlon, nlat, 4, nrp, nthr_seas))
se_rvs_seas1 <- array(as.numeric(NA), c(nlon, nlat, 4, nrp, nthr_seas))
ns_seas1 <- array(as.numeric(NA), c(nlon, nlat, 4, nthr_seas))
rvs_full_seas1 <- array(as.numeric(NA), c(nlon, nlat, nrp, nthr_seas))

fits_seas2 <- list();  length(fits_seas2) <- 4

shapes_seas2 <- array(as.numeric(NA), c(nlon, nlat, 4, nthr))
se_shapes_seas2 <- array(as.numeric(NA), c(nlon, nlat, 4, nthr))
failures_seas2 <- array(as.numeric(NA), c(nlon, nlat, 4, nthr))
rvs_seas2 <- array(as.numeric(NA), c(nlon, nlat, 4, nrp, nthr))
se_rvs_seas2 <- array(as.numeric(NA), c(nlon, nlat, 4, nrp, nthr))
ns_seas2 <- array(as.numeric(NA), c(nlon, nlat, 4, nthr))
rvs_full_seas2 <- array(as.numeric(NA), c(nlon, nlat, nrp, nthr))

emp_quants <- array(as.numeric(NA), c(nlon, nlat, length(emp_qs)))

sum_probs <- function(x, thr_idx, fits) {
    return(sum(sapply(1:4, function(j) {
        if(fits[[j]][[thr_idx]]$info$failure) {
            return(0)
        } else {
            ## Call low-level routines directly to avoid extra overhead.
            mle <- fits[[j]][[thr_idx]]$fit$results$par
            return(extRemes::pevd(x, mle['location'], mle['scale'], mle['shape'],
                                  lower.tail = FALSE, type = "GEV"))
        }}), na.rm = TRUE))
}


objective <- function(x, prob, thr_idx, fits) {
    prob_sum <- sum_probs(x, thr_idx, fits)
    return(abs(prob - prob_sum))
}


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
        se_shapes[focal_lon, focal_lat, wh] <- sapply(fits[wh], function(x) x$se_mle['shape'])        
        rvs[focal_lon, focal_lat, seq_len(nrp), wh] <- sapply(fits[wh], function(x) x$returnValue)
        se_rvs[focal_lon, focal_lat, seq_len(nrp), wh] <- sapply(fits[wh], function(x) x$se_returnValue)

        ## GEV

        ## Shift full last replicate of year 2001 to replace partial initial year (previously discarded).
        blocks <- matrix(c(focal[last], focal[1:(last[1]-1)]), nrow = ndays)
        data <- c(blocks[-leaps,]) # Omit leap days for simplicity.
        for(i in seq_along(blockLens)) {
            blockVals <- matrix(data, nrow = 365*blockLens[i]) 
            mx <- apply(blockVals, 2, max)
            ## Using POT as initial avoids occasional bad optimizations with annual maxima for temperature.
            fits_gev[[i]] <- fit_gev(mx, returnPeriod = returnPeriods/blockLens[i], getParams=TRUE,
                                     initial = as.list(fits[[6]]$mle))
        }

        failures_gev[focal_lon, focal_lat, ] <- sapply(fits_gev, function(x) x$info$failure)
        wh <- !failures_gev[focal_lon, focal_lat,]
        shapes_gev[focal_lon, focal_lat, wh] <- sapply(fits_gev[wh], function(x) x$mle['shape'])
        se_shapes_gev[focal_lon, focal_lat, wh] <- sapply(fits_gev[wh], function(x) x$se_mle['shape'])
        rvs_gev[focal_lon, focal_lat, , wh] <- sapply(fits_gev[wh], function(x) x$returnValue)
        se_rvs_gev[focal_lon, focal_lat, ,wh] <- sapply(fits_gev[wh], function(x) x$se_returnValue)
        

        ## seasonal POT
        yearly <- matrix(blocks[-leaps,], nrow = 365) 
        seas[[1]] <- yearly[c(1:(31+28),(365-31+1):365),]
        seas[[2]] <- yearly[(31+28+1):(31+28+31+30+31),]
        seas[[3]] <- yearly[(31+28+31+30+31+1):(31+28+31+30+31+30+31+31),]
        seas[[4]] <- yearly[(31+28+31+30+31+30+31+31+1):(365-31),]


        ## Approach 1: use same threshold as for full-year analysis,
        ## so that seasonal analysis uses as much data.
        ## There can be many failures due to few observations (or no observations) in certain seasons in locations with seasonality.
        for(j in 1:4) {
            fits_seas1[[j]] <- list(); length(fits_seas1[[j]]) <- nthr_seas
            for(i in seq_len(nthr_seas)) {
                fits_seas1[[j]][[i]] <- fit_pot(seas[[j]][seas[[j]] > thresh[i]], threshold = thresh[i], nBlocks = nBlocks,
                                         returnPeriod = returnPeriods, getParams=TRUE, getFit = TRUE)
                ns_seas1[focal_lon, focal_lat,j,i] <- sum(seas[[j]] > thresh[i])
            }
            failures_seas1[focal_lon, focal_lat,j,] <- sapply(fits_seas1[[j]], function(x) x$info$failure)
            wh <- !failures_seas1[focal_lon, focal_lat,j,]
            if(sum(wh)) {
                shapes_seas1[focal_lon, focal_lat, j, wh] <- sapply(fits_seas1[[j]][wh], function(x) x$mle['shape'])
                se_shapes_seas1[focal_lon, focal_lat, j, wh] <- sapply(fits_seas1[[j]][wh], function(x) x$se_mle['shape'])
                rvs_seas1[focal_lon, focal_lat, j, seq_len(nrp), wh] <- sapply(fits_seas1[[j]][wh], function(x) x$returnValue)
                se_rvs_seas1[focal_lon, focal_lat, j, seq_len(nrp), wh] <- sapply(fits_seas1[[j]][wh], function(x) x$se_returnValue)
            }
        }
        
        ## Now optimize to find the returnValue, x, s.t. sum of exceedance probs over season is desired value.
        maxs <- apply(rvs_seas1[focal_lon, focal_lat, , , ], c(2,3), max, na.rm = TRUE)
        mins <- apply(rvs_seas1[focal_lon, focal_lat, , , ], c(2,3), min, na.rm = TRUE)
        maxs[is.infinite(maxs)] <- NA
        for(pp in seq_along(returnPeriods)) {
            prob <- 1/returnPeriods[pp]
            for(qq in seq_len(nthr_seas)) {
                if(!is.na(maxs[pp,qq])) {
                    if(isTRUE(maxs[pp,qq] == mins[pp,qq])) {  # Valid fit for only one season.
                        rvs_full_seas1[focal_lon, focal_lat, pp, qq] <- maxs[pp,qq]
                    } else {
                        upper <- maxs[pp, qq]
                        ## As upper bound gets large, objective maxes out at the desired probability, so function is
                        ## flat, causing optimization issues. So take steps to find reasonable upper bound.
                        while(sum_probs(upper, qq, fits_seas1) > prob/10 && upper < 200) 
                            upper <- upper + 1
                        if(upper > 100)
                            warning("Reached upper value of 100: ", maxs[pp,qq], " ", pp, " " , qq, " ", focal_lon, " ", focal_lat)
                        if(sum_probs(upper, qq, fits_seas1) == 0)  # Upper bound in flat area; walk back some.
                            while(sum_probs(upper, qq, fits_seas1) < prob/10)
                                upper <- upper - 0.01
                        ## optimResult <- optim(inits[pp,qq], find_annual, prob = prob, thr_idx = qq, fits = fits_seas1, method = 'BFGS')
                        out <- try(optimResult <- optimize(objective, c(maxs[pp,qq], upper), prob, qq, fits_seas1))
                        if(is(out,'try-error')) {
                            warning("optim error: ", pp, " ", qq, " ", focal_lon, " ", focal_lat, " ", maxs[pp,qq])
                            result <- NA
                        } else result <- optimResult$minimum
                        rvs_full_seas1[focal_lon, focal_lat, pp, qq] <- result 
                    }
                }
            }
        }
            

        ## Approach 2: use season-specific thresholds based on full-year quantiles
        ## so that seasonal analysis uses 4x as much data.
        for(j in 1:4) {
            fits_seas2[[j]] <- list(); length(fits_seas2[[j]]) <- nthr
            thresh_seas <- quantile(seas[[j]], qs)
            
            for(i in seq_len(nthr)) {
                fits_seas2[[j]][[i]] <- fit_pot(seas[[j]][seas[[j]] > thresh_seas[i]], threshold = thresh_seas[i], nBlocks = nBlocks,
                                         returnPeriod = returnPeriods, getParams=TRUE, getFit = TRUE)
                ns_seas2[focal_lon, focal_lat,j,i] <- sum(seas[[j]] > thresh_seas[i])
            }
            failures_seas2[focal_lon, focal_lat,j,] <- sapply(fits_seas2[[j]], function(x) x$info$failure)
            wh <- !failures_seas2[focal_lon, focal_lat,j,]
            if(sum(wh)) {
                shapes_seas2[focal_lon, focal_lat, j, wh] <- sapply(fits_seas2[[j]][wh], function(x) x$mle['shape'])
                se_shapes_seas2[focal_lon, focal_lat, j, wh] <- sapply(fits_seas2[[j]][wh], function(x) x$se_mle['shape'])
                rvs_seas2[focal_lon, focal_lat, j, seq_len(nrp), wh] <- sapply(fits_seas2[[j]][wh], function(x) x$returnValue)
                se_rvs_seas2[focal_lon, focal_lat, j, seq_len(nrp), wh] <- sapply(fits_seas2[[j]][wh], function(x) x$se_returnValue)
            }
        }

        ## Now optimize to find the returnValue, x, s.t. sum of exceedance probs over season is desired value.
        maxs <- apply(rvs_seas2[focal_lon, focal_lat, , , ], c(2,3), max, na.rm = TRUE)
        mins <- apply(rvs_seas2[focal_lon, focal_lat, , , ], c(2,3), min, na.rm = TRUE)
        maxs[is.infinite(maxs)] <- NA
        for(pp in seq_along(returnPeriods)) {
            prob <- 1/returnPeriods[pp]
            for(qq in seq_len(nthr_seas)) {
                if(!is.na(maxs[pp,qq])) {
                    if(isTRUE(maxs[pp,qq] == mins[pp,qq])) {  # Valid fit for only one season.
                        rvs_full_seas2[focal_lon, focal_lat, pp, qq] <- maxs[pp,qq]
                    } else {
                        upper <- maxs[pp, qq]
                        ## As upper bound gets large, objective maxes out at the desired probability, so function is
                        ## flat, causing optimization issues. So take steps to find reasonable upper bound.
                        while(sum_probs(upper, qq, fits_seas2) > prob/10 && upper < 200) 
                            upper <- upper + 1
                        if(upper > 100)
                            warning("Reached upper value of 100: ", maxs[pp,qq], " ", pp, " " , qq, " ", focal_lon, " ", focal_lat)
                        if(sum_probs(upper, qq, fits_seas2) == 0)   # Upper bound in flat area; walk back some.
                            while(sum_probs(upper, qq, fits_seas2) < prob/10)
                                upper <- upper - 0.01
                        ## optimResult <- optim(inits[pp,qq], find_annual, prob = prob, thr_idx = qq, fits = fits_seas2, method = 'BFGS')
                        out <- try(optimResult <- optimize(objective, c(maxs[pp,qq], upper), prob, qq, fits_seas2))
                        if(is(out,'try-error')) {
                            warning("optim error: ", pp, " ", qq, " ", focal_lon, " ", focal_lat, " ", maxs[pp,qq])
                            result <- NA
                        } else result <- optimResult$minimum
                        rvs_full_seas2[focal_lon, focal_lat, pp, qq] <- result 
                    }
                }
            }
        }
    }
}

## Naive approach (max over the seasons) as lower bound.
rvs_max_seas1 <- apply(rvs_seas1, c(1,2,4,5), max, na.rm = TRUE)

## Season the max occurs in (for diagnostics).
max_seas1 <- apply(rvs_seas1, c(1,2,4,5), function(x) if(sum(!is.na(x))) which.max(x) else as.numeric(NA))

## Naive approach (max over the seasons) as lower bound.
rvs_max_seas2 <- apply(rvs_seas2, c(1,2,4,5), max, na.rm = TRUE)

## Season the max occurs in (for diagnostics).
max_seas2 <- apply(rvs_seas2, c(1,2,4,5), function(x) if(sum(!is.na(x))) which.max(x) else as.numeric(NA))

## Empirical extremes
max22 <- apply(daily, c(1,2), function(x) max(x[1:(22*365)]))
max100 <- apply(daily, c(1,2), function(x) max(x[1:(100*365)]))
max1000 <- apply(daily, c(1,2), function(x) max(x[1:(1000*365)]))
max10000 <- apply(daily, c(1,2), max)

save(qs, returnPeriods, emp_qs,
    ns, failures, shapes, se_shapes, rvs, se_rvs,
    emp_quants,
    failures_gev, shapes_gev, se_shapes_gev, rvs_gev, se_rvs_gev,
    ns_seas1, failures_seas1, shapes_seas1, se_shapes_seas1, rvs_seas1, se_rvs_seas1, rvs_full_seas1, rvs_max_seas1, max_seas1,
    ns_seas2, failures_seas2, shapes_seas2, se_shapes_seas2, rvs_seas2, se_rvs_seas2, rvs_full_seas2, rvs_max_seas2, max_seas2,
    max22, max100, max1000, max10000,
    file = paste0(var, '_fits.Rda'))


