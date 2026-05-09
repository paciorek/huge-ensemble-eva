
results <- matrix(0,100,3)


for(i in 1:100){
dat <- rnorm(1e7)
mx1 <- matrix(dat, nrow= 1000)
mx2 <- matrix(dat, nrow=10000)
results[i,1]<-fit_gev(apply(mx1,2,max), returnPeriod = 10000)$returnValue
results[i,2]<-fit_gev(apply(mx2,2,max), returnPeriod = 10000)$returnValue
results[i,3]<-fit_gev(apply(mx2,2,max), returnPeriod = 10000/10)$returnValue
}


out1 <- out2 <- array(NA,c(70,30,4))
cnt <- 0
fits1 <- fits2 <- list()
length(fits1) <- length(fits2) <- 2100
for(focal_lon in seq_len(nlon)) {
   for(focal_lat in seq_len(nlat)) {
       cat("Working on latitude", focal_lat, "\n")
       cnt <- cnt+1
        focal <- daily[focal_lon, focal_lat, ]

        quants <- quantile(focal, c(.99, emp_qs, qs))

        ## empirical quantiles
        emp_quants[focal_lon, focal_lat,] <- quants[2:3]

        ## POT
        
        thresh <- quants[4:(length(qs)+3)]  # Thresholds for POT.

        sub <- focal[focal > quants[1]]  # Initial coarse thresholding to reduce computation.
        i <- 6
            fits[[cnt]] <- fit_pot(sub[sub > thresh[i]], threshold = thresh[i], nBlocks = nBlocks, returnPeriod = returnPeriods, getParams=TRUE)
        out1[focal_lon,focal_lat,] <- fits[[cnt]]$returnValue

        ## Shift full last replicate of year 2001 to replace partial initial year (previously discarded).
        blocks <- matrix(c(focal[last], focal[1:(last[1]-1)]), nrow = ndays)
       data <- c(blocks[-leaps,]) # Omit leap days for simplicity.
       i <- 4
            blockVals <- matrix(data, nrow = 365*blockLens[i]) 
            mx <- apply(blockVals, 2, max)
            fits2[[cnt]] <- fit_gev(mx, returnPeriod = returnPeriods/blockLens[i], getParams=TRUE)
       out2[focal_lon,focal_lat,] <- fits2[[cnt]]$returnValue
       
   }
print(focal_lon)
}

out0 <- array(NA,c(70,30,4))
cnt <- 0
fits0 <- list()
length(fits0) <- 2100
for(focal_lon in seq_len(nlon)) {
   for(focal_lat in seq_len(nlat)) {
       cat("Working on latitude", focal_lat, "\n")
       cnt <- cnt+1
        focal <- daily[focal_lon, focal_lat, ]

   
        ## Shift full last replicate of year 2001 to replace partial initial year (previously discarded).
        blocks <- matrix(c(focal[last], focal[1:(last[1]-1)]), nrow = ndays)
       data <- c(blocks[-leaps,]) # Omit leap days for simplicity.
       i <- 4
            blockVals <- matrix(data, nrow = 365*22)
            mx <- apply(blockVals, 2, max)
            fits0[[cnt]] <- fit_gev(mx, returnPeriod = returnPeriods/22, getParams=TRUE)
       out0[focal_lon,focal_lat,] <- fits0[[cnt]]$returnValue
       
   }
print(focal_lon)
}


scales <- array(as.numeric(NA), c(nlon, nlat, nthr))

for(focal_lon in 31:41) { # 25:35) {
    print(focal_lon)
    for(focal_lat in 20:30) { # 13:23) {
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
       scales[focal_lon, focal_lat, wh] <- sapply(fits[wh], function(x) x$mle['scale'])
    }}
data <- daily[31:41,20:30,] # 25:35,13:23,]
save(shapes,scales,data, file = '~/Desktop/tmp8.Rda')

exc <- which(focal > thresh[6])

fit1 <- fit_pot(focal[exc], threshold = thresh[i], nBlocks = nBlocks, returnPeriod = returnPeriods, getParams=TRUE)

sub2 <- climextRemes:::remove_runs(focal[exc],exc)
sub2 <- sub2[!is.na(sub2)]
fit2 <- fit_pot(sub2, threshold = thresh[i], nBlocks = nBlocks, returnPeriod = returnPeriods, getParams=TRUE)



