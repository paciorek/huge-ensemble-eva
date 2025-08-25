var <- commandArgs(trailingOnly=TRUE)[1]
reducer <- ifelse(var == "prec", mean, max)

k2c <- 273.15         # TMP2m is Kelvin
kgpm2ps <- 86400/10   # PRATEsfc is kg/m^2/sec; convert to cm/day

num_days <- c(rep(c(rep(365,3), 366), 5), 365, 365)

ndays <- sum(num_days)
nreps <- 40
nstart <- 12  # number of starting points

total_length <- nreps*ndays*nstart

dm <- c(70,30)
daily <- array(0, c(dm, total_length))

mth_lens <- c(31,28,31,30,31,30,31,31,30,31,30,31)

idx <- 1
for(mth in 1:12) {
    cat("Processing month ", mth, ".\n")
    print(system.time(load(paste0('preds_', mth, '.Rda')))) # 340 sec.

    if(var == 'prec') {
        data <- prec
    } else {
        data <- temp
    }
    rm(temp); rm(prec)
    ntime <- dim(data)[3]
    
    cutlen <- 365*4  
    if(mth > 1)
        cutlen <- cutlen - 4*sum(mth_lens[1:(mth-1)])

    for(i in 1:dm[1]) {
        cat("Processing ", i, ".\n")
        for(j in 1:dm[2]) {
            ## Remove initial (partial) year 2001.
            ## Do units conversion here to avoid making another larger object outside of looping.
            if(var == 'prec') {
                tmp <- data[i, j, (cutlen+1):ntime] * kgpm2ps 
            } else tmp <- data[i, j, (cutlen+1):ntime] - k2c
            if(length(tmp) != nreps*ndays*4)
                stop("incorrect number of time points")
            ## Reduce to daily and concatenate with previous months.
            daily[i,j, idx:(idx+nreps*ndays-1)] <- apply(matrix(tmp, nrow = 4), 2, reducer)
        }
    }
    
    idx <- idx + nreps*ndays
    cat("Next index is ", idx, ".\n")
    save(daily, file = paste0(var, '.Rda')) 
}

