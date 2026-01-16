load("coords.Rda")

nBlocks <- 40*22*12

nlon <- length(lon)
nlat <- length(lat)

ndays <- 365*17+366*5  # leap years in 2004,2008,2012,2016,2020


## Sequence of (daily) quantiles on log scale.
nthr <- 10
qs <- 1-exp(seq(log(.001), log(.00001), length = nthr))
qs <- c(.997262, qs)  # Include quantile giving sample size equal to annual maxima.
nthr <- nthr + 1

## Sequence of block lengths for GEV
blockLens <- c(1, 5, 10, 20, 40, 80, 160)
nlens <- length(blockLens)

returnPeriods <- c(1000,10000,100000,1000000)
nrp <- length(returnPeriods)

## quantiles of daily values corresponding to 1000- and 10000-year return values
emp_qs <- 1-c(1/(1000*365), 1/(10000*365))

