library(stars)
## Loading required package: abind
## Loading required package: sf
## Linking to GEOS 3.10.2, GDAL 3.4.3, PROJ 8.2.0; sf_use_s2() is TRUE
dsn = 'era5-1deg-1940-2022.zarr'

lon <- read_mdim(dsn, variable = 'longitude')[[1]]
lat <- read_mdim(dsn, variable = 'latitude')[[1]]
attributes(lon) <- attributes(lat) <- NULL

lats <- which(lat >= 20 & lat <= 50)
lons <- which(lon >= (360-130) & lon <= (360-60))

offset <- c(lons[1]-1, lats[1]-1, 0)
count <- c(length(lons), length(lats), NA)

prec <- read_mdim(dsn, variable = 'PRATEsfc', offset = offset, count = count)[[1]]
temp <- read_mdim(dsn, variable = 'TMP2m', offset = offset, count = count)[[1]]

class(prec) <- class(temp) <- NULL

k2c <- 273.15         # TMP2m is Kelvin
kgpm2ps <- 86400/10   # PRATEsfc is kg/m^2/sec; convert to cm/day

time <- read_mdim(dsn, variable = 'time')[[1]]
time[1:5]
## 
## Units: [hours since 1940-01-01T12:00:00]
## [1]  0  6 12 18 24
## Since emulator predictions start at 00:00:00, this suggests we omit the first two
## time points to align such that first day is Jan 2, starting at midnight.
## Conveniently since 1940 is a leap year, this still leaves 365 days in the first year.

## 1940-2022
num_days <- c(rep(365,4),         # omit first two time points from Jan 1, 1940
    rep(c(366, rep(365,3)), 20),  # 1944-2019  # year 2000 was a leap year (divisible by 400)
    c(366, 365, 365))             # 2020-2022

ndays <- sum(num_days)

ntime <- dim(temp)[3]

dm <- c(70,30)
temp_era5 <- prec_era5 <- array(0, c(dm, (ntime-2)/4))

for(i in 1:dm[1]) {
    cat("Processing ", i, ".\n")
    for(j in 1:dm[2]) {
        ## Omit first two time points (starting at `3`).
        temp_station <- temp[i,j,3:ntime] - k2c
        prec_station <- prec[i,j,3:ntime] * kgpm2ps
        temp_era5[i,j,] <- apply(matrix(temp_station, nrow = 4), 2, max)
        prec_era5[i,j,] <- apply(matrix(prec_station, nrow = 4), 2, mean)
    }
}

save(prec_era5, temp_era5, file = file.path('..', 'run_analyses', 'era5.Rda'))


