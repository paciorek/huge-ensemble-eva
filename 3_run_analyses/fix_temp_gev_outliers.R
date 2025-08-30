library(climextRemes)
source("init.R")
load('temp_fits.Rda')
source("restrict_contig_US.R")
load('temp.Rda') # 520 sec

tmp1=restrict(rvs[,,3,5])
tmp2=restrict(rvs_gev[,,3])
wh <- which(tmp2-tmp1 > 5 | (is.na(tmp2) & matrix(inside_us, nlon, nlat)))
# [1]  726  876  936  937 1087 1131 1151 1152 1154 1221 1299 1300 1360 1430 1455
# [16] 1548 1566 1568 1569 1570 1598 1618 1619 1621 1639 1642 1657 1690 1691 1704
# [31] 1706 1707 1720 1722 1765 1776 1787 1788 1789 1790 1791 1811 1829 1835 1836
# [46] 1837 1840 1841 1854 1855 1856 1858 1859 1865 1866 1906 1908 1909 1920 1921
# [61] 1922 1925 1926 1970 1971 1972 1973 1974 1976 1978 1991 1992 1995 1996

badlon <- wh%%nlon
badlat <- 1+wh%/%nlon

nt <- dim(daily)[3]
last <- (nt-365+1):nt
## Indices of leap year days (Feb. 29).
leaps <- c(365*3+31+29, 365*6+366+31+29, 365*9+366*2+31+29, 365*12+366*3+31+29,
           365*15+366*4+31+29)


if(FALSE) {  # Check we got the bad ones.
    tmp2 <- rvs_gev[,,3]
    tmp1 <- rvs[,,3,5]
    tmp2[cbind(badlon,badlat)] - tmp1[cbind(badlon,badlat)]
}

rvs_gev_fix = rvs_gev
shapes_gev_fix = shapes_gev
se_rvs_gev_fix = se_rvs_gev
failures_gev_fix = failures_gev

for(i in seq_along(badlon)) {
    focal_lon <- badlon[i]
    focal_lat <- badlat[i]
    focal <- daily[focal_lon, focal_lat, ]

    quants <- quantile(focal, c(.99, emp_qs, qs))
    
    ## empirical
    emp_quants[focal_lon, focal_lat,] <- quants[2:3]
    
    ## POT - fit so can use as starting values.
    
    thresh <- quants[4:(length(qs)+3)]
    
    sub <- focal[focal > quants[1]]
    fit <- fit_pot(sub[sub > thresh[5]], threshold = thresh[5], nBlocks = nBlocks, returnPeriod = returnPeriods, getParams=TRUE)
    ## fit$returnValue
    
    blocks <- matrix(c(focal[last], focal[1:(last[1]-1)]), nrow = 365*17+366*5)
    yearly <- matrix(blocks[-leaps,], nrow = 365)
    
    mx <- apply(yearly, 2, max)
    fit_gev <- fit_gev(mx, returnPeriod = returnPeriods, getParams=TRUE)
    ## fit_gev$returnValue
    
    fit_gev2 <- fit_gev(mx, returnPeriod = returnPeriods, getParams=TRUE,
                        initial = as.list(fit$mle))
    ## fit_gev2$returnValue
    if((!is.null(fit_gev$nllh) && fit_gev$nllh < fit_gev2$nllh) ||
       (!is.null(fit_gev$nllh) && is.null(fit_gev2$nllh))) cat("Found worse LL: ", i, "\n")
    

    failures_gev_fix[focal_lon, focal_lat] <- fit_gev2$info$failure
    if(!fit_gev2$info$failure) {
        shapes_gev_fix[focal_lon, focal_lat] <- fit_gev2$mle['shape']
        rvs_gev_fix[focal_lon, focal_lat, ] <- fit_gev2$returnValue
        se_rvs_gev_fix[focal_lon, focal_lat, ] <- fit_gev2$se_returnValue
    }
 
}

save(failures_gev_fix, shapes_gev_fix, rvs_gev_fix, se_rvs_gev_fix,
     file = 'temp_fits_fixed.Rda')
