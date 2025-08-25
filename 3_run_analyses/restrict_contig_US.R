library(rnaturalearth)  # Version 1.0.1
library(sf)  # Version 1.0.19

load('coords.Rda') # provides lon/lat

# Get US states boundaries
us_states <- ne_states(country = "United States of America", returnclass = "sf")

# Create the contiguous US boundary (excluding Alaska, Hawaii, and territories)
contiguous_us <- us_states[!(us_states$name %in% c("Alaska", "Hawaii")) & 
                          !grepl("Puerto Rico|Virgin Islands|Guam", us_states$name), ]

# Combine all state polygons into a single multipolygon
contiguous_us_combined <- st_union(contiguous_us)

grid_points <- expand.grid(lon = lon, lat = lat)

# Convert to sf object
points_sf <- st_as_sf(grid_points, coords = c("lon", "lat"), crs = 4326)

# Determine which points are inside the contiguous US
inside_us <- st_intersects(points_sf, contiguous_us_combined, sparse = FALSE)

restrict <- function(x) {
    x[!inside_us] <- NA
    return(x)
}
