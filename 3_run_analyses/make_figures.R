plot_dir <- file.path('..', 'figures')

library(ggplot2)
library(gridExtra)
library(fields)
library(dplyr)
library(grid)

source("restrict_contig_US.R")

var <- "prec"

if(var == "prec") {
    units <- "cm"
    varLabel <- "precipitation"
} else {
    units <- "deg. C"
    varLabel <- "temperature"
}

load(paste0(var, "_fits.Rda"))

source("init.R")

qs_round <- round(qs, 5)
labs <-  paste0("thresh = ", qs_round, "%")
n_labs <- paste0("n = ", ns[1,1,])

omi <- c(0,0,.3,0)

par(omi = omi)

# Compare emulator maxes with ERA5 maxes.

load('era5.Rda')

if(var == "prec") {
    era5_max <- apply(prec_era5, c(1,2), max)
} else era5_max <- apply(temp_era5, c(1,2), max)               

if(var == "prec") lim <- c(0,35) else lim <- c(25,60)

# Prepare data for ggplot
data_100 <- data.frame(era5_max = as.vector(restrict(era5_max)), max100 = as.vector(restrict(max100)))
data_1000 <- data.frame(era5_max = as.vector(restrict(era5_max)), max1000 = as.vector(restrict(max1000)))
data_10000 <- data.frame(era5_max = as.vector(restrict(era5_max)), max10000 = as.vector(restrict(max10000)))

# Create individual plots
plot_100 <- ggplot(data_100, aes(x = era5_max, y = max100)) +
    geom_point(size = 0.65) +
    geom_abline(slope = 1, intercept = 0, col = "red") +
    xlab("") + 
    ylab(paste0("max (", units, ") in 100 emulator years")) +
    xlim(lim) + ylim(lim) +
    labs(title = "(a) 100 years", fill = "") +
    theme_minimal()

plot_1000 <- ggplot(data_1000, aes(x = era5_max, y = max1000)) +
    geom_point(size = 0.65) +
    geom_abline(slope = 1, intercept = 0, col = "red") +
    xlab(paste0("max (", units, ") in ERA5 (1940-2022)")) +
    ylab(paste0("max (", units, ") in 1000 emulator years")) +
    xlim(lim) + ylim(lim) +
    labs(title = "(b) 1000 years", fill = "") +
    theme_minimal()

plot_10000 <- ggplot(data_10000, aes(x = era5_max, y = max10000)) +
    geom_point(size = 0.65) +
    geom_abline(slope = 1, intercept = 0, col = "red") +
    xlab("") +
    ylab(paste0("max (", units, ") in 10560 emulator years")) +
    xlim(lim) + ylim(lim) +
    labs(title = "(c) 10560 years", fill = "") +
    theme_minimal()

full_plot <- grid.arrange(plot_100, plot_1000, plot_10000, ncol = 3)

ggsave(filename = file.path(plot_dir, paste0(var, "-compare-era5.pdf")), 
    plot = full_plot, height = 4, width = 10)

plot_100 <- plot_100 + ggtitle("100 emulator years")
plot_1000 <- plot_1000 + ggtitle("1000 emulator years")
plot_10000 <- plot_10000 + ggtitle("10560 emulator years")
full_plot <- grid.arrange(plot_100, plot_1000, plot_10000, ncol = 3)

ggsave(filename = file.path(plot_dir, paste0(var, "-compare-era5.png")), 
    plot = full_plot, height = 1200, width = 3000, units = "px")

# Create maps of extremes climatology.

# Prepare data for ggplot
data_combined <- data.frame(
    lon = rep(lon, length(lat)) - 360,
    lat = rep(lat, each = length(lon)),
    rvs_1000 = as.vector(restrict(rvs[,,1,5])),
    rvs_100000 = as.vector(restrict(rvs[,,3,5])),
    shapes = as.vector(restrict(shapes[,,5]))
)

data_combined <- data_combined |> filter(lon > -130 & lon < -65 & lat > 25 & lat < 50)

# Adjust legend size and position, and increase plot dimensions
plot_rvs_1000 <- ggplot(data_combined, aes(x = lon, y = lat)) +
    geom_tile(aes(fill = rvs_1000)) +
    scale_fill_viridis_c(limits = c(3, 42), na.value = "white") +
    borders("state", colour = "grey", size = 0.25) +
    labs(title = "(a) 1-in-1000 year AEP depth estimates", fill = "") +
    theme_minimal() +
    theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 12), # Reduced text size
        legend.key.size = unit(0.7, "cm"), # Reduced legend key size
        legend.position = "right",
        plot.title = element_text(size = 15) # Adjusted title size
    )

plot_rvs_100000 <- ggplot(data_combined, aes(x = lon, y = lat)) +
    geom_tile(aes(fill = rvs_100000)) +
    scale_fill_viridis_c(limits = c(3, 42),na.value = "white") +
    borders("state", colour = "grey", size = 0.25) +
    labs(title = "(b) 1-in-100000 year AEP depth estimates", fill = "") +
    theme_minimal() +
    theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 12), # Reduced text size
        legend.key.size = unit(0.7, "cm"), # Reduced legend key size
        legend.position = "right",
        plot.title = element_text(size = 15) # Adjusted title size
    )

plot_shapes <- ggplot(data_combined, aes(x = lon, y = lat)) +
    geom_tile(aes(fill = shapes)) +
    scale_fill_gradient2(low = "#D55E00", mid = "white", high = "#0072B2", midpoint = 0, na.value = "white") +
    borders("state", colour = "grey", size = 0.25) +
    labs(title = "(c) shape parameter estimates", fill = "") +
    theme_minimal() +
    theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 12), # Reduced text size
        legend.key.size = unit(0.7, "cm"), # Reduced legend key size
        legend.position = "right",
        plot.title = element_text(size = 15) # Adjusted title size
    )

full_plot <- grid.arrange(plot_rvs_1000, plot_rvs_100000, plot_shapes, ncol = 3)

ggsave(filename = file.path(plot_dir, paste0(var, "-climatology.pdf")), 
    plot = full_plot, height = 3, width = 18)

ggsave(filename = file.path(plot_dir, paste0(var, "-climatology.png")), 
    plot = full_plot, height = 900, width = 5400, units = "px")

## Comparison with empirical quantiles

if(var == "prec") ylim <- c(0,60) else ylim <- c(25,60)

# Prepare data for ggplot
data_combined <- data.frame(
    emp_quants_1 = as.vector(restrict(emp_quants[,,1])),
    emp_quants_2 = as.vector(restrict(emp_quants[,,2])),
    rvs_gev_1 = as.vector(restrict(rvs_gev[,,1])),
    rvs_gev_2 = as.vector(restrict(rvs_gev[,,2]))
)

for(i in c(1,3,5,7)) {
    data_combined[[paste0("rvs_pot_1_", i)]] <- as.vector(restrict(rvs[,,1,i]))
    data_combined[[paste0("rvs_pot_2_", i)]] <- as.vector(restrict(rvs[,,2,i]))
}

# Create individual plots for T=10^3
plot_gev_1 <- ggplot(data_combined, aes(x = emp_quants_1, y = rvs_gev_1)) +
    geom_point(size = 0.65) +
    geom_abline(slope = 1, intercept = 0, col = 'red') +
    labs(x = "", y = "GEV (annual maxima)") +
    ylim(ylim) +
    theme_minimal()

plot_pot_1 <- lapply(c(1,3,5,7), function(i) {
    ggplot(data_combined, aes(x = emp_quants_1, y = .data[[paste0("rvs_pot_1_", i)]]))+
        geom_point(size = 0.65) +
        geom_abline(slope = 1, intercept = 0, col = 'red') +
        labs(x = ifelse(i==3, "empirical quantile-based AEP depth estimates", ""),
            y = paste0("POT, ", n_labs[i])) +
        ylim(ylim) +
        theme_minimal()
})

# Create individual plots for T=10^4
plot_gev_2 <- ggplot(data_combined, aes(x = emp_quants_2, y = rvs_gev_2)) +
    geom_point(size = 0.65) +
    geom_abline(slope = 1, intercept = 0, col = "red") +
    labs(x = "", y = "GEV (annual maxima)") +
    ylim(ylim) +
    theme_minimal()

plot_pot_2 <- lapply(c(1,3,5,7), function(i) {
    ggplot(data_combined, aes(x = emp_quants_2, y = .data[[paste0("rvs_pot_2_", i)]]))+
        geom_point(size = 0.65) +
        geom_abline(slope = 1, intercept = 0, col = "red") +
        labs(x = ifelse(i==3, "empirical quantile-based AEP depth estimates", ""),
            y = paste0("POT, ", n_labs[i])) +
        ylim(ylim) +
        theme_minimal()
})

row1_title <- textGrob("1-in-1000 year AEP depth estimates (cm)", gp = gpar(fontsize = 14, fontface = "bold"))
row2_title <- textGrob("1-in-10000 year AEP depth estimates (cm)", gp = gpar(fontsize = 14, fontface = "bold"))

# Arrange all plots with titles above each row.
full_plot <- grid.arrange(
    row1_title, plot_gev_1, plot_pot_1[[1]], plot_pot_1[[2]], plot_pot_1[[3]], plot_pot_1[[4]],
    row2_title, plot_gev_2, plot_pot_2[[1]], plot_pot_2[[2]], plot_pot_2[[3]], plot_pot_2[[4]],
    ncol = 5,
    layout_matrix = rbind(
        c(1, 1, 1, 1, 1),
        c(2, 3, 4, 5, 6),
        c(7, 7, 7, 7, 7),
        c(8, 9, 10, 11, 12)
    ),
    heights = c(0.5, 4, 0.5, 4) # Adjust the relative heights of rows
)

ggsave(filename = file.path(plot_dir, paste0(var, "-rv-bias.pdf")), plot = full_plot, height = 8, width = 12)

ggsave(filename = file.path(plot_dir, paste0(var, "-rv-bias.png")), plot = full_plot, 
    height = 1600, width = 2400, units = "px")

## Shape parameter and AEP depth stability

# Prepare data for ggplot
data_histograms <- data.frame()
for(i in c(1,3,5,7,9)) {
    values <- as.vector(restrict(shapes[,,i]))
    values <- values[!is.na(values)]
    pct_negative <- 100 * round(mean(values < 0, na.rm = TRUE), 2)
    
    temp_data <- data.frame(
        values = values,
        label = paste(n_labs[i], ", ", pct_negative, "% < 0"),
        index = i
    )
    data_histograms <- rbind(data_histograms, temp_data)
}

# Convert label to factor with levels in the order they appear
data_histograms$label <- factor(data_histograms$label, levels = unique(data_histograms$label))

# Create the histogram plots
lim <- c(-.55, .55)

plot_histograms <- ggplot(data_histograms, aes(x = values)) +
    geom_histogram(bins = 30, fill = "grey", color = "black", alpha = 0.7) +
    facet_wrap(~ label, ncol = 5, scales = "free_y") +
    xlim(lim) + 
    labs(x = "shape parameter estimate",
         y = "frequency") +
    theme_minimal() +
    theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_line(color = "gray95"),
        plot.title = element_text(size = 16, hjust = 0.5),
        strip.text = element_text(size = 10),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 12)
    )

ggsave(filename = file.path(plot_dir, paste0(var, "-shape-hists.pdf")), 
       plot = plot_histograms, height = 3, width = 10)

ggsave(filename = file.path(plot_dir, paste0(var, "-shape-hists.png")), 
       plot = plot_histograms, height = 900, width = 3000, units = "px")




if(var == "prec") {
    ylim <- c(-.5, .45)
} else {
    ylim <- c(-.5, .5)
}

# Prepare data for ggplot
data_stability_shapes <- data.frame(
    shapes_5 = as.vector(restrict(shapes[,,5])),
    shapes_1 = as.vector(restrict(shapes[,,1])),
    shapes_3 = as.vector(restrict(shapes[,,3])),
    shapes_6 = as.vector(restrict(shapes[,,6])),
    shapes_8 = as.vector(restrict(shapes[,,8])),
    shapes_9 = as.vector(restrict(shapes[,,9]))
)

data_stability_rvs <- data.frame(
    rvs_5 = as.vector(restrict(rvs[,,3,5])),
    rvs_1 = as.vector(restrict(rvs[,,3,1])),
    rvs_3 = as.vector(restrict(rvs[,,3,3])),
    rvs_6 = as.vector(restrict(rvs[,,3,6])),
    rvs_8 = as.vector(restrict(rvs[,,3,8])),
    rvs_9 = as.vector(restrict(rvs[,,3,9]))
)

# Create individual plots for shape parameter stability
plot_shapes <- lapply(c(1, 3, 6, 8, 9), function(i) {
    ggplot(data_stability_shapes, aes(x = shapes_5, y = .data[[paste0("shapes_", i)]])) +
        geom_point(size = 0.5) +
        geom_abline(slope = 1, intercept = 0, col = "red") +
        labs(
            x = ifelse(i==6, paste0("shape parameter estimate for ", n_labs[5]), ""),
            y = n_labs[i]
        ) +
        ylim(ylim) +
        theme_minimal()
})

# Create individual plots for AEP depth stability
plot_rvs <- lapply(c(1, 3, 6, 8, 9), function(i) {
    ggplot(data_stability_rvs, aes(x = rvs_5, y = .data[[paste0("rvs_", i)]])) +
        geom_point(size = 0.5) +
        geom_abline(slope = 1, intercept = 0, col = "red") +
        labs(
             x = ifelse(i==6, paste0("AEP depth estimate for ", n_labs[5]), ""),
            y = n_labs[i]
        ) +
        theme_minimal()
})

row1_title <- textGrob("Shape parameter estimates", gp = gpar(fontsize = 14, fontface = "bold"))
row2_title <- textGrob("1-in-100000 year AEP depth estimates (cm)", gp = gpar(fontsize = 14, fontface = "bold"))

# Arrange all plots with titles above each row.
full_plot <- grid.arrange(
    row1_title, plot_shapes[[1]], plot_shapes[[2]], plot_shapes[[3]], plot_shapes[[4]], plot_shapes[[5]],
    row2_title, plot_rvs[[1]], plot_rvs[[2]], plot_rvs[[3]], plot_rvs[[4]], plot_rvs[[5]],
    ncol = 5,
    layout_matrix = rbind(
        c(1, 1, 1, 1, 1),
        c(2, 3, 4, 5, 6),
        c(7, 7, 7, 7, 7),
        c(8, 9, 10, 11, 12)
    ),
    heights = c(0.5, 4, 0.5, 4) # Adjust the relative heights of rows
)

ggsave(filename = file.path(plot_dir, paste0(var, "-stability.pdf")), plot = full_plot, height = 6, width = 10)

ggsave(filename = file.path(plot_dir, paste0(var, "-stability.png")), plot = full_plot, 
    height = 1800, width = 3000, units = "px")

# Alternative version that plots difference from n=499 case.


# Prepare data for ggplot
data_stability_shapes_alt <- data.frame(
    shapes_5 = as.vector(restrict(shapes[,,5])),
    shapes_1_diff = as.vector(restrict(shapes[,,1]) - restrict(shapes[,,5])),
    shapes_3_diff = as.vector(restrict(shapes[,,3]) - restrict(shapes[,,5])),
    shapes_6_diff = as.vector(restrict(shapes[,,6]) - restrict(shapes[,,5])),
    shapes_8_diff = as.vector(restrict(shapes[,,8]) - restrict(shapes[,,5])),
    shapes_9_diff = as.vector(restrict(shapes[,,9]) - restrict(shapes[,,5]))
)

data_stability_rvs_alt <- data.frame(
    rvs_5 = as.vector(restrict(rvs[,,3,5])),
    rvs_1_diff = as.vector(restrict(rvs[,,3,1]) - restrict(rvs[,,3,5])),
    rvs_3_diff = as.vector(restrict(rvs[,,3,3]) - restrict(rvs[,,3,5])),
    rvs_6_diff = as.vector(restrict(rvs[,,3,6]) - restrict(rvs[,,3,5])),
    rvs_8_diff = as.vector(restrict(rvs[,,3,8]) - restrict(rvs[,,3,5])),
    rvs_9_diff = as.vector(restrict(rvs[,,3,9]) - restrict(rvs[,,3,5]))
)

# Create individual plots for shape parameter stability
if(var == "prec") {
    ylim <- c(-.5, .5)
} else {
    ylim <- c(-.5, .5)  # Modify for temperature.
}

plot_shapes <- lapply(c(1, 3, 6, 8, 9), function(i) {
    ggplot(data_stability_shapes_alt, aes(x = shapes_5, y = .data[[paste0("shapes_", i, "_diff")]])) +
        geom_point(size = 0.5) +
        labs(
            x = ifelse(i==6, paste0("shape parameter estimate for ", n_labs[5]), ""),
            y = paste0("difference for ", n_labs[i])
        ) +
        ylim(ylim) +
        theme_minimal()
})

# Create individual plots for AEP depth stability
if(var == "prec") {
    ylim <- c(-20, 20)
} else {
    ylim <- c(-.5, .5)  # Modify for temperature.
}

plot_rvs <- lapply(c(1, 3, 6, 8, 9), function(i) {
    ggplot(data_stability_rvs_alt, aes(x = rvs_5, y = .data[[paste0("rvs_", i, "_diff")]])) +
        geom_point(size = 0.5) +
        labs(
             x = ifelse(i==6, paste0("AEP depth estimate for ", n_labs[5]), ""),
            y = paste0("difference for ", n_labs[i])
        ) +
        ylim(ylim) +
        theme_minimal()
})

row1_title <- textGrob("Shape parameter estimates", gp = gpar(fontsize = 14, fontface = "bold"))
row2_title <- textGrob("1-in-100000 year AEP depth estimates (cm)", gp = gpar(fontsize = 14, fontface = "bold"))

# Arrange all plots with titles above each row.
full_plot <- grid.arrange(
    row1_title, plot_shapes[[1]], plot_shapes[[2]], plot_shapes[[3]], plot_shapes[[4]], plot_shapes[[5]],
    row2_title, plot_rvs[[1]], plot_rvs[[2]], plot_rvs[[3]], plot_rvs[[4]], plot_rvs[[5]],
    ncol = 5,
    layout_matrix = rbind(
        c(1, 1, 1, 1, 1),
        c(2, 3, 4, 5, 6),
        c(7, 7, 7, 7, 7),
        c(8, 9, 10, 11, 12)
    ),
    heights = c(0.5, 4, 0.5, 4) # Adjust the relative heights of rows
)

ggsave(filename = file.path(plot_dir, paste0(var, "-stability-alt.pdf")), plot = full_plot, height = 6, width = 10)

ggsave(filename = file.path(plot_dir, paste0(var, "-stability-alt.png")), plot = full_plot, 
    height = 1800, width = 3000, units = "px")


## Seasonality

if(var == "prec") {
    xlim <- c(0,45)
    ylim <- c(0,60)
} else {
    xlim <- c(25,50)
    ylim <- c(25,50)
}

# Prepare data for ggplot
data_seasonal <- lapply(c(1, 3, 5, 7), function(i) {
    data.frame(
        full_year_rv = as.vector(restrict(rvs[,,3,i])),
        max_seasonal_rv1 = as.vector(restrict(rvs_max_seas1[,,3,i])),
        max_seasonal_rv2 = as.vector(restrict(rvs_max_seas2[,,3,i])),
        label = n_labs[i]
    )
})

# Version 1 of seasonal analysis - use same thresholds as full year analysis

# Create individual plots for each label
plot_seasonal <- lapply(seq_along(data_seasonal), function(i) {
    ggplot(data_seasonal[[i]], aes(x = full_year_rv, y = max_seasonal_rv1)) +
        geom_point(size = 0.75) +
        geom_abline(slope = 1, intercept = 0, col = "red") +
        labs(
            x = paste0("full-year AEP depth (", units, ")"),
            y = ifelse(i==1, paste0("max seasonal AEP depth (", units, ")"), ""),
            title = unique(data_seasonal[[i]]$label)
        ) +
        xlim(xlim) +
        ylim(ylim) +
        theme_minimal()
})

# Arrange all plots in a single row
full_plot <- grid.arrange(
    grobs = plot_seasonal,
    ncol = 4
)


ggsave(filename = file.path(plot_dir, paste0(var, "-rv-seasonal1.pdf")), plot = full_plot, height = 4, width = 9)

ggsave(filename = file.path(plot_dir, paste0(var, "-rv-seasonal1.png")), plot = full_plot, 
height = 1200, width = 2700, units = "px")

# Version 2 of seasonal analysis - use same number of exceedances in each season as in full year

# Create individual plots for each label
plot_seasonal <- lapply(seq_along(data_seasonal), function(i) {
    ggplot(data_seasonal[[i]], aes(x = full_year_rv, y = max_seasonal_rv2)) +
        geom_point(size = 0.75) +
        geom_abline(slope = 1, intercept = 0, col = "red") +
        labs(
            x = paste0("full-year AEP depth (", units, ")"),
            y = ifelse(i==1, paste0("max seasonal AEP depth (", units, ")"), ""),
            title = unique(data_seasonal[[i]]$label)
        ) +
        xlim(xlim) +
        ylim(ylim) +
        theme_minimal()
})

# Arrange all plots in a single row
full_plot <- grid.arrange(
    grobs = plot_seasonal,
    ncol = 4
)

ggsave(filename = file.path(plot_dir, paste0(var, "-rv-seasonal2.pdf")), plot = full_plot, height = 4, width = 9)

ggsave(filename = file.path(plot_dir, paste0(var, "-rv-seasonal2.png")), plot = full_plot, 
height = 1200, width = 2700, units = "px")


## Uncertainty in AEP depth estimates

emp_cv <- se_rvs / rvs

xlim <- c(-.3, .3)
if(var == "prec") ylim <- c(0, .35) else ylim <- c(0,.15)

# Prepare data for ggplot
data_uncertainty <- data.frame(
    shapes_5 = as.vector(restrict(shapes[,,5])),
    emp_cv_10000 = as.vector(restrict(emp_cv[,,2,5])),
    emp_cv_100000 = as.vector(restrict(emp_cv[,,3,5])),
    emp_cv_million = as.vector(restrict(emp_cv[,,4,5]))
)

# Create individual plots for each AEP depth uncertainty
plot_10000 <- ggplot(data_uncertainty, aes(x = shapes_5, y = emp_cv_10000)) +
    geom_point(size = 0.5) +
    xlab("") +
    ylab("AEP depth relative uncertainty") +
    xlim(xlim) +
    ylim(ylim) +
    ggtitle("(a) 1-in-10000 year") +
    theme_minimal()

plot_100000 <- ggplot(data_uncertainty, aes(x = shapes_5, y = emp_cv_100000)) +
    geom_point(size = 0.5) +
    xlab("shape parameter estimate") +
    ylab("") +
    xlim(xlim) +
    ylim(ylim) +
    ggtitle("(b) 1-in-100000 year") +
    theme_minimal()

plot_million <- ggplot(data_uncertainty, aes(x = shapes_5, y = emp_cv_million)) +
    geom_point(size = 0.5) +
    xlab("") +
    ylab("") +
    xlim(xlim) +
    ylim(ylim) +
    ggtitle("(c) 1-in-million year") +
    theme_minimal()

# Arrange all plots in a single row
full_plot <- grid.arrange(plot_10000, plot_100000, plot_million, ncol = 3)

ggsave(filename = file.path(plot_dir, paste0(var, "-rv-uncertainty.pdf")), plot = full_plot, 
    height = 4, width = 8)

ggsave(filename = file.path(plot_dir, paste0(var, "-rv-uncertainty.png")), plot = full_plot, 
    height = 1200, width = 2400, units = "px")

## Temperature results

var <- "temp"

if(var == "prec") {
    units <- "cm"
    varLabel <- "precipitation"
} else {
    units <- "degrees C"
    varLabel <- "temperature"
}

load(paste0(var, "_fits.Rda"))

load('era5.Rda')

if(var == "prec") {
    era5_max <- apply(prec_era5, c(1,2), max)
} else era5_max <- apply(temp_era5, c(1,2), max)               

if(var == "prec") lim <- c(0,35) else lim <- c(25,55)

# Prepare data for ggplot
data_100 <- data.frame(era5_max = as.vector(restrict(era5_max)), max100 = as.vector(restrict(max100)))
data_1000 <- data.frame(era5_max = as.vector(restrict(era5_max)), max1000 = as.vector(restrict(max1000)))
data_10000 <- data.frame(era5_max = as.vector(restrict(era5_max)), max10000 = as.vector(restrict(max10000)))

# Create individual plots
plot_100 <- ggplot(data_100, aes(x = era5_max, y = max100)) +
    geom_point(size = 0.65) +
    geom_abline(slope = 1, intercept = 0, col = "red") +
    xlab("") +
    ylab(paste0("max (", units, ") in 100 emulator years")) +
    xlim(lim) + ylim(lim) +
    labs(title = "(a) 100 years", fill = "") +
    theme_minimal()

plot_1000 <- ggplot(data_1000, aes(x = era5_max, y = max1000)) +
    geom_point(size = 0.65) +
    geom_abline(slope = 1, intercept = 0, col = "red") +
    xlab(paste0("max (", units, ") in ERA5 (1940-2022)")) +
    ylab(paste0("max (", units, ") in 1000 emulator years")) +
    xlim(lim) + ylim(lim) +
    labs(title = "(b) 1000 years", fill = "") +
    theme_minimal()

plot_10000 <- ggplot(data_10000, aes(x = era5_max, y = max10000)) +
    geom_point(size = 0.65) +
    geom_abline(slope = 1, intercept = 0, col = "red") +
    xlab("") +
    ylab(paste0("max (", units, ") in 10560 emulator years")) +
    xlim(lim) + ylim(lim) +
    labs(title = "(c) 10560 years", fill = "") +
    theme_minimal()

full_plot <- grid.arrange(plot_100, plot_1000, plot_10000, ncol = 3)

ggsave(filename = file.path(plot_dir, paste0(var, "-compare-era5.pdf")), 
    plot = full_plot, height = 4, width = 10)

ggsave(filename = file.path(plot_dir, paste0(var, "-compare-era5.png")), 
    plot = full_plot, height = 1200, width = 3000, units = "px")

# Create maps of extremes climatology.

# Prepare data for ggplot
data_combined <- data.frame(
    lon = rep(lon, length(lat)) - 360,
    lat = rep(lat, each = length(lon)),
    rvs_1000 = as.vector(restrict(rvs[,,1,5])),
    rvs_100000 = as.vector(restrict(rvs[,,3,5])),
    shapes = as.vector(restrict(shapes[,,5]))
)

data_combined <- data_combined |> filter(lon > -130 & lon < -65 & lat > 25 & lat < 50)

# Adjust legend size and position, and increase plot dimensions
plot_rvs_1000 <- ggplot(data_combined, aes(x = lon, y = lat)) +
    geom_tile(aes(fill = rvs_1000)) +
    scale_fill_viridis_c(limits = c(26, 52), na.value = "white") +
    borders("state", colour = "grey", size = 0.25) +
    labs(title = "(a) 1-in-1000 year AEP depth estimates", fill = "") +
    theme_minimal() +
    theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 12), # Reduced text size
        legend.key.size = unit(0.7, "cm"), # Reduced legend key size
        legend.position = "right",
        plot.title = element_text(size = 15) # Adjusted title size
    )

plot_rvs_100000 <- ggplot(data_combined, aes(x = lon, y = lat)) +
    geom_tile(aes(fill = rvs_100000)) +
    scale_fill_viridis_c(limits = c(26, 52), na.value = "white") +
    borders("state", colour = "grey", size = 0.25) +
    labs(title = "(b) 1-in-100000 year AEP depth estimates", fill = "") +
    theme_minimal() +
    theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 12), # Reduced text size
        legend.key.size = unit(0.7, "cm"), # Reduced legend key size
        legend.position = "right",
        plot.title = element_text(size = 15) # Adjusted title size
    )

plot_shapes <- ggplot(data_combined, aes(x = lon, y = lat)) +
    geom_tile(aes(fill = shapes)) +
    scale_fill_gradient2(low = "#D55E00", mid = "white", high = "#0072B2", midpoint = 0, na.value = "white") +
    borders("state", colour = "grey", size = 0.25) +
    labs(title = "(c) shape parameter estimates", fill = "") +
    theme_minimal() +
    theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 12), # Reduced text size
        legend.key.size = unit(0.7, "cm"), # Reduced legend key size
        legend.position = "right",
        plot.title = element_text(size = 15) # Adjusted title size
    )

full_plot <- grid.arrange(plot_rvs_1000, plot_rvs_100000, plot_shapes, ncol = 3)

ggsave(filename = file.path(plot_dir, paste0(var, "-climatology.pdf")), 
    plot = full_plot, height = 3, width = 18)

ggsave(filename = file.path(plot_dir, paste0(var, "-climatology.png")), 
    plot = full_plot, height = 900, width = 5400, units = "px")

if(var == "prec") ylim <- c(0,60) else ylim <- c(25,60)


# Prepare data for ggplot
data_combined <- data.frame(
    emp_quants_1 = as.vector(restrict(emp_quants[,,1])),
    emp_quants_2 = as.vector(restrict(emp_quants[,,2])),
    rvs_gev_1 = as.vector(restrict(rvs_gev[,,1])),
    rvs_gev_2 = as.vector(restrict(rvs_gev[,,2]))
)

for(i in c(1,3,5,7)) {
    data_combined[[paste0("rvs_pot_1_", i)]] <- as.vector(restrict(rvs[,,1,i]))
    data_combined[[paste0("rvs_pot_2_", i)]] <- as.vector(restrict(rvs[,,2,i]))
}

plot_gev_1 <- ggplot(data_combined, aes(x = emp_quants_1, y = rvs_gev_1)) +
    geom_point(size = 0.65) +
    geom_abline(slope = 1, intercept = 0, col = 'red') +
    labs(x = "", y = "GEV (annual maxima)") +
    ylim(ylim) +
    theme_minimal()

plot_pot_1 <- lapply(c(1,3,5,7), function(i) {
    ggplot(data_combined, aes(x = emp_quants_1, y = .data[[paste0("rvs_pot_1_", i)]]))+
        geom_point(size = 0.65) +
        geom_abline(slope = 1, intercept = 0, col = 'red') +
        labs(x = ifelse(i==3, "empirical quantile-based AEP depth estimates", ""),
            y = paste0("POT, ", n_labs[i])) +
        ylim(ylim) +
        theme_minimal()
})

# Create individual plots for T=10^4
plot_gev_2 <- ggplot(data_combined, aes(x = emp_quants_2, y = rvs_gev_2)) +
    geom_point(size = 0.65) +
    geom_abline(slope = 1, intercept = 0, col = "red") +
    labs(x = "", y = "GEV (annual maxima)") +
    ylim(ylim) +
    theme_minimal()

plot_pot_2 <- lapply(c(1,3,5,7), function(i) {
    ggplot(data_combined, aes(x = emp_quants_2, y = .data[[paste0("rvs_pot_2_", i)]]))+
        geom_point(size = 0.65) +
        geom_abline(slope = 1, intercept = 0, col = "red") +
        labs(x = ifelse(i==3, "empirical quantile-based AEP depth estimates", ""),
            y = paste0("POT, ", n_labs[i])) +
        ylim(ylim) +
        theme_minimal()
})

row1_title <- textGrob("1-in-1000 year AEP depth estimates (cm)", gp = gpar(fontsize = 14, fontface = "bold"))
row2_title <- textGrob("1-in-10000 year AEP depth estimates (cm)", gp = gpar(fontsize = 14, fontface = "bold"))

# Arrange all plots with titles above each row.
full_plot <- grid.arrange(
    row1_title, plot_gev_1, plot_pot_1[[1]], plot_pot_1[[2]], plot_pot_1[[3]], plot_pot_1[[4]],
    row2_title, plot_gev_2, plot_pot_2[[1]], plot_pot_2[[2]], plot_pot_2[[3]], plot_pot_2[[4]],
    ncol = 5,
    layout_matrix = rbind(
        c(1, 1, 1, 1, 1),
        c(2, 3, 4, 5, 6),
        c(7, 7, 7, 7, 7),
        c(8, 9, 10, 11, 12)
    ),
    heights = c(0.5, 4, 0.5, 4) # Adjust the relative heights of rows
)

ggsave(filename = file.path(plot_dir, paste0(var, "-rv-bias.pdf")), plot = full_plot, height = 8, width = 12)

ggsave(filename = file.path(plot_dir, paste0(var, "-rv-bias.png")), plot = full_plot, 
    height = 1600, width = 2400, units = "px")

## Shape parameter and AEP depth stability

# Prepare data for ggplot
data_histograms <- data.frame()
for(i in c(1,3,5,7,9)) {
    values <- as.vector(restrict(shapes[,,i]))
    values <- values[!is.na(values)]
    pct_negative <- 100 * round(mean(values < 0, na.rm = TRUE), 2)
    
    temp_data <- data.frame(
        values = values,
        label = paste(n_labs[i], ", ", pct_negative, "% < 0"),
        index = i
    )
    data_histograms <- rbind(data_histograms, temp_data)
}

# Convert label to factor with levels in the order they appear
data_histograms$label <- factor(data_histograms$label, levels = unique(data_histograms$label))

# Create the histogram plots
lim <- c(-.55, .55)

plot_histograms <- ggplot(data_histograms, aes(x = values)) +
    geom_histogram(bins = 30, fill = "grey", color = "black", alpha = 0.7) +
    facet_wrap(~ label, ncol = 5, scales = "free_y") +
    xlim(lim) + 
    labs(x = "shape parameter estimate",
         y = "frequency") +
    theme_minimal() +
    theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_line(color = "gray95"),
        plot.title = element_text(size = 16, hjust = 0.5),
        strip.text = element_text(size = 10),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 12)
    )

ggsave(filename = file.path(plot_dir, paste0(var, "-shape-hists.pdf")), 
       plot = plot_histograms, height = 3, width = 10)

ggsave(filename = file.path(plot_dir, paste0(var, "-shape-hists.png")), 
       plot = plot_histograms, height = 900, width = 3000, units = "px")




if(var == "prec") {
    ylim <- c(-.5, .45)
} else {
    ylim <- c(-.5, .5)
}

# Prepare data for ggplot
data_stability_shapes <- data.frame(
    shapes_5 = as.vector(restrict(shapes[,,5])),
    shapes_1 = as.vector(restrict(shapes[,,1])),
    shapes_3 = as.vector(restrict(shapes[,,3])),
    shapes_6 = as.vector(restrict(shapes[,,6])),
    shapes_8 = as.vector(restrict(shapes[,,8])),
    shapes_9 = as.vector(restrict(shapes[,,9]))
)

data_stability_rvs <- data.frame(
    rvs_5 = as.vector(restrict(rvs[,,3,5])),
    rvs_1 = as.vector(restrict(rvs[,,3,1])),
    rvs_3 = as.vector(restrict(rvs[,,3,3])),
    rvs_6 = as.vector(restrict(rvs[,,3,6])),
    rvs_8 = as.vector(restrict(rvs[,,3,8])),
    rvs_9 = as.vector(restrict(rvs[,,3,9]))
)

# Create individual plots for shape parameter stability
plot_shapes <- lapply(c(1, 3, 6, 8, 9), function(i) {
    ggplot(data_stability_shapes, aes(x = shapes_5, y = .data[[paste0("shapes_", i)]])) +
        geom_point(size = 0.5) +
        geom_abline(slope = 1, intercept = 0, col = "red") +
        labs(
            x = ifelse(i==6, paste0("shape parameter estimate for ", n_labs[5]), ""),
            y = n_labs[i]
        ) +
        ylim(ylim) +
        theme_minimal()
})

# Create individual plots for AEP depth stability
plot_rvs <- lapply(c(1, 3, 6, 8, 9), function(i) {
    ggplot(data_stability_rvs, aes(x = rvs_5, y = .data[[paste0("rvs_", i)]])) +
        geom_point(size = 0.5) +
        geom_abline(slope = 1, intercept = 0, col = "red") +
        labs(
             x = ifelse(i==6, paste0("AEP depth estimate for ", n_labs[5]), ""),
            y = n_labs[i]
        ) +
        theme_minimal()
})

row1_title <- textGrob("Shape parameter estimates", gp = gpar(fontsize = 14, fontface = "bold"))
row2_title <- textGrob("1-in-100000 year AEP depth estimates (cm)", gp = gpar(fontsize = 14, fontface = "bold"))

# Arrange all plots with titles above each row.
full_plot <- grid.arrange(
    row1_title, plot_shapes[[1]], plot_shapes[[2]], plot_shapes[[3]], plot_shapes[[4]], plot_shapes[[5]],
    row2_title, plot_rvs[[1]], plot_rvs[[2]], plot_rvs[[3]], plot_rvs[[4]], plot_rvs[[5]],
    ncol = 5,
    layout_matrix = rbind(
        c(1, 1, 1, 1, 1),
        c(2, 3, 4, 5, 6),
        c(7, 7, 7, 7, 7),
        c(8, 9, 10, 11, 12)
    ),
    heights = c(0.5, 4, 0.5, 4) # Adjust the relative heights of rows
)

ggsave(filename = file.path(plot_dir, paste0(var, "-stability.pdf")), plot = full_plot, height = 6, width = 10)

ggsave(filename = file.path(plot_dir, paste0(var, "-stability.png")), plot = full_plot, 
    height = 1800, width = 3000, units = "px")

# Alternative version that plots difference from n=499 case.


# Prepare data for ggplot
data_stability_shapes_alt <- data.frame(
    shapes_5 = as.vector(restrict(shapes[,,5])),
    shapes_1_diff = as.vector(restrict(shapes[,,1]) - restrict(shapes[,,5])),
    shapes_3_diff = as.vector(restrict(shapes[,,3]) - restrict(shapes[,,5])),
    shapes_6_diff = as.vector(restrict(shapes[,,6]) - restrict(shapes[,,5])),
    shapes_8_diff = as.vector(restrict(shapes[,,8]) - restrict(shapes[,,5])),
    shapes_9_diff = as.vector(restrict(shapes[,,9]) - restrict(shapes[,,5]))
)

data_stability_rvs_alt <- data.frame(
    rvs_5 = as.vector(restrict(rvs[,,3,5])),
    rvs_1_diff = as.vector(restrict(rvs[,,3,1]) - restrict(rvs[,,3,5])),
    rvs_3_diff = as.vector(restrict(rvs[,,3,3]) - restrict(rvs[,,3,5])),
    rvs_6_diff = as.vector(restrict(rvs[,,3,6]) - restrict(rvs[,,3,5])),
    rvs_8_diff = as.vector(restrict(rvs[,,3,8]) - restrict(rvs[,,3,5])),
    rvs_9_diff = as.vector(restrict(rvs[,,3,9]) - restrict(rvs[,,3,5]))
)

# Create individual plots for shape parameter stability
if(var == "prec") {
    ylim <- c(-.5, .5)
} else {
    ylim <- c(-.5, .5)  # Modify for temperature.
}

plot_shapes <- lapply(c(1, 3, 6, 8, 9), function(i) {
    ggplot(data_stability_shapes_alt, aes(x = shapes_5, y = .data[[paste0("shapes_", i, "_diff")]])) +
        geom_point(size = 0.5) +
        labs(
            x = ifelse(i==6, paste0("shape parameter estimate for ", n_labs[5]), ""),
            y = paste0("difference for ", n_labs[i])
        ) +
        ylim(ylim) +
        theme_minimal()
})

# Create individual plots for AEP depth stability
if(var == "prec") {
    ylim <- c(-20, 20)
} else {
    ylim <- c(-.5, .5)  # Modify for temperature.
}

plot_rvs <- lapply(c(1, 3, 6, 8, 9), function(i) {
    ggplot(data_stability_rvs_alt, aes(x = rvs_5, y = .data[[paste0("rvs_", i, "_diff")]])) +
        geom_point(size = 0.5) +
        labs(
             x = ifelse(i==6, paste0("AEP depth estimate for ", n_labs[5]), ""),
            y = paste0("difference for ", n_labs[i])
        ) +
        ylim(ylim) +
        theme_minimal()
})

row1_title <- textGrob("Shape parameter estimates", gp = gpar(fontsize = 14, fontface = "bold"))
row2_title <- textGrob("1-in-100000 year AEP depth estimates (cm)", gp = gpar(fontsize = 14, fontface = "bold"))

# Arrange all plots with titles above each row.
full_plot <- grid.arrange(
    row1_title, plot_shapes[[1]], plot_shapes[[2]], plot_shapes[[3]], plot_shapes[[4]], plot_shapes[[5]],
    row2_title, plot_rvs[[1]], plot_rvs[[2]], plot_rvs[[3]], plot_rvs[[4]], plot_rvs[[5]],
    ncol = 5,
    layout_matrix = rbind(
        c(1, 1, 1, 1, 1),
        c(2, 3, 4, 5, 6),
        c(7, 7, 7, 7, 7),
        c(8, 9, 10, 11, 12)
    ),
    heights = c(0.5, 4, 0.5, 4) # Adjust the relative heights of rows
)

ggsave(filename = file.path(plot_dir, paste0(var, "-stability-alt.pdf")), plot = full_plot, height = 6, width = 10)

ggsave(filename = file.path(plot_dir, paste0(var, "-stability-alt.png")), plot = full_plot, 
    height = 1800, width = 3000, units = "px")


## Seasonality

if(var == "prec") {
    xlim <- c(0,45)
    ylim <- c(0,60)
} else {
    xlim <- c(25,50)
    ylim <- c(25,50)
}

# Prepare data for ggplot
data_seasonal <- lapply(c(1, 3, 5, 7), function(i) {
    data.frame(
        full_year_rv = as.vector(restrict(rvs[,,3,i])),
        max_seasonal_rv1 = as.vector(restrict(rvs_max_seas1[,,3,i])),
        max_seasonal_rv2 = as.vector(restrict(rvs_max_seas2[,,3,i])),
        label = n_labs[i]
    )
})

# Create individual plots for each label
plot_seasonal <- lapply(seq_along(data_seasonal), function(i) {
    ggplot(data_seasonal[[i]], aes(x = full_year_rv, y = max_seasonal_rv1)) +
        geom_point(size = 0.75) +
        geom_abline(slope = 1, intercept = 0, col = "red") +
        labs(
            x = paste0("full-year AEP depth (", units, ")"),
            y = ifelse(i==1, paste0("max seasonal AEP depth (", units, ")"), ""),
            title = unique(data_seasonal[[i]]$label)
        ) +
        xlim(xlim) +
        ylim(ylim) +
        theme_minimal()
})

# Arrange all plots in a single row
full_plot <- grid.arrange(
    grobs = plot_seasonal,
    ncol = 4
)


ggsave(filename = file.path(plot_dir, paste0(var, "-rv-seasonal1.pdf")), plot = full_plot, height = 4, width = 9)

ggsave(filename = file.path(plot_dir, paste0(var, "-rv-seasonal1.png")), plot = full_plot, 
        height = 1200, width = 2700, units = "px")

# Create individual plots for each label
plot_seasonal <- lapply(seq_along(data_seasonal), function(i) {
    ggplot(data_seasonal[[i]], aes(x = full_year_rv, y = max_seasonal_rv2)) +
        geom_point(size = 0.75) +
        geom_abline(slope = 1, intercept = 0, col = "red") +
        labs(
            x = paste0("full-year AEP depth (", units, ")"),
            y = ifelse(i==1, paste0("max seasonal AEP depth (", units, ")"), ""),
            title = unique(data_seasonal[[i]]$label)
        ) +
        xlim(xlim) +
        ylim(ylim) +
        theme_minimal()
})

# Arrange all plots in a single row
full_plot <- grid.arrange(
    grobs = plot_seasonal,
    ncol = 4
)


ggsave(filename = file.path(plot_dir, paste0(var, "-rv-seasonal2.pdf")), plot = full_plot, height = 4, width = 9)

ggsave(filename = file.path(plot_dir, paste0(var, "-rv-seasonal2.png")), plot = full_plot, 
        height = 1200, width = 2700, units = "px")

## Uncertainty in AEP depth estimates

emp_cv <- se_rvs / rvs

xlim <- c(-.3, .3)
if(var == "prec") ylim <- c(0, .35) else ylim <- c(0,.15)

# Prepare data for ggplot
data_uncertainty <- data.frame(
    shapes_5 = as.vector(restrict(shapes[,,5])),
    emp_cv_10000 = as.vector(restrict(emp_cv[,,2,5])),
    emp_cv_100000 = as.vector(restrict(emp_cv[,,3,5])),
    emp_cv_million = as.vector(restrict(emp_cv[,,4,5]))
)

# Create individual plots for each AEP depth uncertainty
plot_10000 <- ggplot(data_uncertainty, aes(x = shapes_5, y = emp_cv_10000)) +
    geom_point(size = 0.5) +
    xlab("") +
    ylab("AEP depth relative uncertainty") +
    xlim(xlim) +
    ylim(ylim) +
    ggtitle("(a) 10000-year") +
    theme_minimal()

plot_100000 <- ggplot(data_uncertainty, aes(x = shapes_5, y = emp_cv_100000)) +
    geom_point(size = 0.5) +
    xlab("shape estimate") +
    ylab("") +
    xlim(xlim) +
    ylim(ylim) +
    ggtitle("(b) 100000-year") +
    theme_minimal()

plot_million <- ggplot(data_uncertainty, aes(x = shapes_5, y = emp_cv_million)) +
    geom_point(size = 0.5) +
    xlab("") +
    ylab("") +
    xlim(xlim) +
    ylim(ylim) +
    ggtitle("(c) million-year") +
    theme_minimal()

# Arrange all plots in a single row
full_plot <- grid.arrange(plot_10000, plot_100000, plot_million, ncol = 3)


ggsave(filename = file.path(plot_dir, paste0(var, "-rv-uncertainty.pdf")), plot = full_plot, height = 4, width = 8)

ggsave(filename = file.path(plot_dir, paste0(var, "-rv-uncertainty.png")), plot = full_plot, 
    height = 1200, width = 2400, units = "px")
