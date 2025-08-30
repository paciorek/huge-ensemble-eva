
# Step 3a: Fit EVA models to ensemble output

- Input: `{prec,temp}.Rda` (30 GB and 28 GB respectively)
- Output: `{prec,temp}_fits.Rda` (5 MB)

```bash
## Emulator output (something like 4 hours)
Rscript run_all_fits.R prec
Rscript run_all_fits.R temp
```

The GEV fitting is not very robust -- optimization at many locations fails for temperature. To use better initial values, initialize at the estimates from the POT fits:

```{bash}
Rscript fix_temp_gev_outliers.R
```

- Input: `temp_fits.Rda`
- Output: `temp_fits_fixed.Rda`

This uses R package `climextRemes` v. 0.3.1.

# Step 3b: Fit EVA models to the ERA5 data

- Input: `era5.Rda` (470 MB)
- Output: `{prec,temp}_era5_fits.Rda` (500 KB)


```bash
Rscript run_era5_fits.R prec
Rscript run_era5_fits.R temp
```


This uses R package `climextRemes` v. 0.3.1.

# Step 3c: Make figures

```bash
Rscript make_figures.R
```

Figure files go to the `figures` directory.
