# Step 2a: Process emulator output to subset to contiguous United States (CONUS)

This step subsets to a  box containing the contiguous United States by calling `subset.R` for each of the 12 runs.

Input: 12 `autoregressive_predictions.nc` files produced by the emulator (670 GB each)
Output: `preds_{1..12}.Rda` (15 GB each)

```bash
./subset.sh
```

# Step 2b: Process 6-hourly data to produce daily data for each variable

This step produces daily (total precipitation and maximum temperature) from 6-hourly data, and also combines across the 12 runs. This requires a machine with at least 80 GB memory.

Input: `preds_{1..12}.Rda` (15 GB each)
Output: `{prec,temp}.Rda` (30 GB and 28 GB respectively)


```bash
Rscript process_time.R prec
Rscript process_time.R temp
```

The output files are `{prec,temp}.Rda`.

# Step 2c: Download and process ERA5 data

Some of the analysis involves comparison with ERA5.

The data (i.e., the training data for ACE2-ERA5) are available at `https://console.cloud.google.com/storage/browser/ai2cm-public-requester-pays/2024-11-13-ai2-climate-emulator-v2-amip/data/era5-1deg-1940-2022.zarr`.

Download only the 2m air temperature and precipitation data (plus metadata), using a Google Cloud Platform account.

```bash
acct=<some_account>

mkdir era5-1deg-1940-2022.zarr
cd era5-1deg-1940-2022.zarr
gcloud auth login   # Authorize using a Google account linked to Google Cloud Platform with a payment method.

# Download only the variables needed.
gsutil -u ${acct} -m cp -r   "gs://ai2cm-public-requester-pays/2024-11-13-ai2-climate-emulator-v2-amip/data/era5-1deg-1940-2022.zarr/PRATEsfc"  ${scr}
gsutil -u ${acct} -m cp -r   "gs://ai2cm-public-requester-pays/2024-11-13-ai2-climate-emulator-v2-amip/data/era5-1deg-1940-2022.zarr/TMP2m"  .
gsutil -u ${acct} -m cp -r   "gs://ai2cm-public-requester-pays/2024-11-13-ai2-climate-emulator-v2-amip/data/era5-1deg-1940-2022.zarr/time"   .
gsutil -u ${acct} -m cp -r   "gs://ai2cm-public-requester-pays/2024-11-13-ai2-climate-emulator-v2-amip/data/era5-1deg-1940-2022.zarr/longitude"   .
gsutil -u ${acct} -m cp -r   "gs://ai2cm-public-requester-pays/2024-11-13-ai2-climate-emulator-v2-amip/data/era5-1deg-1940-2022.zarr/latitude"   .
gsutil -u ${acct} -m cp -r   "gs://ai2cm-public-requester-pays/2024-11-13-ai2-climate-emulator-v2-amip/data/era5-1deg-1940-2022.zarr/.zattrs"   .
gsutil -u ${acct} -m cp -r   "gs://ai2cm-public-requester-pays/2024-11-13-ai2-climate-emulator-v2-amip/data/era5-1deg-1940-2022.zarr/.zgroup"   .
gsutil -u ${acct} -m cp -r   "gs://ai2cm-public-requester-pays/2024-11-13-ai2-climate-emulator-v2-amip/data/era5-1deg-1940-2022.zarr/.zmetadata"   .
```

Subset to the contiguous US and produce daily (total precipitation and maximum temperature) from 6-hourly data.

```bash
Rscript subset_era5.R
```


# Software versions

The code was run using R 4.4.0 with packages `ncdf4` v. 1.22 and `stars` v. 0.6.8.
