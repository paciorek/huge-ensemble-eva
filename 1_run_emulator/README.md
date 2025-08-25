This directory contains materials to use the ACE2-ERA5 emulator to produce a large ensemble, used by Chris Paciorek on the UC Berkeley Statistical Computing Facility GPU resources.

ACE2-ERA5 is well-documented in its [GitHub repository](https://github.com/ai2cm/ace), but I provide some additional details here, plus of course, the details of the ensemble that I created. 

## Ensemble details

The ensemble provides 40*22*12=10560 years of 6 hour, 1-degree resolution output, where:
  - runs are initialized from the initial values (in `initial_conditions`) for the first day of each of the 12 months of the year 2001
  - the emulator is run forward for the remainder of 2001 and then for the 21 years 2002-2021, using 2001-2022 forcing variable values (in `forcing_data`)
  - the emulation is repeated 39 times, with the forcing data recycled, so there is a jump in the forcing conditions from 2022-12-31 to 2001-01-01.
  - finally the emulation is continued for one additional year (i.e., 2001 forcings), which is used to replace the (partial) first year (2001 forcing conditions) so that regardless of the start date, we have 40 full replicates of the year 2001 (as well as all the other years of course).

For a regular set of output, one can then omit the output from the first (partial, except for the January 1 start) year of the first repeat. 

The output is put in the `output_22y_40m_${start_date}` directory (one directory for each of the 12 starting months. In particular `autoregressive_predictions.nc` contains the output.

To reduce output size given our analysis is only of temperature and precipitation, the output variables are only 2m air temperature (`TMP2m`) and surface precipitation (`PRATEsfc`).

### Emulator code

The emulator code is in the `fme` Python package. The installation instructions provided in the the ACE2-ERA5 GitHub repository give the following installation steps:

```bash
git clone https://github.com/ai2cm/ace
cd ace
make create_environment
```

I've saved a snapshot of the environment for reproducibility:

```bash
source activate fme
conda env export > fme_saved_january_2025.yml
```

The snapshot has a hard-coded path for the environment location and build hashes that might need to be omitted for the snapshot to be used on other systems.

We used the code on Ubuntu 22.04 (Linux).

### Emulator files

The files were obtained on 2025-01-08 from the [AI2 HuggingFace repository](https://huggingface.co/allenai/ACE2-ERA5/tree/main) from commit `a4ca6cc90f2bf277abdd8c80623936306f144fd0`.

Some of the files are also available on Zenodo (https://zenodo.org/records/10791087, https://zenodo.org/records/14606905)

The model "checkpoint" (i.e., parameter values/weights) are in `ace2_era5_ckpt.tar`.

Initial condition files should be placed in a subdirectory named `initial_conditions`. We used 2001, but some other years are available. One can also download others via the AI2 Google cloud bucket.

Forcing variable files for 2001-2022 should be placed in `forcing_data`. Forcing variable files for other years (not used for this ensemble) should be excluded from the `forcing_data` directory or else the emulator will run for the additional years represented by those forcing variable files.

## Running the emulator

See `run_emulator_one_start.sh` for a simple batch script that takes a single argument, the month for which you want to use the initial conditions file. 

That script then uses the configuration in `inference_config${month}.yaml`. Note that the `experiment_dir` field contains an absolute path for the location of the output files. Change this path and make the relevant directories before running the script.

So creating the full ensemble involves running

```bash
for (( i=1; i<=12; i++ )); do
    run_emulator_one_start.sh ${i}
done
```

The configuration information in `inference_config${month}.yaml` can be understood as follows.

 - The emulator will continue running through all the years in the `forcing_loader: dataset: data_path` directory, so only put data for the years of interest there.
 - The value of `n_forward_steps` should indicate the full number of 6-hour time steps to output, across all of the repeats. The value differs across the 12 months such that the last date in each output file should be for the final "replacement" replication of the year 2001 (see more details below).
 - The value of `forward_steps_in_memory` just relates to computational optimization, giving the number of 6-hour time steps kept in memory. I believe I just kept it at 50 without much experimentation. 

Note that with 22 years (including 5 leap years) there are 8035 days in the time interval and 4*8035=32140 time steps.

For January, `n_forward_steps` is 40*32140 + 365*4 = 1287060
For February, `n_forward_steps` is 40*32140 -31*4 + 365*4 = 1286936
...
For December, `n_forward_steps` is 40*32140 - (365-31)*4 + 365*4 = 1285724

where 365*4 is the additional "replacement" year 2001, and the subtraction accounts for the days missing in the first year 2001 when not starting on January 1.

### Speed and output size

Running the emulator for one of the 12 month starting values (for the 40*22 years) takes ~36 hours on a single A100 GPU. It takes about 5 seconds on an A100 to simulate 50 time steps. (For comparison, it takes about 31 seconds on an A5000.)

For each of the 12 months, the main output file (`autoregressive_predictions.nc`) is ~670 GB.

