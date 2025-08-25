# huge-ensemble-eva
Materials for analysis of huge ML-based climate emulator ensemble output using statistical extreme value analysis

These are the steps involved in reproducing the analysis.

# Step 1: Run emulator

See `1_run_emulator/README.md` for details on running the emulator, including downloading the trained ACE2-ERA5 model (the "checkpoint" containing the trained model weights), the initial conditions files, and the forcing data files. 

# Step 2: Process output

See `2_process_output/README.md` for details on processing the output to reduce to a bounding box around the contiguous United States. 

This also includes details on processing the ERA5 training data for comparison with emulator output.

# Step 3: Do analysis and generate figures

See `3_run_analyses/README.md` for details on fitting the data using EVA (threshold exceedances [POT] and block maxima [GEV]) and making figures.
