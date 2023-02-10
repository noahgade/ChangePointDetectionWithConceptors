# Application Data Description
The data for the application study was obtained from the [CRCNS data shaing website](https://crcns.org/data-sets). We use Section 1 from hc-24 and cite the following works:

Varela, C. and M. A. Wilson (2019). Simultaneous extracellular recordings from midline thamalic neclei, median prefrontal cortex, and ca1 from rats cycling through bouts of sleep and wakefulness. Collaborative Research in Computational Neuroscience.

Varela, C. and M. A. Wilson (2020, June). mpfc spindle cycles organize sparse thalamic activation and recently active ca1 cells during non-rem sleep. *eLife 9*, e48881.

The folder contains:
- crcns_hc-24_data_description (description of the CRCNS dataser hc-24)
- Timestamp_LFP_600Hz.mat (timestamp of LFP measurements)
- HC_LFP_600Hz.mat (LFP measurements taken from the hippocampus)
- PFC_LFP_600Hz.mat (LFP measurements taken from the prefrontal cortex)
- THAL_LFP_600Hz.mat (LFP measurements taken from the midline thalamus)
- SleepTimes.mat (derived data from the authors indicating timestamps of sleep/wake transition)
- HCFilter24.m (code to combine the above .mat files and filter to the delta band at 1-4Hz)
- LFP_Section1_DeltaBand.csv (output .csv file from HCFilter24.m)
- SleepTimes.csv (output .csv file indicating the sleep times from SleepTimes.mat)
