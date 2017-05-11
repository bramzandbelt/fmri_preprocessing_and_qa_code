These directories contain some MATLAB code I wrote for preprocessing and quality assessment of fMRI data. Here's a brief overview of how the code is organized and what it does:

```
├── fmri_multiecho_preprocessing  <- Dir with code for fMRI preprocessing
│   ├── coreg_normalize_smooth.m  <- SPM12 job for coregistration, normalization, and smoothing
│   ├── dicom_import.m            <- SPM12 job for DICOM to NIfTI conversion
│   ├── preproc_mri_stpy.m        <- Preprocesses functional and anatomical MRI data
│   ├── realign_estimate.m        <- SPM12 job for realignment (estimation only)
│   ├── realign_reslice.m         <- SPM12 job for realignment (estimation & reslicing)
│   └── run_preproc_mri_stpy.m    <- Main script, calling preproc_mri_stpy.m
├── quality_check                 <- Dir with code for fMRI quality analysis
│   ├── bztbx_automask.m          <- Creates a rough mask of the brain (based on code from ArtRepair toolbox)
│   ├── bztbx_bpfilt.m            <- Bandpass filters fMRI time series
│   ├── bztbx_createmask.m        <- Creates binary images (masks) of ROIs
│   ├── bztbx_extractts.m         <- Extracts fMRI time series data
│   ├── bztbx_fmrifig.m           <- Makes a figure with quality check results
│   ├── bztbx_multireg.m          <- Makes a text file with multiple regressors (e.g. realigment parameters)
│   ├── bztbx_qa.m                <- Main quality check function
│   ├── bztbx_truncrp.m           <- Truncates realignment parameters
│   └── bztbx_tsnr.m              <- Computes temporal signal to noise ratio of fMRI time series
├── README.md                     <- The file you're reading now
```
