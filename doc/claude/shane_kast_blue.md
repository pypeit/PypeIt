# PypeIt Reduction with Claude: Shane Kast Blue Example

This document describes a workflow using Claude to reduce spectroscopic data from the Shane Kast Blue spectrograph using PypeIt.

## Initial Prompt

Here was the prompt given to Claude:

> In `/home/xavier/Projects/PypeIt/PypeIt-development-suite/RAW_DATA/shane_kast_blue/600_4310_d55` is a set of raw FITS files for the shane_kast_blue spectrograph. Please use the PypeIt code base to generate a PypeIt file. Do your work in `/home/xavier/Projects/PypeIt/PypeIt_Redux/Claude/shane_kast_blue`. If you run python or any scripts, use the "pypeit" environment in conda.

---

## Step 1: PypeIt File Generation

### Summary

Claude successfully generated a PypeIt file for the shane_kast_blue spectrograph data.

### Location

The PypeIt file is located at: `shane_kast_blue_A.pypeit`

### Configuration Details

- **Setup:** A (600/4310 grating with d55 dichroic)
- **Total files:** 25 raw FITS files

### File Classification

The setup identified:

| Frame Type | Files | Description |
|------------|-------|-------------|
| Arc | b1.fits.gz | For wavelength calibration |
| Bias | b14-b23.fits.gz | 10 bias frames |
| Flat | b3-b13.fits.gz | 10 frames (pixel flat, illumination flat, and trace) |
| Science | b27-b28.fits.gz | 2 frames (target: J1217p3905, 1200s exposures) |
| Standard | b24.fits.gz | Feige 66, 30s exposure |

---

## Step 2: Running the Reduction

### Prompt

> Please go ahead and execute the reduction

### Command

```bash
conda run -n pypeit run_pypeit shane_kast_blue_A.pypeit
```

### Result

**Reduction Complete!** The PypeIt reduction finished successfully in 1 minute 38 seconds.

---

## Reduction Output Summary

### Science Products (`Science/`)

#### 1D Extracted Spectra

- `spec1d_b24-Feige66_KASTb_20150520T041246.960.fits` - Standard star (Feige 66)
- `spec1d_b27-J1217p3905_KASTb_20150520T045733.560.fits` - Science target (first exposure)
- `spec1d_b28-J1217p3905_KASTb_20150520T051801.470.fits` - Science target (second exposure)

#### 2D Reduced Spectra

- `spec2d_b24-Feige66_KASTb_20150520T041246.960.fits` (27 MB)
- `spec2d_b27-J1217p3905_KASTb_20150520T045733.560.fits` (27 MB)
- `spec2d_b28-J1217p3905_KASTb_20150520T051801.470.fits` (27 MB)

### Calibration Files (`Calibrations/`)

| File | Description |
|------|-------------|
| `Arc_A_0_DET01.fits` | Arc line image for wavelength calibration |
| `Bias_A_0_DET01.fits` | Master bias frame |
| `Flat_A_0_DET01.fits` | Master flat field |
| `Edges_A_0_DET01.fits.gz` | Slit edge traces |
| `Slits_A_0_DET01.fits.gz` | Slit information |
| `Tilts_A_0_DET01.fits` | Spectral tilt solution |
| `WaveCalib_A_0_DET01.fits` | Wavelength calibration solution |

### Quality Assurance (`QA/`)

- HTML reports for each exposure (b24, b27, b28)
- Master Frame QA report (`MF_A.html`)
- PNG diagnostic plots in the `QA/PNGs/` directory

### Key Features

- Successfully extracted 1 object from each science frame
- Applied heliocentric correction (-22.46 km/s)
- Performed spectral flexure correction (-0.023 pixels)
- All spectra are wavelength calibrated and sky subtracted

You can view the QA HTML files in a browser to inspect the quality of the reduction.

---

## Step 3: Analysis Notebook

### Prompt

> Now generate a Jupyter Notebook describing the outputs and comment on the quality of the reduction.

### Result

Claude created a comprehensive Jupyter notebook that analyzes the PypeIt reduction outputs.

### Notebook Location

`shane_kast_blue_analysis.ipynb`

### Contents

#### 1. Introduction
- Overview of the observation details and configuration
- Summary of calibration frames used

#### 2. Calibration Quality Assessment
- Master bias frame analysis (statistics and visualization)
- Normalized pixel flat examination
- Slit edge detection results
- Wavelength calibration coverage

#### 3. Two-Dimensional Spectra
- Standard star (Feige 66) 2D spectrum display
- Science target (J1217+3905) 2D spectra for both exposures
- Sky model visualizations
- Spatial and spectral profile plots

#### 4. One-Dimensional Extracted Spectra
- Standard star 1D spectrum with S/N analysis
- Science target spectra for both exposures
- Comparison between exposures
- Sky subtraction quality assessment

#### 5. Quality Assessment Summary
- Automated quality metrics calculation
- Detailed assessment of:
  - Calibration quality (bias, flat, wavelength)
  - Science data quality (S/N ratios)
  - Sky subtraction performance
  - Consistency between exposures
- Overall reduction quality rating

#### 6. Conclusions and Next Steps
- Summary of reduction quality
- Recommended next steps (flux calibration, coaddition)
- List of ready-to-use data products

### To Run the Notebook

```bash
cd /home/xavier/Projects/PypeIt/PypeIt_Redux/Claude/shane_kast_blue
conda run -n pypeit jupyter notebook shane_kast_blue_analysis.ipynb
```

The notebook includes detailed visualizations and quantitative assessments that will help you understand the quality of the reduction and identify any potential issues. All analysis is ready to run and will provide comprehensive diagnostics of your Shane Kast Blue spectroscopic data reduction.
