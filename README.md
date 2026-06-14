# In‑situ Soil Organic Carbon Prediction with Vis‑NIR Spectroscopy

[![Python](https://img.shields.io/badge/Python-3.8+-blue)](https://python.org)
[![LightGBM](https://img.shields.io/badge/LightGBM-4.x-green)](https://lightgbm.readthedocs.io/)
[![NGBoost](https://img.shields.io/badge/NGBoost-Probabilistic-orange)](https://github.com/stanfordmlgroup/ngboost)
[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)

**End‑to‑end machine learning pipeline for SOC estimation** | Subterra project | Vis‑NIR (330–2580 nm) | PHENET

## 📂 Repository contents

| File | Description |
|------|-------------|
| `Code_4_absorbance_data_process_visualiazation&model.ipynb` | Main notebook — all steps below |
| `Data_builder.R` | Builds H5 file — compiles experiment and spectra data |
| `lab_data_builder.R` | Alternative builder for lab‑based measurements |
| `H5file_useful_commands.R` | Cheat sheet — useful commands for working with H5 files |
| `Useful_functions.R` | Library of R functions — imported during data building and preprocessing |
| `Extract_data_from_h5_file.R` | Extracts, filters, trims, interpolates, and smooths data from H5 file |
| `RGB_features_extraction.R` | Extracts features from RGB images (supplements spectra analysis) |

## 🧪 Pipeline steps in detail

### 1. Data loading
- Google Drive integration
- Load `AGES_absorbances.csv` into DataFrame
- Round depth to integer

### 2. Spectral filtering & SID outlier removal
- Extract spectral columns (regex `^[SF]\d`)
- Normalize spectra to sum = 1
- Compute SID against mean spectrum
- Remove outliers (SID > 0.02) → **207 spectra removed**

### 3. Savitzky–Golay smoothing
- Optimize window size using lag‑1 autocorrelation of residuals
- **Optimal window:** 11, polyorder = 2

### 4. Target variable corrections (NT_III insertions)
- Load Excel with corrected TC and SIC values (insertions 6,8,9,10)
- Recalculate: **SOC = TC - SIC**

### 5. Feature engineering
- Correlation heatmaps (TC, SOC, Water, SIC, Force, Depth)
- Derived features: `vis_swir_ratio`, `swir_std`

### 6. Depth profiles
- Mean SOC ± SEM across depths
- Facet plots per plot name

### 7. GPS coordinate conversion
- Parse `POINT(Longitude Latitude)` format
- Transform **WGS84 → UTM zone 33N** (pyproj)
- Output: `X_m`, `Y_m`

### 8. Graph‑based & Bayesian smoothing (SOC)
- **Graph smoothing:** kNN on 3D coordinates + depth
- **Bayesian smoothing:** LOWESS global profile + local variance
- Output: `SOC_sm` (smoothed target)

### 9. PCA for plot discrimination
- 50 components on spectra
- ANOVA F‑test → top discriminative PCs: **PC8, PC3**
- Interactive Plotly 3D scatter plot

### 10. EPO (External Parameter Orthogonalization)
- Load wet/dry lab spectra (aligned by Base_ID)
- Compute difference matrix: `D = wet - dry`
- Wilks lambda → **1 component sufficient**
- Apply correction to full spectral data

### 11. Train/test split
- 20% test / 80% train (random state 23)
- **Target:** `SOC_sm`
- **Features:** 
  - PCA of visible spectra (50 PCs)
  - `vis_swir_ratio`
  - `swir_std`
  - PCA of SNV‑detrended SWIR (50 PCs)

### 12. Synthetic data generation (augmentation)
- Extend PC scores beyond original min/max
- Inverse‑transform to spectral space
- Combine synthetic + real data

### 13. Model training

| Model | RMSE | R² | RPD | Note |
|-------|------|----|-----|------|
| **LightGBM** | 0.36 | 0.28 | 1.18 | Baseline |
| **NGBoost** | 0.37 | 0.25 | — | Probabilistic |

**NGBoost uncertainty metrics:**
- 95% prediction interval coverage: **84%**
- CRPS: **0.20** (lower than MAE → well calibrated)

### 14. Model persistence
- Save/load NGBoost model with `joblib`



## About Subterra
In-situ soil properties estimation and digital soil mapping with Vis-NIR spectroscopy.

Sequestering carbon in the soil has a large potential in mitigating climate change. Monitoring soil carbon stocks is however, very labour expensive as the soil samples have to be collected and sent to the lab and analysed. This Use Case aims to analyse soil organic carbon (SOC) stocks without taking soil samples to the lab. Thus, instead of taking physical soil samples we are taking proximal spectral information of the soil at different soil depths using a probe designed for in situ soil sensing. Based on the spectral pattern we try to estimate the soil carbon concentration. 

The Subterra Green which has been developed by our partners S4 Mobile Laboratories is able to insert a metal rod connected with optical fibers to a Visible and NIR spectrometer, down to 90 cm soil depths, taking spectrum (~330-2580nm) every cm. Currently, we are doing soil measurements by collecting and analysing soil samples for ground truthing. Furthermore, using machine learning approach we can teach a model to estimate the SOC content based on a spectral signal.


## Citation

```bibtex
@misc{Baykalov2025scripts,
  Author = {Pavel Baykalov},
  Title = {Subterra-probe-PHENET-project},
  Year = {2025},
  Publisher = {GitHub},
  Journal = {GitHub repository},
  Howpublished = {\url{https://https://github.com/dIcarusb/Subterra-probe-PHENET-project}}
}
```
