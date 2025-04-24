# Bayesian-Optimized Locally Stationary Factor Model

**Author:** Moka Kaleji • Master’s Thesis, University of Bologna • Contact: mohammadkaleji1998@gmail.com

**Affiliation:** Master’s Thesis in Econometrics & Time Series Analysis  

This repository implements hyperparameter tuning via Bayesian optimization for the  
Locally Stationary Factor Model (LSFM), plus out-of-sample forecasting using the  
optimized factors.

---

## Contents

- **lsfmestim_big4.m**  
  - Selects monthly or quarterly data  
  - Splits into training/validation sets  
  - Uses `bayesopt` to choose `(q, h, p)` that minimize MSFE on GDP, Unemp., Inflation, Int. Rate  
  - Produces:  
    - `lsfm_estimation_results.mat` containing `q_opt`, `h_opt`, `p_opt`, `Fhat_opt`, `Lhat_opt`, plus scaling info  
    - Diagnostic plots: variance explained, factor correlation heatmap, rolling correlations  

- **lsfmforecasting_big4.m**  
  - Loads the optimized results  
  - Uses `forecast_LSFM` helper to generate H-step forecasts on the four key indicators  
  - Computes MSFE, RMSE in both normalized and original units  
  - Visualizes: MSFE by horizon, actual vs. forecast, error histograms, factor series, cross-correlograms  

---

## Quickstart

1. **Requirements**  
   - MATLAB R2020b or later  
   - Statistics & Machine Learning Toolbox  
   - Econometrics Toolbox  

2. **Data**  
   Place your processed Excel files here:

3. **Hyperparameter Tuning**  
```matlab
cd path/to/lsfm-bo-factor-model
lsfm_bo_estimation
Optimal: q=..., h=..., p=...
Min CV-MSFE: ...
lsfmforecasting_big4



## See Also
- [Preprocessing-FRED-Quarterly-Dataset](https://github.com/moka-kaleji/Preprocessing-FRED-Quarterly-Dataset)
- [Locally-Stationary-Factor-Model](https://github.com/moka-kaleji/Locally-Stationary-Factor-Model)
- [Dynamic-Factor-Model](https://github.com/moka-kaleji/Dynamic-Factor-Model)
