# Quantitative-Finance
MATLAB scripts for quantitative finance analysis

### Project Overview

This project implements a Heston stochastic volatility model calibration and exotic option pricing workflow. It includes:

### Model Calibration:

Calibrates Heston parameters (κ, θ, η, ρ, V₀) to market volatility surfaces using FFT-based pricing

Minimizes MSE between model and market implied volatilities

Provides parameter standard errors via Hessian-based estimation

### Exotic Option Pricing:

Prices a down-and-in arithmetic Asian average strike call using Monte Carlo simulation

Implements Milstein scheme with truncation to prevent negative variances

Uses antithetic variates for variance reduction 

### Repository Structure

├── main.m                         # Main execution script
├── calibrationSettings.m          # Parameter bounds & optimization settings
├── pricingError.m                 # Calibration loss function
├── CallPricingFFT.m               # FFT-based option pricing
├── priceDownInAsianAvgStrikeCall.m # Monte Carlo pricing for exotic option
├── CharacteristicFunctionLib.m    # Characteristic functions for models
├── hessianest.m                   # Hessian matrix estimator
├── empVolatilitySurfaceData.mat   # Market implied volatility data
└── README.md                      # This documentation

### Requirements

MATLAB (R2020a or newer)

Financial Toolbox (for blsimpv function)

Statistics and Machine Learning Toolbox

### How to Run

1. Clone the repository:
   git clone https://github.com/yourusername/heston-calibration-exotic-pricing.git

2. Open MATLAB and navigate to the project directory
   Run the main script:
   run main.m

   The workflow executes:

      Loads market volatility data (empVolatilitySurfaceData.mat)

      Calibrates Heston parameters to market surface

      Computes standard errors for parameter estimates

      Prices exotic option using calibrated parameters

      Displays calibration results and option price


### Key Features

Robust calibration with parameter constraints

Hessian-based uncertainty quantification

Antithetic variates for Monte Carlo variance reduction

Barrier condition tracking during path simulation

Visual comparison of model vs market volatility surfaces

