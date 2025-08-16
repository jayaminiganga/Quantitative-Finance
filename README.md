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
