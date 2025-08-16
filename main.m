clear all; clc; close all;

%% settings
settings = calibrationSettings();

%% Input data (volatility surfaces)
% struct| data
% Fields| T (time to maturity), K (strike price), r (interest rate), S0
% (underlying price)	
% Load data
load("empVolatilitySurfaceData.mat"); % data
                                                        
% Initial parameters 
parameters0 = settings.parameters0;

% Calibrate model
settings.standardErrors = false;
g = @(parameters) pricingError(data, settings, parameters);
fInitial = g(parameters0);

[parametersFinal, errorFinal, exitFlag] = fminsearch(g, parameters0, settings.calibrOptions);



% Extract final parameters
kappaFinal = parametersFinal(1);
thetaFinal = parametersFinal(2);
etaFinal = parametersFinal(3);
rhoFinal = parametersFinal(4);
V0Final = parametersFinal(5);

% Calculate standard errors for the parameter estimates ((Hessian-based)
n = size(data.IVolSurf, 1) * size(data.IVolSurf, 2); % Number of option contracts
f = g(parametersFinal) * n;

settings.standardErrors = true;
fis = @(x) pricingError(data, settings, x);

% Degrees of freedom (number of parameters)
p = length(parameters0);

% Jacobian matrix
%J = jacobianest(fis, parametersFinal);
% Covariance matrix
%sigma2 = f / (n - p);
%Sigma = sigma2 * inv(J' * J);
% Parameter standard errors
%se = sqrt(diag(Sigma))';



% Numerical Hessian estimation of squared error function
hessFunc = @(params) mean((pricingError(data, settings, params)).^2);

% Compute numerical Hessian at the optimum
H = hessianest(hessFunc, parametersFinal);

% Error variance from in-sample fit
sigma2 = f / (n - p);

% Covariance matrix of parameters
Sigma = sigma2 * inv(H);

% Standard errors and t-values
se = sqrt(diag(Sigma))';
t_values = parametersFinal' ./ se;

% Display results
disp('Calibration Results (Hessian-based):');
disp(['Estimated values: ', num2str(parametersFinal')]);
disp(['t-values: ', num2str(t_values)]);
disp(['Standard errors: ', num2str(se)]);
disp(['In-Sample MSE: ', num2str(errorFinal)]);



%% Price Exotic Option: Down-and-In Asian Call (Arithmetic Avg Strike)
% Model-based Monte Carlo pricing using calibrated parameters

H = 0.85;          % Down barrier
M = 100000;        % Number of simulations
T = 1;             % Maturity time (1 year)
dt = 1/252;        % Daily time steps
S0 = data.S0;      % Current stock price
r = data.r;        % Risk-free rate

% Exotic option pricing

exoticPrice = priceDownInAsianAvgStrikeCall(parametersFinal, H, M, T, dt, S0, r);

disp(['Price of down-and-in arithmetic Asian average strike call: ', num2str(exoticPrice)]);





