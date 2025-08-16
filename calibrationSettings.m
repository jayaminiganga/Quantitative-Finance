function settings = calibrationSettings
% In this function, the settings of the calibration problem are specified

%% GENERAL SETTINGS
% Number of days in a year
settings.tradingDaysInYear = 252;

%% CLOSED-FORM SOLUTION -SETTINGS
settings.n = 13;
settings.model = 'Heston';

%% INITIAL PARAMETER VALUES
%kappa = 8;
%theta = 0.25^2;
%eta = 0.8;
%rho = -0.8;
%V0 = 0.3^2;
%Current Loss: 0.0010126
%Calibration Results (Hessian-based):
%Estimated values: 9.8111    0.014818     0.91109    -0.97919    0.093466
%t-values: 6551.7102      9.5556252      430.18449     -538.16644        49.0389
%Standard errors: 0.0014975   0.0015507   0.0021179   0.0018195   0.0019059
%In-Sample MSE: 0.0010127

%kappa = 4;
%theta = 0.3^2;
%eta = 0.8;
%rho = -0.8;
%V0 = 0.18^2;
%Current Loss: 5.8186e-05
%Calibration Results (Hessian-based):
%Estimated values: 3.7738    0.024098     0.76813    -0.57664    0.014904
%t-values: 0.013282    0.062761    0.020006   -0.073137    0.012332
%Standard errors: 284.1328     0.3839732      38.39552      7.884431       1.20858
%In-Sample MSE: 5.8186e-05

kappa = 6;
theta = 0.3^2;
eta = 0.8;
rho = -0.8;
V0 = 0.18^2;
%Current Loss: 7.6505e-05
%Calibration Results (Hessian-based):
%Estimated values: 8.0281    0.024335       1.361    -0.56784   0.0026683
%t-values: 344.7689      337.1184      361.5268     -356.1226      466.7994
%Standard errors: 0.023286  7.2186e-05   0.0037645   0.0015945  5.7161e-06
%In-Sample MSE: 7.6505e-05
%Price of down-and-in arithmetic Asian average strike call: 0.003679

%kappa = 2;
%theta = 0.04;
%eta = 0.3;
%rho = -0.7;
%V0 = 0.04;
%Current Loss: 5.8186e-05
%Calibration Results (Hessian-based):
%Estimated values: 3.7739    0.024099     0.76814    -0.57664    0.014904
%t-values: 0.013534    0.062767    0.020348   -0.073142    0.012422
%Standard errors: 278.8407     0.3839381      37.75027       7.88391      1.199775
%In-Sample MSE: 5.8186e-05

settings.parameters0 = [kappa; theta; eta; rho; V0];

%% MINUMUM AND MAXIMUM VALUES for the parameters of volatility model
settings.minKappa = 0.5; settings.maxKappa = 12;
settings.minTheta = 0.01^2; settings.maxTheta = 1;
settings.minEta = 0.05^2; settings.maxEta = 2;
settings.minRho = -1; settings.maxRho = 1;
settings.minV0 = 0.01^2; settings.maxV0 = 1;
settings.numberOfVariables = 5;

%% OPTIMIZATION SETTINGS
settings.calibrOptions.MaxFunEvals = 200*6;
settings.calibrOptions.MaxIter = 200*6;
settings.calibrOptions.TolFun = 1e-4;
settings.calibrOptions.TolX = 1e-4;
settings.calibrOptions.Display = 'iter';
settings.calibrOptions.FunValCheck = 'on';

%% DISPLAY SETTINGS, provisional result
settings.displayProvisionalResults = true;

if settings.displayProvisionalResults
    settings.indPlotSurface = 1;
    settings.displayParameters = 1;
    settings.calibrOptions.Display = 'iter';
    settings.calibrOptions.FunValCheck = 'on';
else
    settings.indPlotSurface = 0;
    settings.showParameters = 0;
    settings.calibrOptions.Display = 'off';
    settings.calibrOptions.FunValCheck = 'off';
end


end

