function [loss, modelVolSurf] = pricingError(data, settings, parameters)
    % PRICINGERROR Computes model IV surface and error against market IV surface.
    %
    % Inputs:
    %   - data: Struct with market data (T, K, S0, r, IVolSurf)
    %   - settings: Struct with calibration and model settings
    %   - parameters: Model parameters [kappa, theta, eta, rho, V0]
    %
    % Outputs:
    %   - loss: Mean squared error between market and model implied volatilities
    %   - modelVolSurf: Model implied volatility surface (matrix)
    
    %     % Enforce parameter bounds (only during calibration, not when computing errors)

    if settings.standardErrors == false
        % Apply parameter bounds from settings
        if any([...
            parameters(1) < settings.minKappa, parameters(1) > settings.maxKappa, ...
            parameters(2) < settings.minTheta, parameters(2) > settings.maxTheta, ...
            parameters(3) < settings.minEta,   parameters(3) > settings.maxEta, ...
            parameters(4) < settings.minRho,   parameters(4) > settings.maxRho, ...
            parameters(5) < settings.minV0,    parameters(5) > settings.maxV0 ...
        ])
            loss = 1e6;  % Penalty for constraint violation
            modelVolSurf = NaN;
            return;
        end
    end

    % Reorder parameters to match FFT CallPricing function:
    % {V0, theta, kappa, eta, rho}
    V0    = parameters(5);
    theta = parameters(2);
    kappa = parameters(1);
    eta   = parameters(3);
    rho   = parameters(4);

    % Initialize model volatility surface
    numT = length(data.T);
    numK = length(data.K);
    modelVolSurf = NaN(numT, numK);

    % Compute model option prices and implied vol
    for i = 1:numT
        for j = 1:numK
            T = data.T(i);
            K = data.K(j);
            S0 = data.S0;
            r = data.r;

            try
                % Compute option price using FFT
                price = CallPricingFFT(settings.model, settings.n, S0, K, T, r, 0, ...
                    V0, theta, kappa, eta, rho);

                % Discard non-positive prices (replace with NaN)
                if price <= 0
                    modelVolSurf(i, j) = NaN;
                else
                    % Convert to implied vol using Black-Scholes
                    iv = blsimpv(S0, K, r, T, price, 'Yield', 0, ...
                        'Limit', [1e-4 5], 'Class', {'call'});

                    % Handle conversion failures
                    %if isempty(iv) || isnan(iv)
                        %modelVolSurf(i, j) = NaN;
                    %else
                        %modelVolSurf(i, j) = iv;

                    if isempty(iv) || isnan(iv) || iv == 0
                        modelVolSurf(i, j) = NaN;
                    else
                        modelVolSurf(i, j) = iv;
                    

                    end
                end
            catch
                modelVolSurf(i, j) = NaN;
            end
        end
    end

    % **Only Interpolate Model Volatility Surface (not Market Volatility)**
    % Interpolate missing values in the model surface (linear, fallback = 1e-6)

    nanMask = isnan(modelVolSurf);
    if any(nanMask, 'all')
        [Kmesh, Tmesh] = meshgrid(data.K, data.T);
        modelVolSurf = interp2(Kmesh, Tmesh, modelVolSurf, Kmesh, Tmesh, 'linear', 1e-6);
    end

    % Reference market implied volatility surface
    fixedIVSurf = data.IVolSurf;

    % Compute mean squared error between model and market volatilities
    validMask = ~isnan(fixedIVSurf) & ~isnan(modelVolSurf);
    squaredDiffs = (fixedIVSurf(validMask) - modelVolSurf(validMask)).^2;
    loss = mean(squaredDiffs); % Final calibration objective

    % Visualization during optimization
    if settings.displayProvisionalResults && isfield(settings, 'indPlotSurface') && settings.indPlotSurface
        logK = log(data.K / data.S0); % Log-strike axis
        [logKmesh, Tmesh] = meshgrid(logK, data.T);

        figure(101); clf;

        % Transparent model surface, opaque market surface
        surf(Tmesh, logKmesh, modelVolSurf, 'FaceAlpha', 0.2);
        hold on;
        surf(Tmesh, logKmesh,  fixedIVSurf, 'FaceAlpha', 1);

        % Update the axis labels 
        ylabel('Log-Strike ($\log(K/S_0)$)', 'Interpreter', 'LaTex');
        xlabel ('Maturity time ($T$, in years)', 'Interpreter', 'LaTex');
        zlabel('Implied volatility', 'Interpreter', 'LaTex');
        title('Data (non-transparent), model (transparent) ', 'Interpreter', 'LaTex');


        set(gca, 'XDir', 'reverse');

        view(45, 30);
        grid on;
        colorbar;

        hold off;
        drawnow;
    end

    % Display current parameters and loss during calibration
    if settings.displayProvisionalResults && isfield(settings, 'displayParameters') && settings.displayParameters
        disp('Current Parameters (kappa, theta, eta, rho, V0):');
        disp(parameters');
        disp(['Current Loss: ', num2str(loss)]);
    end
end
