function price = priceDownInAsianAvgStrikeCall(params, H, M, T, dt, S0, r)
    % Extract Heston parameters
    kappa = params(1);
    theta = params(2);
    eta = params(3);
    rho = params(4);
    V0 = params(5);
    
    % Number of time steps
    N = round(T/dt);
    
    % Preallocate arrays
    S = zeros(M, N+1);
    V = zeros(M, N+1);
    S_anti = zeros(M, N+1);
    V_anti = zeros(M, N+1);
    
    % Initialize
    S(:,1) = S0;
    V(:,1) = V0;
    S_anti(:,1) = S0;
    V_anti(:,1) = V0;
    
    % Track barrier hits
    hit_barrier = false(M,1);
    hit_barrier_anti = false(M,1);
    
    % Generate all random numbers upfront for efficiency
    Z1 = randn(M, N);
    e12 = randn(M, N);
    Z2 = rho*Z1 + sqrt(1-rho^2)*e12;
    
    % Antithetic variates
    Z1_anti = -Z1;
    e12_anti = -e12;
    Z2_anti = rho*Z1_anti + sqrt(1-rho^2)*e12_anti;
    
    % Milstein scheme with truncation
    for i = 1:N
        sqrt_dt = sqrt(dt);
        sqrt_V = sqrt(V(:,i));
        sqrt_V_anti = sqrt(V_anti(:,i));
        
        % Main paths
        V(:,i+1) = max(V(:,i) + kappa*(theta - V(:,i))*dt + ...
                   eta*sqrt_V.*Z2(:,i)*sqrt_dt + ...
                   0.25*eta^2*dt*(Z2(:,i).^2 - 1), 0);
        
        S(:,i+1) = S(:,i).*exp((r - 0.5*V(:,i))*dt + sqrt_V.*Z1(:,i)*sqrt_dt);
        
        % Antithetic paths
        V_anti(:,i+1) = max(V_anti(:,i) + kappa*(theta - V_anti(:,i))*dt + ...
                        eta*sqrt_V_anti.*Z2_anti(:,i)*sqrt_dt + ...
                        0.25*eta^2*dt*(Z2_anti(:,i).^2 - 1), 0);
        
        S_anti(:,i+1) = S_anti(:,i).*exp((r - 0.5*V_anti(:,i))*dt + ...
                        sqrt_V_anti.*Z1_anti(:,i)*sqrt_dt);
        
        % Check barrier hits
        hit_barrier = hit_barrier | (S(:,i+1) < H);
        hit_barrier_anti = hit_barrier_anti | (S_anti(:,i+1) < H);
    end
    
    % Calculate payoffs
    %avg_S = mean(S(:,2:end), 2); % Exclude initial price
    avg_S = mean(S, 2);
    payoff = max(S(:,end) - avg_S, 0) .* hit_barrier;
    
    %avg_S_anti = mean(S_anti(:,2:end), 2);
    avg_S_anti = mean(S_anti, 2);
    payoff_anti = max(S_anti(:,end) - avg_S_anti, 0) .* hit_barrier_anti;
    
    % Discount and average with antithetic
    price = mean(0.5 * (exp(-r*T)*payoff + exp(-r*T)*payoff_anti));
end