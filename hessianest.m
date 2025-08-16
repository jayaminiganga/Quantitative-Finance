function H = hessianest(fun, x0)
% Estimate Hessian matrix of a scalar-valued function at x0
% Uses central difference approximation with adaptive step size
%
% Inputs:
%   fun - function handle, returns scalar value
%   x0  - point at which to estimate the Hessian (column vector)
%
% Output:
%   H   - Hessian matrix

    x0 = x0(:); % ensure column vector
    n = length(x0);
    
    % Adaptive step size based on parameter magnitude
    h = max(1e-8, 1e-5 * abs(x0));
    h(h == 0) = 1e-5; % handle zero parameters
    
    H = zeros(n);
    fx = fun(x0); % base evaluation
    
    for i = 1:n
        for j = i:n % exploit symmetry
            % Create step vectors
            hi = h(i);
            hj = h(j);
            
            % Four function evaluations for central difference
            x1 = x0; x1(i) = x1(i) + hi; x1(j) = x1(j) + hj;
            x2 = x0; x2(i) = x2(i) + hi; x2(j) = x2(j) - hj;
            x3 = x0; x3(i) = x3(i) - hi; x3(j) = x3(j) + hj;
            x4 = x0; x4(i) = x4(i) - hi; x4(j) = x4(j) - hj;
            
            f1 = fun(x1);
            f2 = fun(x2);
            f3 = fun(x3);
            f4 = fun(x4);
            
            % Second-order central difference approximation
            H(i,j) = (f1 - f2 - f3 + f4) / (4 * hi * hj);
            H(j,i) = H(i,j); % enforce symmetry
        end
    end
end