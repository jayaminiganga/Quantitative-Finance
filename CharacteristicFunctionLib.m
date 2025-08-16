% This is material illustrating the methods from the book
% Financial Modelling  - Theory, Implementation and Practice with Matlab
% source
% Wiley Finance Series
% ISBN 978-0-470-74489-5
%
% Date: 02.05.2012
%
% Authors:  Joerg Kienitz
%           Daniel Wetterau
%
% Please send comments, suggestions, bugs etc. to neueemail
%
% The authors are not responsible for any loss or damage from applying
% this code
%
% (C) Kienitz, Wetterau


function y = CharacteristicFunctionLib(model,u,lnS,T,r,d,varargin)
%---------------------------------------------------------
% Characteristic Function Library of the following models:
%---------------------------------------------------------
% Black Scholes
% Merton Jump Diffusion
% Heston Stochastic Volatility Model
% Bates Stochastic Volatility / Jump Diffusion Model
% Displaced Heston
% Heston-Hull-White
% Variance Gamma
% VGSIG
% for further details about time-changed model, please refer to 'Stochastic
% volatility for Lévy processes' written by P.Carr,
% et al.
%---------------------------------------------------------
%
optAlfaCalculation = true;

ME1 = MException('VerifyInput:InvalidNrOfArguments',...
    'Invalid number of Input arguments');
ME2 = MException('VerifyInput:InvalidModel',...
    'Undefined Model');

if strcmp(model,'BlackScholes')
    if nargin == 7
        funobj = @cf_bs;
    else 
        throw(ME1)
    end
elseif strcmp(model,'DDHeston')
    if nargin == 13
        funobj = @cf_ddheston;
    else 
        throw(ME1)
    end
elseif strcmp(model,'Heston')
    if nargin == 11
        funobj = @cf_heston;
    else 
        throw(ME1)
    end
elseif strcmp(model,'HestonHullWhite')
   if nargin == 14
       funobj = @cf_hestonhullwhite;
   else
       throw(ME1)
   end

elseif strcmp(model,'NIGSIG')
    if nargin == 13
        funobj = @cf_nigsig;
    else 
        throw(ME1)
    end
elseif strcmp(model,'VGSIG')
    if nargin == 13
        funobj = @cf_vgsig;
    else 
        throw(ME1)
    end
else
    throw(ME2)
end

fval = feval(funobj,u,lnS,T,r,d,varargin{:});

if optAlfaCalculation == true
    y = fval;
else
    y = exp(fval);
end

end


%% Explicit Implementation of the characteristic Functions E[exp(iu*lnS_T)]
%-----------------------------------------------------------------------
   
function y = cf_bs(u,lnS,T,r,d,sigma)
% Black Scholes
    y = 1i*u*(lnS+(r-d-0.5*sigma*sigma)*T) - 0.5*sigma*sigma*u.*u*T;
end


function y = cf_merton(u,lnS,T,r,d,sigma,a,b,lambda)
% Merton Jump Diffusion
    y = cf_bs(u,lnS,T,r,d,sigma) ...
        + cf_jumplognormal(u,a,b,lambda,T);
end

function y = cf_ddheston(u,lnS,T,r,d,V0,theta,kappa,omega,rho,lambda,b)
% Displaced Diffusion Heston
U = 1i*u;
v = 0.5*(lambda*b)^2*U.*(U-1);
theta_star = kappa -rho*omega*lambda*b*U;
gamma = sqrt(theta_star.^2-2*omega^2*v);
Avu = kappa * theta / omega^2 *(2*log(2*gamma./(theta_star ...
    + gamma-exp(-gamma*T).*(theta_star-gamma)))+(theta_star-gamma)*T);
Bvu = 2*v.*(1-exp(-gamma*T))./((theta_star+gamma)...
    .*(1-exp(-gamma*T))+2*gamma.*exp(-gamma*T));

y = Avu + Bvu*V0 + U *(lnS+(r-d)*T);
end


function y = cf_heston(u,lnS,T,r,d,V0,theta,kappa,omega,rho)
% Heston  
alfa = -.5*(u.*u + u*1i);
beta = kappa - rho*omega*u*1i;
omega2 = omega * omega;
gamma = .5 * omega2;

D = sqrt(beta .* beta - 4.0 * alfa .* gamma);

bD = beta - D;
eDt = exp(- D * T);

G = bD ./ (beta + D);
B = (bD ./ omega2) .* ((1.0 - eDt) ./ (1.0 - G .* eDt));
psi = (G .* eDt - 1.0) ./(G - 1.0);
A = ((kappa * theta) / (omega2)) * (bD * T - 2.0 * log(psi));

y = A + B*V0 + 1i*u*(lnS+(r-d)*T);

end

function y = cf_hhw(u,lnS,T,r0,d,V0,theta,kappa,omega,lambda,eta,rho12,rho13,ircurve)
% Heston Hull White with % correlation(variance,rate) = 0
% dr(t) = lambda(r-r(t))dt + eta dW(t); r constant
    
    D1 = sqrt((omega*rho12*1i*u-kappa).^2-omega^2*1i*u.*(1i*u-1));
    g = (kappa-omega*rho12*1i*u-D1)./(kappa-omega*rho12*1i*u+D1);
    
    a = sqrt(theta - .125 * omega^2/kappa);
    b = sqrt(V0) - a;
    ct=.25*omega^2*(1-exp(-kappa))/kappa;
    lambdat=4*kappa*V0*exp(-kappa)/(omega^2*(1-exp(-kappa)));
    d2=4*kappa*theta/omega^2;
    F1 = sqrt(ct*(lambdat-1)+ct*d2+ct*d2/(2*(d2+lambdat)));
    c = -log((F1-a)/b);
    
    I2 = kappa*theta/omega^2*(T*(kappa-omega*rho12*1i*u-D1)-2*log((1-g.*exp(-D1*T))./(1-g)));
    I3 = eta^2*(1i+u).^2/(4*lambda^3)*(3+exp(-2*lambda*T)-4*exp(-lambda*T)-2*lambda*T);
    I4 = -eta*rho13/lambda *(1i*u+u.^2)*(b/c*(1-exp(-c*T))+a*T+a/lambda*(exp(-lambda*T)-1)+b/(c-lambda)*exp(-c*T)*(1-exp(-T*(lambda-c))));
    
    % curve stuff
    date_T = add2date(ircurve.Settle,T);
    Theta = (1-1i*u) * (log(ircurve.getDiscountFactors(date_T))+eta^2/(2*lambda^3)*(T/lambda+2*(exp(-lambda*T)-1)-0.5*(exp(-2*lambda)-1)));
    
    A = I2+I3+I4+Theta;
    BV = (1-exp(-D1*T))./(omega^2.*(1-g.*exp(-D1*T))).*(kappa-omega*rho12*1i*u-D1);
    Br = (1i*u-1)/lambda*(1-exp(-lambda*T));
    
    y = A + 1i*u * (lnS + (r0-d)*T) + BV * V0  + Br * r0;
end

























 
