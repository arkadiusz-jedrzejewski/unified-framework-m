function phi = phi_logit_normal(p, mu, sigma2)
%   Logit-normal preference distribution
%
%   phi = phi_logit_normal(p, mu, sigma2)
%
%   Inputs:
%       p       : discretized preference values
%       mu      : location 
%       sigma2  : squared scale
%
%   Output:
%       phi : normalized logit-normal distribution
    logit = log(p./(1-p));
    phi = exp(-0.5 * (logit - mu).^2 ./ sigma2)./ (p.*(1-p));

    % Normalize
    phi = phi / trapz(p, phi); 
end