function phi = phi_normal(p, mu, sigma2)
%   Normal preference distribution
%
%   phi = phi_normal(p, mu, sigma2)
%
%   Inputs:
%       p       : discretized preference values
%       mu      : mean 
%       sigma2  : variance
%
%   Output:
%       phi : normalized normal distribution

    phi = exp(-0.5 * (p - mu).^2 ./ sigma2);

    % Normalize
    phi = phi / trapz(p, phi); 
end