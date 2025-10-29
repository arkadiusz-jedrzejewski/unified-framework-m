function phi = phi_mixture_normals(p, mu1, mu2, sigma2)
%   Mixture of two normal distributions
%
%   phi = phi_mixture_normals(p, mu1, mu2, sigma)
%
%   Inputs:
%       p       : discretized preference values
%       mu1     : mean of the first dist
%       mu2     : mean of the second dist
%       sigma2  : variance
%
%   Output:
%       phi : normalized mixture of normal distributions

    phi = exp(-0.5 * (p - mu1).^2 ./ sigma2) + exp(-0.5 * (p - mu2).^2 ./ sigma2); 

    % Normalize
    phi = phi / trapz(p, phi);  
end