function phi = phi_beta(p, alpha, beta)
%   Beta distribution
%
%   phi = phi_beta(p, alpha, beta)
%
%   Inputs:
%       p       : discretized preference values
%       alpha   : shape parameter
%       beta    : shape parameter
%
%   Output:
%       phi : normalized mixture of normal distributions

    phi = gamma(alpha+beta)/(gamma(alpha)*gamma(beta))*p.^(alpha-1).*(1-p).^(beta-1);
    
    % Normalize
    phi = phi / trapz(p, phi);  
end