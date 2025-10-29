function phi = phi_uniform(p)
%   Uniform distribution
%
%   phi = phi_uniform(p)
%
%   Inputs:
%       p       : discretized preference values
%
%   Output:
%       phi : normalized uniform distribution

    phi = ones(size(p));

    % Normalize
    phi = phi / trapz(p, phi);
end