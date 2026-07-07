function phi = phi_Bernoulli(p_mean)
%   Bernoulli distribution
%
%   phi = phi_Bernoulli(p_mean)
%
%   Inputs:
%       p_mean       : mean
%
%   Output:
%       phi : Bernoulli distribution

phi = [1-p_mean;p_mean];

end