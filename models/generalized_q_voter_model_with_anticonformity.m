function model = generalized_q_voter_model_with_anticonformity(q_a, q_c)
%   Constructs the transition probabilites for the generalized q-voter
%   model with anticonformity.
%
%   model = generalized_q_voter_model_with_anticonformity(q_a, q_c)
%   returns a struct containing function handles for:
%
%       X_BA(a) - transition probability from B to A under anticonformity
%       X_AB(a) - transition probability from A to B under anticonformity
%       Y_BA(a) - transition probability from B to A under conformity
%       Y_AB(a) - transition probability from A to B under conformity
%
%   and their derivatives:
%
%       dX_BA(a), dX_AB(a), dY_BA(a), dY_AB(a)
%
%   Inputs:
%       q_a : nonlinearity parameter (exponent) controlling anticonformity
%       q_c : nonlinearity parameter (exponent) controlling conformity
%
%   Output: 
%       model : struct with 8 function handles

    % Transition probabilites
    model.X_BA = @(a) (1 - a) ^ q_a;
    model.X_AB = @(a) a ^ q_a;
    model.Y_BA = @(a) a ^ q_c;
    model.Y_AB = @(a) (1 - a) ^ q_c;

    % Derivatives 
    model.dX_BA = @(a) - q_a * (1 - a) ^ (q_a - 1);
    model.dX_AB = @(a) q_a * a ^ (q_a - 1);
    model.dY_BA = @(a) q_c * a ^ (q_c - 1);
    model.dY_AB = @(a) q_c * (1 - a) ^ (q_c - 1);
end