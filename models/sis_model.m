function model = sis_model(beta, gamma)
%   Constructs the transition probabilites for the sis model.
%
%   model = sis_model(beta, gamma)
%   returns a struct containing function handles for:
%
%       X_BA(a) - transition probability from B to A under infection
%       X_AB(a) - transition probability from A to B under infection
%       Y_BA(a) - transition probability from B to A under recovery
%       Y_AB(a) - transition probability from A to B under recovery
%
%   and their derivatives:
%
%       dX_BA(a), dX_AB(a), dY_BA(a), dY_AB(a)
%
%   Inputs:
%       beta  : parameter controlling infection
%       gamma : parameter controlling recovery
%
%   Output: 
%       model : struct with 8 function handles

    % Transition probabilites
    model.X_BA = @(a) beta * a;
    model.X_AB = @(a) 0;
    model.Y_BA = @(a) 0;
    model.Y_AB = @(a) gamma;

    % Derivatives 
    model.dX_BA = @(a) beta;
    model.dX_AB = @(a) 0;
    model.dY_BA = @(a) 0;
    model.dY_AB = @(a) 0;
end