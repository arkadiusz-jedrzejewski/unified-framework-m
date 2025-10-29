function model = dynamical_system_model_decision_making(q, m)
%   Constructs the transition probabilites for the dynamical system model
%   of decision-making.
%
%   model = dynamical_system_model_decision_making(q, m)
%   returns a struct containing function handles for:
%
%       X_BA(a) - transition probability from B to A under individual
%           learning
%       X_AB(a) - transition probability from A to B under individual
%           learning
%       Y_BA(a) - transition probability from B to A under social learning
%       Y_AB(a) - transition probability from A to B under social learning
%
%   and their derivatives:
%
%       dX_BA(a), dX_AB(a), dY_BA(a), dY_AB(a)
%
%   Inputs:
%       q : parameter (exponent) controlling conformity (social learning)
%       m : relative merit of options
%
%   Output: 
%       model : struct with 8 function handles
    
    % Social learning function
    S = @(a) 0.5 * (2 * a) ^ q * (0 <= a  && a < 0.5) + (1 - 0.5 * (2 * (1 - a)) ^ q) * ( 0.5 <= a && a <= 1);
    
    % Derivative of the social learning function
    dS = @(a) q * (2 * a) ^ (q - 1) * (0 <= a  && a < 0.5) +  q * (2 * (1 - a)) ^ (q  - 1) * ( 0.5 <= a && a <= 1);
    
    % Transition probabilites
    model.X_BA = @(a) m;
    model.X_AB = @(a) 1 - m;
    model.Y_BA = @(a) S(a);
    model.Y_AB = @(a) S(1-a);
    
    % Derivatives 
    model.dX_BA = @(a) 0;
    model.dX_AB = @(a) 0;
    model.dY_BA = @(a) dS(a);
    model.dY_AB = @(a) -dS(1-a);
end