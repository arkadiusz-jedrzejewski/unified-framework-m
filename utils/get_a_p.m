function a_ps = get_a_p(model, ps, a)
%   Computes a_p for each preference value in ps and given a from the
%   equation for fixed-points under quenched dynamics
%
%   aps = get_a_p(model, ps, a)
%
%   Inputs:
%       model   : struct containing transition probabilites:
%                   model.X_BA(a), model.X_AB(a), model.Y_BA(a),
%                   model.Y_AB(a) 
%       ps      : verctor of preferences
%       a       : the fraction of agents in state A   
%
%   Outputs:
%       a_ps  : values of a_p for given a corresponding to each element in
%               ps

    % Equation for fixed-points
    a_ps = (model.Y_BA(a) - ps * (model.Y_BA(a) - model.X_BA(a))) ./ ...
            (model.Y_BA(a) + model.Y_AB(a) + ps * (model.X_BA(a) + model.X_AB(a) - model.Y_BA(a) - model.Y_AB(a)));
end