function [t_vals, a_vals] = solve_ode_annealed(T, p_mean, model, a_0)
%   Solves a differntial equation for a model under annealed dynamics with
%   mean preference p_mean
%
%   [t_vals, a_vals] = solve_ode_annealed(T, p_mean, model, a_0)
%   integrates the differential equation describing the time
%   evolution of the fraction of agents in state A. 
% 
%   Inputs:
%       T       : total simulation time
%       p_mean  : mean preference
%       model   : struct containing transition probabilites:
%                   model.X_BA(a), model.X_AB(a), model.Y_BA(a),
%                   model.Y_AB(a) 
%       a_0     : initial fractions of agents in state A
%  
%   Outputs:
%       t_vals  : time points returned by the ODE solver
%       a_vals  : total traction of agents in state A as a function of
%                   time, a(t) 

    % Transition probabilites
    P_BA = @(p, a) p .* model.X_BA(a) + (1 - p) .* model.Y_BA(a); 
    P_AB = @(p, a) p .* model.X_AB(a) + (1 - p) .* model.Y_AB(a);

    % Define ODE function
    odefun = @(t, a) da_dt(t, a, p_mean, P_BA, P_AB);
    
    % Solve for a(t)
    tspan = [0, T];
    [t_vals, a_vals] = ode89(odefun, tspan, a_0);
end

% Derivative da/dt
function da = da_dt(~, a, p, P_BA, P_AB)
    % Evaluate transition probabilites
    P_BA_values = P_BA(p, a);
    P_AB_values = P_AB(p, a);

    % Compute da/dt
    da = P_BA_values .* (1 - a) - P_AB_values .* a;
end
