function [t_vals, a_vals, a_p_vals] = solve_ode_quenched_continuous_phi(T, p, phi, model, a_p0)
%   Solves differntial equations for a model under quenched dynamics with a
%   continuous preference distribution phi(p), discretized for numerical
%   integration
%
%   [t_vals, a_vals, a_p_vals] = solve_ode_quenched_continuous_phi(T, p, phi, model, a_p0) 
%   integrates the system of differential equations describing the time
%   evolution of the fraction of agents in state A with preference p, a_p,
%   for each preference value p. The total fraction of agents in state A,
%   a, is then computed as a weighted average over the preference
%   distribution phi(p).
% 
%   Inputs:
%       T       : total simulation time
%       p       : discretized preference values
%       phi     : preference distribution
%       model   : struct containing transition probabilites:
%                   model.X_BA(a), model.X_AB(a), model.Y_BA(a),
%                   model.Y_AB(a) 
%       a_p0    : initial fractions of agents in state A for each p
%  
%   Outputs:
%       t_vals  : time points returned by the ODE solver
%       a_vals  : total traction of agents in state A as a function of
%                   time, a(t) 
%       a_p_vals: fractions of agents in state A for each p as a function
%                   of time, a_p(t) 

    % Transition probabilites
    P_BA = @(p, a) p .* model.X_BA(a) + (1 - p) .* model.Y_BA(a); 
    P_AB = @(p, a) p .* model.X_AB(a) + (1 - p) .* model.Y_AB(a);

    % Define ODE function
    odefun = @(t, a_p) da_dt(t, a_p, p, phi, P_BA, P_AB);
    
    % Solve the system
    tspan = [0, T];
    [t_vals, a_p_vals] = ode89(odefun, tspan, a_p0);
    
    % Compute a(t) from a_p(t)
    a_vals = trapz(p, a_p_vals .* phi', 2);
end

% Derivative da_p/dt
function da_p = da_dt(~, a_p, p, phi, P_BA, P_AB)
    % Compute the total fraction of agents in state A
    a = trapz(p, a_p .*phi); 

    % Evaluate transition probabilites for all values of preferences
    P_BA_values = P_BA(p, a); 
    P_AB_values = P_AB(p, a);  

    % Compute da_p/dt
    da_p = P_BA_values .* (1 - a_p) - P_AB_values .* a_p; 
end
