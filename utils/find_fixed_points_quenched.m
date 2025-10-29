function [a_vals, p_vals] = find_fixed_points_quenched(model, ps, phi, is_continuous)
%   Finds fixed-points of the quenched dynamics
%
%   Inputs:
%       model           : struct containing transition probabilites:
%                           model.X_BA(a), model.X_AB(a), model.Y_BA(a),
%                           model.Y_AB(a) 
%       ps              : discretized preference vaues
%       phi             : discretized preference distribution
%       is_continuous   : logical flag
%
%   Outputs:
%       a_vals  : fixed-points, stationary values of a
%       p_vals  : corresponding preference values
%
%   Notes:
%       - the interval [0, 1] is divided into subintervals to detect sign
%         changes on its borders in the equation for fixed-points
%       - roots in each subinterval are refined using fzero

    % Number of subintervals in [0, 1]
    num_a_sub = 107;    
    a_grid = linspace(0, 1, num_a_sub + 1);

    % Initialize outputs
    a_vals = [];
    p_vals = [];

    % Mean preference
    if is_continuous
        p_mean = trapz(ps,ps.*phi);
    else
        p_mean = ps'*phi;
    end
    
    % Loop over subintervals
    for i = 1:num_a_sub
        a_lower = a_grid(i);
        a_upper = a_grid(i + 1);
        
        % Check for sign change
        if sign(fixed_point_eq(a_lower, model, ps, phi, is_continuous)) ~= ...
           sign(fixed_point_eq(a_upper, model, ps, phi, is_continuous))
            
            % Refine the fixed point
            root = fzero(@(a) fixed_point_eq(a, model, ps, phi, is_continuous), [a_lower, a_upper]);
            
            % Store the fixed point and the corresponding p value
            a_vals = [a_vals; root];
            p_vals = [p_vals; p_mean];
        end
    end
end

function val = fixed_point_eq(a, model, ps, phi, is_continuous)
    % Compute a_p from the equation for fixed-points
    aps = get_a_p(model, ps, a);

    % Evaluate self-consistency equation
    if is_continuous
        val = a - trapz(ps,aps.*phi);
    else
        val = a - aps'*phi;
    end
end