function stability_flag = determine_stability(a, model, ps, phi, is_continuous)
%   Checks stability of flixed points
%
%   Inputs:
%       a               : fraction of agents in state A
%       model           : struct containing transition probabilites:
%                           model.X_BA(a), model.X_AB(a), model.Y_BA(a),
%                           model.Y_AB(a) 
%       ps              : discretized preference vaues
%       phi             : discretized preference distribution  
%       is_continuous   : logical flag
%
%   Outputs:
%       stability_flag  : integer indicating:
%                           -1 : stable
%                           1  : unstable
%                           0  : undefined

    % Compute statinary values of a_p
    a_ps = get_a_p(model, ps, a);

    % Precompute transition probabilites
    X_BA_a = model.X_BA(a);
    Y_BA_a = model.Y_BA(a);
    X_AB_a = model.X_AB(a);
    Y_AB_a = model.Y_AB(a);
    dX_BA_a = model.dX_BA(a);
    dY_BA_a = model.dY_BA(a);
    dX_AB_a = model.dX_AB(a);
    dY_AB_a = model.dY_AB(a);
    
    % Compute Jacobian matrix 
    off_diag = ((1 - a_ps) .* (ps .* dX_BA_a + (1 - ps) .* dY_BA_a)) .* phi' ...
                - (a_ps .* (ps .* dX_AB_a + (1 - ps) .* dY_AB_a)) .* phi';
    if is_continuous
        off_diag = off_diag*(ps(2)-ps(1));
    end

    diag_terms = ps .* X_BA_a + (1 - ps) .* Y_BA_a + ps .* X_AB_a + (1 - ps) .* Y_AB_a;

    jacobian_matrix = off_diag;
    jacobian_matrix(1:length(ps) + 1:end) = jacobian_matrix(1:length(ps) + 1:end) - diag_terms';

    % Compute eigenvalues
    eigenvalues = eig(jacobian_matrix);
    real_parts = real(eigenvalues);

    % Determine stability
    if all(real_parts < 0)
        stability_flag = -1;
    elseif any(real_parts > 0)
        stability_flag = 1;
    else
        stability_flag = 0;
    end
end