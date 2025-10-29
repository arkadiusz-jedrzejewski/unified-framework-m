function model = noisy_threshold_q_voter_model(q, q_0)
%   Constructs the transition probabilites for the noisy thershold q-voter
%   model.
%
%   model = noisy_threshold_q_voter_model(q, q_0)
%   returns a struct containing function handles for:
%
%       X_BA(a) - transition probability from B to A under independence
%       X_AB(a) - transition probability from A to B under independence
%           learning
%       Y_BA(a) - transition probability from B to A under conformity
%       Y_AB(a) - transition probability from A to B under conformity
%
%   and their derivatives:
%
%       dX_BA(a), dX_AB(a), dY_BA(a), dY_AB(a)
%
%   Inputs:
%       q  : number of agnets in the infuance panel
%       q0 : threshold
%
%   Output: 
%       model : struct with 8 function handles
    
    % Transition probabilites for independence
    model.X_BA = @(a) 1/2;
    model.X_AB = @(a) 1/2;
    
    % Transition probabilites for conformity
    model.Y_BA = @(a) 0;
    model.Y_AB = @(a) 0;
    for i = q_0:q
        model.Y_BA = @(a) model.Y_BA(a) + nchoosek(q, i) * (a^i) * (1 - a)^(q - i);
    end
    for i = q_0:q
        model.Y_AB = @(a) model.Y_AB(a) + nchoosek(q, i) * (a^(q - i)) * (1 - a)^i;
    end
    
    % Derivatives
    model.dX_BA = @(a) 0;
    model.dX_AB = @(a) 0;
    model.dY_BA = @(a) 0;
    model.dY_AB = @(a) 0;
    for i = q_0:q
        model.dY_BA = @(a) model.dY_BA(a) + nchoosek(q, i) * ( i * a^(i - 1) * (1 - a)^(q - i) - (q - i) * a^i * (1 - a)^(q - i - 1));
    end
    for i = q_0:q
        model.dY_AB = @(a) model.dY_AB(a) + nchoosek(q, i) * ((q - i) * a^(q - i - 1) * (1 - a)^i - i * a^(q - i) * (1 - a)^(i - 1));
    end
end
