function model = ising_model_with_competing_temperatures_Metropolis(T1,T2,Q1,Q2)
%   Constructs the transition probabilites for the Ising model with
%   competing temperatures and Metropolis dynamics
%
%   model = ising_model_with_competing_temperatures_Metropolis(T1,T2,Q1,Q2)
%   returns a struct containing function handles for:
%
%       X_BA(a) - transition probability from B to A under T1
%       X_AB(a) - transition probability from A to B under T1
%           learning
%       Y_BA(a) - transition probability from B to A under T2
%       Y_AB(a) - transition probability from A to B under T2
%
%   and their derivatives:
%
%       dX_BA(a), dX_AB(a), dY_BA(a), dY_AB(a)
%
%   Inputs:
%       T1 : temperature of the first heat reservoir
%       T2 : temperature of the second heat reservoir
%       Q1 : number of interacting spins in the first heat reservoir
%       Q2 : number of interacting spins in the second heat reservoir
%   Output: 
%       model : struct with 8 function handles
   
    % revers temperatures
    beta1 = 1/T1;
    beta2 = 1/T2;

    % Initialization probabilites
    model.X_BA = @(a) 0;
    model.X_AB = @(a) 0;
    model.Y_BA = @(a) 0;
    model.Y_AB = @(a) 0;

    % Initialization derivatives
    model.dX_BA = @(a) 0;
    model.dX_AB = @(a) 0;
    model.dY_BA = @(a) 0;
    model.dY_AB = @(a) 0;

    % T1
    for i = 0:Q1
        model.X_BA = @(a) model.X_BA(a) + nchoosek(Q1, i) * (a^i) * (1 - a)^(Q1 - i) * min(1,exp(2*beta1*(2*i-Q1)));
        model.X_AB = @(a) model.X_AB(a) + nchoosek(Q1, i) * (a^i) * (1 - a)^(Q1 - i) * min(1,exp(-2*beta1*(2*i-Q1)));
        model.dX_BA = @(a) model.dX_BA(a) + nchoosek(Q1, i) * (a^(i - 1)*i*(1 - a)^(Q1 - i) - a^i*(Q1 - i)*(1 - a)^(Q1 - i - 1)) * min(1,exp(2*beta1*(2*i-Q1)));
        model.dX_AB = @(a) model.dX_AB(a) + nchoosek(Q1, i) * (a^(i - 1)*i*(1 - a)^(Q1 - i) - a^i*(Q1 - i)*(1 - a)^(Q1 - i - 1)) * min(1,exp(-2*beta1*(2*i-Q1)));
    end

    % T2
    for i = 0:Q2
        model.Y_BA = @(a) model.Y_BA(a) + nchoosek(Q2, i) * (a^i) * (1 - a)^(Q2 - i) * min(1,exp(2*beta2*(2*i-Q2)));
        model.Y_AB = @(a) model.Y_AB(a) + nchoosek(Q2, i) * (a^i) * (1 - a)^(Q2 - i) * min(1,exp(-2*beta2*(2*i-Q2)));
        model.dY_BA = @(a) model.dY_BA(a) + nchoosek(Q2, i) * (a^(i - 1)*i*(1 - a)^(Q2 - i) - a^i*(Q2 - i)*(1 - a)^(Q2 - i - 1)) * min(1,exp(2*beta2*(2*i-Q2)));
        model.dY_AB = @(a) model.dY_AB(a) + nchoosek(Q2, i) * (a^(i - 1)*i*(1 - a)^(Q2 - i) - a^i*(Q2 - i)*(1 - a)^(Q2 - i - 1)) * min(1,exp(-2*beta2*(2*i-Q2)));
    end
end
