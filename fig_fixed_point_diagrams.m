% Fixed-point diagrams for models with diffrent distributions
% ------------------------------------------------------------------
% This script plots fixed-point diagrams for two models under quenched
% dynamics:
%   (1) the noisy threshold q-voter model (violating the balancing
%       condition) and
%   (2) the dynamical system model of decision-making (satisfying the
%       balancing condition).
% 
% For each model, the script plots fixed-point diagrams for three diffrent
% distributions of preferences:
%   (1) sliding one-point distribution,
%   (2) extending uniform distribution,
%   (3) Bernoulli distribution.
%
% Generated data was used for Figs. 4 and 5
% ------------------------------------------------------------------

clc; clear; close all;
addpath('models');
addpath('distributions')
addpath('utils')

%% Models
% ------------------------------------------------------------------

% Model that does not satisfy the balancing condition
%model = noisy_threshold_q_voter_model(10, 9);

% Model that satisfies the balancing condition
model = dynamical_system_model_decision_making(1.5, 0.6);

%% One-point distribution
% ------------------------------------------------------------------
mean_ps = 0:0.0009:0.4;
[a_sol, p_sol, stability] = deal([]);

for mean_p = mean_ps
    ps = mean_p; 
    phi = 1;  

    % Find fixed points
    [a_vals, p_vals] = find_fixed_points_quenched(model, ps, phi, false);
    a_sol = [a_sol; a_vals];
    p_sol = [p_sol; p_vals];

    % Check stability
    for a = a_vals'
        stability = [stability; determine_stability(a, model, ps, phi, false)];
    end
end

% Plot diagram
subplot(1,3,1)
plot_fixed_points(p_sol, a_sol, stability)
title('One-point distribution');

%% Uniform distribution
% ------------------------------------------------------------------
mean_ps = 0.001:0.0009:0.4;
[a_sol, p_sol, stability] = deal([]);
N = 500; % Number of discretized points for phi
% NOTE: accuracy of the solutions depends on N (especially stability 
% clasification). Increase N when higher precision is required.

for mean_p = mean_ps
    %display(mean_p)
    ps = linspace(0, mean_p*2, N)'; 
    phi = phi_uniform(ps);
    
    % Find fixed points
    [a_vals, p_vals] = find_fixed_points_quenched(model, ps, phi, true);
    a_sol = [a_sol; a_vals];
    p_sol = [p_sol; p_vals];

    % Check stability
    for a = a_vals'
        stability = [stability; determine_stability(a, model, ps, phi, true)];
    end
end

% Plot diagram
subplot(1,3,2)
plot_fixed_points(p_sol, a_sol, stability)
title('Uniform distribution');

%% Bernoulli distribution
% ------------------------------------------------------------------
mean_ps = 0:0.0009:0.4;
[a_sol, p_sol, stability] = deal([]);

for mean_p = mean_ps
    ps = [0, 1]'; 
    phi = [1-mean_p, mean_p]';  

    % Find fixed points
    [a_vals, p_vals] = find_fixed_points_quenched(model, ps, phi, false);
    a_sol = [a_sol; a_vals];
    p_sol = [p_sol; p_vals];

    % Check stability
    for a = a_vals'
        stability = [stability; determine_stability(a, model, ps, phi, false)];
    end
end

% Plot diagram
subplot(1,3,3)
plot_fixed_points(p_sol, a_sol, stability)
title('Bernoulli distribution');

%% Helper functions
% ------------------------------------------------------------------
function plot_fixed_points(p_sol, a_sol, stability)
    % Sort points
    [a_sol,idx]=sort(a_sol);
    p_sol = p_sol(idx);
    stability = stability(idx);

    % Plot stable (black) and unstable (red) points
    mask_stable = stability < 0;
    mask_unstable = stability >= 0;
    
    hold on;
    plot(p_sol(mask_stable), a_sol(mask_stable), 'k.');
    plot(p_sol(mask_unstable), a_sol(mask_unstable), 'r.');
    xlabel('p'); 
    ylabel('a^*');
    ylim([0, 1]);
end