% Fixed-point diagram for model with different distributions
% ------------------------------------------------------------------
% This script plots fixed-point diagrams for SIS models under quenched
% dynamics for three different distributions of preferences:
%   (1) one-point distribution,
%   (2) uniform distribution,
%   (3) normal distribution.
%
% Generated data was used for Fig. 7
% ------------------------------------------------------------------

clc; clear; close all;
addpath('models');
addpath('distributions')
addpath('utils')
figure

%% Model
% ------------------------------------------------------------------

% SIS model, which does not satisfy the balancing condition
model_builder = @sis_model;
betas = 0:0.01:1;
gamma = 0.2;

%% One-point distribution
% ------------------------------------------------------------------
ps = 0.5; 
phi = 1; 
is_continuous = false;

[a_sol, p_sol, beta_sol, stability] = ...
    scan_var(betas, model_builder, gamma, ps, phi, is_continuous);

% Plot diagram
subplot(1,3,1)
plot_fixed_points(beta_sol, a_sol, stability)
title('One-point distribution');

%% Uniform distribution
% ------------------------------------------------------------------
N=500;
mean_p = 0.5;
ps = linspace(0.0001,0.9999, N)'; 
phi = phi_uniform(ps);
is_continuous = true;

[a_sol, p_sol, beta_sol, stability] = ...
    scan_var(betas, model_builder, gamma, ps, phi, is_continuous);

% Plot diagram
subplot(1,3,2)
plot_fixed_points(beta_sol, a_sol, stability)
title('Unifrom distribution');

%% Normal distribution
% ------------------------------------------------------------------
ps=linspace(0.0001,0.9999,N)';
phi = phi_normal(ps, mean_p, 1/40);
is_continuous = true;

[a_sol, p_sol, beta_sol, stability] = ...
    scan_var(betas, model_builder, gamma, ps, phi, is_continuous);

% Plot diagram
subplot(1,3,3)
plot_fixed_points(beta_sol, a_sol, stability)
title('Normal distribution');

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
    xlabel('\beta'); 
    ylabel('a^*');
    ylim([0, 1]);
end

function [a_sol, p_sol, beta_sol, stability] = ...
    scan_var(betas, model_builder, gamma, ps, phi, is_continuous)

    % Preallocate as empty
    [a_sol, p_sol, stability, beta_sol] = deal([]);

    for beta = betas
        disp(beta)

        % Build model using the passed function handle
        model = model_builder(beta, gamma);

        % Find fixed points
        [a_vals, p_vals] = find_fixed_points_quenched(model, ps, phi, is_continuous);

        a_sol = [a_sol; a_vals];
        p_sol = [p_sol; p_vals];
        beta_sol = [beta_sol; beta * ones(size(p_vals))];

        % Check stability
        for a = a_vals'
            stability = [stability; ...
                determine_stability(a, model, ps, phi, is_continuous)];
        end
    end
end