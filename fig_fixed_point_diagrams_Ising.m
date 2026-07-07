% Fixed-point diagram for spin models with different distributions
% ------------------------------------------------------------------
% This script plots fixed-point diagrams for two spin models with 
% competing temperatures under quenched dynamics and:
%   (1) Glauber dynamics and
%   (2) Metropilis dynamics.
% 
% For each model, the script plots fixed-point diagrams for three different
% distributions of preferences:
%   (1) one-point distribution,
%   (2) logit-normal distribution,
%   (3) Bernoulli distribution.
%
% Generated data was used for Fig. 6.
% ------------------------------------------------------------------

clc; clear; close all;
addpath('models');
addpath('distributions')
addpath('utils')

%% Models
% ------------------------------------------------------------------

% Parameters
T2 = 1;
Q2 = 3;
Q1 = 4;

%% One-point distribution
% ------------------------------------------------------------------
ps = 0.5; 
phi = 1; 
is_continuous = false;

model_builder = @ising_model_with_competing_temperatures_Metropolis;
T1s = 2.0:0.001:4;
[a_sol, p_sol, t_sol, stability] = ...
    scan_T1(T1s, model_builder, T2, Q1, Q2, ps, phi, is_continuous);

% Plot diagram
subplot(2,3,1)
plot_fixed_points(t_sol, a_sol, stability)
title('One-point distribution');

model_builder = @ising_model_with_competing_temperatures_Glauber;
T1s = 2.0:0.005:6;
[a_sol, p_sol, t_sol, stability] = ...
    scan_T1(T1s, model_builder, T2, Q1, Q2, ps, phi, is_continuous);

% Plot diagram
subplot(2,3,4)
plot_fixed_points(t_sol, a_sol, stability)

%% Logit-normal distribution
% ------------------------------------------------------------------
N=500;
mean_p = 0.5;
ps=linspace(1e-4,1-1e-4,N)';
phi=phi_logit_normal(ps,0,2^2);
is_continuous = true;

model_builder = @ising_model_with_competing_temperatures_Metropolis;
T1s = 2.0:0.001:4;
[a_sol, p_sol, t_sol, stability] = ...
    scan_T1(T1s, model_builder, T2, Q1, Q2, ps, phi, is_continuous);

% Plot diagram
subplot(2,3,2)
plot_fixed_points(t_sol, a_sol, stability)
title('Logit-normal distribution');

model_builder = @ising_model_with_competing_temperatures_Glauber;
T1s = 2.0:0.005:6;
[a_sol, p_sol, t_sol, stability] = ...
    scan_T1(T1s, model_builder, T2, Q1, Q2, ps, phi, is_continuous);

% Plot diagram
subplot(2,3,5)
plot_fixed_points(t_sol, a_sol, stability)

%% Bernoulli distribution
% ------------------------------------------------------------------
mean_p = 0.5;
ps = [0, 1]'; 
phi = [1-mean_p, mean_p]';
is_continuous = false;

model_builder = @ising_model_with_competing_temperatures_Metropolis;
T1s = 2.0:0.005:4;
[a_sol, p_sol, t_sol, stability] = ...
    scan_T1(T1s, model_builder, T2, Q1, Q2, ps, phi, is_continuous);

% Plot diagram
subplot(2,3,3)
plot_fixed_points(t_sol, a_sol, stability)
title('Bernoulli distribution');

model_builder = @ising_model_with_competing_temperatures_Glauber;
T1s = 2.0:0.005:6;
[a_sol, p_sol, t_sol, stability] = ...
    scan_T1(T1s, model_builder, T2, Q1, Q2, ps, phi, is_continuous);

% Plot diagram
subplot(2,3,6)
plot_fixed_points(t_sol, a_sol, stability)

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
    xlabel('T_1'); 
    ylabel('a^*');
    ylim([0, 1]);
end

function [a_sol, p_sol, t_sol, stability] = ...
    scan_T1(T1s, model_builder, T2, Q1, Q2, ps, phi, is_continuous)

    % Preallocate as empty
    [a_sol, p_sol, stability, t_sol] = deal([]);

    for T1 = T1s
        disp(T1)

        % Build model using the passed function handle
        model = model_builder(T1, T2, Q1, Q2);

        % Find fixed points
        [a_vals, p_vals] = find_fixed_points_quenched(model, ps, phi, is_continuous);

        a_sol = [a_sol; a_vals];
        p_sol = [p_sol; p_vals];
        t_sol = [t_sol; T1 * ones(size(p_vals))];

        % Check stability
        for a = a_vals'
            stability = [stability; ...
                determine_stability(a, model, ps, phi, is_continuous)];
        end
    end
end
