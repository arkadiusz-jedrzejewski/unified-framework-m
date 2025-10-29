% Time evolution for models under annealed and quenched dynamics
% ------------------------------------------------------------------
% This script compares time evolution of two generalized q-voter models
% with anticonformity:
%   (1) one satisfying the balancing condition and
%   (2) one violating it.
% 
% For each model, the script computes both the anneald and quenched 
% dynamics using diffrent preference distributions.
% 
% Generated data was used for Fig. 3.
% ==================================================================

clc; clear; close all;
addpath('models');
addpath('distributions');
addpath('utils')

%% Parameters
% ------------------------------------------------------------------
T = 400;                % Time horizon
                        % Parameters for quenched dynamics
N = 500;                % Discretization of p, p in [0, 1]
p = linspace(0, 1, N)'; 
a_p0 = 1 * ones(N, 1);  % Initial a_p(t=0)
                        % Parameters for annealed dynamics
p_mean = 0.5;           % Mean preference
a_0 = 1;                % Initial value of a(t=0)

%% Distributions
% ------------------------------------------------------------------
phi_uniform = phi_uniform(p);
phi_normal = phi_normal(p, 0.5, 1/800);  
phi_mixture_normals = phi_mixture_normals(p, 0, 1, 1/800);
phi_beta = phi_beta(p, 2, 2);

% Dispay means for verification
fprintf('Mean uniform: %.3f\n', trapz(p,p.*phi_uniform))
fprintf('Mean normal: %.3f\n', trapz(p,p.*phi_normal))
fprintf('Mean mixture: %.3f\n',trapz(p,p.*phi_mixture_normals))
fprintf('Mean beta: %.3f\n', trapz(p,p.*phi_beta))

%% Models
% ------------------------------------------------------------------
model_balanced = generalized_q_voter_model_with_anticonformity(6,6);
model_unbalanced = generalized_q_voter_model_with_anticonformity(6,2);

%% Solve ODEs for Quenched dynamics 
% ------------------------------------------------------------------
% Uniform distribution
[t_vals, a_vals] = solve_ode_quenched_continuous_phi(T, p, phi_uniform, model_balanced, a_p0);

subplot(4,3,[3,6])
plot(t_vals, a_vals, 'b', 'LineWidth', 2);
ylabel('a');
title('Quenched dynamics');
ylim([0 1])
hold on;

[t_vals, a_vals] = solve_ode_quenched_continuous_phi(T, p, phi_uniform, model_unbalanced, a_p0);

subplot(4,3,[9,12])
plot(t_vals, a_vals, 'b', 'LineWidth', 2);
xlabel('t'); ylabel('a');
ylim([0 1])
hold on;

% Mixture of two normal distributions
[t_vals, a_vals] = solve_ode_quenched_continuous_phi(T, p, phi_mixture_normals, model_balanced, a_p0);

subplot(4,3,[3,6])
plot(t_vals, a_vals, 'r', 'LineWidth', 2);

[t_vals, a_vals] = solve_ode_quenched_continuous_phi(T, p, phi_mixture_normals, model_unbalanced, a_p0);

subplot(4,3,[9,12])
plot(t_vals, a_vals, 'r', 'LineWidth', 2);

% Normal distibution
[t_vals, a_vals] = solve_ode_quenched_continuous_phi(T, p, phi_normal, model_balanced, a_p0);

subplot(4,3,[3,6])
plot(t_vals, a_vals, 'g', 'LineWidth', 2);

[t_vals, a_vals] = solve_ode_quenched_continuous_phi(T, p, phi_normal, model_unbalanced, a_p0);

subplot(4,3,[9,12])
plot(t_vals, a_vals, 'g', 'LineWidth', 2);

% Beta distibution
[t_vals, a_vals] = solve_ode_quenched_continuous_phi(T, p, phi_beta, model_balanced, a_p0);

subplot(4,3,[3,6])
plot(t_vals, a_vals, 'm', 'LineWidth', 2);
hold off

[t_vals, a_vals, a_p_vals] = solve_ode_quenched_continuous_phi(T, p, phi_beta, model_unbalanced, a_p0);

subplot(4,3,[9,12])
plot(t_vals, a_vals, 'm', 'LineWidth', 2);
hold off

%% Solve ODEs for annealed dynamics 
% ------------------------------------------------------------------
[t_vals, a_vals] = solve_ode_annealed(T, p_mean, model_balanced, a_0);

subplot(4,3,[2,5])
plot(t_vals, a_vals, 'k', 'LineWidth', 2);
ylabel('a');
title('Annealed dynamics');
ylim([0 1])

[t_vals, a_vals] = solve_ode_annealed(T, p_mean, model_unbalanced, a_0);

subplot(4,3,[8,11])
plot(t_vals, a_vals, 'k', 'LineWidth', 2);
xlabel('t'); 
ylabel('a');
ylim([0 1])

%% Plot distributions
% ------------------------------------------------------------------
subplot(4,3,1)
plot(p, phi_uniform, 'b', 'LineWidth', 2);
xlabel('p'); 
ylabel('\phi(p)');
title('Preference distribution');

subplot(4,3,4)
plot(p, phi_normal, 'g', 'LineWidth', 2);
xlabel('p'); 
ylabel('\phi(p)');
title('Preference distribution');

subplot(4,3,7)
plot(p, phi_mixture_normals, 'r', 'LineWidth', 2);
xlabel('p'); 
ylabel('\phi(p)');
title('Preference distribution');

subplot(4,3,10)
plot(p, phi_beta, 'm', 'LineWidth', 2);
xlabel('p'); 
ylabel('\phi(p)');
title('Preference distribution');