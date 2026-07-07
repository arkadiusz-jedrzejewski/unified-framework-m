% Time evolution for models under annealed and quenched dynamics
% ------------------------------------------------------------------
% This script compares time evolution of 3 generalized q-voter models
% with anticonformity:
%   (1) one satisfying the balancing condition and
%   (2) two violating it.
% 
% For each model, the script computes both the annealed and quenched 
% dynamics using different preference distributions.
% 
% Generated data was used for Fig. 5.
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
phi_degenerate = 1;   % blue
phi_normal = phi_normal(p, 0.5, 1/40);  % green
phi_uniform = phi_uniform(p); % red
p_Bernoulli = [0;1];
phi_Bernoulli = phi_Bernoulli(p_mean); % purple

% Display means for verification
fprintf('Mean uniform: %.3f\n', trapz(p,p.*phi_uniform))
fprintf('Mean normal: %.3f\n', trapz(p,p.*phi_normal))
fprintf('Mean degenerate: %.3f\n',trapz(p,p.*phi_degenerate))
fprintf('Mean Bernoulli: %.3f\n', p_Bernoulli'*phi_Bernoulli)

%% Models
% ------------------------------------------------------------------
model_balanced = generalized_q_voter_model_with_anticonformity(6,6,1);
model_unbalanced1 = generalized_q_voter_model_with_anticonformity(6,2,1);
model_unbalanced2 = generalized_q_voter_model_with_anticonformity(6,6,4);

%% Solve ODEs for Quenched dynamics 
% ------------------------------------------------------------------
% Degenerate distribution
[t_vals, a_vals] = solve_ode_quenched_discrete_phi(T, p_mean, phi_degenerate, model_balanced, a_0);

subplot(3,2,2)
plot(t_vals, a_vals, 'b', 'LineWidth', 2);
ylabel('a');
title('Quenched dynamics');
ylim([0 1])
hold on;

[t_vals, a_vals] = solve_ode_quenched_discrete_phi(T, p_mean, phi_degenerate, model_unbalanced1, a_0);

subplot(3,2,4)
plot(t_vals, a_vals, 'b', 'LineWidth', 2);
xlabel('t'); ylabel('a');
ylim([0 1])
hold on;

[t_vals, a_vals] = solve_ode_quenched_discrete_phi(T, p_mean, phi_degenerate, model_unbalanced2, a_0);

subplot(3,2,6)
plot(t_vals, a_vals, 'b', 'LineWidth', 2);
xlabel('t'); ylabel('a');
ylim([0 1])
hold on;

% Uniform distribution
[t_vals, a_vals] = solve_ode_quenched_continuous_phi(T, p, phi_uniform, model_balanced, a_p0);

subplot(3,2,2)
plot(t_vals, a_vals, 'r', 'LineWidth', 2);

[t_vals, a_vals] = solve_ode_quenched_continuous_phi(T, p, phi_uniform, model_unbalanced1, a_p0);

subplot(3,2,4)
plot(t_vals, a_vals, 'r', 'LineWidth', 2);

[t_vals, a_vals] = solve_ode_quenched_continuous_phi(T, p, phi_uniform, model_unbalanced2, a_p0);

subplot(3,2,6)
plot(t_vals, a_vals, 'r', 'LineWidth', 2);

% Normal distribution
[t_vals, a_vals] = solve_ode_quenched_continuous_phi(T, p, phi_normal, model_balanced, a_p0);

subplot(3,2,2)
plot(t_vals, a_vals, 'g', 'LineWidth', 2);

[t_vals, a_vals] = solve_ode_quenched_continuous_phi(T, p, phi_normal, model_unbalanced1, a_p0);

subplot(3,2,4)
plot(t_vals, a_vals, 'g', 'LineWidth', 2);

[t_vals, a_vals] = solve_ode_quenched_continuous_phi(T, p, phi_normal, model_unbalanced2, a_p0);

subplot(3,2,6)
plot(t_vals, a_vals, 'g', 'LineWidth', 2);

% Bernoulli distribution
a_p0 = [1 1]';
[t_vals, a_vals] = solve_ode_quenched_discrete_phi(T, p_Bernoulli, phi_Bernoulli, model_balanced, a_p0);

subplot(3,2,2)
plot(t_vals, a_vals, 'm', 'LineWidth', 2);
hold off

[t_vals, a_vals, a_p_vals] = solve_ode_quenched_discrete_phi(T, p_Bernoulli, phi_Bernoulli, model_unbalanced1, a_p0);

subplot(3,2,4)
plot(t_vals, a_vals, 'm', 'LineWidth', 2);
hold off

[t_vals, a_vals, a_p_vals] = solve_ode_quenched_discrete_phi(T, p_Bernoulli, phi_Bernoulli, model_unbalanced2, a_p0);

subplot(3,2,6)
plot(t_vals, a_vals, 'm', 'LineWidth', 2);
hold off

%% Solve ODEs for annealed dynamics 
% ------------------------------------------------------------------
[t_vals, a_vals] = solve_ode_annealed(T, p_mean, model_balanced, a_0);

subplot(3,2,1)
plot(t_vals, a_vals, 'k', 'LineWidth', 2);
ylabel('a');
title('Annealed dynamics');
ylim([0 1])

[t_vals, a_vals] = solve_ode_annealed(T, p_mean, model_unbalanced1, a_0);

subplot(3,2,3)
plot(t_vals, a_vals, 'k', 'LineWidth', 2);
xlabel('t'); 
ylabel('a');
ylim([0 1])

[t_vals, a_vals] = solve_ode_annealed(T, p_mean, model_unbalanced2, a_0);

subplot(3,2,5)
plot(t_vals, a_vals, 'k', 'LineWidth', 2);
xlabel('t'); 
ylabel('a');
ylim([0 1])