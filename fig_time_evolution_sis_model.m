% Time evolution for SIS models with different preference distributions
% ------------------------------------------------------------------
% This script compares time evolution of SIS models with different 
% preference distributions.
% 
% Generated data was used for Fig. 7.
% ==================================================================

clc; clear; close all;
addpath('models');
addpath('distributions');
addpath('utils')

%% Parameters
% ------------------------------------------------------------------
T = [0:200];            % Time horizon
                        % Parameters for quenched dynamics
                        % Discretization of p, p in [0, 1]
N = 500;                
p = linspace(0.0001, 0.9999, N)'; 
                        % Initial a_p(t=0)
a_p0 = 0.01 * ones(N, 1);  
                        
p_mean = 0.5;           % Mean preference
a_0 = 0.01;             % Initial value of a(t=0)

%% Distributions
% ------------------------------------------------------------------
phi_uniform = phi_uniform(p);
phi_normal = phi_normal(p, 0.5, 1/40);  
phi_mixture_normals = 1;

% Display means for verification
fprintf('Mean uniform: %.3f\n', trapz(p,p.*phi_uniform))
fprintf('Mean normal: %.3f\n', trapz(p,p.*phi_normal))
fprintf('Mean degenerate: %.3f\n',p_mean)


%% Models
% ------------------------------------------------------------------
beta1 = 0.3;
beta2 = 0.6;
gamma = 0.2;
model_beta1 = sis_model(beta1,gamma);
model_beta2 = sis_model(beta2,gamma);

%% Solve ODEs for Quenched dynamics 
% ------------------------------------------------------------------
% Uniform distribution
[t_vals, a_vals] = solve_ode_quenched_continuous_phi(T, p, phi_uniform, model_beta2, a_p0);

subplot(3,3,[3,6,9])
plot(t_vals, a_vals, 'r', 'LineWidth', 2);
ylabel('a');
title('Quenched dynamics');
ylim([0 1])
hold on;

[t_vals, a_vals] = solve_ode_quenched_continuous_phi(T, p, phi_uniform, model_beta1, a_p0);

subplot(3,3,[2,5,8])
plot(t_vals, a_vals, 'r', 'LineWidth', 2);
xlabel('t'); ylabel('a');
ylim([0 1])
hold on;

% Normal distribution
[t_vals, a_vals] = solve_ode_quenched_continuous_phi(T, p, phi_normal, model_beta2, a_p0);

subplot(3,3,[3,6,9])
plot(t_vals, a_vals, 'g', 'LineWidth', 2);

[t_vals, a_vals] = solve_ode_quenched_continuous_phi(T, p, phi_normal, model_beta1, a_p0);

subplot(3,3,[2,5,8])
plot(t_vals, a_vals, 'g', 'LineWidth', 2);


%% Solve ODEs for annealed dynamics / one-point distribution
% ------------------------------------------------------------------
[t_vals, a_vals] = solve_ode_annealed(T, p_mean, model_beta2, a_0);

subplot(3,3,[3,6,9])
plot(t_vals, a_vals, 'b', 'LineWidth', 2);
xlabel('t');
ylabel('a');
title(['\beta=' num2str(beta2)]);
ylim([0 1])

[t_vals, a_vals] = solve_ode_annealed(T, p_mean, model_beta1, a_0);

subplot(3,3,[2,5,8])
plot(t_vals, a_vals, 'b', 'LineWidth', 2);
title(['\beta=' num2str(beta1)]);
xlabel('t'); 
ylabel('a');
ylim([0 1])

%% Plot distributions
% ------------------------------------------------------------------
subplot(3,3,1)
plot(p, phi_uniform, 'b', 'LineWidth', 2);
xlabel('p'); 
ylabel('\phi(p)');
title('Preference distribution');

subplot(3,3,4)
plot(p, phi_normal, 'g', 'LineWidth', 2);
xlabel('p'); 
ylabel('\phi(p)');
title('Preference distribution');

subplot(3,3,7)
plot(p_mean, phi_mixture_normals, '.r', 'LineWidth', 2);
xlabel('p'); 
ylabel('\phi(p)');
title('Preference distribution');

