%% ORBIT DETERMINATION OF A MARS ORBITER FROM EARTH-BASED RADIOMETRIC TRACKING
% Batch weighted least-squares estimation of the state of a Mars orbiter
% (plus GM of Mars and atmospheric density) from range and range-rate
% measured by a rotating ground station on the Earth, 2.1 AU away.
%
%   1. State and parameter estimation with a priori information
%   2. Information content of range vs range-rate
%   3. Covariance propagation: altitude and its 1-sigma over 3 sidereal days
%   4. Drag coefficient added to the solve-for parameters
%
% Simplified 2D problem: planets fixed, orbit in the Earth equatorial plane,
% no light-time correction, forces = Mars monopole + exponential-atmosphere drag.
% Author: Girolamo Castronovo - Sapienza University of Rome, 2024

clear; clc; close all;
addpath('functions');

SAVE_FIGURES = false;          % true -> export PNGs to /figures

%% Model parameters (km, s, kg)
p.R_m  = 3389;                 % Mars radius [km]
p.R_e  = 6378;                 % Earth radius [km]
p.AU   = 149597871;            % astronomical unit [km]
p.t_sd = 86164;                % sidereal day [s]
p.z_0  = 35;                   % reference altitude of the density model [km]
p.H    = 14.7;                 % density scale height [km]
p.area = 37.5e-6;              % cross-section [km^2]
p.m_sc = 2000;                 % spacecraft mass [kg]
p.C_d  = 2;                    % drag coefficient (Tasks 1-3)

p.d_earth        = 2.1*p.AU;   % Mars-Earth distance [km]
p.station_angle0 = deg2rad(200);  % station longitude at t0 in the Mars-centred frame

p.sigma_range = 5e-4;          % range accuracy [km]
p.sigma_rate  = 4e-8;          % range-rate accuracy [km/s]

n_iter = 5;

%% Observations: [t (s), range (km), range-rate (km/s)]
data = load(fullfile('data', 'observables.txt'));
obs.t     = data(:, 1);
obs.range = data(:, 2);
obs.rate  = data(:, 3);
fprintf('%d epochs, %.1f h of tracking\n', numel(obs.t), obs.t(end)/3600);

%% A priori state and uncertainties: [x, y, vx, vy, GM, rho_0]
X0     = [2506.16, 2505.29, -2.53608, 2.61703, 42830, 0.0];
sigma0 = [0.3,     0.3,     5e-5,     5e-5,    3,     0.03];
names  = {'x [km]', 'y [km]', 'vx [km/s]', 'vy [km/s]', 'GM [km^3/s^2]', 'rho_0'};

%% TASK 1 - Batch weighted least squares
fprintf('\nTask 1: batch least squares\n');
est = batch_least_squares(X0, sigma0, obs, p, n_iter);
sigma_post = sqrt(diag(est.P))';

disp(table(est.X', (est.X - X0)', sigma_post', ...
     'RowNames', names, 'VariableNames', {'Estimate', 'Correction', 'Sigma'}));

fig = plot_residuals(obs.t, est, [1, 2, n_iter], p, 'Residuals - Task 1');
save_figure(fig, 'residuals', SAVE_FIGURES);

%% TASK 2 - Which observable carries more information?
% Covariance of a range-only and of a range-rate-only solution (no a priori),
% built from the H matrix of the last iteration.
H_range = est.H(1:2:end, :);
H_rate  = est.H(2:2:end, :);
P_range = inv(H_range'*H_range / p.sigma_range^2);
P_rate  = inv(H_rate'*H_rate   / p.sigma_rate^2);

sig_range = sqrt(diag(P_range));
sig_rate  = sqrt(diag(P_rate));
fprintf('\nTask 2: 1-sigma of range-only vs range-rate-only solutions\n');
disp(table(sig_range, sig_rate, sig_range./sig_rate, 'RowNames', names, ...
     'VariableNames', {'Range only', 'Range-rate only', 'Ratio'}));

% The residuals of Task 1 give an empirical noise level: check that the
% conclusion does not depend on the declared accuracies
rms_range = sqrt(mean(est.res_range(:, end).^2));
rms_rate  = sqrt(mean(est.res_rate(:, end).^2));
ratio_emp = sqrt(diag(inv(H_range'*H_range/rms_range^2))) ./ sqrt(diag(inv(H_rate'*H_rate/rms_rate^2)));
fprintf('Using the post-fit RMS as noise level, range-only sigmas are still %.1f-%.1f times larger.\n', ...
        min(ratio_emp), max(ratio_emp));

%% TASK 3 - Altitude and 1-sigma over 3 sidereal days (covariance mapping)
t3   = linspace(0, 3*p.t_sd, 20000)';
opts = odeset('RelTol', 5e-14, 'AbsTol', 5e-14);
[~, S] = ode113(@(t, s) mars_dynamics_stm(t, s, p, 6), t3, [est.X'; reshape(eye(6), [], 1)], opts);

alt       = zeros(numel(t3), 1);
sigma_alt = zeros(numel(t3), 1);
for i = 1:numel(t3)
    r   = norm(S(i, 1:2));
    Phi = reshape(S(i, 7:end), 6, 6);
    P_t = Phi*est.P*Phi';                      % covariance mapped to time t
    g   = [S(i, 1:2)/r, 0, 0, 0, 0];           % gradient of r = |(x, y)| (R_m is exact)
    alt(i)       = r - p.R_m;
    sigma_alt(i) = sqrt(g*P_t*g');
end
fprintf('\nTask 3: altitude %.0f-%.0f km, 1-sigma %.2f-%.2f mm\n', ...
        min(alt), max(alt), 1e6*min(sigma_alt), 1e6*max(sigma_alt));

fig = figure('Name', 'Altitude and uncertainty');
tiledlayout(2, 1, 'TileSpacing', 'compact');
nexttile; plot(t3/3600, alt, 'Color', [0 0.6 0.8]); grid on;
ylabel('altitude [km]'); title('Altitude over 3 sidereal days');
nexttile; plot(t3/3600, sigma_alt*1e6, 'Color', [0.85 0.33 0.10]); grid on;
xlabel('time [h]'); ylabel('1\sigma altitude [mm]'); title('Mapped altitude uncertainty');
save_figure(fig, 'altitude_uncertainty', SAVE_FIGURES);

%% TASK 4 - Drag coefficient as additional solve-for parameter (a priori 1 +/- 0.5)
fprintf('\nTask 4: C_D added to the state\n');
X0_4     = [X0, 1];
sigma0_4 = [sigma0, 0.5];
names_4  = [names, {'C_D'}];

est4 = batch_least_squares(X0_4, sigma0_4, obs, p, n_iter);
sigma_post_4 = sqrt(diag(est4.P))';
disp(table(est4.X', sigma_post_4', (sigma_post_4./sigma0_4)', 'RowNames', names_4, ...
     'VariableNames', {'Estimate', 'Sigma', 'Sigma_over_apriori'}));

corr_rho_cd = est4.P(6,7)/(sigma_post_4(6)*sigma_post_4(7));
fprintf('Correlation rho_0 - C_D: %.3f\n', corr_rho_cd);
fprintf('rho_0*C_D: %.5e (Task 1, C_D = %g) vs %.5e (Task 4)\n', ...
        est.X(6)*p.C_d, p.C_d, est4.X(6)*est4.X(7));
% Drag only depends on the product rho_0*C_D: the data constrain the product,
% not the two parameters separately (C_D stays at its a priori uncertainty).

fig = plot_residuals(obs.t, est4, [1, n_iter], p, 'Residuals - Task 4');
save_figure(fig, 'residuals_task4', SAVE_FIGURES);


%% ===================== Local functions =====================
function save_figure(fig, name, save_flag)
if save_flag
    if ~exist('figures', 'dir'), mkdir('figures'); end
    exportgraphics(fig, fullfile('figures', [name '.png']), 'Resolution', 200);
end
end
