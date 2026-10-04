%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Simulate 100 trials with the same parameters, estimate μ and σ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;


% ------------------------------------------------------------------------
neuron_parameters;   % Load neuron and network parameters.
% ------------------------------------------------------------------------
N_input = 100;

%% Simulation parameters
NUM_TRIALS     = 100;     % Number of trials
input_mean     = 0.01;    % DC drive amplitude
RAMP_THRESHOLD = 25;      % Hz – population mean FR threshold (θ)
alpha          = 1.03;    % Unused parameter
t_total        = 1000;    % ms per trial
% ------------------------------------------------------------------------

T_steps  = t_total / dt;
time_vec = (0:T_steps-1)*dt;
dt_window = 50;
% For μ, σ estimation
all_dX      = [];         % collect all ΔX before the hit across trials
all_dX_dt   = [];         % Firing-rate increments per unit time
X0_all      = nan(NUM_TRIALS,1);
theta_all   = nan(NUM_TRIALS,1);

% For cross-time stats
ramp_store  = nan(NUM_TRIALS, T_steps);
crossTimes  = nan(NUM_TRIALS, 1);

% Small margin to stay away from boundaries when estimating μ, σ
eps_low   = 0.5;          % Hz above 0 (reflecting bound)
eps_high  = 0.5;          % Hz below θ

fprintf('Running %d trials …\n', NUM_TRIALS);

for i = 1:NUM_TRIALS
i
    % ---------- re-load / re-initialize per-trial state -----------------
    neuron_parameters;

    % Initialize external input weights.
    mean_wff = 0.001/10;
    W_ff     = (rand(N,N_input) < sparsity);
    W_ff     = W_ff/mean(W_ff,"all")*mean_wff;

    v_ex            = v_rest * ones(N,T_steps);
    g_input_to_ex   = zeros(N,T_steps);
    g_ex_to_ex      = zeros(N,T_steps);
    is_refractory_ex= zeros(N,T_steps);
    is_spike_ex     = zeros(N,T_steps);
    r_ex            = zeros(N,T_steps);
    s_ex            = zeros(N,T_steps);
    s_input         = zeros(N_input,T_steps);

    % ----------------------------- TIME LOOP -----------------------------
    for t = 2:T_steps
        % Poisson external spikes
        is_spike_input = rand(N_input,1) < 0.005;
        s_input(:, t) = s_input(:, t-1) - (s_input(:, t-1)*dt/20) + ...
                        rho * is_spike_input .* (1 - s_input(:, t-1));
        % s_input(:, t) = 0.02;

        %---------------- Refractory handling ----------------------------%
        is_refrac = is_refractory_ex(:,t)==1;
        v_ex(is_refrac,t) = v_rest;

        %---------------- Spike detection --------------------------------%
        is_spike           = v_ex(:,t-1) >= v_th;
        is_spike_ex(:,t)   = is_spike;
        v_ex(is_spike,t)   = v_hold;
        if t < T_steps - t_refractory
            is_refractory_ex(is_spike,t+1:t+t_refractory) = 1;
        end

        %---------------- Synaptic gating --------------------------------%
        s_ex(:,t) = s_ex(:,t-1) - (dt/tau_se)*s_ex(:,t-1) ...
                    + rho*is_spike.*(1 - s_ex(:,t-1));

        %---------------- Low‑pass firing rate ---------------------------%
        r_ex(:,t) = r_ex(:,t-1) + (is_spike/dt - r_ex(:,t-1))*(dt/tau_w);

        %---------------- Membrane update (non‑refrac / non‑spike) -------%
        free = ~is_refrac & ~is_spike;
        v_ex(free,t) = v_ex(free,t-1) + ...
            ( g_L*(E_l - v_ex(free,t-1)) + ...
              (g_input_to_ex(free,t-1)+g_ex_to_ex(free,t-1)).* ...
              (E_e - v_ex(free,t-1)) + ...
              norm_noise*randn(sum(free),1) ) * dt / C_m;

        %---------------- Conductances -----------------------------------%
        g_input_to_ex(:,t) = W_ff * s_input(:,t);
        g_ex_to_ex(:,t)    = W    * s_ex(:,t);
    end
    % --------------------------- end time loop ---------------------------

    % Decision variable: population mean FR (Hz)
    mean_FR          = 1000*mean(r_ex(:,1:T_steps),1);   % Hz
    ramp_store(i,:)  = mean_FR;
    X                = mean_FR;                       % column

    % Initial value
    X0_all(i) = X(1);

    % First‑passage to θ
    idx = find(X > RAMP_THRESHOLD, 1, 'first');
    if ~isempty(idx)
        crossTimes(i) = idx * dt;                        % ms
        theta_all(i)  = X(idx);                          % empirical value at hit
    else
        theta_all(i)  = NaN;
    end

    % --------------- Collect increments for μ, σ estimation -------------
    % Use only pre-hit samples; if never hit, use first 500 points
    if ~isempty(idx)
        idx_use = 1:(idx-1);
    else
        idx_use = 1:500;
    end

    % keep away from the bounds
    interior = X(idx_use) > eps_low & (RAMP_THRESHOLD - X(idx_use)) > eps_high;

    idx_use  = idx_use(interior);
    if numel(idx_use) >= 2
        data = X(idx_use);
        len = length(data);
        bin_length = floor(len/dt_window);
        data = reshape(data(1:dt_window*bin_length), [dt_window, bin_length ]);
        data = mean(data);
        dX  = diff(data);
        all_dX = [all_dX; dX(:)];
       
    end
end

%% ---------------- μ & σ estimation from increments -----------------------

mu_hat     = mean(all_dX) / dt_window;
sigma_hat  = sqrt( var(all_dX - mu_hat*dt_window) / dt_window );

fprintf('\nEstimated parameters (from increments):\n');
fprintf('  mu     = %.4f (units of X/s)\n', mu_hat);
fprintf('  sigma  = %.4f (units of X/s^{1/2})\n', sigma_hat);
fprintf('  X0     = %.4f ± %.4f  (mean ± SD over trials)\n', ...
        mean(X0_all,'omitnan'), std(X0_all,'omitnan'));
fprintf('  theta  = %.4f ± %.4f  (empirical at hit)\n', ...
        mean(theta_all,'omitnan'), std(theta_all,'omitnan'));

% ------------------- quick diagnostics ----------------------------------
%% ---------------- μ & σ estimation from increments -----------------------
% all_dX was computed from 50-ms coarse-grained population-rate bins.
% Therefore mu_hat is in Hz/ms and sigma_hat is in Hz/sqrt(ms).

mu_hat     = mean(all_dX) / dt_window;
sigma_hat  = sqrt(var(all_dX - mu_hat * dt_window) / dt_window);

fprintf('\nEstimated effective DDM parameters:\n');
fprintf('  mu     = %.5f Hz/ms\n', mu_hat);
fprintf('  sigma  = %.5f Hz/sqrt(ms)\n', sigma_hat);
fprintf('  X0     = %.4f ± %.4f Hz\n', ...
        mean(X0_all,'omitnan'), std(X0_all,'omitnan'));
fprintf('  theta  = %.4f ± %.4f Hz\n', ...
        mean(theta_all,'omitnan'), std(theta_all,'omitnan'));

%% ------------------------------------------------------------------------
% Reproduce DDM vs RNN comparison panel
% -------------------------------------------------------------------------

n_plot = 10;                         % number of example trajectories shown
theta  = RAMP_THRESHOLD;             % decision boundary
x0     = 0;                           % start near zero, as in the figure

ddm_trials = nan(n_plot, T_steps);

rng(2);

for k = 1:n_plot

    X = nan(1, T_steps);
    X(1) = x0;

    for t = 2:T_steps

        dX = mu_hat * dt + sigma_hat * sqrt(dt) * randn;
        X(t) = X(t-1) + dX;

        % Reflecting lower boundary at zero
        if X(t) < 0
            X(t) = -X(t);
        end

        % Absorbing upper decision boundary
        if X(t) >= theta
            X(t) = theta;
            X(t+1:end) = NaN;
            break;
        end
    end

    ddm_trials(k, :) = X;
end

%% Prepare RNN example trajectories
rnn_trials = ramp_store(1:n_plot, :);

for k = 1:n_plot
    idx_hit = find(rnn_trials(k, :) >= theta, 1, 'first');

    if ~isempty(idx_hit)
        rnn_trials(k, idx_hit) = theta;
        rnn_trials(k, idx_hit+1:end) = NaN;
    end
end

%% Mean drift line
t_hit_line = theta / mu_hat;
t_line = 0:dt:min(t_total, t_hit_line);
mean_drift_line = mu_hat * t_line;

%% Plot
figDDM = figure('Color', 'w', ...
    'Units', 'centimeters', ...
    'Position', [3 3 18 7]);

% ---------------- DDM panel ----------------
ax1 = subplot(1,2,1);
hold(ax1, 'on');

for k = 1:n_plot
    plot(ax1, time_vec, ddm_trials(k, :), ...
        'Color', [0.65 0.65 0.65], ...
        'LineWidth', 0.9);
end

plot(ax1, t_line, mean_drift_line, ...
    'k', ...
    'LineWidth', 2.2);

yline(ax1, theta, '--', ...
    'Color', [0.35 0.35 0.35], ...
    'LineWidth', 1.5);

xlabel(ax1, 'Time (ms)');
ylabel(ax1, 'Decision variable X');
title(ax1, 'DDM', ...
    'FontWeight', 'bold');

text(ax1, 430, theta + 2.0, 'Decision boundary', ...
    'FontSize', 9);

text(ax1, 430, 2.0, 'Reflecting boundary', ...
    'FontSize', 9);

legend(ax1, {'DDM trials', 'Mean drift'}, ...
    'Location', 'southeast', ...
    'Box', 'on');

xlim(ax1, [0 t_total]);
ylim(ax1, [0 theta + 5]);

set(ax1, ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'FontSize', 9, ...
    'LineWidth', 0.9, ...
    'Layer', 'top');

% ---------------- RNN panel ----------------
ax2 = subplot(1,2,2);
hold(ax2, 'on');

for k = 1:n_plot
    plot(ax2, time_vec, rnn_trials(k, :), ...
        'Color', [0.62 0.84 0.84], ...
        'LineWidth', 0.9);
end

plot(ax2, t_line, mean_drift_line, ...
    'Color', [0.20 0.45 0.58], ...
    'LineWidth', 2.2);

yline(ax2, theta, '--', ...
    'Color', [0.35 0.35 0.35], ...
    'LineWidth', 1.5);

xlabel(ax2, 'Time (ms)');
ylabel(ax2, 'Firing rate (Hz)');
title(ax2, 'RNN', ...
    'FontWeight', 'bold');

legend(ax2, {'RNN trials', 'Mean drift'}, ...
    'Location', 'southeast', ...
    'Box', 'on');

xlim(ax2, [0 t_total]);
ylim(ax2, [0 theta + 5]);

set(ax2, ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'FontSize', 9, ...
    'LineWidth', 0.9, ...
    'Layer', 'top');

%% Save
print(figDDM, '-dpdf', '-r600', 'Fig4d_DDM_vs_RNN.pdf');
print(figDDM, '-dsvg', '-r600', 'Fig4d_DDM_vs_RNN.svg');
print(figDDM, '-dpng', '-r600', 'Fig4d_DDM_vs_RNN.png');