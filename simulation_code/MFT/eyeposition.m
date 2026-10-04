close all; clear all;
rng(1)% to ensure the reproducibility

%% Structure Parameters
neural_paprameters_eyeposotion; % load parameters

% num_trials = 3;
num_trials = 1;

%% Time Parameters
t_total = 20000;             % total time of a trial (ms)
T_steps = t_total / dt;      % total number of time steps
timeline = (1:T_steps) * dt;

%% Input Parameters
r_input_bkg = 0.004;         % background input

%% Firing-rate estimation from spikes
rate_win_ms = 200;                          % moving window size in ms
win_steps = max(1, round(rate_win_ms / dt)); % window size in time steps

%% Initialize Variables to Store Results
r_ex_all = zeros(num_trials, N, T_steps);    % firing rate from spikes, Hz
s_ex_all = zeros(num_trials, N, T_steps);    % synaptic activation
v_ex_all = zeros(num_trials, N, T_steps);    % membrane potential
spike_trials = cell(num_trials, 1);          % spike matrices
s_input_all = zeros(num_trials, N, T_steps); % input profile

for i = 1:num_trials

    v_ex = v_rest * ones(N, T_steps);
    r_input = zeros(N, T_steps);
    s_input = zeros(N, T_steps);

    g_input_to_ex = zeros(N, T_steps);
    g_ex_to_ex = zeros(N, T_steps);

    is_refractory_ex = zeros(N, T_steps);
    is_spike_ex = zeros(N, T_steps);

    s_ex = zeros(N, T_steps);

    %% Time Step Loop
    for t = 2:T_steps

        %% Input profile
        s_input(:, t) = r_input_bkg;

        if t > 0 && t < 50
            s_input(:, t) = 0.18/5 + r_input_bkg;
        end

        if t > 1000 && t < 1050
            s_input(:, t) = 0.05/5 + r_input_bkg;
        end

        if t > 5000 && t < 5050
            s_input(:, t) = 0.2/5 + r_input_bkg;
        end

        if t > 7000 && t < 7050
            s_input(:, t) = 0.25/5 + r_input_bkg;
        end

        %% Excitatory Neuron Dynamics

        % Refractory neurons
        is_refractory = is_refractory_ex(:, t) == 1;
        v_ex(is_refractory, t) = v_rest;

        % Spikes
        is_spike = v_ex(:, t - 1) >= v_th;
        is_spike_ex(:, t) = is_spike;

        % Reset after spike
        v_ex(is_spike, t) = v_hold;

        if t < T_steps - t_refractory
            is_refractory_ex(is_spike, t+1 : t + t_refractory) = 1;
        end

        % Synaptic activation
        s_ex(:, t) = s_ex(:, t - 1) - ...
            (s_ex(:, t - 1) * dt ./ tau_se) + ...
            rho * is_spike .* (1 - s_ex(:, t - 1));

        % Membrane potential update
        not_refractory = ~is_refractory & ~is_spike;

        v_ex(not_refractory, t) = v_ex(not_refractory, t - 1) + ...
            ((g_L * (E_l - v_ex(not_refractory, t - 1)) + ...
            (g_input_to_ex(not_refractory, t - 1) + ...
             g_ex_to_ex(not_refractory, t - 1)) .* ...
            (E_e - v_ex(not_refractory, t - 1)) + ...
            norm_noise * randn(sum(not_refractory), 1)) * dt / C_m);

        % Conductances
        g_input_to_ex(:, t) = W_ff * s_input(:, t);
        g_ex_to_ex(:, t) = W * s_ex(:, t);
    end

    %% Convert spikes to firing rate using moving window
    % is_spike_ex is N x T_steps.
    % movmean gives spikes/bin; multiply by 1000/dt to convert to Hz.
    rate_from_spikes = movmean(double(is_spike_ex), win_steps, 2) * (1000 / dt);

    %% Store trial results
    s_ex_all(i, :, :) = s_ex;
    r_ex_all(i, :, :) = rate_from_spikes;  % Hz
    v_ex_all(i, :, :) = v_ex;
    spike_trials{i} = is_spike_ex;
    s_input_all(i, :, :) = s_input;
end

%% ------------------------------------------------------------------------
% Plot input
% -------------------------------------------------------------------------
figure(1); clf; hold on

color = [0 0 0] ;

input_mean = squeeze(mean(s_input_all(1, :, :), 2));

plot(timeline, smooth(input_mean, 20), ...
    'Color', [color, 0.5], ...
    'LineWidth', 1);

xlabel('Time (ms)');
ylabel('Input');
title('Input');

axis([4000 t_total 0 r_input_bkg * 15]);



%% ------------------------------------------------------------------------
%% ------------------------------------------------------------------------
% plot  10 neurons
% Trial-averaged firing rates
% Single-neuron traces are smoothed with a 100-ms moving window
% -------------------------------------------------------------------------
neuron_ids = 1:10:100;
n_neurons = numel(neuron_ids);

single_smooth_ms = 200;
single_smooth_steps = max(1, round(single_smooth_ms / dt));

figFiring = figure('Color', 'w', ...
    'Units', 'centimeters', ...
    'Position', [3 3 18 9]);

ax = axes(figFiring);
hold(ax, 'on');

cmap = lines(n_neurons);

% Trial-averaged single-neuron firing rates
for k = 1:n_neurons
    thisRate = squeeze(mean(r_ex_all(:, neuron_ids(k), :), 1));
    thisRate = movmean(thisRate, single_smooth_steps);
    
    plot(ax, timeline, thisRate, ...
        'Color', [.7 .7 .7], ...
        'LineWidth', 0.8);
end

% Population mean firing rate, averaged across trials and neurons
popMean = squeeze(mean(mean(r_ex_all, 1), 2));

plot(ax, timeline, popMean, ...
    'k', ...
    'LineWidth', 2, ...
    'DisplayName', 'Population mean');

xlim(ax, [4000 10000]);
ylim(ax, 'tight');

xlabel(ax, 'Time (ms)');
ylabel(ax, 'Firing rate (Hz)');

set(ax, ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'FontSize', 8, ...
    'LineWidth', 0.6, ...
    'Layer', 'top');

legend(ax, {'Single neurons', 'Population mean'}, ...
    'Location', 'northeastoutside', ...
    'Box', 'off', ...
    'FontSize', 7);

% print(figFiring, '-dpdf', '-r600', 'Fig_FiringRates.pdf');
%% ------------------------------------------------------------------------
% Raster plot – first 50 neurons
% Uses one example trial
% -------------------------------------------------------------------------
neurons2plot = 1:50;
trial2plot = 1;

figRaster = figure('Color', 'w', ...
    'Units', 'centimeters', ...
    'Position', [3 3 18 6]);

ax = axes(figRaster);
hold(ax, 'on');

spmat = spike_trials{trial2plot};

% Highlight stimulation epochs
patch(ax, [0 100 100 0], ...
    [0 0 numel(neurons2plot)+1 numel(neurons2plot)+1], ...
    [0.85 0.85 0.85], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.5);

patch(ax, [2000 2100 2100 2000], ...
    [0 0 numel(neurons2plot)+1 numel(neurons2plot)+1], ...
    [0.85 0.85 0.85], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.5);

% Plot spikes
for ii = 1:numel(neurons2plot)
    n = neurons2plot(ii);
    ts = timeline(spmat(n, :) == 1) - 5000;

    plot(ax, ts, ii * ones(size(ts)), ...
        '.', ...
        'MarkerSize', 4, ...
        'Color', 'k');
end

set(ax, ...
    'YDir', 'reverse', ...
    'YTick', 1:numel(neurons2plot), ...
    'YTickLabel', neurons2plot, ...
    'XLim', [-1000 5000], ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'FontSize', 8, ...
    'LineWidth', 0.6, ...
    'Layer', 'top');

xlabel(ax, 'Time (ms)');
ylabel(ax, 'Neuron #');
title(ax, 'Spike raster');

% print(figRaster, '-dpdf', '-r600', 'Fig_Raster.pdf');
