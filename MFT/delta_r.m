
%% ------------------------------------------------------------------------
% Fig. 4c: Response is proportional to input
% Brief pulse amplitude sweep for the oculomotor-like integrator.
% Delta r is measured as the population firing rate in the 600-700 ms
% window immediately after pulse offset.
% -------------------------------------------------------------------------

%% Structure Parameters
neural_paprameters_eyeposotion;   % load parameters

%% Simulation parameters
num_trials = 10;

input_amps = 0.01:0.005:0.05;     % pulse amplitude per input channel
n_amp = numel(input_amps);

r_input_bkg = 0.004;              % background input

t_total = 1000;                   % ms
T_steps = t_total / dt;
timeline = (0:T_steps-1) * dt;

pulse_on_ms  = 500;
pulse_off_ms = 600;

post_on_ms  = 600;
post_off_ms = 700;

pulse_idx = timeline >= pulse_on_ms & timeline < pulse_off_ms;
post_idx  = timeline >= post_on_ms  & timeline < post_off_ms;

%% Firing-rate estimation for optional trace inspection
rate_win_ms = 200;
win_steps = max(1, round(rate_win_ms / dt));

%% Store results
delta_r_trials = nan(n_amp, num_trials);
pop_rate_example = nan(n_amp, T_steps);

%% ------------------------------------------------------------------------
% Amplitude sweep
% -------------------------------------------------------------------------
for ia = 1:n_amp

    I_ext = input_amps(ia);

    fprintf('Running input amplitude %.3f (%d/%d)\n', I_ext, ia, n_amp);

    for itrial = 1:num_trials

        %% Initialize variables
        v_ex = v_rest * ones(N, T_steps);

        s_input = r_input_bkg * ones(N, T_steps);
        s_input(:, pulse_idx) = r_input_bkg + I_ext;

        g_input_to_ex = zeros(N, T_steps);
        g_ex_to_ex = zeros(N, T_steps);

        is_refractory_ex = zeros(N, T_steps);
        is_spike_ex = zeros(N, T_steps);

        s_ex = zeros(N, T_steps);

        % Initial conductances
        g_input_to_ex(:, 1) = W_ff * s_input(:, 1);
        g_ex_to_ex(:, 1) = W * s_ex(:, 1);

        %% Time loop
        for t = 2:T_steps

            % Refractory neurons
            is_refractory = is_refractory_ex(:, t) == 1;
            v_ex(is_refractory, t) = v_rest;

            % Spikes
            is_spike = v_ex(:, t - 1) >= v_th;
            is_spike_ex(:, t) = is_spike;

            % Reset after spike
            v_ex(is_spike, t) = v_hold;

            if t < T_steps - t_refractory
                is_refractory_ex(is_spike, t+1 : t+t_refractory) = 1;
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

        %% Population response in the post-pulse window
        post_window_sec = sum(post_idx) * dt / 1000;

        post_rate_Hz = sum(is_spike_ex(:, post_idx), 'all') / ...
            (N * post_window_sec);

        delta_r_trials(ia, itrial) = post_rate_Hz;

        %% Optional example population trace
        if itrial == 1
            rate_from_spikes = movmean(double(is_spike_ex), win_steps, 2) * ...
                (1000 / dt);

            pop_rate_example(ia, :) = mean(rate_from_spikes, 1);
        end
    end
end

%% Summary statistics
delta_r_mean = mean(delta_r_trials, 2, 'omitnan');
delta_r_sd   = std(delta_r_trials, 0, 2, 'omitnan');

%% ------------------------------------------------------------------------
% Plot Fig. 4c
% -------------------------------------------------------------------------
fig = figure('Color', 'w', ...
    'Units', 'centimeters', ...
    'Position', [4 4 8 6]);

ax = axes(fig);
hold(ax, 'on');

line_col = [0.28 0.47 0.65];
mark_col = [0.65 0.78 0.90];

errorbar(ax, input_amps, delta_r_mean, delta_r_sd, ...
    '-o', ...
    'Color', line_col, ...
    'MarkerFaceColor', mark_col, ...
    'MarkerEdgeColor', line_col, ...
    'LineWidth', 1.6, ...
    'MarkerSize', 5, ...
    'CapSize', 5);

xlabel(ax, 'Input pulse amplitude');
ylabel(ax, '\Delta r (Hz)');
title(ax, 'Response is proportional to input', ...
    'FontWeight', 'bold');

xlim(ax, [0.01 0.05]);
ylim(ax, [0 30]);

xticks(ax, [0.01 0.02 0.03 0.04 0.05]);
yticks(ax, [0 10 20 30]);

set(ax, ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'FontSize', 9, ...
    'LineWidth', 0.8, ...
    'Layer', 'top');

