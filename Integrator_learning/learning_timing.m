%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decision network coupled to a recurrent integrator ("timer") that learns
% to shape its ramping activity through reward-modulated Hebbian plasticity.
%
% Model:
%   - Input layer (two stimuli)
%   - Two-choice decision circuit (excitatory & inhibitory populations)
%   - Recurrent integrator (timer) network
%
% Once the timer's mean activity crosses threshold, the crossing is detected
% online and, after a fixed delay, the decision-to-timer input is reduced for
% the remainder of the trial.
%
% Trials run as three phases: test (no weight update) -> learning
% (weight update on) -> test (no weight update).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
neuron_parameters_timing;

%% GLOBAL PARAMETERS ------------------------------------------------------
 % t_target = 700;            % Reward target time
 % w_ff_mean = 0.00003;       % Initial feedforward weight to timer

% 
 t_target = 1500;            % parameter pakage for short to long example
 w_ff_mean = 0.000055;       




threshould = 32;            % Threshold firing rate in Hz

% Decision-to-timer input gain before and after the delayed reduction
decision_gain_before = 1.0;     % Gain before reduction
decision_gain_after  = 0.2;     % Gain after reduction
decision_reduce_delay_ms = 50;  % Delay after threshold crossing before reduction (ms)

%% STRUCTURE PARAMETERS ---------------------------------------------------
N_integrator  = 100;        % Integrator neurons
npp           = 50;         % Neurons per decision population
n_exc_total   = 2*npp;      % Total excitatory decision neurons
n_inh_total   = 2*npp;      % Total inhibitory decision neurons
N_input       = 2*npp;      % Input neurons

% Trial schedule: test (no learning) -> learning -> test (no learning)
num_test_before     = 2;   % Test trials before learning (no weight update)
num_learning_trials = 100;  % Learning trials (weight update on)
num_test_after      = 2;   % Test trials after learning (no weight update)
num_trials          = num_test_before + num_learning_trials + num_test_after;

T         = 3000 / dt;
t_total   = 3001;
t_steps   = t_total / dt;

%% DECISION-NET WEIGHTS ---------------------------------------------------
rho       = 0.5;
wij_mean  = 0.002;          % E to E
m_mean    = 0.0015;         % I to E
p_mean    = 0.0015;         % E to I
k_rec     = 0.00;           % I to I

W_ji = zeros(n_exc_total, n_exc_total);
M_ki = zeros(n_exc_total, n_inh_total);
P_ik = zeros(n_inh_total, n_exc_total);
K_jk = zeros(n_inh_total, n_inh_total);

W_ji(1:npp,1:npp)                 = wij_mean*(rand(npp)<rho);
W_ji(npp+1:2*npp,npp+1:2*npp)     = wij_mean*(rand(npp)<rho);

M_ki(1:npp,npp+1:2*npp)           = m_mean *(rand(npp)<rho);
M_ki(npp+1:2*npp,1:npp)           = m_mean *(rand(npp)<rho);

P_ik(1:npp,1:npp)                 = p_mean *(rand(npp)<rho);
P_ik(npp+1:2*npp,npp+1:2*npp)     = p_mean *(rand(npp)<rho);

K_jk(1:npp,1:npp)                 = k_rec  *(rand(npp)<rho);
K_jk(npp+1:2*npp,npp+1:2*npp)     = k_rec  *(rand(npp)<rho);

rec_identity_decision = zeros(n_exc_total,n_exc_total);
rec_identity_decision(1:npp,1:npp)             = W_ji(1:npp,1:npp)>0;
rec_identity_decision(npp+1:2*npp,npp+1:2*npp) = W_ji(npp+1:2*npp,npp+1:2*npp)>0;

%% INTEGRATOR-NET WEIGHTS -------------------------------------------------
W_ff = w_ff_mean*rand(N_integrator,n_exc_total) .* ...
    (rand(N_integrator,n_exc_total)<rho);

W_ff(:,1:npp) = W_ff(:,1:npp) * ...
    w_ff_mean/4/mean(W_ff(:,1:npp),"all");

W_ff(:,npp+1:end) = W_ff(:,npp+1:end) * ...
    w_ff_mean/4/mean(W_ff(:,npp+1:end),"all");

W_ff_id = W_ff > 0;

w_rec_mean  = 0.000156;
W_rec       = w_rec_mean*rand(N_integrator,N_integrator) .* ...
    (rand(N_integrator,N_integrator)<rho);

rec_identity_timer = W_rec > 0;

w_ff_decision_mean = 0.001;
W_ff_decision = zeros(n_exc_total,N_input);

W_ff_decision(1:npp,1:N_input/2) = ...
    w_ff_decision_mean*rand(npp,N_input/2);

W_ff_decision(npp+1:n_exc_total,N_input/2+1:N_input) = ...
    w_ff_decision_mean*rand(npp,N_input/2);

%% REWARD PARAMETERS ------------------------------------------------------
delay_time  = 25/dt;
t_reward    = t_target + delay_time;

rew_vect = zeros(t_total,1);
del_time = round(delay_time/2);

for m = t_reward
    rew_vect(m-del_time:m+del_time) = 1;
    rew_vect(m+del_time+1)          = 2;
end

%% BOOK-KEEPING ARRAYS ----------------------------------------------------
sc_R_it               = zeros(N_integrator,t_steps,num_trials+1);
plot_R_it             = zeros(1,t_steps,num_trials+1);
rec_vect              = zeros(1,num_trials+1);

trial_type            = zeros(num_trials,1);
trial_type(1:2:end)   = 1;
trial_type(2:2:end)   = 2;

choice1_trials        = find(trial_type==1);
choice2_trials        = find(trial_type==2);

first_choice1_trial   = choice1_trials(1);
first_choice2_trial   = choice2_trials(1);
last_choice1_trial    = choice1_trials(end);
last_choice2_trial    = choice2_trials(end);

mean_activity_inh_group1  = zeros(num_trials,t_steps);
mean_activity_inh_group2  = zeros(num_trials,t_steps);
mean_T_p_timer_all_trials = zeros(num_trials,t_steps);
mean_T_d_timer_all_trials = zeros(num_trials,t_steps);

decision_threshold_cross_time_all = NaN(num_trials,1);
decision_input_reduce_time_all    = NaN(num_trials,1);

first_trial_choice1_input    = [];
first_trial_choice1_decision = [];
first_trial_choice1_r_timer  = [];

first_trial_choice2_input    = [];
first_trial_choice2_decision = [];
first_trial_choice2_r_timer  = [];

last_trial_choice1_input     = [];
last_trial_choice1_decision  = [];
last_trial_choice1_r_timer   = [];

last_trial_choice2_input     = [];
last_trial_choice2_r_ex      = [];
last_trial_choice2_r_timer   = [];

time_activity_cross = zeros(num_trials,1);
time_trace_cross    = zeros(num_trials,1);
meanwin             = zeros(num_trials,1);

%% MAIN PROGRAM -----------------------------------------------------------
for l = 1:num_trials

    choice = trial_type(l);

    % Learning (weight update) is enabled only during the middle phase.
    learning_on = (l > num_test_before) && ...
        (l <= num_test_before + num_learning_trials);

    %% Initialization -----------------------------------------------------
    trace_refractory = 0;

    % Excitatory neurons in decision network
    v_ex = v_rest * ones(n_exc_total, t_steps);
    r_ex = zeros(n_exc_total, t_steps);
    s_ex = zeros(n_exc_total, t_steps);
    s_input_decision = zeros(n_exc_total, t_steps);
    is_refractory_ex = zeros(n_exc_total, t_steps);
    is_spike_ex = zeros(n_exc_total, t_steps);

    % Inhibitory neurons in decision network
    v_inh = v_rest * ones(n_inh_total, t_steps);
    r_inh = zeros(n_inh_total, t_steps);
    s_inh = zeros(n_inh_total, t_steps);
    is_refractory_inh = zeros(n_inh_total, t_steps);
    is_spike_inh = zeros(n_inh_total, t_steps);

    % Timer network neurons
    v_timer = zeros(N_integrator, t_steps) + v_rest;
    r_timer = zeros(N_integrator, t_steps);
    s_timer = zeros(N_integrator, t_steps);
    g_decision_to_timer = zeros(N_integrator, t_steps);
    g_timer_to_timer = zeros(N_integrator, t_steps);
    is_refractory_timer = zeros(N_integrator, t_steps);
    is_spike_timer = zeros(N_integrator, t_steps);

    % Eligibility traces
    T_p_decision_to_timer = zeros(N_integrator, n_exc_total);
    T_d_decision_to_timer = zeros(N_integrator, n_exc_total);

    T_p_input_to_decision = zeros(n_exc_total, N_input);
    T_d_input_to_decision = zeros(n_exc_total, N_input);

    mean_T_p_timer = zeros(1, t_steps);
    mean_T_d_timer = zeros(1, t_steps);

    % Threshold-gated delayed reduction of decision-to-timer input
    decision_input_reduced = false;
    decision_threshold_cross_time = NaN;
    decision_input_reduce_time = NaN;
    decision_input_gain = decision_gain_before;

    %% Time Step Loop -----------------------------------------------------
    for t = 2:t_steps

        %% Input Layer Dynamics ------------------------------------------
        if choice == 1 && t < 500
            s_input_decision(1:npp, t) = 0.01*rand(npp, 1);
            s_input_decision(npp+1:end, t) = 0.0001*rand(npp, 1);
        end

        if choice == 2 && t < 500
            s_input_decision(npp+1:end, t) = 0.01*rand(npp, 1);
            s_input_decision(1:npp, t) = 0.0001*rand(npp, 1);
        end

        %% Decision Network Dynamics -------------------------------------

        % Excitatory neurons
        is_refractory = is_refractory_ex(:, t) == 1;
        v_ex(is_refractory, t) = v_rest;

        is_spike = (v_ex(:, t-1) >= v_th);
        is_spike_ex(:, t) = is_spike;

        not_spike = (~is_spike & ~is_refractory);

        v_ex(is_spike, t) = v_hold;

        if t < t_steps - t_refractory
            is_refractory_ex(is_spike, t+1:t+t_refractory) = 1;
        end

        s_ex(:, t) = s_ex(:, t-1) - ...
            (s_ex(:, t-1) * dt / tau_se_decision_input) + ...
            rho_syn * is_spike .* (1 - s_ex(:, t-1));

        v_ex(not_spike, t) = v_ex(not_spike, t-1) + ...
            ((randn(sum(not_spike),1) * norm_noise) + ...
            g_L * (E_l - v_ex(not_spike, t-1)) + ...
            (W_ji(not_spike, :) * s_ex(:, t-1)) .* ...
            (E_e - v_ex(not_spike, t-1)) + ...
            (M_ki(not_spike, :) * s_inh(:, t-1)) .* ...
            (E_i - v_ex(not_spike, t-1)) + ...
            W_ff_decision(not_spike, :) * s_input_decision(:, t-1) .* ...
            (E_e - v_ex(not_spike, t-1))) * (dt / C_m);

        r_ex(:, t) = r_ex(:, t-1) + ...
            (is_spike / dt - r_ex(:, t-1)) * (dt / tau_w);

        %% Inhibitory neurons --------------------------------------------

        is_refractory = is_refractory_inh(:, t) == 1;
        v_inh(is_refractory, t) = v_rest;

        is_spike = (v_inh(:, t-1) >= v_th_i);
        is_spike_inh(:, t) = is_spike;

        not_spike = (~is_spike & ~is_refractory);

        v_inh(is_spike, t) = v_hold;

        if t < t_steps - t_refractory
            is_refractory_inh(is_spike, t+1:t+t_refractory) = 1;
        end

        s_inh(:, t) = s_inh(:, t-1) - ...
            (s_inh(:, t-1) * dt / tau_si) + ...
            rho_syn * is_spike .* (1 - s_inh(:, t-1));

        v_inh(not_spike, t) = v_inh(not_spike, t-1) + ...
            ((randn(sum(not_spike),1) * norm_noise) + ...
            g_L * (E_l - v_inh(not_spike, t-1)) + ...
            (P_ik(not_spike, :) * s_ex(:, t-1)) .* ...
            (E_e - v_inh(not_spike, t-1)) + ...
            (K_jk(not_spike, :) * s_inh(:, t-1)) .* ...
            (E_i - v_inh(not_spike, t-1))) * (dt / C_m);

        r_inh(:, t) = r_inh(:, t-1) + ...
            (is_spike / dt - r_inh(:, t-1)) * (dt / tau_w);

        %% Timer Network Dynamics ----------------------------------------

        % Decision-to-timer input (gain reduced after the delayed crossing)
        g_decision_to_timer(:, t) = decision_input_gain * ...
            (W_ff * s_ex(:, t - 1));

        % Recurrent timer input
        g_timer_to_timer(:, t) = W_rec * s_timer(:, t - 1);

        % Total timer conductance
        g_total_timer = g_decision_to_timer(:, t) + g_timer_to_timer(:, t);

        % Timer membrane dynamics
        is_refractory = is_refractory_timer(:, t) == 1;
        v_timer(is_refractory, t) = v_rest;

        is_spike = (v_timer(:, t - 1) >= v_th);
        is_spike_timer(:, t) = is_spike;

        not_spike = (~is_spike & ~is_refractory);

        v_timer(is_spike, t) = v_hold;

        if t < t_steps - t_refractory
            is_refractory_timer(is_spike, t + 1:t + t_refractory) = 1;
        end

        s_timer(:, t) = s_timer(:, t - 1) - ...
            (s_timer(:, t - 1) * dt / tau_se_timer) + ...
            rho_syn * is_spike .* (1 - s_timer(:, t - 1));

        v_timer(not_spike, t) = v_timer(not_spike, t - 1) + ...
            ((randn(sum(not_spike), 1) * norm_noise) + ...
            g_L * (E_l - v_timer(not_spike, t - 1)) + ...
            g_total_timer(not_spike) .* ...
            (E_e - v_timer(not_spike, t - 1))) * (dt / C_m);

        r_timer(:, t) = r_timer(:, t - 1) + ...
            (is_spike / dt - r_timer(:, t - 1)) * (dt / tau_w);

        % Detect the threshold crossing, then reduce the input after the delay
        mean_activity_timer_now = mean(1000 * r_timer(:, t));
        current_time_ms = t * dt;

        if isnan(decision_threshold_cross_time) && mean_activity_timer_now >= threshould
            decision_threshold_cross_time = current_time_ms;
            decision_input_reduce_time = decision_threshold_cross_time + decision_reduce_delay_ms;
        end

        if ~decision_input_reduced && ~isnan(decision_input_reduce_time) && ...
                current_time_ms >= decision_input_reduce_time
            decision_input_reduced = true;
            decision_input_gain = decision_gain_after;
        end

        %% Learning -------------------------------------------------------
        if t > D

            r_pre = r_ex(:, t - D) .* (r_ex(:, t - D) > 0.01);
            r_post = r_timer(:, t - 1) .* (r_timer(:, t - 1) > 0.01);

            H_d = eta_d * (r_pre * r_post')';
            H_p = eta_p * (r_pre * r_post')';

            del_T_p = (-T_p_decision_to_timer + ...
                H_p .* (T_max_p - T_p_decision_to_timer)) * ...
                (dt / tau_p);

            del_T_d = (-T_d_decision_to_timer + ...
                H_d .* (T_max_d - T_d_decision_to_timer)) * ...
                (dt / tau_d);

            T_p_decision_to_timer = T_p_decision_to_timer + del_T_p;
            T_d_decision_to_timer = T_d_decision_to_timer + del_T_d;

            T_p_decision_to_timer(T_p_decision_to_timer > T_max_p) = T_max_p;
            T_d_decision_to_timer(T_d_decision_to_timer > T_max_d) = T_max_d;

            % Weight update only during the learning phase
            if rew_vect(t) == 1 && learning_on
                del_W_decision_to_timer = eta_ff * ...
                    (T_p_decision_to_timer - T_d_decision_to_timer) * ...
                    (2 * dt / delay_time);

                W_ff = W_ff + W_ff_id .* del_W_decision_to_timer;
                W_ff = W_ff .* (W_ff > 0) .* (1 - isnan(W_ff));
            end

            mean_T_p_timer(t) = mean(T_p_decision_to_timer(:));
            mean_T_d_timer(t) = mean(T_d_decision_to_timer(:));
        end
    end

    %% Store trial-level results -----------------------------------------
    decision_threshold_cross_time_all(l) = decision_threshold_cross_time;
    decision_input_reduce_time_all(l) = decision_input_reduce_time;

    mean_T_p_timer_all_trials(l, :) = mean_T_p_timer;
    mean_T_d_timer_all_trials(l, :) = mean_T_d_timer;

    mean_activity_timer = mean(1000 * r_timer(:, 1:t_total), 1);

    % Threshold crossing time (not the delayed reduction time)
    time_activity_cross(l) = decision_threshold_cross_time;

    %% Find trace crossing time ------------------------------------------
    difference_traces = mean_T_p_timer - mean_T_d_timer;
    
    crossing_indices_trace = find( ...
        difference_traces(2:end) .* difference_traces(1:end-1) < 0);

    if ~isempty(crossing_indices_trace)
        crossing_index_trace = crossing_indices_trace(1);
        time_trace_cross(l) = crossing_index_trace * dt;
    else
        time_trace_cross(l) = NaN;
    end

    %% Decision making ----------------------------------------------------
    decision_window = (t_steps - 500):t_steps;

    avg_activity_group1 = mean(mean(r_ex(1:npp, decision_window), 2));
    avg_activity_group2 = mean(mean(r_ex(npp+1:2*npp, decision_window), 2));

    if avg_activity_group1 > avg_activity_group2
        network_choice = 1;
    else
        network_choice = 2;
    end

    H = zeros(n_exc_total, n_exc_total);

    %% Store weights and activity ----------------------------------------
    meanwin(l) = mean(W_ff, "all");

    sc_R_it(:, :, l) = is_spike_timer;
    plot_R_it(:, :, l) = mean(r_timer, 1);

    mean_activity_inh_group1(l, :) = mean(r_inh(1:npp, :), 1);
    mean_activity_inh_group2(l, :) = mean(r_inh(npp+1:end, :), 1);

    %% Store first and last trial examples -------------------------------
    if l == first_choice1_trial
        first_trial_choice1_input = s_input_decision;
        first_trial_choice1_decision = r_ex;
        first_trial_choice1_r_timer = r_timer;
    end

    if l == first_choice2_trial
        first_trial_choice2_input = s_input_decision;
        first_trial_choice2_decision = r_ex;
        first_trial_choice2_r_timer = r_timer;
    end

    if l == last_choice1_trial
        last_trial_choice1_input = s_input_decision;
        last_trial_choice1_decision = r_ex;
        last_trial_choice1_r_timer = r_timer;
    end

    if l == last_choice2_trial
        last_trial_choice2_input = s_input_decision;
        last_trial_choice2_r_ex = r_ex;
        last_trial_choice2_r_timer = r_timer;
    end

    %% Online plotting during learning -----------------------------------
    figure(1);
    hold on;

    plot_every = 10;

    if rem(l,plot_every) == 0 || rem(l,plot_every) == 1 || ...
            l == num_trials || l == num_trials - 1

        mean_activity_timer_plot = mean(1000 * r_timer(:, 1:t_total), 1);

        % Clamp the colour index to the number of preset colours
        ci1 = min(round(l/plot_every)+1, size(color_lines_choice1,1));
        ci2 = min(round(l/plot_every)+1, size(color_lines_choice2,1));

        if choice == 1
            plot(0:dt:t_total - dt, mean_activity_timer_plot, ...
                'LineWidth', 1, ...
                'Color', color_lines_choice1(ci1,:));
        else
            plot(0:dt:t_total - dt, mean_activity_timer_plot, ...
                'LineWidth', 1, ...
                'Color', color_lines_choice2(ci2,:));
        end
    end

    xlabel('Time (ms)');
    ylabel('Mean Activity (Hz)');
    title(sprintf('Mean Activity of Integrator - Trial %d', l));

    end

%% Identify the before/after test-phase trials per choice type ------------
before_trials = (1:num_test_before)';
after_trials  = ((num_trials - num_test_after + 1):num_trials)';

before_choice1 = before_trials(trial_type(before_trials) == 1);
before_choice2 = before_trials(trial_type(before_trials) == 2);
after_choice1  = after_trials(trial_type(after_trials) == 1);
after_choice2  = after_trials(trial_type(after_trials) == 2);

time_axis = 0:dt:t_total - dt;

%% Plot mean LTP/LTD eligibility traces for the test phases ----------------
% Averaged over the no-learning trials of each choice type, before vs. after
% the learning phase.
figure(1);
clf;

window_size_samples = 100;

% Choice 1 eligibility traces
Tp_before_c1 = mean(mean_T_p_timer_all_trials(before_choice1, :), 1);
Td_before_c1 = mean(mean_T_d_timer_all_trials(before_choice1, :), 1);
Tp_after_c1  = mean(mean_T_p_timer_all_trials(after_choice1, :), 1);
Td_after_c1  = mean(mean_T_d_timer_all_trials(after_choice1, :), 1);

% Choice 2 eligibility traces
Tp_before_c2 = mean(mean_T_p_timer_all_trials(before_choice2, :), 1);
Td_before_c2 = mean(mean_T_d_timer_all_trials(before_choice2, :), 1);
Tp_after_c2  = mean(mean_T_p_timer_all_trials(after_choice2, :), 1);
Td_after_c2  = mean(mean_T_d_timer_all_trials(after_choice2, :), 1);

Tp_before = (Tp_before_c1+Tp_before_c2)/2;
Td_before = (Td_before_c1+Td_before_c2)/2;
Tp_after = (Tp_after_c1+Tp_after_c2)/2;
Td_after = (Td_after_c1+Td_after_c2)/2;

hold on;
plot(time_axis, Tp_before, 'Color', [line_colors(2,:),0.5], 'LineWidth', 1.5);
plot(time_axis, Td_before, 'Color', [line_colors(4,:),0.5], 'LineWidth', 1.5);
plot(time_axis, Tp_after,  'Color', line_colors(2,:),       'LineWidth', 1.5);
plot(time_axis, Td_after,  'Color', line_colors(4,:),       'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Mean eligibility trace');
title('Mean eligibility trace)');
xlim([0 2500]);
legend({'LTP before','LTD before','LTP after','LTD after'}, 'Location','best');
box off;



%% Plot mean timer ramp for the test phases -------------------------------
% Average ramp across the no-learning trials of each choice type, before vs.
% after the learning phase.
figure(2);
clf;
hold on;

ramp_before_c1 = 1000 * squeeze(mean(plot_R_it(1, :, before_choice1), 3));
ramp_before_c2 = 1000 * squeeze(mean(plot_R_it(1, :, before_choice2), 3));
ramp_after_c1  = 1000 * squeeze(mean(plot_R_it(1, :, after_choice1), 3));
ramp_after_c2  = 1000 * squeeze(mean(plot_R_it(1, :, after_choice2), 3));

plot(time_axis, smooth(ramp_before_c1, window_size_samples), ...
    'Color', [line_colors(1,:),0.5], 'LineWidth', 3);
plot(time_axis, smooth(ramp_before_c2, window_size_samples), ...
    'Color', [line_colors(3,:),0.5], 'LineWidth', 3);
plot(time_axis, smooth(ramp_after_c1, window_size_samples), ...
    'Color', line_colors(1,:), 'LineWidth', 3);
plot(time_axis, smooth(ramp_after_c2, window_size_samples), ...
    'Color', line_colors(3,:), 'LineWidth', 3);

xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
xlim([0 2500]);
legend({'Before Choice 1','Before Choice 2', ...
        'After Choice 1','After Choice 2'}, 'Location','best');
box off;

%% Plot mean feedforward weight -------------------------------------------
figure(3);
clf;

plot(meanwin, ...
    'Color', line_colors(1,:), ...
    'LineWidth', 1.5);

xlabel('Trial Number');
ylabel('Mean W_{ff}');
title('Evolution of Mean Weight');
box off;

%% Plot example activity in each layer ------------------------------------
figure(4);
clf;

% Input layer
subplot(3,1,1);
maxY = -inf;
hold on;

for neuron = 1:3
    smoothed_input_choice1_last = smoothdata( ...
        last_trial_choice1_input(neuron, :) * 1000, ...
        'movmean', window_size_samples);

    plot(0:dt:t_total - dt, smoothed_input_choice1_last, '-', ...
        'Color', line_colors(1,:), ...
        'LineWidth', 2);

    maxY = max(maxY, max(smoothed_input_choice1_last));
end

for neuron = 1:3
    smoothed_input_choice2_last = smoothdata( ...
        last_trial_choice2_input(neuron, :) * 1000, ...
        'movmean', window_size_samples);

    plot(0:dt:t_total - dt, smoothed_input_choice2_last, '-', ...
        'Color', line_colors(3,:), ...
        'LineWidth', 2);

    maxY = max(maxY, max(smoothed_input_choice2_last));
end

xline(500, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 2);
ylabel('Input rate (Hz)');
xlim([0 1500]);

if maxY > 0
    ylim([0, maxY * 1.2]);
end

box off;

% Decision layer
subplot(3,1,2);
maxY = -inf;
hold on;

for neuron = 1:3
    smoothed_decision_choice1_last = smoothdata( ...
        last_trial_choice1_decision(neuron, :) * 1000, ...
        'movmean', window_size_samples);

    plot(0:dt:t_total - dt, smoothed_decision_choice1_last, '-', ...
        'Color', line_colors(1,:), ...
        'LineWidth', 2);

    maxY = max(maxY, max(smoothed_decision_choice1_last));
end

for neuron = 1:3
    smoothed_decision_choice2_last = smoothdata( ...
        last_trial_choice2_r_ex(neuron, :) * 1000, ...
        'movmean', window_size_samples);

    plot(0:dt:t_total - dt, smoothed_decision_choice2_last, '-', ...
        'Color', line_colors(3,:), ...
        'LineWidth', 2);

    maxY = max(maxY, max(smoothed_decision_choice2_last));
end

xline(500, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 2);
ylabel('Decision rate (Hz)');
xlim([0 1500]);

if maxY > 0
    ylim([0, maxY * 1.2]);
end

box off;

% Timer layer
subplot(3,1,3);
maxY = -inf;
hold on;

for neuron = 1:3
    smoothed_timer_choice1_last = smoothdata( ...
        last_trial_choice1_r_timer(neuron, :) * 1000, ...
        'movmean', window_size_samples);

    plot(0:dt:t_total - dt, smoothed_timer_choice1_last, '-', ...
        'Color', line_colors(1,:), ...
        'LineWidth', 2);

    maxY = max(maxY, max(smoothed_timer_choice1_last));
end

for neuron = 1:3
    smoothed_timer_choice2_last = smoothdata( ...
        last_trial_choice2_r_timer(neuron, :) * 1000, ...
        'movmean', window_size_samples);

    plot(0:dt:t_total - dt, smoothed_timer_choice2_last, '-', ...
        'Color', line_colors(3,:), ...
        'LineWidth', 2);

    maxY = max(maxY, max(smoothed_timer_choice2_last));
end

xline(500, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Integrator rate (Hz)');
xlim([0 1500]);

if maxY > 0
    ylim([0, maxY * 1.2]);
end

box off;

%% Plot trace crossing time vs activity threshold time --------------------
figure(5);
clf;

plot(time_trace_cross(1:2:end), time_activity_cross(1:2:end), 'o', ...
    'Color', line_colors(1,:), ...
    'MarkerFaceColor', line_colors(1,:));

hold on;

plot(time_trace_cross(2:2:end), time_activity_cross(2:2:end), 'o', ...
    'Color', line_colors(3,:), ...
    'MarkerFaceColor', line_colors(3,:));

xlabel('Trace crossing time (ms)');
ylabel('Threshold crossing time (ms)');
title('Trace Crossing Time vs Timer Threshold Time');
box off;



%% Plot binned timer mean activity every 10 trials ------------------------
figure(7);
clf;
hold on;

bin_size = 10;
num_bins = ceil(num_trials / bin_size);

% Gradually darker red
light_red = [1.00, 0.75, 0.75];
dark_red  = [0.55, 0.00, 0.00];

red_colors = zeros(num_bins, 3);

for b = 1:num_bins
    if num_bins == 1
        red_colors(b,:) = dark_red;
    else
        alpha = (b - 1) / (num_bins - 1);
        red_colors(b,:) = light_red * (1 - alpha) + dark_red * alpha;
    end
end

for b = 1:num_bins
    trial_start = (b - 1) * bin_size + 1;
    trial_end   = min(b * bin_size, num_trials);
    trial_idx   = trial_start:trial_end;

    % plot_R_it stores mean(r_timer,1), so multiply by 1000 to convert to Hz.
    binned_mean_activity = 1000 * squeeze(mean(plot_R_it(1,:,trial_idx), 3));

    plot(time_axis, binned_mean_activity, ...
        'Color', red_colors(b,:), ...
        'LineWidth', 2.5, ...
        'DisplayName', sprintf('Trials %d-%d', trial_start, trial_end));
end

yline(threshould, '--', ...
    '', ...
    'LineWidth', 1.5, ...
    'Color', [0.3 0.3 0.3]);

xlabel('Time (ms)');
ylabel('Mean activity (Hz)');
title('Binned Mean Activity');
xlim([0 2500]);
legend('Location', 'best');
box off;