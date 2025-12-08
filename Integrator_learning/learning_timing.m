%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:  Simulate a two-module spiking network (Decision ↔ Timer).
%           The script runs multiple learning trials and visualizes activity 
%           and trace dynamics across trials.
%
% Summary:
%   • Decision network: two excitatory pools (size = npp each) with
%     mutual inhibition via two inhibitory pools; receives brief Poisson
%     drive selecting Choice 1 or Choice 2 each trial.
%   • Timer(integrator) network: recurrent excitatory population (N_timer) driven by
%     decision E neurons; forms a ramp/accumulator
%   • Learning: pre/post rate products form LTP/LTD eligibility traces
%     (T_p, T_d), learning happended at the time of reward (t_stim = 1500;)


%
% Dependencies:
%   • `slanCL.m` (or `slanCL/` directory) for color palettes

%
% Author:    Bolu Feng
% Date:      2025-11-06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Clear workspace and close all figures
clear all
close all


% addpath slanCL
% line_colors = slanCL(617);
% color_lines_choice1 = slanCL(277);
% color_lines_choice2 = slanCL(526);

% if don't use slanCL color palettes
nNeeded = 10;
% Default MATLAB palette (classic "lines" scheme)
line_colors = lines(nNeeded);  % replaces: slanCL(617)
color_lines_choice1 = repmat(line_colors(1,:), nNeeded, 1);  % replaces: slanCL(277)
color_lines_choice2 = repmat(line_colors(2,:), nNeeded, 1);  % replaces: slanCL(526)


%% Structure Parameters
% Number of neurons in populations
N_timer = 100; % Number of neurons in timer network (integrator network)
npp = 50;     % Neurons per population in decision network
n_exc_total = 2 * npp;    % Total number of excitatory neurons in decision network
n_inh_total = 2 * npp;    % Total number of inhibitory neurons in decision network
N_input = 2*npp; % Number of input neurons

num_trials = 200;             % Number of learning trials

%% Time Parameters
dt = 1;                      % Time step size (ms)
T = 3000 / dt;               % Duration of stimulus presentation (time steps)
delta = 600 / dt;            % Time between stimuli (time steps)
t_total = 3001;              % Total time of a trial (ms)
t_steps = t_total / dt;      % Total number of time steps
D = 10 / dt;                 % Intrinsic delay (time steps)
tau_w = 50;                  % Time window for firing rate dynamics (ms)

%% Input Parameters
eG = 0.5;                    % Input gain for excitatory neurons
iG = 0.5;                    % Input gain for inhibitory neurons

% Input layer to decision neurons
W_in = cell(2,1);
W_in{1} = 0.002*rand(n_exc_total,1);
W_in{2} = 0.002*rand(n_exc_total,1);
% Stimulus 1 excites group 1
W_in{1}(1:npp) = eG; % Excitatory input to group 1
% Stimulus 2 excites group 2
W_in{2}(npp+1:2*npp) = eG; % Excitatory input to group 2

% Input layer to timer network remains unchanged
p_r = 0.02 * dt;             % Poisson rate for input spikes (40 Hz)

%% Weight Matrix Parameters
% Decision network connectivity
rho = 0.5;                % Sparsity of connections
w_rec = 0.002;         % Recurrent excitatory weight
i_mean = 0.0015;         % Inhibitory weight from inhibitory to excitatory neurons
p_mean = 0.0015;         % Excitatory to inhibitory weight
k_rec = 0.00;    % Recurrent inhibitory weight (inhibitory to inhibitory)

% Initialize Weight Matrices for Decision Network
W_ji = zeros(n_exc_total, n_exc_total);       % Excitatory to excitatory weight matrix
M_ki = zeros(n_exc_total, n_inh_total);       % Inhibitory to excitatory weight matrix
P_ik = zeros(n_inh_total, n_exc_total);       % Excitatory to inhibitory weight matrix
K_jk = zeros(n_inh_total, n_inh_total);       % Inhibitory to inhibitory weight matrix

% Randomly initialize weights with sparsity rho
% Excitatory to excitatory
W_ji(1:npp, 1:npp) = w_rec * (rand(npp) < rho);
W_ji(npp+1:2*npp, npp+1:2*npp) = w_rec * (rand(npp) < rho);

% Inhibitory to excitatory (mutual inhibition between groups)
M_ki(1:npp, npp+1:2*npp) = i_mean * (rand(npp) < rho);      % Inh group 2 to Exc group 1
M_ki(npp+1:2*npp, 1:npp) = i_mean * (rand(npp) < rho);      % Inh group 1 to Exc group 2

% Excitatory to inhibitory
P_ik(1:npp, 1:npp) = p_mean * (rand(npp) < rho);            % Exc group 1 to Inh group 1
P_ik(npp+1:2*npp, npp+1:2*npp) = p_mean * (rand(npp) < rho);% Exc group 2 to Inh group 2

% Inhibitory to inhibitory (recurrent inhibitory connections)
K_jk(1:npp, 1:npp) = k_rec * (rand(npp) < rho);        % Inh group 1 recurrent connections
K_jk(npp+1:2*npp, npp+1:2*npp) = k_rec * (rand(npp) < rho); % Inh group 2 recurrent connections

% Identify recurrent synapses
rec_identity_decision = zeros(n_exc_total, n_exc_total);
rec_identity_decision(1:npp, 1:npp) = W_ji(1:npp, 1:npp) > 0;
rec_identity_decision(npp+1:2*npp, npp+1:2*npp) = W_ji(npp+1:2*npp, npp+1:2*npp) > 0;

% Decision neurons to timer network
W_decision_to_timer = 0.00014 *rand(N_timer, n_exc_total).* (rand(N_timer, n_exc_total) < rho);
W_decision_to_timer_id = W_decision_to_timer > 0;

% Recurrent connections in timer network
w_timer_mean = 0.000135;            % Recurrent connections in timer network
W_timer = w_timer_mean * rand(N_timer, N_timer) .* (rand(N_timer, N_timer) < rho);
rec_identity_timer = W_timer > 0;

w_input_to_decision_mean = 0.001;  
W_input_to_decision = zeros(n_exc_total, N_input);
W_input_to_decision(1:npp, 1:N_input/2) = w_input_to_decision_mean*rand(npp, N_input/2);
W_input_to_decision(npp + 1: n_exc_total, N_input/2 + 1: N_input) = w_input_to_decision_mean*rand(npp, N_input/2);
% %% Membrane Dynamics Parameters
rho_syn = 1 / 7;            % Fractional change of synaptic activation
tau_se_decision_input = 10; % Time constant for excitatory synaptic activation (ms) for decision and input net
tau_se_timer = 80;          % Time constant for excitatory synaptic activation (ms) for timer net
tau_si = 10;                % Time constant for inhibitory synaptic activation (ms)
norm_noise = 0.1;           % Normalization factor for noise
C_m = 0.2;                  % Membrane capacitance (µF/cm²)
g_L = 0.01;                 % Leak conductance (mS/cm²)
E_i = -70;                  % Inhibitory reversal potential (mV)
E_l = -60;                  % Leak reversal potential (mV)
E_e = -5;                   % Excitatory reversal potential (mV)
v_th = -55;                 % Threshold potential for excitatory neurons (mV)
v_th_i = -50;               % Threshold potential for inhibitory neurons (mV)
v_rest = -60;               % Resting potential (mV)
v_hold = -61;               % Reset potential after spike (mV)
t_refractory = 2 / dt;      % Refractory period duration (time steps)

%% Learning Parameters
tau_p = 1500;                
    tau_d = 1000;                
    T_max_p = 0.032;             
    T_max_d = 0.038;             
    eta_p = 1.5e4;%30 * 350*1.3*30*5;        
    eta_d = 1.5e3;%3 * 350*1.3*30/20;         
    eta_ff = 0.00005;      
% eta_ff = 0.0002;

% eta_ff = 0.00006;             % Learning rate for recurrent connections

%% Reward Parameters
t_stim = 1500;           % Times when reward are presented (time steps)
delay_time = 25 / dt;        % Reward delay (time steps)
t_reward = t_stim(1:end) + delay_time; % Times when reward is given
rew_vect = zeros(t_total, 1); % Reward vector initialization
del_time = round(delay_time / 2); % Half of reward delay

% Set up reward vector with reward window
for m = t_reward
    rew_vect(m - del_time : m + del_time) = 1; % Reward window
    rew_vect(m + del_time + 1) = 2;            % End of reward window
end

% Initialize variables for recording data
sc_R_it = zeros(N_timer, t_steps, num_trials + 1); % Spike counts
rec_vect = zeros(1, num_trials + 1); % Mean recurrent weights over trials

%% Trial Types
trial_type = zeros(num_trials, 1);
trial_type(1:2:num_trials) = 1;  % Stimulus type 1
trial_type(2:2:num_trials) = 2;  % Stimulus type 2

% Find indices of trials for choice 1 and choice 2
choice1_trials = find(trial_type == 1);
choice2_trials = find(trial_type == 2);

% First trial of choice 1 and choice 2
first_choice1_trial = choice1_trials(1);
first_choice2_trial = choice2_trials(1);

% Last trial of choice 1 and choice 2
last_choice1_trial = choice1_trials(end);
last_choice2_trial = choice2_trials(end);

%% Variables to Store Inhibitory Neuron Activity
mean_activity_inh_group1 = zeros(num_trials, t_steps);
mean_activity_inh_group2 = zeros(num_trials, t_steps);

%% Initialize storage for mean eligibility traces before the trial loop
mean_T_p_timer_all_trials = zeros(num_trials, t_steps);
mean_T_d_timer_all_trials = zeros(num_trials, t_steps);

% Initialize variables to store neuron dynamics for first and last trials
first_trial_choice1_input = [];
first_trial_choice1_decision = [];
first_trial_choice1_r_timer = [];

first_trial_choice2_input = [];
first_trial_choice2_decision = [];
first_trial_choice2_r_timer = [];

last_trial_choice1_input = [];
last_trial_choice1_decision = [];
last_trial_choice1_r_timer = [];

last_trial_choice2_input = [];
last_trial_choice2_r_ex = [];
last_trial_choice2_r_timer = [];

%% Initialize arrays to store crossing times
time_activity_cross = zeros(num_trials,1);
time_trace_cross = zeros(num_trials,1);

%% Main Program
for l = 1:num_trials
    % Determine current stimulus
    choice = trial_type(l);
    
    %% Initialization
    trace_refractory = 0; % Initialize trace refractory period
    
    % Input layer
    s_input = zeros(N_input, t_steps); % Synaptic activations for input
    input_spike = poissrnd(p_r, N_input, t_steps); % Poisson input spike trains
    
    % Decision network neurons
    % Excitatory neurons
    v_ex = v_rest * ones(n_exc_total, t_steps);
    r_ex = zeros(n_exc_total, t_steps);
    s_ex = zeros(n_exc_total, t_steps);
    r_input_decision = zeros(n_exc_total, t_steps);
    is_refractory_ex = zeros(n_exc_total, t_steps);
    is_spike_ex = zeros(n_exc_total, t_steps);
    s_input_decision = zeros(n_exc_total, t_steps);
    
    % Inhibitory neurons
    v_inh = v_rest * ones(n_inh_total, t_steps);
    r_inh = zeros(n_inh_total, t_steps);
    s_inh = zeros(n_inh_total, t_steps);
    is_refractory_inh = zeros(n_inh_total, t_steps);
    is_spike_inh = zeros(n_inh_total, t_steps);
    
    % Timer network neurons
    v_timer = zeros(N_timer, t_steps) + v_rest; % Membrane potentials for timer neurons
    r_timer = zeros(N_timer, t_steps); % Firing rates for timer neurons
    s_timer = zeros(N_timer, t_steps); % Synaptic activations for timer neurons
    g_decision_to_timer = zeros(N_timer, t_steps); % Conductances from decision neurons to timer neurons
    g_timer_to_timer = zeros(N_timer, t_steps); % Recurrent excitatory conductance in timer network
    is_refractory_timer = zeros(N_timer, t_steps); % Refractory periods for timer neurons
    is_spike_timer = zeros(N_timer, t_steps); % Spikes in timer neurons
    
    % Variables to accumulate Hebbian terms
    T_p_decision_to_timer = zeros(N_timer, n_exc_total); % Eligibility traces for LTP (decision to timer)
    T_d_decision_to_timer = zeros(N_timer, n_exc_total); % Eligibility traces for LTD (decision to timer)
    % T_p_input_to_decision = zeros( n_exc_total, N_input); % Eligibility traces for LTP (decision to timer)
    % T_d_input_to_decision = zeros( n_exc_total, N_input); % Eligibility traces for LTD (decision to timer)
    
    % Initialize mean eligibility traces before the time loop
    mean_T_p_timer = zeros(1, t_steps);
    mean_T_d_timer = zeros(1, t_steps);
    
    %% Time Step Loop
    for t = 2:t_steps
        %% Input Layer Dynamics
        % Generate Poisson input spike trains for decision network
        input_spike_decision = poissrnd(p_r, n_exc_total, 1);
        if t > 1 && t < 500
            is_spike_input_decision = (input_spike_decision == 1);
        else
            is_spike_input_decision = 0;
        end
        s_input_decision(:, t) = s_input_decision(:, t-1) - (s_input_decision(:, t-1) * dt / tau_se_decision_input) + ...
                                rho_syn * W_in{choice} .* is_spike_input_decision .* (1 - s_input_decision(:, t-1));
        r_input_decision(:, t) = r_input_decision(:, t-1) + (W_in{choice}.*is_spike_input_decision / dt - r_input_decision(:, t-1)) * (dt / tau_w);
        
        %% Decision Network Dynamics
        % Excitatory Neuron Dynamics
        is_refractory = is_refractory_ex(:, t) == 1;
        v_ex(is_refractory, t) = v_rest;

        % Check for spikes
        is_spike = (v_ex(:, t-1) >= v_th);
        is_spike_ex(:, t) = is_spike;

        % Update refractory status
        not_spike = (~is_spike & ~is_refractory);
        v_ex(is_spike, t) = v_hold;
        if t < t_steps - t_refractory
            is_refractory_ex(is_spike, t+1:t+t_refractory) = 1;
        end

        % Update synaptic activation
        s_ex(:, t) = s_ex(:, t-1) - (s_ex(:, t-1) * dt / tau_se_decision_input) + ...
                     rho_syn * is_spike .* (1 - s_ex(:, t-1));

        % Update membrane potential
        v_ex(not_spike, t) = v_ex(not_spike, t-1) + ...
            ((randn(sum(not_spike),1) * norm_noise) + ...
            g_L * (E_l - v_ex(not_spike, t-1)) + ...
            (W_ji(not_spike, :) * s_ex(:, t-1)) .* (E_e - v_ex(not_spike, t-1)) + ...
            (M_ki(not_spike, :) * s_inh(:, t-1)) .* (E_i - v_ex(not_spike, t-1)) + ...
            W_input_to_decision(not_spike, :) * s_input_decision(:, t-1) .* (E_e - v_ex(not_spike, t-1))) * (dt / C_m);

        % Update firing rates
        r_ex(:, t) = r_ex(:, t-1) + (is_spike / dt - r_ex(:, t-1)) * (dt / tau_w);

        %% Inhibitory Neuron Dynamics
        is_refractory = is_refractory_inh(:, t) == 1;
        v_inh(is_refractory, t) = v_rest;

        % Check for spikes
        is_spike = (v_inh(:, t-1) >= v_th_i);
        is_spike_inh(:, t) = is_spike;

        % Update refractory status
        not_spike = (~is_spike & ~is_refractory);
        v_inh(is_spike, t) = v_hold;
        if t < t_steps - t_refractory
            is_refractory_inh(is_spike, t+1:t+t_refractory) = 1;
        end

        % Update synaptic activation
        s_inh(:, t) = s_inh(:, t-1) - (s_inh(:, t-1) * dt / tau_si) + ...
                      rho_syn * is_spike .* (1 - s_inh(:, t-1));

        % Update membrane potential (including recurrent inhibitory connections)
        v_inh(not_spike, t) = v_inh(not_spike, t-1) + ...
            ((randn(sum(not_spike),1) * norm_noise) + ...
            g_L * (E_l - v_inh(not_spike, t-1)) + ...
            (P_ik(not_spike, :) * s_ex(:, t-1)) .* (E_e - v_inh(not_spike, t-1)) + ...
            (K_jk(not_spike, :) * s_inh(:, t-1)) .* (E_i - v_inh(not_spike, t-1))) * (dt / C_m);

        % Update firing rates
        r_inh(:, t) = r_inh(:, t-1) + (is_spike / dt - r_inh(:, t-1)) * (dt / tau_w);

        %% Timer Network Dynamics
        % Compute input conductance from decision neurons
        g_decision_to_timer(:, t) = W_decision_to_timer * s_ex(:, t - 1);
        % Compute recurrent conductance within timer network
        g_timer_to_timer(:, t) = W_timer * s_timer(:, t - 1);
        % Total conductance to timer neurons
        g_total_timer = g_decision_to_timer(:, t) + g_timer_to_timer(:, t);

        % Update membrane potentials for timer neurons
        is_refractory = is_refractory_timer(:, t) == 1;
        v_timer(is_refractory, t) = v_rest; % Reset membrane potential for refractory neurons
        is_spike = (v_timer(:, t - 1) >= v_th); % Check if neuron spikes
        is_spike_timer(:, t) = is_spike;
        not_spike = (~is_spike & ~is_refractory); % Non-spiking and non-refractory neurons
        v_timer(is_spike, t) = v_hold; % Set membrane potential to reset after spike
        if t < t_steps - t_refractory
            is_refractory_timer(is_spike, t + 1:t + t_refractory) = 1;
        end
        % Update synaptic activations for timer neurons
        s_timer(:, t) = s_timer(:, t - 1) - (s_timer(:, t - 1) * dt / tau_se_timer) + rho_syn * is_spike .* (1 - s_timer(:, t - 1));
        % Update membrane potentials for non-spiking neurons
        v_timer(not_spike, t) = v_timer(not_spike, t - 1) + ...
            ((randn(sum(not_spike), 1) * norm_noise) + ...
            g_L * (E_l - v_timer(not_spike, t - 1)) + ... % Noise and decay term
            g_total_timer(not_spike) .* ...
            (E_e - v_timer(not_spike, t - 1))) * (dt / C_m); % Input and recurrent terms
        % Update firing rates for timer neurons
        r_timer(:, t) = r_timer(:, t - 1) + (is_spike / dt - r_timer(:, t - 1)) * (dt / tau_w);

         %% Learning
        if t > D
            %% Learning for Decision to Timer Neurons Connections
            r_pre = r_ex(:, t - D) .* (r_ex(:, t - D) > 0.01);
            r_post = r_timer(:, t - 1) .* (r_timer(:, t - 1) > 0.01);
            H_d = eta_d * (r_pre * r_post')'; % LTD term
            H_p = eta_p * (r_pre * r_post')'; % LTP term
            % Update eligibility traces
            del_T_p = (-T_p_decision_to_timer + H_p .* (T_max_p - T_p_decision_to_timer)) * (dt / tau_p); % LTP trace
            del_T_d = (-T_d_decision_to_timer + H_d .* (T_max_d - T_d_decision_to_timer)) * (dt / tau_d); % LTD trace
            T_p_decision_to_timer = T_p_decision_to_timer + del_T_p;
            T_d_decision_to_timer = T_d_decision_to_timer + del_T_d;
            T_p_decision_to_timer(T_p_decision_to_timer > T_max_p) = T_max_p;
            T_d_decision_to_timer(T_d_decision_to_timer > T_max_d) = T_max_d;
            % Update synaptic weights during the reward window
            % if rew_vect(t) == 1 && l > 2 && l < num_trials -1 && choice == 1% Reward window
            %     del_W_decision_to_timer = eta_ff * (T_p_decision_to_timer - T_d_decision_to_timer) * (2 * dt / delay_time); % Weight update
            %     W_decision_to_timer = W_decision_to_timer + W_decision_to_timer_id .* del_W_decision_to_timer;
            %     W_decision_to_timer = W_decision_to_timer .* (W_decision_to_timer > 0) .* (1 - isnan(W_decision_to_timer)); % Ensure weights remain positive and finite
            % end

            % if t > 800 && t < 830 && l > 2 && l < num_trials -1 && choice == 1% Reward window
            %     del_W_decision_to_timer = eta_ff * (T_p_decision_to_timer - T_d_decision_to_timer) * (2 * dt / delay_time); % Weight update
            %     W_decision_to_timer = W_decision_to_timer + W_decision_to_timer_id .* del_W_decision_to_timer;
            %     W_decision_to_timer = W_decision_to_timer .* (W_decision_to_timer > 0) .* (1 - isnan(W_decision_to_timer)); % Ensure weights remain positive and finite
            % end
            % if t > 1200 && t < 1230 && l > 2 && l < num_trials -1 && choice == 2% Reward window
            %     del_W_decision_to_timer = eta_ff * (T_p_decision_to_timer - T_d_decision_to_timer) * (2 * dt / delay_time); % Weight update
            %     W_decision_to_timer = W_decision_to_timer + W_decision_to_timer_id .* del_W_decision_to_timer;
            %     W_decision_to_timer = W_decision_to_timer .* (W_decision_to_timer > 0) .* (1 - isnan(W_decision_to_timer)); % Ensure weights remain positive and finite
            % end
            if t > t_stim && t < t_stim + 30 && l > 2 && l < num_trials -1 %% reward window
                del_W_decision_to_timer = eta_ff * (T_p_decision_to_timer - T_d_decision_to_timer) * (2 * dt / delay_time); % Weight update
                W_decision_to_timer = W_decision_to_timer + W_decision_to_timer_id .* del_W_decision_to_timer;
                W_decision_to_timer = W_decision_to_timer .* (W_decision_to_timer > 0) .* (1 - isnan(W_decision_to_timer)); % Ensure weights remain positive and finite
            end

   
            % Record Mean Eligibility Traces
            mean_T_p_timer(t) = mean(T_p_decision_to_timer(:));
            mean_T_d_timer(t) = mean(T_d_decision_to_timer(:));
        end
    end % End of time step loop

    %% Store mean eligibility traces after the trial
    mean_T_p_timer_all_trials(l, :) = mean_T_p_timer;
    mean_T_d_timer_all_trials(l, :) = mean_T_d_timer;

    %% Compute mean activity of the timer neurons
    mean_activity_timer = mean(1000 * r_timer(:, 1:t_total), 1); % Mean firing rate of all timer neurons

    %% Find the first time when mean_activity_timer crosses 20 Hz
    crossing_index_activity = find(mean_activity_timer > 28, 1, 'first');
    if ~isempty(crossing_index_activity)
        time_activity_cross(l) = crossing_index_activity * dt; % Time in ms
    else
        time_activity_cross(l) = NaN;
    end

    %% Find the time when mean_T_p_timer and mean_T_d_timer cross
    difference_traces = mean_T_p_timer - mean_T_d_timer;
    crossing_indices_trace = find(difference_traces(2:end) .* difference_traces(1:end-1) <= 0);
    if ~isempty(crossing_indices_trace)
        % Take the first crossing
        crossing_index_trace = crossing_indices_trace(end);
        time_trace_cross(l) = crossing_index_trace * dt; % Time in ms
    else
        time_trace_cross(l) = NaN;
    end

    %% Decision Making and Reward/Punishment
    % Compute average firing rates over decision window (e.g., last 500 ms)
    decision_window = (t_steps - 500):t_steps;
    avg_activity_group1 = mean(mean(r_ex(1:npp, decision_window), 2));
    avg_activity_group2 = mean(mean(r_ex(npp+1:2*npp, decision_window), 2));

    % Determine network's decision based on higher activity
    if avg_activity_group1 > avg_activity_group2
        network_choice = 1;
    else
        network_choice = 2;
    end

    % Reset Hebbian terms
    H = zeros(n_exc_total, n_exc_total);

    %% Update Weights and Plotting
    meanwin(l) = mean(W_decision_to_timer, "all");
    sc_R_it(:, :, l) = is_spike_timer; % Record spike counts for plotting
    plot_R_it(:, :, l) = mean(r_timer, 1); % Compute mean firing rates for each population

    %% Store Mean Inhibitory Activity for Each Group
    mean_activity_inh_group1(l, :) = mean(r_inh(1:npp, :), 1);
    mean_activity_inh_group2(l, :) = mean(r_inh(npp+1:end, :), 1);

    %% Store neuron dynamics for specified trials
    if l == first_choice1_trial
        % Store data for first trial of choice 1
        first_trial_choice1_input = r_input_decision;
        first_trial_choice1_decision = r_ex;
        first_trial_choice1_r_timer = r_timer;
    end
    if l == first_choice2_trial
        % Store data for first trial of choice 2
        first_trial_choice2_input = r_input_decision;
        first_trial_choice2_decision = r_ex;
        first_trial_choice2_r_timer = r_timer;
    end
    if l == last_choice1_trial
        % Store data for last trial of choice 1
        last_trial_choice1_input = r_input_decision;
        last_trial_choice1_decision = r_ex;
        last_trial_choice1_r_timer = r_timer;
    end
    if l == last_choice2_trial
        % Store data for last trial of choice 2
        last_trial_choice2_input = r_input_decision;
        last_trial_choice2_r_ex = r_ex;
        last_trial_choice2_r_timer = r_timer;
    end

    %% Custom Plotting Section
    figure(1);
    % Plot mean activity of the timer neurons
    % subplot(3, 1, 1);
    hold on;
    plot_every = 25;
    if rem(l,plot_every) == 0 || rem(l,plot_every) == 1 || l == num_trials || l == num_trials -1
    mean_activity_timer_plot = mean(1000 * r_timer(:, 1:t_total), 1); % Mean firing rate of all timer neurons
    if choice == 1
        plot(0:dt:t_total - dt, mean_activity_timer_plot, 'LineWidth', 1, 'Color', color_lines_choice1(round(l/plot_every)+1,:));
    else
        plot(0:dt:t_total - dt, mean_activity_timer_plot, 'LineWidth', 1, 'Color', color_lines_choice2(round(l/plot_every)+1,:));
    end
    end
    xlabel('Time (ms)');
    ylabel('Mean Activity (Hz)');
    title(sprintf('Mean Activity of Timer Neurons - Trial %d', l));

    % Plot mean activity of the decision neurons (groups separately)
    figure(2);
    subplot(2, 1, 1);
    hold on;
    mean_activity_decision_group1 = mean(1000 * r_ex(1:npp, 1:t_total), 1); % Mean firing rate of decision group 1
    mean_activity_decision_group2 = mean(1000 * r_ex(npp+1:2*npp, 1:t_total), 1); % Mean firing rate of decision group 2

    plot(0:dt:t_total - dt, mean_activity_decision_group1, 'LineWidth', 2, 'DisplayName', 'Excitatory Group 1');
    plot(0:dt:t_total - dt, mean_activity_decision_group2, 'LineWidth', 2, 'DisplayName', 'Excitatory Group 2');
    xlabel('Time (ms)');
    ylabel('Mean Activity (Hz)');
    title(sprintf('Mean Activity of Excitatory Neurons - Trial %d', l));
    % legend();

    % Plot mean activity of the inhibitory neurons (groups separately)
    subplot(2, 1, 2);
    hold on;
    for ip = 1:5
        plot(0:dt:t_total - dt, 1000 * r_ex(ip, 1:t_total), 'Color',line_colors(1,:) ,'LineWidth', 2, 'DisplayName', 'Excitatory Group 1');
        plot(0:dt:t_total - dt, 1000 * r_ex(ip+npp, 1:t_total), 'Color',line_colors(3,:)','LineWidth', 2, 'DisplayName', 'Excitatory Group 2');
    end
    % mean_activity_inh1 = mean(1000 * r_inh(1:npp, 1:t_total), 1); % Mean firing rate of inhibitory group 1
    % mean_activity_inh2 = mean(1000 * r_inh(npp+1:2*npp, 1:t_total), 1); % Mean firing rate of inhibitory group 2
    % 
    % plot(0:dt:t_total - dt, mean_activity_inh1, 'LineWidth', 2, 'DisplayName', 'Inhibitory Group 1');
    % plot(0:dt:t_total - dt, mean_activity_inh2, 'LineWidth', 2, 'DisplayName', 'Inhibitory Group 2');
    xlabel('Time (ms)');
    ylabel('Mean Activity (Hz)');
    title(sprintf('Mean Activity of Inhibitory Neurons - Trial %d', l));
    % legend();
end % End of main loop over trials

% Define window size for smoothing (100 ms)
window_size_ms = 100; % 100 milliseconds
window_size_samples = round(window_size_ms / dt); % Convert to number of samples

% Ensure window size is at least 1
window_size_samples = max(window_size_samples, 1);

% % Display the figure
% figure(2);
% plot(meanwin);
% xlabel('Trial Number');
% ylabel('Mean Recurrent Weight in Decision Network');
% title('Evolution of Mean Recurrent Weight in Decision Network over Trials');
%%
% Plot the mean LTP and LTD eligibility traces for each trial
figure(3);clf
hold on;
% Smooth the LTP and LTD traces
smoothed_mean_T_p_timer_all_trials = smoothdata(mean_T_p_timer_all_trials, 2, 'movmean', window_size_samples);
smoothed_mean_T_d_timer_all_trials = smoothdata(mean_T_d_timer_all_trials, 2, 'movmean', window_size_samples);

% Plot Mean LTP and LTD Eligibility Traces
plot(0:dt:t_total - dt, mean_T_p_timer_all_trials(1, :), 'Color', [line_colors(2,:),0.5], 'LineWidth', 1.5);

plot(0:dt:t_total - dt, mean_T_d_timer_all_trials(1, :), 'Color', [line_colors(4,:),0.5],  'LineWidth', 1.5);
plot(0:dt:t_total - dt, mean_T_p_timer_all_trials(100, :), 'Color', line_colors(2,:), 'LineWidth', 1.5);
plot(0:dt:t_total - dt, mean_T_d_timer_all_trials(100, :), 'Color', line_colors(4,:), 'LineWidth', 1.5);

xlabel('Time (ms)');
ylabel('Mean LTP Eligibility Trace');
title('Mean LTP and LTD Eligibility Traces over Time for Each Trial');
xlim([0 2500]); % Set x-axis limits to 0-2s
% legend();
%%
figure(4);clf
hold on;

% 
% Plot Mean LTP and LTD Eligibility Traces
plot(0:dt:t_total - dt, 1000*smooth(plot_R_it(1,:,1),window_size_samples), 'Color', [line_colors(1,:),0.5], 'LineWidth', 3);
hold on;
plot(0:dt:t_total - dt, 1000*smooth(plot_R_it(1,:,2),window_size_samples), 'Color', [line_colors(3,:),0.5],  'LineWidth', 3);
xlabel('Time (ms)');
ylabel('Mean LTP Eligibility Trace');
title('Mean LTP and LTD Eligibility Traces over Time for Each Trial');
xlim([0 2500]); % Set x-axis limits to 0-2s

% Plot Mean LTP and LTD Eligibility Traces
plot(0:dt:t_total - dt, 1000*smooth(plot_R_it(1,:,199),window_size_samples), 'Color', [line_colors(1,:)], 'LineWidth', 3);
hold on;
plot(0:dt:t_total - dt, 1000*smooth(plot_R_it(1,:,200),window_size_samples), 'Color', [line_colors(3,:)],  'LineWidth', 3);
xlabel('Time (ms)');
ylabel('Mean LTP Eligibility Trace');
title('Mean LTP and LTD Eligibility Traces over Time for Each Trial');
xlim([0 2500]); % Set x-axis limits to 0-2s
%%

figure(8);clf
hold on;
% Smooth the LTP and LTD traces
smoothed_mean_T_p_timer_all_trials = smoothdata(mean_T_p_timer_all_trials, 2, 'movmean', window_size_samples);
smoothed_mean_T_d_timer_all_trials = smoothdata(mean_T_d_timer_all_trials, 2, 'movmean', window_size_samples);

% Plot Mean LTP and LTD Eligibility Traces
plot(0:dt:t_total - dt, smoothed_mean_T_p_timer_all_trials(1, :), 'Color', [.7 .7 .7], 'LineWidth', 2);

plot(0:dt:t_total - dt, smoothed_mean_T_d_timer_all_trials(1, :), '--','Color', [.7 .7 .7],  'LineWidth', 2);
plot(0:dt:t_total - dt, smoothed_mean_T_p_timer_all_trials(num_trials -1, :), 'Color', [.2 .2 .2], 'LineWidth', 2);
plot(0:dt:t_total - dt, smoothed_mean_T_d_timer_all_trials(num_trials -1, :), '--','Color', [.2 .2 .2], 'LineWidth', 2);

xlabel('Time (ms)');
ylabel('Mean LTP Eligibility Trace');
title('Mean LTP and LTD Eligibility Traces over Time for Each Trial');
xlim([0 2500]); % Set x-axis limits to 0-2s
% hold on;
%%
% x_shade = [1450, 1550, 1550, 1450];
% y_limits = ylim; % Get current y-axis limits for this subplot
% y_shade = [y_limits(1), y_limits(1), y_limits(2), y_limits(2)];
% patch(x_shade, y_shade, 'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% hold off;

% Define window size for smoothing (100 ms)
window_size_ms = 100; % 100 milliseconds
window_size_samples = round(window_size_ms / dt); % Convert to number of samples

% Ensure window size is at least 1
window_size_samples = max(window_size_samples, 1);

% Display the figure for Mean Recurrent Weight
figure(6);
plot(meanwin, 'Color', line_colors(1,:), 'LineWidth', 1.5);
xlabel('Trial Number');
ylabel('Mean Recurrent Weight in Decision Network');
title('Evolution of Mean Recurrent Weight in Decision Network over Trials');
% grid on; % Optional: Add grid for better readability

% Plot the mean LTP and LTD eligibility traces for each trial
figure(3);

% Define window size for smoothing (100 ms)
window_size_ms = 50; % 100 milliseconds
window_size_samples = round(window_size_ms / dt); % Convert to number of samples

% Ensure window size is at least 1
window_size_samples = max(window_size_samples, 1);

% === Plotting for First Trials ===
figure;

% Subplot 1: Input Layer
subplot(3,1,1);
maxY = -inf; % Initialize maxY
% Input layer
% Plot neuron 1-3 in input layer for Choice 1 (blue)
for neuron = 1:3
    smoothed_input_choice1 = smoothdata(first_trial_choice1_input(neuron, :) * 1000, 'movmean', window_size_samples);
    plot(0:dt:t_total - dt, smoothed_input_choice1,'-', 'Color', line_colors(1,:),'LineWidth', 2);
    hold on;
    maxY = max(maxY, max(smoothed_input_choice1));
end
% Plot neuron 1-3 in input layer for Choice 2 (red)
for neuron = 1:3
    smoothed_input_choice2 = smoothdata(first_trial_choice2_input(neuron, :) * 1000, 'movmean', window_size_samples);
    plot(0:dt:t_total - dt, smoothed_input_choice2, '-', 'Color', line_colors(3,:),'LineWidth', 2);
    hold on;
    maxY = max(maxY, max(smoothed_input_choice2));
end

% Add vertical gray lines
xline(500, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 2);

% Remove title
% title('Input Layer Neuron Dynamics - First Trials');

% Remove x-label
% xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
xlim([0 1500]); % Set x-axis limits to 0-1.5s (1500 ms)
ylim([0, maxY * 1.2]); % Adjust y-axis limits based on data
% legend('Choice 1', 'Choice 2');
% grid on; % Optional
hold off;

% Subplot 2: Decision Network
subplot(3,1,2);
maxY = -inf; % Initialize maxY
% Decision network
% Plot firing rates of neuron 1-3 in decision net for Choice 1 (blue)
for neuron = 1:3
    smoothed_decision_choice1 = smoothdata(first_trial_choice1_decision(neuron, :) * 1000, 'movmean', window_size_samples);
    plot(0:dt:t_total - dt, smoothed_decision_choice1, '-', 'Color', line_colors(1,:), 'LineWidth', 2, 'DisplayName', 'Choice 1 Decision');
    hold on;
    maxY = max(maxY, max(smoothed_decision_choice1));
end
% Plot firing rates of neuron 1-3 in decision net for Choice 2 (red)
for neuron = 1:3
    smoothed_decision_choice2 = smoothdata(first_trial_choice2_decision(neuron, :) * 1000, 'movmean', window_size_samples);
    plot(0:dt:t_total - dt, smoothed_decision_choice2, '-', 'Color', line_colors(3,:), 'LineWidth', 2, 'DisplayName', 'Choice 2 Decision');
    hold on;
    maxY = max(maxY, max(smoothed_decision_choice2));
end

% Add vertical gray lines
xline(500, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 2);

% Remove title
% title('Decision Network Neuron Dynamics - First Trials');

% Remove x-label
% xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
xlim([0 1500]);
ylim([0, maxY * 1.2]); % Adjust y-axis limits based on data
% grid on; % Optional
hold off;

% Subplot 3: Timer Network
subplot(3,1,3);
maxY = -inf; % Initialize maxY
% Timer network
% Plot firing rates of neuron 1-3 in timer net for Choice 1 (blue)
for neuron = 1:3
    smoothed_timer_choice1 = smoothdata(first_trial_choice1_r_timer(neuron, :) * 1000, 'movmean', window_size_samples);
    plot(0:dt:t_total - dt, smoothed_timer_choice1, '-', 'Color', line_colors(1,:), 'LineWidth', 2, 'DisplayName', 'Choice 1 Timer');
    hold on;
    maxY = max(maxY, max(smoothed_timer_choice1));
end
% Plot firing rates of neuron 1-3 in timer net for Choice 2 (red)
for neuron = 1:3
    smoothed_timer_choice2 = smoothdata(first_trial_choice2_r_timer(neuron, :) * 1000, 'movmean', window_size_samples);
    plot(0:dt:t_total - dt, smoothed_timer_choice2, '-', 'Color', line_colors(3,:), 'LineWidth', 2, 'DisplayName', 'Choice 2 Timer');
    hold on;
    maxY = max(maxY, max(smoothed_timer_choice2));
end

% Add vertical gray lines
xline(500, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 2);

% Remove title
% title('Timer Network Neuron Dynamics - First Trials');

% Keep x-label only for the bottom plot
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
xlim([0 1500]);
ylim([0, maxY * 1.2]); % Adjust y-axis limits based on data
% legend('Choice 1', 'Choice 2');
% grid on; % Optional
hold off;

% === Plotting for Last Trials ===
figure;

% Subplot 1: Input Layer
subplot(3,1,1);
maxY = -inf; % Initialize maxY
% Input layer
% Plot neuron 1-3 in input layer for Choice 1 (blue)
for neuron = 1:3
    smoothed_input_choice1_last = smoothdata(last_trial_choice1_input(neuron, :) * 1000, 'movmean', window_size_samples);
    plot(0:dt:t_total - dt, smoothed_input_choice1_last, '-', 'Color', line_colors(1,:), 'LineWidth', 2, 'DisplayName', 'Choice 1 Input');
    hold on;
    maxY = max(maxY, max(smoothed_input_choice1_last));
end
% Plot neuron 1-3 in input layer for Choice 2 (red)
for neuron = 1:3
    smoothed_input_choice2_last = smoothdata(last_trial_choice2_input(neuron, :) * 1000, 'movmean', window_size_samples);
    plot(0:dt:t_total - dt, smoothed_input_choice2_last, '-', 'Color', line_colors(3,:), 'LineWidth', 2, 'DisplayName', 'Choice 2 Input');
    hold on;
    maxY = max(maxY, max(smoothed_input_choice2_last));
end

% Add vertical gray lines
xline(500, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 2);

% Remove title
% title('Input Layer Neuron Dynamics - Last Trials');

% Remove x-label
% xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
xlim([0 1500]);
ylim([0, maxY * 1.2]); % Adjust y-axis limits based on data
% legend('Choice 1', 'Choice 2');
% grid on; % Optional
hold off;

% Subplot 2: Decision Network
subplot(3,1,2);
maxY = -inf; % Initialize maxY
% Decision network
% Plot firing rates of neuron 1-3 in decision net for Choice 1 (blue)
for neuron = 1:3
    smoothed_decision_choice1_last = smoothdata(last_trial_choice1_decision(neuron, :) * 1000, 'movmean', window_size_samples);
    plot(0:dt:t_total - dt, smoothed_decision_choice1_last, '-', 'Color', line_colors(1,:), 'LineWidth', 2, 'DisplayName', 'Choice 1 Decision');
    hold on;
    maxY = max(maxY, max(smoothed_decision_choice1_last));
end
% Plot firing rates of neuron 1-3 in decision net for Choice 2 (red)
for neuron = 1:3
    smoothed_decision_choice2_last = smoothdata(last_trial_choice2_r_ex(neuron, :) * 1000, 'movmean', window_size_samples);
    plot(0:dt:t_total - dt, smoothed_decision_choice2_last, '-', 'Color', line_colors(3,:), 'LineWidth', 2, 'DisplayName', 'Choice 2 Decision');
    hold on;
    maxY = max(maxY, max(smoothed_decision_choice2_last));
end

% Add vertical gray lines
xline(500, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 2);

% Remove title
% title('Decision Network Neuron Dynamics - Last Trials');

% Remove x-label
% xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
xlim([0 1500]);
ylim([0, maxY * 1.2]); % Adjust y-axis limits based on data
% grid on; % Optional
hold off;

% Subplot 3: Timer Network
subplot(3,1,3);
maxY = -inf; % Initialize maxY
% Timer network
% Plot firing rates of neuron 1-3 in timer net for Choice 1 (blue)
for neuron = 1:3
    smoothed_timer_choice1_last = smoothdata(last_trial_choice1_r_timer(neuron, :) * 1000, 'movmean', window_size_samples);
    plot(0:dt:t_total - dt, smoothed_timer_choice1_last, '-', 'Color', line_colors(1,:), 'LineWidth', 2, 'DisplayName', 'Choice 1 Timer');
    hold on;
    maxY = max(maxY, max(smoothed_timer_choice1_last));
end
% Plot firing rates of neuron 1-3 in timer net for Choice 2 (red)
for neuron = 1:3
    smoothed_timer_choice2_last = smoothdata(last_trial_choice2_r_timer(neuron, :) * 1000, 'movmean', window_size_samples);
    plot(0:dt:t_total - dt, smoothed_timer_choice2_last, '-', 'Color', line_colors(3,:), 'LineWidth', 2, 'DisplayName', 'Choice 2 Timer');
    hold on;
    maxY = max(maxY, max(smoothed_timer_choice2_last));
end

% Add vertical gray lines
xline(500, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 2);

% Remove title
% title('Timer Network Neuron Dynamics - Last Trials');

% Keep x-label only for the bottom plot
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
xlim([0 1500]);
ylim([0, maxY * 1.2]); % Adjust y-axis limits based on data
% legend('Choice 1', 'Choice 2');
% grid on; % Optional
hold off;

%% Plotting the Crossing Times
figure;

plot(time_trace_cross(1:2:end), time_activity_cross(1:2:end), 'o', 'Color', line_colors(1,:), 'MarkerFaceColor', line_colors(1,:));
hold on;
plot(time_trace_cross(2:2:end), time_activity_cross(2:2:end), 'o', 'Color', line_colors(3,:), 'MarkerFaceColor', line_colors(3,:));
hold on;
% plot(400:2500,400:2500,'--','Color',[0.5 0.5 0.5],'LineWidth',2)
xlabel('Trace crossing time (ms)');
ylabel('Threshold crossing time (ms)');
title('Time when Mean Timer Activity crosses 20 Hz per Trial');
% grid on;
% axis([300 2500 300 2500])


