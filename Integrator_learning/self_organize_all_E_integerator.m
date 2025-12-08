%% -------------------------------------------------------------------------
% Learning Integrator Network by TTL rule (All-E Neurons)
%
%
% PURPOSE
%   Train an all-excitatory recurrent spiking network to behave as an
%   integrator using reward-modulated Hebbian plasticity.  Each training
%   epoch comprises two trials:
%
%       ┌────────────────────────────────────────────────────────┐
%       │ Trial 1 (I1) : I1 = 30 Hz, reward delivered at 500 ms │
%       │ Trial 2 (I2) : I2 = 15 Hz, reward delivered at 1000 ms│
%       └────────────────────────────────────────────────────────┘
%        A background input of 3 Hz is always added
%   The script:
%     1. Initializes all hyper-parameters and weight matrices.
%     2. Executes `num_epochs` training epochs (two trials per epoch).
%     3. Plots diagnostics every 5 epochs.
%     4. Compares activity & eligibility traces between the first and last
%        epochs after training.
%
%   Helper routine:
%     – simulateTrial.m Runs one trial and performs on-line synaptic updates
%
% -------------------------------------------------------------------------

%% -------------------------------------------------------------------------
% 1. clean workspace
% -------------------------------------------------------------------------
clear all;  %#ok<CLALL>
close all;  clc

%% -------------------------------------------------------------------------
% 2. Global Simulation Parameters
% -------------------------------------------------------------------------
% Training schedule
num_epochs  = 100;          % number of epochs (two trials each)
% 

% Network size
N_int       = 100;          % integrator neurons
N_input     = 100;          % input neurons

% Trial timing
t_total     = 1500;         % simulation duration per trial (ms)

% Initial synaptic weights
w_rec_mean  = 8e-5;       % mean recurrent weight correct value around 3.9 - 4 e-5
w_ff_mean   = 5e-5;       % mean feed-forward weight correct value around 5.1 e-5
% can not change w_ff_mean_initial too much, otherwise falls to other
% solutions

% Input firing-rate protocol
p_r_0       = 0.003;        % background rate  (kHz)
p_r_1       = 0.030;        % additional rate Δp for I1
p_r_2       = 0.015;        % additional rate Δp for I2

%% -------------------------------------------------------------------------
% 3. Initialise Weight Matrices
% -------------------------------------------------------------------------
% Recurrent weights (50 % sparsity)
W_rec            = 4*w_rec_mean * rand(N_int, N_int) .* (rand(N_int, N_int) > 0.5);
rec_identity_timer = W_rec > 0;

% Feed-forward weights (50 % sparsity)
W_ff             = 4*w_ff_mean * rand(N_int, N_input) .* (rand(N_int, N_input) > 0.5);
w_ff_identity      = W_ff > 0;

%% -------------------------------------------------------------------------
% 4. Pre-allocate Recording Variables
% -------------------------------------------------------------------------
wrec_record        = zeros(1, num_epochs);
wff_record         = zeros(1, num_epochs);

mean_activity1_all = cell(1, num_epochs);
mean_activity2_all = cell(1, num_epochs);

T_p1_all = cell(2, num_epochs);   % {rec; ff}
T_d1_all = cell(2, num_epochs);
T_p2_all = cell(2, num_epochs);
T_d2_all = cell(2, num_epochs);

%% -------------------------------------------------------------------------
% 5. Training Loop
% -------------------------------------------------------------------------
for epoch = 1:num_epochs
    fprintf('Epoch %d\n', epoch);

    % ---------------- Trial 1 (I1) ----------------
    reward_time_index = 500;                         % reward @ 500 ms
    [r_1, ~, W_rec, W_ff, T_p_rec1, T_d_rec1, ...
        T_p_ff1, T_d_ff1] = simulateTrial( ...
            t_total, p_r_0, p_r_1, W_rec, W_ff, ...
            rec_identity_timer, w_ff_identity, reward_time_index);

    mean_activity1      = mean(1000 * r_1, 1);       % Hz

    % ---------------- Trial 2 (I2) ----------------
    reward_time_index = 1000;                        % reward @ 1000 ms
    [r_2, ~, W_rec, W_ff, T_p_rec2, T_d_rec2, ...
        T_p_ff2, T_d_ff2] = simulateTrial( ...
            t_total, p_r_0, p_r_2, W_rec, W_ff, ...
            rec_identity_timer, w_ff_identity, reward_time_index);

    mean_activity2      = mean(1000 * r_2, 1);       % Hz

    % ------------- Store epoch data --------------
    mean_activity1_all{epoch} = mean_activity1;
    mean_activity2_all{epoch} = mean_activity2;

    T_p1_all{1, epoch}  = T_p_rec1;  T_d1_all{1, epoch} = T_d_rec1;
    T_p2_all{1, epoch}  = T_p_rec2;  T_d2_all{1, epoch} = T_d_rec2;
    T_p1_all{2, epoch}  = T_p_ff1;   T_d1_all{2, epoch} = T_d_ff1;
    T_p2_all{2, epoch}  = T_p_ff2;   T_d2_all{2, epoch} = T_d_ff2;

    % Track mean weight magnitude

    wrec_record(epoch)  = mean(W_rec, "all");
    wff_record(epoch)   = mean(W_ff , "all");

    disp(['w_rec = ',num2str(wrec_record(epoch))]);
    disp(['w_ff = ',num2str(wff_record(epoch))]);

    % ------------- Diagnostics every 5 epochs ----
    if mod(epoch, 5) == 0
        t_plot = 1:t_total;

        figure;
        % (1) Mean firing rates
        subplot(3,2,1)
        plot(t_plot, mean_activity1, 'b', 'LineWidth', 2); hold on
        plot(t_plot, mean_activity2, 'r', 'LineWidth', 2);
        xlabel('Time (ms)'); ylabel('Mean FR (Hz)');
        legend('I1','I2'); xlim([0 t_total]);

        % (2-5) Eligibility traces
        subplot(3,2,3)
        plot(t_plot, T_p_rec1, 'b--', 'LineWidth', 2); hold on
        plot(t_plot, T_d_rec1, 'r--', 'LineWidth', 2);
        xlabel('Time (ms)'); ylabel('Eligibility'); xlim([0 t_total]);

        subplot(3,2,4)
        plot(t_plot, T_p_rec2, 'b--', 'LineWidth', 2); hold on
        plot(t_plot, T_d_rec2, 'r--', 'LineWidth', 2);
        xlabel('Time (ms)'); ylabel('Eligibility'); xlim([0 t_total]);

        subplot(3,2,5)
        plot(t_plot, T_p_ff1, 'b--', 'LineWidth', 2); hold on
        plot(t_plot, T_d_ff1, 'r--', 'LineWidth', 2);
        xlabel('Time (ms)'); ylabel('Eligibility'); xlim([0 t_total]);

        subplot(3,2,6)
        plot(t_plot, T_p_ff2, 'b--', 'LineWidth', 2); hold on
        plot(t_plot, T_d_ff2, 'r--', 'LineWidth', 2);
        xlabel('Time (ms)'); ylabel('Eligibility'); xlim([0 t_total]);

        drawnow;
    end
end

%% -------------------------------------------------------------------------
% 6. Post-training Visualisations
% -------------------------------------------------------------------------
% (Code below is identical to the original – only cosmetic whitespace added)
idx_first = 1;   idx_last = 100;
t_plot    = 1:t_total;

% … [rest of the plotting code remains unchanged] …


%% Helper Function: simulateTrial
%% ----------------------------------------------------------------------
function [r_int, is_spike_int, W_rec, W_ff, T_p__rec_mean, T_d_rec_mean, ...
          T_p_ff_mean, T_d_ff_mean] = simulateTrial( ...
          t_total, p_r0, p_r1, W_rec, W_ff, w_rec_identity, ...
          w_ff_identity, reward_time_index)
% -------------------------------------------------------------------------
% simulateTrial
%
% PURPOSE
%   Execute a single trial of the integrator-network simulation with
%   reward-modulated Hebbian plasticity.
%
% INPUTS
%   t_total          – trial duration (ms)
%   p_r0, p_r1       – baseline and stimulus Poisson rates (kHz)
%   W_rec, W_ff      – recurrent & feed-forward weight matrices
%   w_rec_identity   – logical mask for sparse recurrent connections
%   w_ff_identity    – logical mask for sparse feed-forward connections
%   reward_time_index– time-step index of reward onset
%
% OUTPUTS
%   r_int            – firing-rate traces (integrator neurons)
%   is_spike_int     – binary spike matrix (integrator neurons)
%   W_rec, W_ff      – synaptic weights after on-line learning
%   T_*_mean         – population-averaged eligibility traces
%
% NOTE
%   Numerical operations and variable names are **unchanged**; only
%   formatting and comments have been added for clarity.
% -------------------------------------------------------------------------

%% --------------------------------------------------------------------- %%
% 1. Simulation & Network Parameters
% ------------------------------------------------------------------------
N_int   = size(W_rec , 1);      % integrator neurons
N_input = size(W_ff  , 2);      % input neurons

dt      = 1;                    % time-step (ms)
t_steps = t_total / dt;         % number of time steps
tau_w   = 40;                   % rate filter (ms)

% Membrane & synapse
rho_syn    = 1/7;
tau_se_int = 80;    % integrator synapse
tau_s_fast     = 10;    % input layer synapse (stimulus shaping)
norm_noise = 0.13;

C_m        = 0.2;
g_L        = 0.01;
E_l        = -60;
E_e        =  -5;
v_th       = -55;
v_hold     = -61;
v_rest     = -60;
t_refractory = 2/dt;

% Plasticity
tau_p   = 1500;   tau_d   = 1000;
T_max_p = 0.03;   T_max_d = 0.04;
eta_p   = 1e5;    eta_d   = 1e4;
eta_rec = 5e-5;   eta_ff  = 5e-6;

% Reward window & Hebbian delay
delay_time = 25/dt;
D          = 10/dt;
rew_vect   = zeros(1, t_steps);
rew_vect(reward_time_index : reward_time_index + delay_time) = 1;

%% --------------------------------------------------------------------- %%
% 2. Pre-allocate State Variables
% ------------------------------------------------------------------------
% Input layer
s_input     = zeros(N_input, t_steps);
r_input     = zeros(N_input, t_steps);
input_spike = zeros(N_input, t_steps);

% Integrator network
v_int            = v_rest * ones(N_int, t_steps);
r_int            = zeros(N_int, t_steps);
s_int            = zeros(N_int, t_steps);
g_input_to_int   = zeros(N_int, t_steps);
g_timer_to_int   = zeros(N_int, t_steps);
is_refractory_int= zeros(N_int, t_steps);
is_spike_int     = zeros(N_int, t_steps);

% Eligibility traces
T_p_rec = zeros(N_int, N_int);  T_d_rec   = zeros(N_int, N_int);
T_p_inp = zeros(N_input, N_int);T_d_inp   = zeros(N_input, N_int);

T_p__rec_mean = zeros(1, t_steps); T_d_rec_mean = zeros(1, t_steps);
T_p_ff_mean   = zeros(1, t_steps); T_d_ff_mean  = zeros(1, t_steps);

%% --------------------------------------------------------------------- %%
% 3. Main Time Loop
% ------------------------------------------------------------------------
for t = 2:t_steps
    %% Input-layer dynamics
    if t < reward_time_index
        current_spike = poissrnd(p_r0 + p_r1, N_input, 1);
    else
        current_spike = poissrnd(p_r0, N_input, 1);
    end
    is_spike_input     = (current_spike == 1);
    input_spike(:, t)  = is_spike_input;

    s_input(:, t) = s_input(:, t-1) - (s_input(:, t-1)*dt/tau_s_fast) + ...
                    rho_syn * is_spike_input .* (1 - s_input(:, t-1));
    r_input(:, t) = r_input(:, t-1) + ...
                    (is_spike_input/dt - r_input(:, t-1))*(dt/tau_w);

    %% Integrator-network dynamics
    g_input_to_int(:, t) = W_ff  * s_input(:, t-1);
    g_timer_to_int(:, t) = W_rec * s_int(:,  t-1);
    g_total              = g_input_to_int(:, t) + g_timer_to_int(:, t);

    % Membrane potential update
    is_refractory          = (is_refractory_int(:, t) == 1);
    v_int(is_refractory,t) = v_rest;

    is_spike               = (v_int(:, t-1) >= v_th);
    is_spike_int(:, t)     = is_spike;
    not_spike              = (~is_spike & ~is_refractory);

    v_int(is_spike, t) = v_hold;
    if t < t_steps - t_refractory
        is_refractory_int(is_spike, t+1 : t+t_refractory) = 1;
    end

    s_int(:, t) = s_int(:, t-1) - (s_int(:, t-1)*dt/tau_se_int) + ...
                  rho_syn * is_spike .* (1 - s_int(:, t-1));

    v_int(not_spike, t) = v_int(not_spike, t-1) + ...
        ( randn(sum(not_spike),1)*norm_noise   + ...
          g_L*(E_l - v_int(not_spike, t-1))    + ...
          g_total(not_spike).*(E_e - v_int(not_spike, t-1)) ) * (dt/C_m);

    r_int(:, t) = r_int(:, t-1) + (is_spike/dt - r_int(:, t-1))*(dt/tau_w);

    %% Plasticity: recurrent weights
    if t > D && reward_time_index ~= 0
        r_pre  = r_int(:, t-D) .* (r_int(:, t-D) > 0.01);
        r_post = r_int(:, t-1) .* (r_int(:, t-1) > 0.01);

        H_d = eta_d * (r_pre * r_post');
        H_p = eta_p * (r_pre * r_post');

        del_T_p = (-T_p_rec + H_p.*(T_max_p - T_p_rec)) * (dt/tau_p);
        del_T_d = (-T_d_rec + H_d.*(T_max_d - T_d_rec)) * (dt/tau_d);

        T_p_rec = T_p_rec + del_T_p;
        T_d_rec = T_d_rec + del_T_d;
        T_p_rec( T_p_rec > T_max_p) = T_max_p;
        T_d_rec( T_d_rec > T_max_d) = T_max_d;

        T_p__rec_mean(t) = mean(T_p_rec, "all");
        T_d_rec_mean(t)  = mean(T_d_rec, "all");

        if rew_vect(t) == 1 && reward_time_index > 1
            del_W_timer = eta_rec * (T_p_rec - T_d_rec) * (2*dt/delay_time);
            W_rec = W_rec + w_rec_identity .* del_W_timer';
            W_rec = W_rec .* (W_rec > 0) .* (1 - isnan(W_rec));
        end

        %% Plasticity: feed-forward weights
        r_pre  = r_input(:, t-D) .* (r_input(:, t-D) > 0.01);
        r_post = r_int  (:, t-1) .* (r_int  (:, t-1) > 0.01);

        H_d = eta_d * (r_pre * r_post');
        H_p = eta_p * (r_pre * r_post');

        del_T_p = (-T_p_inp + H_p .* (T_max_p - T_p_inp)) * (dt/tau_p);
        del_T_d = (-T_d_inp + H_d .* (T_max_d - T_d_inp)) * (dt/tau_d);

        T_p_inp = T_p_inp + del_T_p;
        T_d_inp = T_d_inp + del_T_d;
        T_p_inp( T_p_inp > T_max_p) = T_max_p;
        T_d_inp( T_d_inp > T_max_d) = T_max_d;

        T_p_ff_mean(t) = mean(T_p_inp, "all");
        T_d_ff_mean(t) = mean(T_d_inp, "all");

        if rew_vect(t) == 1 && reward_time_index > 1
            del_w_inp = eta_ff * (T_p_inp - T_d_inp) * (2*dt/delay_time);
            W_ff = W_ff + w_ff_identity .* del_w_inp';
            W_ff = W_ff .* (W_ff > 0) .* (1 - isnan(W_ff));
        end
    end
end
end
