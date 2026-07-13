%% -------------------------------------------------------------------------
%   Post-training analysis, need to run self_organize_all_E_integrator.m
%   first
%   
%   mean ± SEM of LI vs. <W_rec> and save all per-trial mean traces
% -------------------------------------------------------------------------
Wrec_f_mean0 = 3.95e-5;                         %  reference (Wrec0 for labeling)

% Keep trained weights as the base, then rescale mean recurrent strength
W_rec_trained = W_rec;
W_ff_trained  = W_ff;
orig_mean_rec = mean(W_rec_trained, "all");

% (2) sweep target means
nCond             = 21;
Wrec_target_means = linspace(0.7, 1.3, nCond) * orig_mean_rec;
xvals_ratio       = Wrec_target_means / orig_mean_rec;   % normalized x-axis
rate = orig_mean_rec/Wrec_f_mean0;
% Inputs to test
p_r1_I1   = p_r_1;
p_r1_I2   = p_r_2;
p_r1_IMID = 0.5*(p_r_1 + p_r_2);

% (3) repeats
nRepeats = 20;
t_total  = 2000;

% Disable learning during evaluation by zeroing masks; keep stimulus on
reward_time_index_test = 2000;          % stimulus lasts ~full trial
zero_mask_rec = zeros(size(W_rec_trained));
zero_mask_ff  = zeros(size(W_ff_trained));

% -------------------------------------------------------------------------
% Storage for per-trial population-mean FR traces (Hz) for every condition
% Dimensions: [nCond, nRepeats, t_total]
% -------------------------------------------------------------------------
meanFR_I1   = zeros(nCond, nRepeats, t_total);
meanFR_I2   = zeros(nCond, nRepeats, t_total);
meanFR_IMID = zeros(nCond, nRepeats, t_total);

% -------------------------------------------------------------------------
% Storage for per-trial LI (R^2) for every condition
% Dimensions: [nCond, nRepeats]
% -------------------------------------------------------------------------
LI_I1   = nan(nCond, nRepeats);
LI_I2   = nan(nCond, nRepeats);
LI_IMID = nan(nCond, nRepeats);

% (Optional) also keep R^2 from trial-averaged curves if you still want it
R2_I1_mean   = zeros(1, nCond);
R2_I2_mean   = zeros(1, nCond);
R2_IMID_mean = zeros(1, nCond);

% Utility for robust ramp-fit index selection
get_fit_idx = @(y) ( ...
    function_handle_select_fit_idx(y) ...
);

for ii = 1:nCond
    fprintf('Condition %d / %d\n', ii, nCond);

    % Rescale recurrent weights to target mean (preserve sparsity pattern)
    scale_factor   = Wrec_target_means(ii) / orig_mean_rec;
    W_rec_testbase = W_rec_trained * scale_factor;

    % Accumulators for trial-averaged curves (for optional R^2 on averages)
    yacc_I1   = zeros(1, t_total);
    yacc_I2   = zeros(1, t_total);
    yacc_IMID = zeros(1, t_total);

    % =========================
    % -------- I1 trials ------
    % =========================
    for rr = 1:nRepeats
        [r_eval, ~, ~, ~, ~, ~, ~, ~] = simulateTrial( ...
            t_total, p_r_0, p_r1_I1, ...
            W_rec_testbase, W_ff_trained, ...
            zero_mask_rec, zero_mask_ff, ...
            reward_time_index_test);

        y_cur = mean(1000*r_eval, 1);                % Hz, pop. mean FR
        meanFR_I1(ii, rr, :) = y_cur;                % save per-trial trace
        yacc_I1 = yacc_I1 + y_cur;                   % for avg curve

        % ---- per-trial LI (R^2) ----
        fit_idx = get_fit_idx(y_cur);
        x   = fit_idx(:);
        y   = y_cur(fit_idx).';
        X   = [x, ones(size(x))];
        bet = X \ y;                   % OLS with intercept
        yhat= X*bet;
        SSE = sum((y - yhat).^2);
        SST = sum((y - mean(y)).^2 + eps);
        LI_I1(ii, rr) = max(0, 1 - SSE/SST);   % clamp at [0,1]
    end
    y_I1 = yacc_I1 / nRepeats;   % average curve (optional visualization)

    % =========================
    % -------- I2 trials ------
    % =========================
    for rr = 1:nRepeats
        [r_eval, ~, ~, ~, ~, ~, ~, ~] = simulateTrial( ...
            t_total, p_r_0, p_r1_I2, ...
            W_rec_testbase, W_ff_trained, ...
            zero_mask_rec, zero_mask_ff, ...
            reward_time_index_test);

        y_cur = mean(1000*r_eval, 1);
        meanFR_I2(ii, rr, :) = y_cur;
        yacc_I2 = yacc_I2 + y_cur;

        fit_idx = get_fit_idx(y_cur);
        x   = fit_idx(:);
        y   = y_cur(fit_idx).';
        X   = [x, ones(size(x))];
        bet = X \ y;
        yhat= X*bet;
        SSE = sum((y - yhat).^2);
        SST = sum((y - mean(y)).^2 + eps);
        LI_I2(ii, rr) = max(0, 1 - SSE/SST);
    end
    y_I2 = yacc_I2 / nRepeats;

    % =======================================
    % -------- (I1+I2)/2 (mid) trials -------
    % =======================================
    for rr = 1:nRepeats
        [r_eval, ~, ~, ~, ~, ~, ~, ~] = simulateTrial( ...
            t_total, p_r_0, p_r1_IMID, ...
            W_rec_testbase, W_ff_trained, ...
            zero_mask_rec, zero_mask_ff, ...
            reward_time_index_test);

        y_cur = mean(1000*r_eval, 1);
        meanFR_IMID(ii, rr, :) = y_cur;
        yacc_IMID = yacc_IMID + y_cur;

        fit_idx = get_fit_idx(y_cur);
        x   = fit_idx(:);
        y   = y_cur(fit_idx).';
        X   = [x, ones(size(x))];
        bet = X \ y;
        yhat= X*bet;
        SSE = sum((y - yhat).^2);
        SST = sum((y - mean(y)).^2 + eps);
        LI_IMID(ii, rr) = max(0, 1 - SSE/SST);
    end
    y_IMID = yacc_IMID / nRepeats;

    % -------- Optional: R^2 on trial-averaged curves --------
    fit_idx = get_fit_idx(y_I1); x = fit_idx(:); y = y_I1(fit_idx).'; X=[x,ones(size(x))];
    bet = X\y; yhat = X*bet; SSE = sum((y-yhat).^2); SST = sum((y-mean(y)).^2 + eps);
    R2_I1_mean(ii) = max(0, 1 - SSE/SST);

    fit_idx = get_fit_idx(y_I2); x = fit_idx(:); y = y_I2(fit_idx).'; X=[x,ones(size(x))];
    bet = X\y; yhat = X*bet; SSE = sum((y-yhat).^2); SST = sum((y-mean(y)).^2 + eps);
    R2_I2_mean(ii) = max(0, 1 - SSE/SST);

    fit_idx = get_fit_idx(y_IMID); x = fit_idx(:); y = y_IMID(fit_idx).'; X=[x,ones(size(x))];
    bet = X\y; yhat = X*bet; SSE = sum((y-yhat).^2); SST = sum((y-mean(y)).^2 + eps);
    R2_IMID_mean(ii) = max(0, 1 - SSE/SST);
end

%% (5) Plot mean ± SD of LI vs. <W_rec>
LI_I1_mean   = mean(LI_I1,   2);  LI_I1_sem   = std(LI_I1,   0, 2)/sqrt(nRepeats);
LI_I2_mean   = mean(LI_I2,   2);  LI_I2_sem   = std(LI_I2,   0, 2)/sqrt(nRepeats);
LI_IMID_mean = mean(LI_IMID, 2);  LI_IMID_sem = std(LI_IMID, 0, 2)/sqrt(nRepeats);

figure; hold on; box on; grid on;
errorbar(xvals_ratio/ 0.9762, LI_I1_mean,   LI_I1_sem,   '-o', 'LineWidth', 1.8, 'MarkerSize', 5);
errorbar(xvals_ratio/ 0.9762, LI_I2_mean,   LI_I2_sem,   '-s', 'LineWidth', 1.8, 'MarkerSize', 5);
errorbar(xvals_ratio/ 0.9762, LI_IMID_mean, LI_IMID_sem, '-^', 'LineWidth', 1.8, 'MarkerSize', 5);
xlabel('\langle W_{rec} \rangle / \langle W_{rec} \rangle_{trained}');
ylabel('LI (R^2), mean \pm SD');
legend({'I_1','I_2','(I_1+I_2)/2'}, 'Location','best');
title('Linearity (per-trial LI) vs. mean recurrent weight');

% (Optional) If you still want the R^2 from trial-averaged curves:
% figure; hold on; box on; grid on;
% plot(xvals_ratio, R2_I1_mean,   '-o', 'LineWidth', 2, 'MarkerSize', 5);
% plot(xvals_ratio, R2_I2_mean,   '-s', 'LineWidth', 2, 'MarkerSize', 5);
% plot(xvals_ratio, R2_IMID_mean, '-^', 'LineWidth', 2, 'MarkerSize', 5);
% xlabel('\langle W_{rec} \rangle / \langle W_{rec} \rangle_{trained}');
% ylabel('R^2 of linear fit (avg curve)');
% legend({'I_1','I_2','(I_1+I_2)/2'}, 'Location','best');
% title('Linearity (avg-curve R^2) vs. mean recurrent weight');

%% Save all per-trial mean FR traces and LI arrays


%% ----------------------------------------------------------------------
% Helper for robust window selection (10%→70% of late fixed point).
% Implemented as a local function body via a function-handle wrapper above.
% This avoids code duplication while keeping everything in this file.
function fit_idx = function_handle_select_fit_idx(y)
    T = numel(y);
    tail_start = max(1, T-100);
    fp = mean(y(tail_start:T));
    if ~isfinite(fp) || fp <= 0
        fp = max(1, mean(y(max(1, round(0.8*T)):T))); % fallback
    end
    fp = max(fp, 40);
    i10 = find(y >= 0.0*fp, 1, 'first');
    i70 = find(y >= 0.70*fp, 1, 'first');
    if isempty(i10), i10 = max(1, round(0.2*T)); end
    if isempty(i70), i70 = min(T-10, i10 + 300); end
    
    if i70 <= i10 + 5
        i70 = min(T-10, i10 + 100);
    end
    fit_idx = i10:i70;
end

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
% NOTE: Learning is effectively disabled during evaluation when identity
%       masks (w_rec_identity, w_ff_identity) are zeros.
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
tau_se_int     = 80;            % integrator synapse
tau_s_fast     = 10;            % input layer synapse (stimulus shaping)
norm_noise     = 0.13;

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
t_steps_i  = t_steps;
rew_vect   = zeros(1, t_steps_i);
if reward_time_index > 0
    idx_end = min(t_steps_i, reward_time_index + delay_time);
    rew_vect(reward_time_index:idx_end) = 1;
end

%% --------------------------------------------------------------------- %%
% 2. Pre-allocate State Variables
% ------------------------------------------------------------------------
% Input layer
s_input     = zeros(N_input, t_steps_i);
r_input     = zeros(N_input, t_steps_i);
input_spike = zeros(N_input, t_steps_i);

% Integrator network
v_int            = v_rest * ones(N_int, t_steps_i);
r_int            = zeros(N_int, t_steps_i);
s_int            = zeros(N_int, t_steps_i);
g_input_to_int   = zeros(N_int, t_steps_i);
g_timer_to_int   = zeros(N_int, t_steps_i);
is_refractory_int= zeros(N_int, t_steps_i);
is_spike_int     = zeros(N_int, t_steps_i);

% Eligibility traces
T_p_rec = zeros(N_int, N_int);  T_d_rec   = zeros(N_int, N_int);
T_p_inp = zeros(N_input, N_int);T_d_inp   = zeros(N_input, N_int);

T_p__rec_mean = zeros(1, t_steps_i); T_d_rec_mean = zeros(1, t_steps_i);
T_p_ff_mean   = zeros(1, t_steps_i); T_d_ff_mean  = zeros(1, t_steps_i);

%% --------------------------------------------------------------------- %%
% 3. Main Time Loop
% ------------------------------------------------------------------------
for t = 2:t_steps_i
    %% Input-layer dynamics
    if t < reward_time_index || reward_time_index == 0
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
    if t < t_steps_i - t_refractory
        is_refractory_int(is_spike, t+1 : t+t_refractory) = 1;
    end

    s_int(:, t) = s_int(:, t-1) - (s_int(:, t-1)*dt/tau_se_int) + ...
                  rho_syn * is_spike .* (1 - s_int(:, t-1));

    v_int(not_spike, t) = v_int(not_spike, t-1) + ...
        ( randn(sum(not_spike),1)*norm_noise   + ...
          g_L*(E_l - v_int(not_spike, t-1))    + ...
          g_total(not_spike).*(E_e - v_int(not_spike, t-1)) ) * (dt/C_m);

    r_int(:, t) = r_int(:, t-1) + (is_spike/dt - r_int(:, t-1))*(dt/tau_w);

    %% Plasticity (disabled if masks are zeros)
    if t > D && reward_time_index ~= 0
        % Recurrent
        r_pre  = r_int(:, t-D) .* (r_int(:, t-D) > 0.01);
        r_post = r_int(:, t-1) .* (r_int(:, t-1) > 0.01);
        H_d = eta_d * (r_pre * r_post');   H_p = eta_p * (r_pre * r_post');
        del_T_p = (-T_p_rec + H_p.*(T_max_p - T_p_rec)) * (dt/tau_p);
        del_T_d = (-T_d_rec + H_d.*(T_max_d - T_d_rec)) * (dt/tau_d);
        T_p_rec = T_p_rec + del_T_p;       T_d_rec = T_d_rec + del_T_d;
        T_p_rec( T_p_rec > T_max_p) = T_max_p;
        T_d_rec( T_d_rec > T_max_d) = T_max_d;
        T_p__rec_mean(t) = mean(T_p_rec, "all");
        T_d_rec_mean(t)  = mean(T_d_rec, "all");
        if rew_vect(t) == 1 && reward_time_index > 1
            del_W_timer = eta_rec * (T_p_rec - T_d_rec) * (2*dt/delay_time);
            W_rec = W_rec + w_rec_identity .* del_W_timer';
            W_rec = W_rec .* (W_rec > 0) .* (1 - isnan(W_rec));
        end

        % Feed-forward
        r_pre  = r_input(:, t-D) .* (r_input(:, t-D) > 0.01);
        r_post = r_int  (:, t-1) .* (r_int  (:, t-1) > 0.01);
        H_d = eta_d * (r_pre * r_post');   H_p = eta_p * (r_pre * r_post');
        del_T_p = (-T_p_inp + H_p .* (T_max_p - T_p_inp)) * (dt/tau_p);
        del_T_d = (-T_d_inp + H_d .* (T_max_d - T_d_inp)) * (dt/tau_d);
        T_p_inp = T_p_inp + del_T_p;       T_d_inp = T_d_inp + del_T_d;
        T_p_inp( T_p_inp > T_max_p) = T_max_p;
        T_d_inp( T_d_inp > T_max_d) = T_max_d;
        T_p_ff_mean(t) = mean(T_p_inp, "all");
        T_d_ff_mean(t) = mean(T_d_inp, "all");
    end
end
end
