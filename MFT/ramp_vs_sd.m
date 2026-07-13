%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decision-to-timer weight sweep
%
% This script simulates a decision network coupled to a timer network.
% Across weight conditions, it measures the first time at which the timer
% population activity crosses a fixed threshold. It then plots the standard
% deviation of crossing time against the mean crossing time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;  close all;  clc;

%% Simulation setup
n_cond     = 20;
num_trials = 20;

mean_cross = nan(1, n_cond);
sd_cross   = nan(1, n_cond);
timer_curve_store = cell(1, n_cond);

fprintf('Running %d weight conditions × %d trials …\n', n_cond, num_trials);

%% Model parameters
N_timer   = 100;
npp       = 50;
n_exc_total = 2*npp;
n_inh_total = 2*npp;
N_input   = 100;

dt        = 1;
T         = 4000/dt;
delta     = 600/dt;
t_total   = 3001;
t_steps   = t_total/dt;
D         = 10/dt;
tau_w     = 50;

eG = 0.4;
iG = 0.5;

%% Stimulus input to decision network
W_in = cell(2,1);
W_in{1} = 0.01*rand(n_exc_total,1);
W_in{2} = 0.01*rand(n_exc_total,1);
W_in{1}(1:npp)             = eG;
W_in{2}(npp+1:2*npp)       = eG;

p_r  = 0.02*dt;

%% Decision-network connectivity
rho     = 0.5;
w_rec   = 0.002;
i_mean  = 0.0015;
p_mean  = 0.0015;
k_rec   = 0.00;

W_ji = zeros(n_exc_total);
M_ki = zeros(n_exc_total, n_inh_total);
P_ik = zeros(n_inh_total, n_exc_total);
K_jk = zeros(n_inh_total);

W_ji(1:npp,1:npp)                 = w_rec*(rand(npp)<rho);
W_ji(npp+1:end,npp+1:end)         = w_rec*(rand(npp)<rho);

M_ki(1:npp,npp+1:end)             = i_mean*(rand(npp)<rho);
M_ki(npp+1:end,1:npp)             = i_mean*(rand(npp)<rho);

P_ik(1:npp,1:npp)                 = p_mean*(rand(npp)<rho);
P_ik(npp+1:end,npp+1:end)         = p_mean*(rand(npp)<rho);

K_jk(1:npp,1:npp)                 = k_rec*(rand(npp)<rho);
K_jk(npp+1:end,npp+1:end)         = k_rec*(rand(npp)<rho);

rec_identity_decision = zeros(n_exc_total);
rec_identity_decision(1:npp,1:npp)                 = W_ji(1:npp,1:npp)>0;
rec_identity_decision(npp+1:end,npp+1:end)         = W_ji(npp+1:end,npp+1:end)>0;

%% Timer-network connectivity
w_timer_mean = 0.000158/4;  
W_timer = w_timer_mean*rand(N_timer,N_timer).*(rand(N_timer,N_timer)<rho);
W_timer = W_timer./mean(W_timer,"all")*w_timer_mean;

rec_identity_timer = W_timer>0;

%% Input-to-decision weights
w_input_to_decision_mean = 0.001;
W_input_to_decision = zeros(n_exc_total,N_input);
W_input_to_decision(1:npp,1:N_input/2)                     = w_input_to_decision_mean*rand(npp,N_input/2);
W_input_to_decision(npp+1:end,N_input/2+1:end)             = w_input_to_decision_mean*rand(npp,N_input/2);

%% Membrane and synaptic parameters
rho_syn  = 1/7;
tau_se_decision_input = 10;
tau_se_timer = 80;
tau_si  = 10;
norm_noise = 0.13;
C_m      = 0.2;
g_L      = 0.01;
E_i      = -70;
E_l      = -60;
E_e      = -5;
v_th     = -55;
v_th_i   = -50;
v_rest   = -60;
v_hold   = -61;
t_refractory = 2/dt;

%% Timer threshold and initial decision-to-timer weights
ramp_threshold = 32;
 w0 = 0.00008;

W_decision_to_timer      = w0*rand(N_timer,n_exc_total).*(rand(N_timer,n_exc_total)<rho);
W_decision_to_timer = W_decision_to_timer./mean(W_decision_to_timer,"all")*w0*rho*0.5;
W_decision_to_timer_id   = W_decision_to_timer>0;  %#ok<NASGU>

%% Weight sweep
for cond = 1:n_cond
   W_decision_to_timer = 0.95*W_decision_to_timer;

    crossTimes      = nan(num_trials,1);
    timer_FR_trials = nan(num_trials, t_steps);

    %% Trial loop
    for l = 1:num_trials

        %% State variables
        s_input_decision = zeros(n_exc_total, t_steps);
        r_input_decision = zeros(n_exc_total, t_steps);

        v_ex           = v_rest  * ones(n_exc_total, t_steps);
        r_ex           = zeros(n_exc_total, t_steps);
        s_ex           = zeros(n_exc_total, t_steps);
        is_refractory_ex = zeros(n_exc_total, t_steps);
        is_spike_ex      = zeros(n_exc_total, t_steps);

        v_inh          = v_rest * ones(n_inh_total, t_steps);
        r_inh          = zeros(n_inh_total, t_steps);
        s_inh          = zeros(n_inh_total, t_steps);
        is_refractory_inh = zeros(n_inh_total, t_steps);
        is_spike_inh      = zeros(n_inh_total, t_steps);

        v_timer        = v_rest * ones(N_timer, t_steps);
        r_timer        = zeros(N_timer, t_steps);
        s_timer        = zeros(N_timer, t_steps);
        g_decision_to_timer = zeros(N_timer, t_steps);
        g_timer_to_timer    = zeros(N_timer, t_steps);
        is_refractory_timer = zeros(N_timer, t_steps);
        is_spike_timer      = zeros(N_timer, t_steps);

        choice = 1;

        %% Time loop
        for t = 2:t_steps

            %% Input layer
            input_spike_decision = poissrnd(p_r, n_exc_total, 1);
            if t>1 && t<500
                is_spike_input_decision = (input_spike_decision==1);
            else
                is_spike_input_decision = 0;
            end
            s_input_decision(:,t) = s_input_decision(:,t-1) ...
                   - (s_input_decision(:,t-1)*dt/tau_se_decision_input) ...
                   + rho_syn * W_in{choice} .* is_spike_input_decision .* (1-s_input_decision(:,t-1));
            r_input_decision(:,t) = r_input_decision(:,t-1) ...
                   + (W_in{choice}.*is_spike_input_decision/dt - r_input_decision(:,t-1))*(dt/tau_w);

            %% Decision excitatory population
            is_refrac = is_refractory_ex(:,t)==1;
            v_ex(is_refrac,t) = v_rest;
            is_spike = (v_ex(:,t-1) >= v_th);
            is_spike_ex(:,t) = is_spike;
            not_spike = ~is_spike & ~is_refrac;
            v_ex(is_spike,t) = v_hold;
            if t < t_steps - t_refractory
                is_refractory_ex(is_spike, t+1:t+t_refractory) = 1;
            end
            s_ex(:,t) = s_ex(:,t-1) - (s_ex(:,t-1)*dt/tau_se_decision_input) ...
                       + rho_syn * is_spike .* (1-s_ex(:,t-1));
            v_ex(not_spike,t) = v_ex(not_spike,t-1) + ...
               ((randn(sum(not_spike),1)*norm_noise) ...
               + g_L*(E_l - v_ex(not_spike,t-1)) ...
               + (W_ji(not_spike,:) * s_ex(:,t-1)) .* (E_e - v_ex(not_spike,t-1)) ...
               + (M_ki(not_spike,:) * s_inh(:,t-1)) .* (E_i - v_ex(not_spike,t-1)) ...
               + W_input_to_decision(not_spike,:) * s_input_decision(:,t-1) .* (E_e - v_ex(not_spike,t-1)) ...
               )*(dt/C_m);
            r_ex(:,t) = r_ex(:,t-1) + (is_spike/dt - r_ex(:,t-1))*(dt/tau_w);

            %% Decision inhibitory population
            is_refrac = is_refractory_inh(:,t)==1;
            v_inh(is_refrac,t) = v_rest;
            is_spike = (v_inh(:,t-1) >= v_th_i);
            is_spike_inh(:,t) = is_spike;
            not_spike = ~is_spike & ~is_refrac;
            v_inh(is_spike,t) = v_hold;
            if t < t_steps - t_refractory
                is_refractory_inh(is_spike, t+1:t+t_refractory) = 1;
            end
            s_inh(:,t) = s_inh(:,t-1) - (s_inh(:,t-1)*dt/tau_si) ...
                       + rho_syn * is_spike .* (1-s_inh(:,t-1));
            v_inh(not_spike,t) = v_inh(not_spike,t-1) + ...
               ((randn(sum(not_spike),1)*norm_noise) ...
               + g_L*(E_l - v_inh(not_spike,t-1)) ...
               + (P_ik(not_spike,:) * s_ex(:,t-1)) .* (E_e - v_inh(not_spike,t-1)) ...
               + (K_jk(not_spike,:) * s_inh(:,t-1)) .* (E_i - v_inh(not_spike,t-1)) ...
               )*(dt/C_m);
            r_inh(:,t) = r_inh(:,t-1) + (is_spike/dt - r_inh(:,t-1))*(dt/tau_w);

            %% Timer network
            g_decision_to_timer(:,t) = W_decision_to_timer * s_ex(:,t-1);
            g_timer_to_timer(:,t)    = W_timer          * s_timer(:,t-1);
            g_total_timer = g_decision_to_timer(:,t) + g_timer_to_timer(:,t);

            is_refrac = is_refractory_timer(:,t)==1;
            v_timer(is_refrac,t) = v_rest;
            is_spike = (v_timer(:,t-1) >= v_th);
            is_spike_timer(:,t) = is_spike;
            not_spike = ~is_spike & ~is_refrac;
            v_timer(is_spike,t) = v_hold;
            if t < t_steps - t_refractory
                is_refractory_timer(is_spike,t+1:t+t_refractory) = 1;
            end
            s_timer(:,t) = s_timer(:,t-1) - (s_timer(:,t-1)*dt/tau_se_timer) ...
                         + rho_syn*is_spike .* (1-s_timer(:,t-1));
            v_timer(not_spike,t) = v_timer(not_spike,t-1) + ...
               ((randn(sum(not_spike),1)*norm_noise) ...
               + g_L*(E_l - v_timer(not_spike,t-1)) ...
               + g_total_timer(not_spike) .* (E_e - v_timer(not_spike,t-1)) ...
               )*(dt/C_m);
            r_timer(:,t) = r_timer(:,t-1) + (is_spike/dt - r_timer(:,t-1))*(dt/tau_w);

        end

        %% First threshold-crossing time
        mean_activity_timer = mean(1000*r_timer(:,1:t_total),1);
        idx_cross = find(mean_activity_timer > ramp_threshold, 1, 'first');
        if ~isempty(idx_cross)
            crossTimes(l) = idx_cross*dt;
        end
        timer_FR_trials(l,:) = mean_activity_timer;

    end

    %% Crossing-time statistics
    mean_cross(cond) = mean(crossTimes,'omitnan');
    sd_cross(cond)   = std (crossTimes,'omitnan');
    timer_curve_store{cond} = timer_FR_trials;

    fprintf('  cond %02d  w=%.5f   mean=%.0f ms  SD=%.0f ms\n', ...
            cond, w0, mean_cross(cond), sd_cross(cond));
end

%% Summary plot with zero-intercept linear fit
figure('Name','mean_vs_sd','Color','w');
plot(mean_cross, sd_cross, 'o','LineWidth',2);
hold on;

fit_idx = isfinite(mean_cross) & isfinite(sd_cross);
x_fit = mean_cross(fit_idx);
y_fit = sd_cross(fit_idx);

k = (x_fit * y_fit') / (x_fit * x_fit');

x_line = linspace(min(x_fit), max(x_fit), 200);
y_line = k * x_line;

plot(x_line, y_line, 'r-', 'LineWidth', 2);

xlabel('Mean ramp-cross time  (ms)');
ylabel('SD of ramp-cross time (ms)');
title('Timing variability vs mean across W_{decision→timer}');
legend({'Simulation', sprintf('Fit: y = %.3gx', k)}, 'Location', 'best');
grid on;