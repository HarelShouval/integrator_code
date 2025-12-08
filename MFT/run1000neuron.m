% clear all       
% close all        

%% Structure Parameters
% IO_function;                          % run IO_function n
% load IO_baserun.mat
%%
neuron_parameters;

N = 1000;                   % Neurons per population

w_rec = 0.000132/10;             % Recurrent connections in layer 5
w_rec_mean = 0.000132*1/4/10;
W = w_rec * rand(N, N) .* (rand(N, N) < sparsity);  % Excitatory to excitatory weight matrix
W = W.*w_rec_mean/mean(W,"all"); % make sure the mean is correct

num_trials = 5;                % trials to run

%% Time Parameters
t_total = 10000;               % trial length (ms)
T_steps = t_total / dt;        % time steps per trial

%% Input Parameters  (currently unused)
% W_input = (rand(N,1) < sparsity);   % binary input weights
% W_in   = W_in * 0.2;                % scale input weights

%% Preallocate result arrays
r_ex_all = zeros(num_trials, N, T_steps);
s_ex_all = zeros(num_trials, N, T_steps);

mean_wff = 0.0016*1/4*0.67/10;
W_in = (rand(N,N) < sparsity);   % new external weights

W_in = W_in/mean(W_in,"all")*mean_wff;

for i = 1:num_trials
    %% Trial‑level variables
    v_ex            = v_rest * ones(N, T_steps);   % membrane voltage
    input           = zeros(N, T_steps);           % external drive
    g_input_to_ex   = zeros(N, T_steps);           % input conductance
    g_ex_to_ex      = zeros(N, T_steps);           % recurrent conductance
    is_refractory_ex= zeros(N, T_steps);           % refractory flag
    is_spike_ex     = zeros(N, T_steps);           % spike flag

    r_ex = zeros(N, T_steps);      % firing rate
    s_ex = zeros(N, T_steps);      % synaptic activation
    ds1  = zeros(N, T_steps);      % decay term
    ds2  = zeros(N, T_steps);      % growth term

    %% External input settings
    input_spike = zeros(N, T_steps);
    input_mean  = 0.01;            % baseline rate

   

    %% Time loop
    for t = 2:T_steps
        %% External drive (constant)
        input(:, t) = input_mean;

        %% Refractory handling
        is_refractory       = is_refractory_ex(:, t) == 1;
        v_ex(is_refractory, t) = v_rest;            % hold at rest

        %% Spike detection
        is_spike        = v_ex(:, t-1) >= v_th;
        is_spike_ex(:,t)= is_spike;
        v_ex(is_spike, t) = v_hold;                 % after spike

        %% Mark refractory period
        if t < T_steps - t_refractory
            is_refractory_ex(is_spike, t+1:t+t_refractory) = 1;
        end

        %% Synaptic update
        s_ex(:, t) = s_ex(:, t-1) - (s_ex(:, t-1) * dt / tau_se) + ...
                     rho * is_spike .* (1 - s_ex(:, t-1));
        ds1(:, t)  =  s_ex(:, t-1) * dt / tau_se;           % decay part
        ds2(:, t)  =  rho * is_spike .* (1 - s_ex(:, t-1)); % growth part

        %% Rate update
        r_ex(:, t) = r_ex(:, t-1) + (is_spike / dt - r_ex(:, t-1)) * (dt / tau_w);

        %% Voltage update for active cells
        not_ref = ~is_refractory & ~is_spike;
        v_ex(not_ref, t) = v_ex(not_ref, t-1) + ...
            (g_L * (E_l - v_ex(not_ref, t-1)) + ...            % leak
            (g_input_to_ex(not_ref, t-1) + g_ex_to_ex(not_ref, t-1)) .* ...
            (E_e - v_ex(not_ref, t-1)) + ...                   % synaptic
            norm_noise * randn(sum(not_ref),1)) * dt / C_m;   % noise

        %% Conductances
        g_input_to_ex(:, t) = W_in * input(:, t);
        g_ex_to_ex(:, t)    = W * s_ex(:, t);
    end  % time loop

    %% Save results for this trial
    s_ex_all(i,:,:) = s_ex;
    r_ex_all(i,:,:) = r_ex;
end  % trial loop

%% Plotting and analysis (main code unchanged below) ----------------------


%% --- Figures -----------------------------------------------------------
line_colors = slanCL(210);

figure(8); clf; hold on
plot(timeline(1:10000), squeeze(mean(s_ex_all, 2)), 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
plot(timeline, S, 'Color', line_colors(8,:), 'LineWidth', 4);
xlabel('Time (ms)'); ylabel('S'); axis([0 1000 0 0.5]);
legend('Simulations','','','','','MFT', 'Location','northwest','Box','off');

figure(7); clf; hold on
for k = 1:5
    plot(timeline(1:10000), smooth(1000*squeeze(r_ex_all(2,k+10,1:10000)),200), 'LineWidth', 2);
end
plot(timeline(1:10000), 1000*squeeze(mean(r_ex_all(2,:,1:10000),2)), 'k','LineWidth',2);
xlabel('Time (ms)'); ylabel('Firing rate (Hz)'); axis([0 1000 0 50]);
legend('','','','','','Mean firing rate','Location','northwest','Box','off');

figure(6); clf; hold on
plot(timeline(1:10000), 1000*squeeze(mean(r_ex_all,2)), 'Color',[line_colors(8,:) 0.5],'LineWidth',1);
plot(timeline, 1000*R, 'Color', line_colors(8,:), 'LineStyle','--', 'LineWidth',4);
xlabel('Time (ms)'); ylabel('Firing rate (Hz)');

figure(9); clf; yyaxis left;
plot(timeline, S, 'Color',[0.3059 0.4745 0.6549],'LineWidth',4);
ylabel('Synaptic activation'); ylim([0 0.5]); xlim([0 1000]);
yyaxis right;
plot(timeline, 1000*R,'Color',[0.8824 0.3412 0.3490],'LineWidth',4);
ylabel('Firing rate (Hz)'); axis([0 1000 0 65]); xlabel('Time (ms)');
legend('MFT S','MFT f','Location','northwest');

%% Phase‑plane data
plot_s        = 0.0001:0.0001:1;
plot_ds       = zeros(size(plot_s));
plot_decay    = zeros(size(plot_s));
plot_ds_theory= zeros(size(plot_s));

for k = 1:length(plot_s)
    R_tmp        = 1.04*f(plot_s(k)*W_MEAN + I_ext, beta_fit);
     % R_tmp  = 1e-3*fr_fun(plot_s(k)*W_MEAN + I_ext);
    plot_ds(k)   = R_tmp * rho * (1 - plot_s(k));
    plot_decay(k)= plot_s(k) / tau_se;

    R_tmp        = ft(plot_s(k)*W_MEAN + I_ext, beta_fit);
    plot_ds_theory(k) = R_tmp * rho * (1 - plot_s(k));
end
%%
figure(4); clf; hold on; color = line_colors(6,:);
plot(plot_s, plot_decay, 'Color',[0.2 0.2 0.2],'LineWidth',3);
% scatter(plot_s2, plot_ds2, 20, color,'filled','MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.6);
plot(plot_s, plot_ds, 'Color', color,'LineWidth',3);
axis([0 0.5 0 6e-3]); xlabel('S'); ylabel('dS/dt');
legend('Decay','Source (data)','Source (sim)','Location','northwest','Box','off');

%% Time‑resolved CV and Fano (first 1 s)
window_ms      = 100;                              % window size
num_windows    = 1000 / window_ms;                 % 10 windows
samples_per_win= window_ms / dt;

cv_time   = zeros(num_trials, num_windows);
fano_time = zeros(num_windows, N);

time_axis = (0:num_windows-1) * window_ms;

for w = 1:num_windows
    t0 = (w-1)*samples_per_win + 1;
    t1 = w*samples_per_win;

    for tr = 1:num_trials
        r_seg = squeeze(r_ex_all(tr,:,t0:t1));      % [N x samples]
        cv_time(tr, w) = mean(std(r_seg,0,1));      % across neurons
    end

    spk_cnt = squeeze(sum(r_ex_all(:,:,t0:t1),3));  % [trials x N]
    fano_time(w,:) = var(spk_cnt,0,1) ./ (mean(spk_cnt,1) + eps);
end

mean_cv   = mean(cv_time,1);
mean_fano = mean(fano_time,2)';

figure(10); clf;
subplot(2,1,1); plot(time_axis, 1000*mean_cv,'-o','Color',[0.6118 0.4588 0.3725],'LineWidth',2);
ylabel('SD across neurons (Hz)'); set(gca,'XTick',[]);
subplot(2,1,2); plot(time_axis, mean_fano,'-o','Color',[0.6902 0.4784 0.6314],'LineWidth',2);
xlabel('Time (ms)'); ylabel('Fano Factor');

%% ------------------------------------------------ Functions ------------
function f_val = f(g_e, beta)
    % Piecewise firing‑rate approximation
    C_m = 0.2; g_L = 0.01; E_L = -60; E_E = -5; V_theta = -55; V_rest = -60; gc = 2.5e-3;
    if g_e < gc
        f_val = (beta(1)*g_e^4 + beta(2)*g_e^3 + beta(3)*g_e^2 + beta(4)*g_e + beta(5)) * 1e-3;
    else
        num = g_e*E_E + g_L*E_L - (g_e+g_L)*V_rest;
        den = g_e*E_E + g_L*E_L - (g_e+g_L)*V_theta;
        if den<=0 || num<=0
            f_val = 0;
        else
            phi = (C_m/(g_e+g_L))*log(num/den);
            f_val = (phi>0 && isreal(phi)) * 1/(4+phi);
        end
    end
end




function f_val = ft(g_e, beta)
    % Threshold branch only (for phase‑plane)
    C_m = 0.2; g_L = 0.01; E_L = -60; E_E = -5; V_theta = -55; V_rest = -60;
    num = g_e*E_E + g_L*E_L - (g_e+g_L)*V_rest;
    den = g_e*E_E + g_L*E_L - (g_e+g_L)*V_theta;
    if den<=0 || num<=0
        f_val = 0;
    else
        phi = (C_m/(g_e+g_L))*log(num/den);
        f_val = (phi>0 && isreal(phi)) * 1/(4+phi);
    end
end
