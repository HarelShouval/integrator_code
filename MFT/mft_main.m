%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Simulations vs. Mean-Field Comparison
%
%  • Runs “num_trials” independent simulations 
%    - Each trial lasts 10 s (dt defined in neuron_parameters.m).
%    - External drive is fixed (input_mean = 0.01) 
%  • After the loop, computes a mean-field approximation and plots:
%        1. Ensemble-averaged  s(t)  vs. mean-field  S(t)
%        2. Ensemble-averaged  r(t)  vs. mean-field  R(t)
%        3. Phase-plane dS/dt (source–sink curves)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close all;clear all;
% IO_function;  %  I/O 

%% ------------------------------------------------------------------------
neuron_parameters;            % loads dt, N, W, etc.
num_trials = 10;            
alpha = 1.03; % a small empirical number to correct MFT
%% Time parameters
t_total  = 10000;             % ms per trial
T_steps  = t_total / dt;

%% Pre-allocate output arrays
r_ex_all = zeros(num_trials, N, T_steps);
s_ex_all = zeros(num_trials, N, T_steps);
spk_all  = zeros(num_trials, N, T_steps);  % raster

%=========================================================================%
%                               TRIAL LOOP                                %
%=========================================================================%
for i = 1:num_trials
    %-------------------- state variables (per trial) --------------------%
    v_ex            = v_rest * ones(N,T_steps);
    input           = zeros(N,T_steps);          % external drive
    g_input_to_ex   = zeros(N,T_steps);
    g_ex_to_ex      = zeros(N,T_steps);
    is_refractory_ex= zeros(N,T_steps);
    is_spike_ex     = zeros(N,T_steps);
    r_ex            = zeros(N,T_steps);
    s_ex            = zeros(N,T_steps);

    input_mean = 0.01;                          % constant drive

    %----------------------------- TIME LOOP -----------------------------%
    for t = 2:T_steps

        % External input (DC)
        input(:,t) = input_mean;

        % Refractory handling
        is_refractory = is_refractory_ex(:,t)==1;
        v_ex(is_refractory,t) = v_rest;

        % Spike detection
        is_spike           = v_ex(:,t-1) >= v_th;
        is_spike_ex(:,t)   = is_spike;
        v_ex(is_spike,t)   = v_hold;
        if t < T_steps - t_refractory
            is_refractory_ex(is_spike,t+1:t+t_refractory) = 1;
        end

        % Synaptic gating
        s_ex(:,t) = s_ex(:,t-1) - (dt/tau_se)*s_ex(:,t-1) ...
                    + rho*is_spike.*(1 - s_ex(:,t-1));
        

        % Low-pass firing rate
        r_ex(:,t) = r_ex(:,t-1) + (is_spike/dt - r_ex(:,t-1))*(dt/tau_w);
        % s_ex(1:5,t) = 0;
        % r_ex(1:5,t) = 0;

        % Membrane update (non-refractory / non-spiking)
        free = ~is_refractory & ~is_spike;
        v_ex(free,t) = v_ex(free,t-1) + ...
            ( g_L*(E_l - v_ex(free,t-1)) + ...
              (g_input_to_ex(free,t-1)+g_ex_to_ex(free,t-1)).* ...
              (E_e - v_ex(free,t-1)) + ...
              norm_noise*randn(sum(free),1) ) * dt / C_m;

        % Conductances
        g_input_to_ex(:,t) = W_ff * input(:,t);
        g_ex_to_ex(:,t)    = W    * s_ex(:,t);

        % Save spikes
        spk_all(i,:,t) = is_spike;
    end %--------------------------- end time loop -----------------------%

    % Store trial results
    s_ex_all(i,:,:) = s_ex;
    r_ex_all(i,:,:) = r_ex;
end %============================= end trial loop =========================%

%% --------------------- Mean-field curve ---------------------------
runtime   = 10000;                % ms for MF sim
timeline  = 0:dt:runtime;
runlen    = numel(timeline);

S = zeros(1,runlen);              % MF synaptic activation
R = zeros(1,runlen);              % MF rate trace

W_MEAN = mean(W,'all') * N;
I_ext  = mean(W_ff,'all') * input_mean * N;

for it = 1:runlen-1
    ge        = S(it)*W_MEAN + I_ext;
    R(it)     = alpha * f(ge,p) * tanh(it/50);
    S(it+1)   = S(it) + dt*( R(it)*rho*(1-S(it)) - S(it)/tau_se );
end

%% -------------------------- Plotting -----------------------------------
line_colors = slanCL(210);

figure(8); clf; hold on
plot(timeline(1:10000), squeeze(mean(s_ex_all,2)), ...
     'Color',[.5 .5 .5],'LineWidth',2);
plot(timeline, S, 'Color', line_colors(8,:), 'LineWidth',4);
xlabel('Time (ms)'); ylabel('S'); xlim([0 1000]);
legend('Simulations','MFT','Box','off','Location','northwest');

figure(9); hold on;
subplot(2,1,1); hold on
plot(timeline(1:10000), squeeze(mean(s_ex_all,2)), ...
     'Color',[.5 .5 .5],'LineWidth',1);
plot(timeline, S, 'Color', line_colors(8,:), 'LineWidth',2);
ylabel('S'); xlim([0 2000]);

subplot(2,1,2); hold on
plot(timeline(1:10000), 1000*squeeze(mean(r_ex_all,2)), ...
     'Color',[.5 .5 .5],'LineWidth',1);
plot(timeline, 1000*R, 'Color', line_colors(8,:), 'LineWidth',2);
xlabel('Time (ms)'); ylabel('Firing rate (Hz)'); xlim([0 2000]);



figure(10); hold on;
subplot(2,1,1); hold on
plot(timeline(1:10000), squeeze(mean(s_ex_all,2)), ...
     'Color',[.5 .5 .5],'LineWidth',1);

plot(timeline(1:10000), mean(squeeze(mean(s_ex_all,2))), ...
     'Color',[.5 .5 .5]*0,'LineWidth',2);
% plot(timeline, S, 'Color', line_colors(8,:), 'LineWidth',2);
ylabel('S'); xlim([0 2000]);

subplot(2,1,2); hold on
plot(timeline(1:10000), 1000*squeeze(mean(r_ex_all,2)), ...
     'Color',[.5 .5 .5],'LineWidth',1);
plot(timeline(1:10000), 1000*mean(squeeze(mean(r_ex_all,2))), ...
     'Color',[.5 .5 .5]*0,'LineWidth',2);
% plot(timeline, 1000*R, 'Color', line_colors(8,:), 'LineWidth',2);
xlabel('Time (ms)'); ylabel('Firing rate (Hz)'); xlim([0 2000]);

%% ------------------ MFT source and sink curves ---------------------------------
plot_s = 1e-4:1e-4:1;
plot_ds = zeros(size(plot_s));
plot_decay = plot_s/tau_se;

for k = 1:numel(plot_s)
    R_tmp      = f(plot_s(k)*W_MEAN + I_ext, p);
    plot_ds(k) = alpha*R_tmp*rho*(1-plot_s(k));
end

figure(4); clf; hold on
plot(plot_s, plot_decay,'Color',[.2 .2 .2],'LineWidth',3);
plot(plot_s, plot_ds,'Color', line_colors(6,:),'LineWidth',3);
xlabel('S'); ylabel('dS/dt'); axis([0 0.5 0 6e-3]);
legend('Decay','Source','Box','off','Location','northwest');

%% ----------------------------------------- I-O Functions ------------

function f_val = f(g_e, beta)
    %  I-0 function fitted by 4th order polynoimail
    C_m = 0.2; g_L = 0.01; E_L = -60; E_E = -5; V_theta = -55; V_rest = -60; gc = 3e-3;
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
    %  Analytical I-0 function 
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
