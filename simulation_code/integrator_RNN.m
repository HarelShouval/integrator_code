dt = 0.1;                    % Time step size (ms)/
% dt = 1;                    % Time step size (ms)/
N = 100;                     % Neurons per population
tau_w = 40;                  % Time window for firing rate dynamics (ms)
%% Membrane Dynamics Parameters
rho = 1 / 7;                 % Fractional change of synaptic activation with input spikes
tau_se = 80;
norm_noise = 0.13;           % sd for gaussion noise (per sqrt(ms))
norm_noise = norm_noise / sqrt(dt);   % Scale noise amplitude with the time step.

C_m = 0.2;                   % Membrane capacitance (µF/cm²)
g_L = 0.01;                  % Leak conductance (mS/cm²)
E_l = -60;                   % Leak reversal potential (mV)
E_e = -5;                    % Excitatory reversal potential (mV)
v_hold =  -61; 
v_th = -55;                  % Threshold potential for excitatory neurons (mV)
v_rest = -60;                % Resting potential (mV)
t_refractory = round((4-2*dt)/dt); % Refractory period duration (time steps)

%% Weight Matrix Parameters
sparsity = 0.5;            
w_rec = 0.000156;            % Wrec
w_rec_mean = w_rec*1/4*1;
W = w_rec * rand(N, N) .* (rand(N, N) < sparsity); 
W(1:N+1:end) = 0;       
if w_rec~= 0
    W = W.*w_rec_mean/mean(W,"all"); % make sure the mean is correct
else
    W = W.*0;
end

mean_wff = 0.00005*3
W_ff = (rand(N,N) < sparsity);   % new external weights
W_ff = W_ff/mean(W_ff,"all")*mean_wff;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Simulations (no MFT): plot 5 example trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_trials = 5;               % Number of example trials
%% Time parameters
t_total  = 1000;             % ms per trial
T_steps  = round(t_total / dt);  % robust to any dt

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

        % Synaptic gating (Euler)
        s_ex(:,t) = s_ex(:,t-1) - (dt/tau_se)*s_ex(:,t-1) ...
                    + rho*is_spike.*(1 - s_ex(:,t-1));
        
        % Low-pass firing rate (Euler)
        r_ex(:,t) = r_ex(:,t-1) + (is_spike/dt - r_ex(:,t-1))*(dt/tau_w);

        % Membrane update (non-refractory / non-spiking)
        free = ~is_refractory & ~is_spike;

        noiseterm = norm_noise*randn(sum(free),1);
        noiseterm(1:10) = 3*norm_noise*randn(10,1);

        v_ex(free,t) = v_ex(free,t-1) + ...
            ( g_L*(E_l - v_ex(free,t-1)) + ...
              (g_input_to_ex(free,t-1)+g_ex_to_ex(free,t-1)).* ...
              (E_e - v_ex(free,t-1)) + ...
              noiseterm ) * dt / C_m;  % dt-robust via scaled norm_noise above

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

%% -------------------------- Plotting: 5 example trials ------------------
figure(1); clf;

subplot(2,1,1); hold on;
for i = 1:num_trials
    plot((0:T_steps-1)*dt, squeeze(mean(s_ex_all(i,:,:),2)), 'LineWidth', 1.2);
end
ylabel('Mean S');
title('Five example trials');
xlim([0 t_total]);
box off;

subplot(2,1,2); hold on;
for i = 1:num_trials
    plot((0:T_steps-1)*dt, 1000*squeeze(mean(r_ex_all(i,:,:),2)), 'LineWidth', 1.2);
end
xlabel('Time (ms)'); ylabel('Mean firing rate (Hz)');
xlim([0 t_total]);
box off;


