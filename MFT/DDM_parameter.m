%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Simulate 100 trials with the same parameters, estimate μ and σ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;


% ------------------------------------------------------------------------
neuron_parameters;   % <— your file that defines dt, N, W, etc.
% ------------------------------------------------------------------------
N_input = 100;

%% ---------------- user‑set simulation parameters -----------------------
NUM_TRIALS     = 100;     % <- 100 trials
input_mean     = 0.01;    % DC drive amplitude (if you keep it)
RAMP_THRESHOLD = 25;      % Hz – population mean FR threshold (θ)
alpha          = 1.03;    % (kept from your code; unused here)
t_total        = 3000;    % ms per trial
% ------------------------------------------------------------------------

T_steps  = t_total / dt;
time_vec = (0:T_steps-1)*dt;
dt_window = 50;
% For μ, σ estimation
all_dX      = [];         % collect all ΔX before the hit across trials
all_dX_dt   = [];         % (optional) ΔX/dt if you want it directly
X0_all      = nan(NUM_TRIALS,1);
theta_all   = nan(NUM_TRIALS,1);

% For cross-time stats
ramp_store  = nan(NUM_TRIALS, T_steps);
crossTimes  = nan(NUM_TRIALS, 1);

% Small margin to stay away from boundaries when estimating μ, σ
eps_low   = 0.5;          % Hz above 0 (reflecting bound)
eps_high  = 0.5;          % Hz below θ

fprintf('Running %d trials …\n', NUM_TRIALS);

for i = 1:NUM_TRIALS
i
    % ---------- re-load / re-initialize per-trial state -----------------
    neuron_parameters;

    % (re)define your external weights if you need them unchanged per trial
    mean_wff = 0.001/10;
    W_ff     = (rand(N,N_input) < sparsity);
    W_ff     = W_ff/mean(W_ff,"all")*mean_wff;

    v_ex            = v_rest * ones(N,T_steps);
    g_input_to_ex   = zeros(N,T_steps);
    g_ex_to_ex      = zeros(N,T_steps);
    is_refractory_ex= zeros(N,T_steps);
    is_spike_ex     = zeros(N,T_steps);
    r_ex            = zeros(N,T_steps);
    s_ex            = zeros(N,T_steps);
    s_input         = zeros(N_input,T_steps);

    % ----------------------------- TIME LOOP -----------------------------
    for t = 2:T_steps
        % Example: Poisson external spikes (replace with your constant drive if you want)
        is_spike_input = rand(N_input,1) < 0.005;
        s_input(:, t) = s_input(:, t-1) - (s_input(:, t-1)*dt/20) + ...
                        rho * is_spike_input .* (1 - s_input(:, t-1));
        % s_input(:, t) = 0.02;

        %---------------- Refractory handling ----------------------------%
        is_refrac = is_refractory_ex(:,t)==1;
        v_ex(is_refrac,t) = v_rest;

        %---------------- Spike detection --------------------------------%
        is_spike           = v_ex(:,t-1) >= v_th;
        is_spike_ex(:,t)   = is_spike;
        v_ex(is_spike,t)   = v_hold;
        if t < T_steps - t_refractory
            is_refractory_ex(is_spike,t+1:t+t_refractory) = 1;
        end

        %---------------- Synaptic gating --------------------------------%
        s_ex(:,t) = s_ex(:,t-1) - (dt/tau_se)*s_ex(:,t-1) ...
                    + rho*is_spike.*(1 - s_ex(:,t-1));

        %---------------- Low‑pass firing rate ---------------------------%
        r_ex(:,t) = r_ex(:,t-1) + (is_spike/dt - r_ex(:,t-1))*(dt/tau_w);

        %---------------- Membrane update (non‑refrac / non‑spike) -------%
        free = ~is_refrac & ~is_spike;
        v_ex(free,t) = v_ex(free,t-1) + ...
            ( g_L*(E_l - v_ex(free,t-1)) + ...
              (g_input_to_ex(free,t-1)+g_ex_to_ex(free,t-1)).* ...
              (E_e - v_ex(free,t-1)) + ...
              norm_noise*randn(sum(free),1) ) * dt / C_m;

        %---------------- Conductances -----------------------------------%
        g_input_to_ex(:,t) = W_ff * s_input(:,t);
        g_ex_to_ex(:,t)    = W    * s_ex(:,t);
    end
    % --------------------------- end time loop ---------------------------

    % Decision variable: population mean FR (Hz)
    mean_FR          = 1000*mean(r_ex(:,1:T_steps),1);   % Hz
    ramp_store(i,:)  = mean_FR;
    X                = mean_FR;                       % column

    % Initial value
    X0_all(i) = X(1);

    % First‑passage to θ
    idx = find(X > RAMP_THRESHOLD, 1, 'first');
    if ~isempty(idx)
        crossTimes(i) = idx * dt;                        % ms
        theta_all(i)  = X(idx);                          % empirical value at hit
    else
        theta_all(i)  = NaN;
    end

    % --------------- Collect increments for μ, σ estimation -------------
    % Use only pre-hit samples; if never hit, use first 500 points
    if ~isempty(idx)
        idx_use = 1:(idx-1);
    else
        idx_use = 1:500;
    end

    % keep away from the bounds
    interior = X(idx_use) > eps_low & (RAMP_THRESHOLD - X(idx_use)) > eps_high;

    idx_use  = idx_use(interior);
    if numel(idx_use) >= 2
        data = X(idx_use);
        len = length(data);
        bin_length = floor(len/dt_window);
        data = reshape(data(1:dt_window*bin_length), [dt_window, bin_length ]);
        data = mean(data);
        dX  = diff(data);
        all_dX = [all_dX; dX(:)];
       
    end
end

%% ---------------- μ & σ estimation from increments -----------------------

mu_hat     = mean(all_dX) / dt_window;
sigma_hat  = sqrt( var(all_dX - mu_hat*dt_window) / dt_window );

fprintf('\nEstimated parameters (from increments):\n');
fprintf('  mu     = %.4f (units of X/s)\n', mu_hat);
fprintf('  sigma  = %.4f (units of X/s^{1/2})\n', sigma_hat);
fprintf('  X0     = %.4f ± %.4f  (mean ± SD over trials)\n', ...
        mean(X0_all,'omitnan'), std(X0_all,'omitnan'));
fprintf('  theta  = %.4f ± %.4f  (empirical at hit)\n', ...
        mean(theta_all,'omitnan'), std(theta_all,'omitnan'));

% ------------------- quick diagnostics ----------------------------------
% 1) Normality of residual increments
residuals = (all_dX - mu_hat*dt) / sqrt(dt);
figure('Color','w'); histogram(residuals, 40, 'Normalization','pdf');
title('Residual increments (normalized)'); xlabel('z'); ylabel('pdf');

% 2) Autocorrelation of residuals (should be ~white)
try
    figure('Color','w'); autocorr(residuals, 30);
    title('ACF of residuals');
catch
    % If you don't have Econometrics Toolbox, just skip
end

% 3) Cross-time summary
figure('Color','w'); histogram(crossTimes, 30);
xlabel('First‑passage time (ms)'); ylabel('Count');
title('Distribution of ramp‑cross times');

plot(ramp_store(1:10,:)');