%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Simulations + Ramp‑Cross Statistics with a Constant 0–500 ms Input
%
%  ── What this script does ───────────────────────────────────────────────
%  • Runs NUM_TRIALS independent network simulations (10 s each).
%  • Each trial receives a constant external drive (input_mean) for the
%      first STIM_DURATION ms, then 0 afterward.          [same as code 1]
%  • After every trial the population‑mean firing‑rate trace is examined:
%        – The first time it exceeds RAMP_THRESHOLD (34 Hz) is recorded.
%  • Prints mean ± SD of those crossing times and plots the single‑trial
%    ramp traces (grey) together with the threshold (red dashed line).
%  • All of the original mean‑field calculations and phase‑plane plots
%    from code 2 are preserved so you can still compare to MFT.
%
%  ── Requirements ────────────────────────────────────────────────────────
%  • neuron_parameters_taus20.m  (same as before; defines dt, N, W, etc.)
%  • slanCL.m  (optional – for nicer colours; comment out if not present)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;

% ------------------------------------------------------------------------
neuron_parameters;      % loads dt, N, W, rho, tau_se, etc.
% ------------------------------------------------------------------------

%% ---------------- user‑set simulation parameters -----------------------
NUM_TRIALS     = 1000;     % trials (code 1 used 20)
input_mean     = 0.01;   % amplitude of DC drive (same as old code 2)
STIM_DURATION  = 500;    % ms – constant drive from 0…500 ms
RAMP_THRESHOLD = 34;     % Hz – population mean FR threshold
alpha          = 1.03;   % empirical factor for mean‑field
% ------------------------------------------------------------------------

%% Time parameters
t_total  = 4000;                     % ms per trial (10 s)
T_steps  = t_total / dt;
time_vec = (0:T_steps-1)*dt;          % convenient time axis


n_cond = 1;
%% Pre‑allocate output arrays
% r_ex_all   = zeros(NUM_TRIALS, N, T_steps);   % low‑pass FRs (excitatory)
% s_ex_all   = zeros(NUM_TRIALS, N, T_steps);   % synaptic gating
% spk_all    = zeros(NUM_TRIALS, N, T_steps);   % raster
ramp_store = nan(NUM_TRIALS, T_steps);        % population mean FR
crossTimes = nan(NUM_TRIALS, 1);              % 1st threshold crossing (ms)
s_input = zeros(N, T_steps);


%=========================================================================%
%                               TRIAL LOOP                                %
%=========================================================================%
fprintf('Running %d trials …\n', NUM_TRIALS);
for ic = 1:n_cond
for i = 1:NUM_TRIALS
    %---------------- state variables (per trial) ------------------------%
    neuron_parameters;
    mean_wff = 0.0001*(0.8 + 2/20*1.5);
    W_ff = (rand(N,N) < sparsity);   % new external weights
    W_ff = W_ff/mean(W_ff,"all")*mean_wff;


    v_ex            = v_rest * ones(N,T_steps);
    input           = zeros(N,T_steps);          
    g_input_to_ex   = zeros(N,T_steps);
    g_ex_to_ex      = zeros(N,T_steps);
    is_refractory_ex= zeros(N,T_steps);
    is_spike_ex     = zeros(N,T_steps);
    r_ex            = zeros(N,T_steps);
    s_ex            = zeros(N,T_steps);

    %----------------------------- TIME LOOP -----------------------------%
    for t = 2:T_steps
        % ----------- Constant drive for 0–STIM_DURATION ms --------------
       

        is_spike_input = rand(N,1) < 0.004;
   

         
        s_input(:, t) = s_input(:, t-1) - (s_input(:, t-1)*dt/20) + ...
                         rho * is_spike_input .* (1 - s_input(:, t-1));
        % s_input(:, t) = input_mean;
      
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

        %---------------- Save spikes (optional raster) ------------------%
        % spk_all(i,:,t) = is_spike;
      
    end %--------------------------- end time loop -----------------------%

    %---------------- After trial: ramp‑cross statistics -----------------%
    mean_FR          = 1000*mean(r_ex(:,1:T_steps),1);   % Hz
    ramp_store(i,:)  = mean_FR;                          % keep full trace
    idx              = find(mean_FR > RAMP_THRESHOLD,1,'first');
    if ~isempty(idx)
        crossTimes(i) = idx*dt;                          % ms
    end

    %---------------- Save trial results (original arrays) --------------%

    % fprintf('  Trial %02d  cross = %.0f ms\n', i, crossTimes(i));
end %============================= end trial loop =========================%

%% -------------------- Summary statistics & quick plot -------------------
mean_cross(ic) = mean(crossTimes,'omitnan');
sd_cross(ic)   = std (crossTimes,'omitnan');
fprintf('\nSummary: mean = %.0f ms   SD = %.0f ms  (n = %d trials)\n', ...
        mean_cross(ic), sd_cross(ic), NUM_TRIALS);


end

figure(1);clf;
plot(mean_cross,sd_cross)

%% --------------------- Mean‑field curve (unchanged) --------------------
%% ----------------------  summary plot  ---------------------------------
figure(1); clf; hold on;


fit_idx = 2: min(19, numel(mean_cross));
x_fit_pts = mean_cross(fit_idx)-50;
y_fit_pts = sd_cross(fit_idx);

% 1) scatter the data as dots
scatter(x_fit_pts, y_fit_pts, 60, 'filled', ...
        'MarkerFaceColor',[0 .4 .8], 'MarkerEdgeColor','k');

xlabel('Mean ramp‑cross time  (ms)');
ylabel('SD of ramp‑cross time (ms)');
title('Timing variability vs mean across W_{decision→timer}');
% grid on;

% -- choose the points to fit: indices 3…18 (or as many as exist)


% -- slope for y = kx (least‑squares with intercept = 0)
% k = (x_fit_pts * y_fit_pts') / (x_fit_pts * x_fit_pts');
% 
% % -- line to plot
% x_line = linspace(min(x_fit_pts), max(x_fit_pts), 100);
% y_line = k * x_line;
% plot(x_line, y_line, 'r-', 'LineWidth', 2);
% 



x_pts     = mean_cross(fit_idx);
y_pts     = sd_cross(fit_idx);

% 3) least‑squares slope for y = k·x^{3/2}
x_pow     = x_pts .^ 1.5;                 % x^{3/2}
k         = (x_pow * y_pts') / (x_pow * x_pow');   % closed‑form LS solution

% 4) plot the fitted curve
x_line = linspace(min(x_pts), max(x_pts), 200);
y_line = k * x_line .^ 1.5;
plot(x_line, y_line, 'r-', 'LineWidth', 2);

%%
% -------------------------------------------------------------------------
%  Fit crossTimes with an inverse‑Gaussian (Wald) and plot
% -------------------------------------------------------------------------
%   Requires the Statistics & Machine Learning Toolbox (fitdist).
%   crossTimes :  Nx1 (or 1xN) vector of first‑passage times (ms)
% -------------------------------------------------------------------------

% Example: crossTimes already in workspace
% crossTimes = [...];        % your data

% 1) Fit inverse‑Gaussian  p(T | μ, λ)
pd = fitdist(crossTimes(:), 'InverseGaussian');   % returns μ (mu) and λ (lambda)
mu_hat  = pd.mu;          % mean of the IG
lam_hat = pd.lambda;      % shape parameter

fprintf('Inverse‑Gaussian fit:  μ = %.1f ms   λ = %.2e\n', mu_hat, lam_hat);

% 2) Plot histogram (PDF normalised)
figure('Color','w');  hold on
nbins = max(10, round(sqrt(numel(crossTimes))));  % Freedman–Diaconis rule-ish
histogram(crossTimes, 20, 'Normalization', 'pdf', ...
          'FaceColor', [0.25 0.45 0.75], 'EdgeColor', 'none');

% 3) Overlay fitted PDF
t_grid = linspace(min(crossTimes), max(crossTimes), 400);
pdf_fit = pdf(pd, t_grid);     % inverse‑Gaussian pdf evaluated at t_grid
plot(t_grid, pdf_fit, 'r', 'LineWidth', 2);

% 4) Labels
xlabel('First‑passage time  (ms)', 'FontSize', 11);
ylabel('Probability density',      'FontSize', 11);
title('Inverse‑Gaussian fit to ramp first‑passage times');
legend({'Data (histogram)', 'Fitted inverse‑Gaussian'}, 'Location','northeast');

box off
