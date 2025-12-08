%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Single-Neuron gE–rate Curve (W_rec = 0)
%
%  * Purpose
%    – Set recurrent weights to zero so that N identical neurons run in 
%      parallel with identical feed-forward drive; each trial therefore 
%      mimics one “isolated” cell driven by a constant excitatory conductance.
%    – Sweep the steady input level from -5 pS to +35 pS (1000 steps).
%    – Simulate 5 s of activity per trial, take the last 100 ms, and measure:
%         • mean excitatory conductance  g_E
%         • mean firing rate             r
%    – Scatter-plot the measured points and overlay the analytic f–I curve 
%      given by helper function f(…).
%
%  * Output
%      Figure 11: “Firing rate vs. Conductance”
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set W_rec = 0 to get I-0 curve. 
clear all; close all;

neuron_parameters;                 % neuron & network constants
W = W.*0;                          % wipe recurrent weights 
W_ff = 0.001 *ones(N,N); % all Wff are the same, when  W_rec = 0, equivelent to 
% run 100 identical single neurons at the same time.

%% Colormap
line_colors = [ 0.3059    0.4745    0.6549
    0.9490    0.5569    0.1686
    0.8824    0.3412    0.3490
    0.4627    0.7176    0.6980
    0.3490    0.6314    0.3098
    0.9294    0.7882    0.2824
    0.6902    0.4784    0.6314
    1.0000    0.6157    0.6549
    0.6118    0.4588    0.3725
    0.7294    0.6902    0.6745];

num_trials = 1000;


%% Metrics to record (per trial)
mean_ge_over_trials = zeros(num_trials,1);
mean_fr_over_trials = zeros(num_trials,1);
T_steps = 5000;
input           = zeros(N,T_steps);

%% Network connectivity




%% Trial loop
for i = 1:num_trials

    i
    % neuron_parameters;% regenerate W


    input(:,:) = ((i-1)/num_trials)*0.04 - 0.005;          
    % sweep input level. we allow negative (inhibitory) input to make sure
    % the neuron can fire at low fring rate.

    % State variables
    v_ex            = v_rest * ones(N,T_steps);
    
    g_input_to_ex   = zeros(N,T_steps);
    g_ex_to_ex      = zeros(N,T_steps);
    is_refractory_ex= zeros(N,T_steps);
    is_spike_ex     = zeros(N,T_steps);
    r_ex            = zeros(N,T_steps);
    s_ex            = zeros(N,T_steps);
    ds1             = zeros(N,T_steps);
    ds2             = zeros(N,T_steps);

    %% Time integration
    for t = 2:T_steps
       

        %% Refractory handling
        is_refractory = (is_refractory_ex(:,t)==1);
        v_ex(is_refractory,t) = v_rest;

        %% Spike detection
        is_spike             = v_ex(:,t-1) >= v_th;
        is_spike_ex(:,t)     = is_spike;
        v_ex(is_spike,t)     = v_hold;
        if t < T_steps - t_refractory
            is_refractory_ex(is_spike,t+1:t+t_refractory) = 1;
        end

        %% Synaptic gating (s dynamics)
        s_ex(:,t)  = s_ex(:,t-1) - (s_ex(:,t-1)*dt/tau_se) + rho*is_spike.*(1 - s_ex(:,t-1));
      
        %% Firing‑rate trace (low‑pass of spikes)
        r_ex(:,t)  = r_ex(:,t-1) + (is_spike/dt - r_ex(:,t-1)) * (dt/tau_w);

        %% Membrane potential update (non‑refractory & non‑spiking)
        not_ref = ~is_refractory & ~is_spike;
        v_ex(not_ref,t) = v_ex(not_ref,t-1) + ...
            (( g_L*(E_l - v_ex(not_ref,t-1)) + ...               % leak
               (g_input_to_ex(not_ref,t-1) + g_ex_to_ex(not_ref,t-1)).* ...
               (E_e - v_ex(not_ref,t-1)) + ...                   % synaptic
               norm_noise*randn(sum(not_ref),1) ) * dt / C_m );  % noise

        %% Conductances
        g_input_to_ex(:,t) = W_ff * input(:,t);
        g_ex_to_ex(:,t)    = W * s_ex(:,t);
    end

    %% Average over last 100 ms
    idx = T_steps-99:T_steps;
    ge_mat  = g_input_to_ex(:,idx) + g_ex_to_ex(:,idx);
    fr_mat  = 1000 * r_ex(:,idx);          % spikes/ms -> Hz

    mean_ge_over_trials(i) = mean(ge_mat,"all");
    mean_fr_over_trials(i) = mean(fr_mat,"all");

end

idx = mean_ge_over_trials >= 0;
mean_ge_over_trials = mean_ge_over_trials(idx);
mean_fr_over_trials = mean_fr_over_trials(idx);

%% Scatter plot (simulation)
figure(11);clf;hold on;
scatter(mean_ge_over_trials, mean_fr_over_trials, 20, "o", ...
        "MarkerEdgeColor", line_colors(1,:), "MarkerFaceColor", line_colors(1,:), ...
        "MarkerEdgeAlpha",0.5, "MarkerFaceAlpha",0.5);
xlabel('Mean g_e (last 100 ms)'); ylabel('Mean firing rate (Hz)');
title('Firing rate vs. Conductance across trials'); grid on;

%Theoretical curve using f(...)
ge_theory = linspace(0, max(mean_ge_over_trials)*1.2, 200);
fr_theory = zeros(size(ge_theory));
for k = 1:length(ge_theory)
    fr_theory(k) = 1000 * f(0, ge_theory(k), g_L, v_th, v_rest, C_m, E_l, E_e, ...
                            t_refractory, rho, tau_se, 0);
end


plot(ge_theory, fr_theory, 'Color',[0.5 0.5 0.5], 'LineWidth',2);
legend({'Simulation','Theory'}, 'Location','NorthWest');

%% Polynomial fit (order 4)
[p,~] = polyfit(mean_ge_over_trials, mean_fr_over_trials, 4);
ge_fit  = linspace(0, max(mean_ge_over_trials), 200);
fr_fit  = polyval(p, ge_fit);




%% Plot polynomial fit (order 4)
plot(ge_fit, fr_fit, 'LineWidth', 2, 'Color', line_colors(2,:));
legend({'Simulation','Theory','Poly fit (deg 4)'}, 'Location','NorthWest');

%% Goodness-of-fit metrics at observed x
x = mean_ge_over_trials(:);
y = mean_fr_over_trials(:);
yhat = polyval(p, x);

resid = y - yhat;
SSE = sum(resid.^2);
SST = sum( (y - mean(y)).^2 );
n   = numel(y);
k   = 5;                         % degree 4 polynomial => 5 coefficients

R2     = 1 - SSE/SST;
R2_adj = 1 - (SSE/(n - k)) / (SST/(n - 1));
MSE    = SSE / n;
RMSE   = sqrt(MSE);

% Annotate main figure
annotation_text = sprintf('R^2 = %.4f (adj %.4f)\\nRMSE = %.3f Hz', R2, R2_adj, RMSE);
text(0.05, 0.95, annotation_text, 'Units','normalized', ...
     'VerticalAlignment','top', 'BackgroundColor',[1 1 1 0.85], 'EdgeColor',[0.8 0.8 0.8]);

%% Residuals plot
figure(12); clf; hold on;
plot(x, resid, '.', 'MarkerSize', 10);
yline(0,'--');
grid on; box on;
xlabel('Mean g_e (last 100 ms)');
ylabel('Residual (Hz)');
title('Residuals of polynomial fit');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Theoritical f-I function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = f(S, I_ext, g_L, V_th, V_rest, C_m, E_L, E_E, ...
              t_ref, rho, tau_se, W)
    s_inf = (rho * S) / (rho * S + 1 / tau_se);
    g_inf = W * s_inf + I_ext;
    num   = g_inf*E_E + g_L*E_L - (g_inf + g_L)*V_rest;
    den   = g_inf*E_E + g_L*E_L - (g_inf + g_L)*V_th;
    if den<=0 || num<=0, r = 0; return; end
    phi = C_m/(g_inf + g_L) * log(num/den);
    if phi<0 || ~isreal(phi), r = 0; else, r = 1/(t_ref + 2 + phi); end 
    % actual t_ref = t_ref + 2 ms ( 2ms t_ref + dt before spike + dt to rest
    % V
 
 
end
