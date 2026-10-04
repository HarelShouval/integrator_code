
dt = 1;                      % Time step size (ms)
N = 100;                   % Neurons per population
tau_w = 40;                % Time window for firing rate dynamics (ms)
%% Membrane Dynamics Parameters
rho = 1 / 7;                 % Fractional change of synaptic activation with input spikes
tau_se = 80;
norm_noise = 0.13;           % sd for gaussion noise


C_m = 0.2;                   % Membrane capacitance (µF/cm²)
g_L = 0.01;                  % Leak conductance (mS/cm²)
E_l = -60;                   % Leak reversal potential (mV)
E_e = -5;                    % Excitatory reversal potential (mV)
v_hold =  -61; 
v_th = -55;                  % Threshold potential for excitatory neurons (mV)
v_rest = -60;                % Resting potential (mV)
t_refractory = 2 / dt;       % Refractory period duration (time steps)

%% Weight Matrix Parameters
sparsity = 0.5;            
w_rec = 0.000158;              % Wrec
w_rec_mean = w_rec*1/4*1;
W = w_rec * rand(N, N) .* (rand(N, N) < sparsity); 
W(1:N+1:end) = 0;       
if w_rec~= 0
    W = W.*w_rec_mean/mean(W,"all"); % make sure the mean is correct
else
    W = W.*0;
end


mean_wff = 0.00007*3*1.0;
W_ff = (rand(N,N) < sparsity);   % new external weights
W_ff = W_ff/mean(W_ff,"all")*mean_wff;

% % %% Weight Matrix Parameters  – Gaussian distribution, same mean
% sparsity    = 1;                    % connection probability
% w_rec       = 0.000158;               % reference scale
% w_rec_mean  = w_rec * 1/4;            % desired mean of W (incl. zeros)
% 
% % --- Connection mask ----------------------------------------------------
% mask = rand(N,N) < sparsity;          % connectivity pattern
% mask(1:N+1:end) = 0;                  % remove self-connections
% 
% % --- Draw Gaussian weights  --------------------------------------------
% mu_w = w_rec_mean;                    % mean of non-zero entries before rescale
% sigma_w = mu_w/4;                     % Standard deviation
% 
% W = mu_w + sigma_w * randn(N,N);      % normal N(μ, σ²)
% W(W<0) = 0;                     % ensure excitatory (>0)
% W = W .* mask;                        % apply sparsity mask
% 
% % --- Rescale so ⟨W⟩ (including zeros) = w_rec_mean ----------------------
% W = W * (w_rec_mean / mean(W,"all"));
% 
% 
% %% Weight Matrix Parameters  – Gamma distribution with heavy tail
% sparsity    = 0.2;                    % connection probability
% w_rec       = 0.000158;               % reference scale
% w_rec_mean  = w_rec * 1/4;            % desired mean of W (incl. zeros)
% 
% % --- Connection mask ----------------------------------------------------
% mask = rand(N,N) < sparsity;          % sparsity pattern
% mask(1:N+1:end) = 0;                  % remove self-connections
% 
% % --- Draw heavy-tailed Gamma weights -----------------------------------
% %   shape < 1  →  heavy right tail
% gamma_k     = 0.5;                    % shape  (heavier tail as k → 0)
% gamma_theta = 1;                      % scale  (arbitrary – we rescale later)
% 
% W = gamrnd(gamma_k, gamma_theta, N, N);
% W = W .* mask;                        % apply sparsity mask  (zeros elsewhere)
% 
% % --- Rescale so ⟨W⟩ (including zeros) = w_rec_mean ----------------------
% W = W * (w_rec_mean / mean(W,"all")); % final mean exactly w_rec_mean
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



