
dt = 1;                      % Time step size (ms)
N = 100;                   % Neurons per population
tau_w = 40;                % Time window for firing rate dynamics (ms)
%% Membrane Dynamics Parameters
rho = 1 / 7;                 % Fractional change of synaptic activation with input spikes
tau_se = 80;                 % Time constant for excitatory synaptic activation (ms)
tau_si = 10;
norm_noise = 0.13;           % sd for gaussion noise
 % norm_noise = 0.08;   

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
% w_rec = 1.21e-04   ;
w_rec_mean = w_rec*1/4;
W = w_rec * rand(N, N) .* (rand(N, N) < sparsity); 
W(1:N+1:end) = 0;       
if w_rec~= 0
    W = W.*w_rec_mean/mean(W,"all"); % make sure the mean is correct
else
    W = W.*0;
end


mean_wff = 0.0001;
W_ff = (rand(N,N) < sparsity);   % new external weights
W_ff = W_ff/mean(W_ff,"all")*mean_wff;

