SELF-ORGANIZED NEURAL INTEGRATORS IN NOISY SPIKING
NETWORKS (MATLAB)
================================================

This repository contains MATLAB code for simulating and analyzing recurrent spiking networks that exhibit integrator-like dynamics through biologically inspired learning rules (TTL) and mean-field theory (MFT) approximations.

------------------------------------------------
Directory Structure
------------------------------------------------
.
├── Integrator_learning/             
│   ├── learning_timing.m	# learning cortical timing
│   ├── LI_vs_wrec.m		# analyze liberality vs W_rec
│   └── self_organize_all_E_integrator.m 	# self-organize integrator
├── MFT/                             
│   ├── DDM_parameter.m		# estimate equivelent DDM parameters from integrator simulations
│   ├── integrator_RNN.m	# integrator model (demonstration, no learning, for any time step)
│   ├── IO_function.m		# helper function of mft_main.m, extracting IO curve from simulations
│   ├── mft_main.m		# run MFT analyze	
│   ├── neuron_parameters.m	# parameters for MFT analyze (100 neurons)
│   ├── neuron_parameters_1000.m	# parameters for MFT analyze (1000 neurons)
│   ├── neuron_parameters_taus20.m	# parameters for MFT analyze	(faster synaptic constant)
│   ├── ramp_vs_sd.m		# demonstrate Wieber's law (linear scaling)
│   ├── run1000neuron.m		# main function to run N = 1000 neurons



------------------------------------------------
How to Use Each Script
------------------------------------------------

[Integrator_learning folder]

1. self_organize_all_E_integrator.m
   Trains an excitatory-only network to become a linear integrator. (Figure 2)
   Usage:
       cd Integrator_learning
       self_organize_all_E_integrator


2. LI_vs_wrec.m (Figure 2)
   Evaluates linearity vs mean recurrent weight after learning.
   Usage:
       cd Integrator_learning
       LI_vs_wrec
   Scales W_rec from 0.7×–1.3× baseline (W_rec0) and computes the Linearity Index (R²).

3. learning_timing.m (Figure 3)
   Simulates a two-choice timing task combining decision and timer networks.
   Usage:
       cd Integrator_learning
       learning_timing
   Builds decision and timer networks, applies reward-modulated TTL rule, and plots timing and weight evolution.

[MFT folder]

4. mft_main.m (Figure 1)
   Main script for mean-field analysis and source–sink plotting. It calls IO_function.m and plot I/O curve too.
	To test a faster synaptic constant, change neuron_parameters to neuron_parameters_tau20 in the code.
   Usage:
       cd MFT
       mft_main


5. integrator_RNN.m 
   Demonstrates the integrator mechanism in recurrent network (N = 100), with any time step dt.
   Usage:
       cd MFT
       integrator_RNN

6. DDM_parameter.m (Figure 5)
   Extracts drift-diffusion model (DDM) parameters from noisy integrator trials.
   Usage:
       cd MFT
       DDM_parameter

7. ramp_vs_sd.m (Figure 3)
   Demonstrate linear scaling between variability and ramp time. 
   Usage:
       cd MFT
       ramp_vs_sd

8. run1000neuron.m (Figure S)
   Runs a large-scale simulation (N=1000 neurons) to test robustness.
   Usage:
       cd MFT
       run1000neuron

9. slow_ramping.m
    Runs a slow ramp simulation (N=100 neurons)
    Usage:
        cd MFT
        super_long




- May replace missing color map functions (e.g., slanCL) with MATLAB defaults.


------------------------------------------------
Contact
------------------------------------------------
For questions or updates, contact [bolu.feng@rice.edu].
