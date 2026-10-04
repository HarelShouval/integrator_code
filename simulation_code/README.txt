SELF-ORGANIZED NEURAL INTEGRATORS IN NOISY SPIKING NETWORKS
MATLAB CODE
============================================================

This repository contains MATLAB code for simulating noisy spiking-network
integrators, mean-field theory (MFT), two-trace learning (TTL), timing,
oculomotor persistence, and drift--diffusion dynamics.

------------------------------------------------------------
Directory Structure
------------------------------------------------------------

All simulation scripts and parameter files are in simulation_code/.
Open this folder in MATLAB before running the scripts below.

------------------------------------------------------------
How to Use the Scripts
------------------------------------------------------------

[Learning simulations]

1. self_organize_all_E_integerator.m
   Trains an excitatory network to become an approximate linear integrator
   using the TTL rule.

   Reproduces the main analyses in Figure 2c 2d 2f

   Usage:
       self_organize_all_E_integerator


2. LI_vs_wrec.m
   Evaluates the Linearity Index as a function of mean recurrent weight.

   Reproduces Figure 2e. Error bars are standard errors.

   Usage:
       LI_vs_wrec


3. learning_timing.m
   Simulates the two-stage decision--timing model and TTL-dependent
   adaptation of ramp speed.

   Reproduces Figure 3f 3h. Example parameters are included for both
   short-to-long and long-to-short delay switching.

   Usage:
       learning_timing


[Mean-field theory and network simulations]

4. mft_main.m
   Performs the mean-field analysis and plots the I/O relation and
   source--sink dynamics.

   Reproduces Figure 1 and related supplementary analyses.
   The script calls IO_function.m for IO function.

   To test different conditions, load the corresponding parameter file:
       neuron_parameters.m          Figure 1 related
       neuron_parameters_1000.m     N = 1000 Figure S1
       neuron_parameters_taus20.m   tau_s = 20 ms Figure S2

   Usage:
       mft_main


5. integrator_RNN.m
   Demonstrates integrator-like dynamics in a recurrent spiking network
   without learning.

   Usage:
       integrator_RNN


6. eyeposition.m and delta_r.m
   Simulate oculomotor-like persistent activity and the relationship
   between input strength and firing-rate change.

   These scripts reproduce the simulation panels in Figure 4b--c.

   Usage:
       eyeposition
       delta_r


7. DDM_parameter.m
   Estimates effective drift--diffusion parameters from integrator-network
   simulations.

   Reproduces the drift--diffusion analysis in Figure 4d.

   Usage:
       DDM_parameter


8. ramp_vs_sd.m
   Examines the linear relationship between mean ramp time and its
   trial-to-trial variability.

   Reproduces Figure 3h.

   Usage:
       ramp_vs_sd





------------------------------------------------------------
Contact
------------------------------------------------------------

For questions, contact:
bolu.feng@rice.edu