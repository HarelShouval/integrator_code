Experimental figures from single-trial PSTHs
MATLAB R2023a; Statistics and Machine Learning Toolbox (prctile).

Data download (Zenodo):
  https://zenodo.org/records/23089809
The single-trial PSTH data file is not included in the GitHub repository.
Download All_trial_PSTHs.h5 from the Zenodo record (extract the archive
if necessary), then place it at:
  exp_code/All_trial_PSTHs.h5
Keep the filename unchanged before running the scripts below.

In MATLAB, open this folder and run:
  reproduce_fig3_and_s5

Or run separately:
  run_fig3_bcd_from_psth
  run_figS5_exp_from_psth

Input: All_trial_PSTHs.h5
ALM only: 1290 neurons, 2806 trials, 8 sessions, 2 male mice.
Raw spike data and previously generated figures are not required.

Outputs:
  output_from_psth/       Fig. 3b-d: FIG, SVG, PDF, PNG, MAT, CSV
  output_S5_from_psth/    S5 experimental panels: FIG, SVG, PNG, MAT, CSV

Fig. 3c reference line: 32 Hz; 76 ramp-classified neurons.
S5 displays five bins: first 100/150 post-switch trials. The archive
retains all trials. Simulation panels are not included in this package.

To open a saved figure visibly:
  openfig('path/to/figure.fig','new','visible')
