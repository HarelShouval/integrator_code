Experimental figures from single-trial PSTHs
MATLAB R2023a; Statistics and Machine Learning Toolbox (prctile).

Data download:
  https://zenodo.org/records/23089809
Download All_trial_PSTHs.h5 (extract the archive if necessary) and place
it in this folder. The PSTH file is not included in the GitHub repository.

In MATLAB, open this folder and run:
  run_all

To reproduce individual figures:
  figure3     Figure 3b-d
  figureS5    Figure S5 experimental panels

The three scripts include their helper functions; no additional code files
are required. The CSV files provide session, neuron, and trial metadata;
the analyses read the required metadata directly from the HDF5 file.

Outputs:
  output_from_psth/       Figure 3b-d: FIG, SVG, PDF, PNG, MAT, CSV
  output_S5_from_psth/    Figure S5: FIG, SVG, PNG, MAT, CSV

The data include 1290 ALM neurons, 2806 trials, 8 sessions, and 2 male mice.
Figure 3c uses a 32-Hz reference line and 76 ramp-classified neurons.
Figure S5 displays five bins spanning the first 100 or 150 post-switch
trials, depending on switch direction. All trials remain in the data file.

To replot saved Figure S5 results:
  figureS5('plot')

To open a saved figure:
  openfig('path/to/figure.fig','new','visible')
