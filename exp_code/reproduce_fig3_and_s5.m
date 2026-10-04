% Read All_trial_PSTHs.h5 and export experimental figures.
projectDir=fileparts(mfilename('fullpath'));
run(fullfile(projectDir,'run_fig3_bcd_from_psth.m'));
run_figS5_exp_from_psth;
