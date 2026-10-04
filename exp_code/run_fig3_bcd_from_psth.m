% Read All_trial_PSTHs.h5 and export experimental figures.
projectDir=fileparts(mfilename('fullpath'));
addpath(projectDir);
cfg=fig3_config(projectDir);
cfg=fig3_load_psth_settings(cfg);
cfg.referenceRate=32; % Figure 3c display reference; does not affect neuron selection.
cfg.previousFigures=''; % No dependency on any previously generated figures.
cfg.outputDir=fullfile(projectDir,'output_from_psth');
[analysis,figures]=fig3_pipeline(cfg);
fprintf('\nPSTH-only reconstruction complete: %s\n',cfg.outputDir);
