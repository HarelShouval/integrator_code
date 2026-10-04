function cfg = fig3_config(projectDir)
cfg.dataFile = '';
cfg.inputSource = 'psth';
cfg.psthFile = fullfile(projectDir,'All_trial_PSTHs.h5');
cfg.sessions = {};
cfg.exampleNeuronIDs = [703 1090 780]; % Example neuron IDs in concatenated session order
cfg.outputDir = fullfile(projectDir,'output');
cfg.histTime = -1:0.001:5;
cfg.crop = 201:5801;               % Time support: -0.8 to 4.8 s
cfg.time = -0.8:0.001:4.8;        % Mode-fit time grid, including floating-point endpoints
cfg.rateWindow = 200;             % 200-ms boxcar rate estimate
cfg.dt = 0.001;
cfg.displayWindow = 100;          % Smooth each trial before calculating the mean and SD
cfg.goTimes = [3.3 2.3];
cfg.sampleOffset = 1.3;
cfg.selectionPercentile = 90;     % Percentile threshold for pooled absolute ramp-mode coefficients
cfg.referenceRate = 32;           % Display reference line for Fig. 3c
cfg.xLimits = [-0.8 4.5];
cfg.Rcolor = [0.25 0.50 0.72];
cfg.Lcolor = [0.90 0.30 0.33];
cfg.sizeB = [22 15]; cfg.sizeC = [22 8]; cfg.sizeD = [23 12];
cfg.fontSize = 16; cfg.lineWidth = 2; cfg.axisWidth = 1.5;
cfg.visible = 'on';
cfg.exportFigures = true;
cfg.saveProcessedData = true;     % Contains no spike timestamps or raw electrophysiology
cfg.previousFigures = '';
end
