function cfg=fig3_load_psth_settings(cfg)
assert(isfile(cfg.psthFile),'PSTH file not found. Run export_all_trial_psth.m first.');
assert(h5readatt(cfg.psthFile,'/','complete')==1,'PSTH export is incomplete');
assert(strcmp(char(h5readatt(cfg.psthFile,'/','schema_version')),'1.0'),'Unsupported PSTH schema');
cfg.inputSource='psth';
cfg.dataFile=''; % Use the public PSTH file as the data source.
cfg.sessions=strsplit(char(h5readatt(cfg.psthFile,'/','session_order')),',');
cfg.time=double(h5read(cfg.psthFile,'/time_s'))';
cfg.exampleNeuronIDs=double(h5read(cfg.psthFile,'/analysis/example_neuron_ids'))';
cfg.goTimes=double(h5read(cfg.psthFile,'/analysis/go_times_s'))';
cfg.selectionPercentile=double(h5read(cfg.psthFile,'/analysis/selection_percentile'));
cfg.referenceRate=double(h5read(cfg.psthFile,'/analysis/reference_rate_Hz'));
cfg.displayWindow=double(h5read(cfg.psthFile,'/analysis/display_window_samples'));
cfg.sampleOffset=double(h5read(cfg.psthFile,'/analysis/sample_offset_s'));
cfg.rateWindow=double(h5read(cfg.psthFile,'/preprocessing/rate_window_samples'));
cfg.dt=double(h5read(cfg.psthFile,'/preprocessing/dt_s'));
end
