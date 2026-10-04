function S=fig3_load_psth_session(cfg,k)
S.file=cfg.psthFile;S.group=['/sessions/' cfg.sessions{k}];
S.unitIndex=double(h5read(S.file,[S.group '/source_unit_index']))';
S.globalNeuronIDs=double(h5read(S.file,[S.group '/global_neuron_id']))';
S.nUnits=numel(S.unitIndex);S.nTime=numel(cfg.time);
S.delay=double(h5read(S.file,[S.group '/delay_duration_s']))';
S.sampleDuration=double(h5read(S.file,[S.group '/sample_duration_s']))';
S.outcomeFlags=h5read(S.file,[S.group '/outcome_flags']);
S.originalTrialIDs=double(h5read(S.file,[S.group '/trial_id']))';
S.nTotalTrials=numel(S.originalTrialIDs);
assert(isequal(S.originalTrialIDs,1:S.nTotalTrials),'Trial order mismatch');
assert(isequal(size(S.outcomeFlags),[6 S.nTotalTrials]));
info=h5info(S.file,[S.group '/psth_Hz']);
assert(isequal(info.Dataspace.Size,[S.nTime S.nUnits S.nTotalTrials]),'PSTH dimensions mismatch');
code=[1 4 2 5];masks={S.delay>1.5,S.delay<1.5};
modeMask=false(S.nTotalTrials,1);plotMask=modeMask;
for d=1:2,for q=1:4
    S.trialIDs{d,q}=find(masks{d} & S.outcomeFlags(code(q),:)==1);
    S.nTrials(d,q)=numel(S.trialIDs{d,q});
    assert(S.nTrials(d,q)>0,'Empty PSTH condition');
    modeMask(S.trialIDs{d,q})=true;
    if q<=2,plotMask(S.trialIDs{d,q})=true;end
end,end
assert(isequal(modeMask,logical(h5read(S.file,[S.group '/used_for_mode_fitting']))));
assert(isequal(plotMask,logical(h5read(S.file,[S.group '/used_for_plotted_curves']))));
end
