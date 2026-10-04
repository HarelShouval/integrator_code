function [A,figures] = fig3_pipeline(cfg)
fromPSTH=isfield(cfg,'inputSource') && strcmp(cfg.inputSource,'psth');
if fromPSTH
    assert(isfile(cfg.psthFile),'PSTH file not found: %s',cfg.psthFile);
else
    assert(isfile(cfg.dataFile),'Data file not found: %s',cfg.dataFile);
end
if ~isfolder(cfg.outputDir),mkdir(cfg.outputDir);end
ns=numel(cfg.sessions);nt=numel(cfg.time);nb=numel(cfg.histTime);
A.cfg=cfg;A.time=cfg.time;A.conditionNames={'correct_right','correct_left','error_right','error_left'};
A.meanRate={zeros(0,nt,4),zeros(0,nt,4)};
A.sessions=struct([]);offset=0;
for k=1:ns
    fprintf('Condition means %d/%d: %s\n',k,ns,cfg.sessions{k});S=fig3_load_session(cfg,k);
    if fromPSTH,assert(isequal(S.globalNeuronIDs,offset+(1:S.nUnits)),'Global neuron order changed');end
    A.sessions(k).name=cfg.sessions{k};A.sessions(k).mouse=extractBefore(string(cfg.sessions{k}),'_');
    A.sessions(k).firstNeuron=offset+1;A.sessions(k).lastNeuron=offset+S.nUnits;
    A.sessions(k).sourceUnitIndices=S.unitIndex;A.sessions(k).nTrials=S.nTrials;
    A.sessions(k).trialIDs=S.trialIDs;
    for d=1:2
        means=zeros(S.nUnits,nt,4);
        for q=1:4
            ids=S.trialIDs{d,q};
            if fromPSTH
                total=zeros(S.nUnits,nt);
                for tr=ids,total=total+fig3_read_trial_psth(S,tr);end
                means(:,:,q)=total/numel(ids);
            else
                counts=zeros(S.nUnits,nb);
                for j=1:S.nUnits
                    cells=S.spikes(ids,S.unitIndex(j));
                    spikes=cellfun(@(v)v(:),cells,'UniformOutput',false);
                    spikes=vertcat(spikes{:});assert(all(isfinite(spikes)),'Nonfinite spike times');
                    counts(j,:)=hist(spikes,cfg.histTime)/numel(ids);
                end
                means(:,:,q)=fig3_rate(counts,cfg,false);
            end
        end
        A.meanRate{d}=cat(1,A.meanRate{d},means);
    end
    offset=offset+S.nUnits;
end
A.nNeurons=offset;A.nSessions=ns;A.nMice=numel(unique([A.sessions.mouse]));
assert(all(cfg.exampleNeuronIDs<=offset),'Example neuron index outside selected sessions');
for d=1:2
    [A.basis{d},A.modeDetails{d}]=fig3_modes(A.meanRate{d},cfg.time,cfg.goTimes(d));
end
coeff=[abs(A.basis{1}(:,7));abs(A.basis{2}(:,7))];
A.rampThreshold=prctile(coeff,cfg.selectionPercentile);
A.rampNeuronIDs=find(abs(A.basis{1}(:,7))>=A.rampThreshold & abs(A.basis{2}(:,7))>=A.rampThreshold);
assert(~isempty(A.rampNeuronIDs),'No ramp-classified neurons');
fprintf('Selected %d/%d ramping neurons, %d sessions, %d mice.\n',numel(A.rampNeuronIDs),offset,ns,A.nMice);
for d=1:2,for q=1:2,A.panelC.mean{d,q}=mean(A.meanRate{d}(A.rampNeuronIDs,:,q),1);end,end
A.panelC.referenceRate=cfg.referenceRate;
A.panelD.mean=zeros(nt,2,2,2); A.panelD.variance=zeros(nt,2,2,2);
A.panelD.sessionTrials=cell(ns,2,2);
A.panelB.examples=struct([]);
for k=1:ns
    fprintf('All-trial variability %d/%d: %s\n',k,ns,cfg.sessions{k});S=fig3_load_session(cfg,k);
    ix=A.sessions(k).firstNeuron:A.sessions(k).lastNeuron;
    ex=find(ismember(cfg.exampleNeuronIDs,ix));localIDs=cfg.exampleNeuronIDs(ex)-ix(1)+1;
    for d=1:2
        weights=A.basis{d}(ix,[7 2]);weights(:,2)=-weights(:,2); % Sign convention for the displayed choice mode
        for q=1:2
            ids=S.trialIDs{d,q};n=numel(ids);assert(n>=2,'SD requires two trials');
            projected=zeros(n,nt,2);exampleRates=zeros(n,nt,numel(ex));
            for tr=1:n
                if fromPSTH
                    rates=fig3_read_trial_psth(S,ids(tr));
                    Z=smoothdata(weights'*rates,2,'movmean',cfg.displayWindow);
                else
                    counts=zeros(S.nUnits,nb);
                    for j=1:S.nUnits,counts(j,:)=hist(S.spikes{ids(tr),S.unitIndex(j)},cfg.histTime);end
                    Z=fig3_rate(weights'*counts,cfg,true);
                end
                projected(tr,:,:)=reshape(Z',[1 nt 2]);
                if ~isempty(ex)
                    if fromPSTH
                        E=smoothdata(rates(localIDs,:),2,'movmean',cfg.displayWindow);
                    else
                        E=fig3_rate(counts(localIDs,:),cfg,true);
                    end
                    exampleRates(tr,:,:)=reshape(E',[1 nt numel(ex)]);
                end
            end
            A.panelD.sessionTrials{k,d,q}=single(projected); % Processed rates, not spike times
            sessionMean=reshape(mean(projected,1),nt,2);
            sessionVar=reshape(var(projected,0,1),nt,2);
            A.panelD.mean(:,:,d,q)=A.panelD.mean(:,:,d,q)+sessionMean;
            A.panelD.variance(:,:,d,q)=A.panelD.variance(:,:,d,q)+sessionVar;
            reference=fig3_rate_from_means(A.meanRate{d}(ix,:,q),weights,cfg);
            assert(max(abs(sessionMean-reference),[],'all')<1e-8,'All-trial means are inconsistent');
            for ee=1:numel(ex)
                e=ex(ee);Y=exampleRates(:,:,ee);
                A.panelB.examples(e).globalNeuronID=cfg.exampleNeuronIDs(e);
                A.panelB.examples(e).session=cfg.sessions{k};
                A.panelB.examples(e).sourceUnitIndex=S.unitIndex(localIDs(ee));
                A.panelB.examples(e).nTrials(d,q)=n;
                A.panelB.examples(e).trialIDs{d,q}=ids;
                A.panelB.examples(e).rates{d,q}=single(Y);
                A.panelB.examples(e).mean{d,q}=mean(Y,1);
                A.panelB.examples(e).sd{d,q}=std(Y,0,1);
            end
        end
    end
end
A.panelD.sd=sqrt(A.panelD.variance);
A.panelD.definition=['Estimated SD of a pseudopopulation projection = sqrt(sum of within-session ' ...
    'sample variances). All condition-matched trials are used. Preserves within-session covariance; ' ...
    'assumes independence across sessions. Not SEM, confidence intervals, or between-mouse variability.'];
fig3_selfcheck(A);
figures=fig3_plot(A,cfg);
fig3_export(A,cfg,figures);
end
function Y=fig3_rate_from_means(M,W,cfg)
Y=smoothdata((W'*M)',1,'movmean',cfg.displayWindow);
end
