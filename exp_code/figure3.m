function [analysis, figures] = figure3
% Reproduce Figure 3b-d from the Zenodo single-trial PSTHs.
projectDir = fileparts(mfilename('fullpath'));
cfg = default_settings(projectDir);
cfg = load_data_settings(cfg);
cfg.referenceRate = 32;
cfg.outputDir = fullfile(projectDir,'output_from_psth');
[analysis, figures] = analyze_trials(cfg);
end

function cfg = default_settings(projectDir)
cfg.psthFile = fullfile(projectDir,'All_trial_PSTHs.h5');
cfg.sessions = {};
cfg.exampleNeuronIDs = [703 1090 780]; % Example neuron IDs in concatenated session order
cfg.outputDir = fullfile(projectDir,'output');
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
end

function cfg=load_data_settings(cfg)
assert(isfile(cfg.psthFile),'Download All_trial_PSTHs.h5 from https://zenodo.org/records/23089809 and place it in this folder.');
assert(h5readatt(cfg.psthFile,'/','complete')==1,'PSTH export is incomplete');
assert(strcmp(char(h5readatt(cfg.psthFile,'/','schema_version')),'1.0'),'Unsupported PSTH schema');
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

function S=load_session(cfg,k)
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

function Y=read_trial(S,trial)
Y=double(h5read(S.file,[S.group '/psth_Hz'],[1 1 trial],[S.nTime S.nUnits 1]))';
assert(all(isfinite(Y),'all') && all(Y>=0,'all'),'Invalid PSTH values');
end

function [Q,details] = fit_modes(M,t,goTime)
cr=M(:,:,1); cl=M(:,:,2); er=M(:,:,3); el=M(:,:,4);
avg=(cr+cl)/2; A=zeros(size(cr,1),8);
A(:,1)=meanwindow((cr+er-cl-el)/2,t,0,1.3);
A(:,2)=meanwindow((cr+el-cl-er)/2,t,1.3,goTime);
A(:,3)=meanwindow((cr+cl-er-el)/2,t,goTime,goTime+1.3);
A(:,4)=meanwindow(cr-cl,t,.2,.4);
A(:,5)=meanwindow(cr-cl,t,goTime-.3,goTime-.1);
A(:,6)=meanwindow(cr-cl,t,goTime+.1,goTime+.3);
A(:,7)=meanwindow(avg,t,goTime-.3,goTime-.1)-meanwindow(avg,t,-.3,-.1);
A(:,8)=meanwindow(avg,t,goTime,goTime+.1)-meanwindow(avg,t,goTime-.1,goTime);
lengths=sqrt(sum(A.^2,1));
assert(all(lengths>eps),'A task-defined mode is zero');
A=A./lengths;
assert(rank(A)==8,'Task-defined modes are dependent; original mode indices would be ambiguous');
Q=zeros(size(A));
for i=1:8
    u=A(:,i);v=zeros(size(u));
    for j=1:i-1,v=v-dot(u,Q(:,j))*Q(:,j);end
    w=u+v;Q(:,i)=w/norm(w);
end
assert(norm(Q'*Q-eye(8),'fro')<1e-8,'Mode orthogonality check failed');
details.rawNormalizedModes=A;
details.order={'stimulus','choice','outcome','early_sample','late_delay','early_response','ramp','go'};
details.fit='All available trials, descriptive in-sample modes; no held-out validation';
end
function y=meanwindow(X,t,a,b)
sel=t>a & t<b;assert(any(sel),'Empty mode window');y=mean(X(:,sel),2);
end

function [A,figures] = analyze_trials(cfg)
assert(isfile(cfg.psthFile),'PSTH file not found: %s',cfg.psthFile);
if ~isfolder(cfg.outputDir),mkdir(cfg.outputDir);end
ns=numel(cfg.sessions);nt=numel(cfg.time);
A.cfg=cfg;A.time=cfg.time;A.conditionNames={'correct_right','correct_left','error_right','error_left'};
A.meanRate={zeros(0,nt,4),zeros(0,nt,4)};
A.sessions=struct([]);offset=0;
for k=1:ns
    fprintf('Condition means %d/%d: %s\n',k,ns,cfg.sessions{k});S=load_session(cfg,k);
    assert(isequal(S.globalNeuronIDs,offset+(1:S.nUnits)),'Global neuron order changed');
    A.sessions(k).name=cfg.sessions{k};A.sessions(k).mouse=extractBefore(string(cfg.sessions{k}),'_');
    A.sessions(k).firstNeuron=offset+1;A.sessions(k).lastNeuron=offset+S.nUnits;
    A.sessions(k).sourceUnitIndices=S.unitIndex;A.sessions(k).nTrials=S.nTrials;
    A.sessions(k).trialIDs=S.trialIDs;
    for d=1:2
        means=zeros(S.nUnits,nt,4);
        for q=1:4
            ids=S.trialIDs{d,q};
            total=zeros(S.nUnits,nt);
            for tr=ids,total=total+read_trial(S,tr);end
            means(:,:,q)=total/numel(ids);
        end
        A.meanRate{d}=cat(1,A.meanRate{d},means);
    end
    offset=offset+S.nUnits;
end
A.nNeurons=offset;A.nSessions=ns;A.nMice=numel(unique([A.sessions.mouse]));
assert(all(cfg.exampleNeuronIDs<=offset),'Example neuron index outside selected sessions');
for d=1:2
    [A.basis{d},A.modeDetails{d}]=fit_modes(A.meanRate{d},cfg.time,cfg.goTimes(d));
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
    fprintf('All-trial variability %d/%d: %s\n',k,ns,cfg.sessions{k});S=load_session(cfg,k);
    ix=A.sessions(k).firstNeuron:A.sessions(k).lastNeuron;
    ex=find(ismember(cfg.exampleNeuronIDs,ix));localIDs=cfg.exampleNeuronIDs(ex)-ix(1)+1;
    for d=1:2
        weights=A.basis{d}(ix,[7 2]);weights(:,2)=-weights(:,2); % Sign convention for the displayed choice mode
        for q=1:2
            ids=S.trialIDs{d,q};n=numel(ids);assert(n>=2,'SD requires two trials');
            projected=zeros(n,nt,2);exampleRates=zeros(n,nt,numel(ex));
            for tr=1:n
                rates=read_trial(S,ids(tr));
                Z=smoothdata(weights'*rates,2,'movmean',cfg.displayWindow);
                projected(tr,:,:)=reshape(Z',[1 nt 2]);
                if ~isempty(ex)
                    E=smoothdata(rates(localIDs,:),2,'movmean',cfg.displayWindow);
                    exampleRates(tr,:,:)=reshape(E',[1 nt numel(ex)]);
                end
            end
            A.panelD.sessionTrials{k,d,q}=single(projected); % Processed rates, not spike times
            sessionMean=reshape(mean(projected,1),nt,2);
            sessionVar=reshape(var(projected,0,1),nt,2);
            A.panelD.mean(:,:,d,q)=A.panelD.mean(:,:,d,q)+sessionMean;
            A.panelD.variance(:,:,d,q)=A.panelD.variance(:,:,d,q)+sessionVar;
            reference=project_condition_mean(A.meanRate{d}(ix,:,q),weights,cfg);
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
check_results(A);
figures=plot_panels(A,cfg);
export_results(A,cfg,figures);
end
function Y=project_condition_mean(M,W,cfg)
Y=smoothdata((W'*M)',1,'movmean',cfg.displayWindow);
end

function F=plot_panels(A,cfg)
t=A.time(:); colors={cfg.Rcolor,cfg.Lcolor};
F.b=newfig('Fig. 3b: all trials',cfg.sizeB,cfg);
for e=1:3
    E=A.panelB.examples(e);
    for d=1:2
        ax=axes(F.b,'Position',[.17+(d-1)*.43 .735-(e-1)*.30 .38 .23]);
        hold(ax,'on'); hh=gobjects(1,2);
        for q=1:2,hh(q)=band(ax,t,E.mean{d,q},E.sd{d,q},colors{q},cfg);end
        yl=limits(E.mean{d,1},E.mean{d,2},E.sd{d,1},E.sd{d,2});yl(1)=0;
        format(ax,yl,cfg); marks(ax,d,yl,e==1,cfg);
        if d==1
            ylabel(ax,'Firing rate (Hz)');
            text(ax,-.34,.5,sprintf('N%d',e),'Units','normalized','FontName','Arial', ...
                'FontSize',cfg.fontSize,'FontWeight','bold','HorizontalAlignment','right','Clipping','off');
        end
        if e==3,xlabel(ax,'Time (s)');end
        if e==1 && d==2, key(ax,hh,cfg);end
    end
end
F.c=newfig('Fig. 3c: ramp-classified neurons',cfg.sizeC,cfg);
tl=tiledlayout(F.c,1,2,'TileSpacing','compact','Padding','compact');
for d=1:2
    ax=nexttile(tl);hold(ax,'on');hh=gobjects(1,2);
    for q=1:2,hh(q)=plot(ax,t,A.panelC.mean{d,q},'Color',colors{q},'LineWidth',cfg.lineWidth);end
    format(ax,[0 50],cfg); marks(ax,d,[0 50],true,cfg);
    yline(ax,cfg.referenceRate,'--','Color',[.6 .6 .6],'LineWidth',cfg.axisWidth,'HandleVisibility','off');
    xlabel(ax,'Time (s)');if d==1,ylabel(ax,'Firing rate (Hz)');else,key(ax,hh,cfg);end
end
title(tl,sprintf('Average PSTHs of ramping neurons (%d/%d)',numel(A.rampNeuronIDs),A.nNeurons), ...
    'FontName','Arial','FontWeight','normal','FontSize',cfg.fontSize);
F.d=newfig('Fig. 3d: session-based trial variability',cfg.sizeD,cfg);
tl=tiledlayout(F.d,2,2,'TileSpacing','compact','Padding','compact');
modeNames={'Ramp mode','Choice mode'};
for m=1:2
    for d=1:2
        ax=nexttile(tl);hold(ax,'on');hh=gobjects(1,2);
        for q=1:2,hh(q)=band(ax,t,A.panelD.mean(:,m,d,q),A.panelD.sd(:,m,d,q),colors{q},cfg);end
        yl=limits(A.panelD.mean(:,m,d,1),A.panelD.mean(:,m,d,2),A.panelD.sd(:,m,d,1),A.panelD.sd(:,m,d,2));
        format(ax,yl,cfg);marks(ax,d,yl,m==1,cfg);
        if d==1,ylabel(ax,{modeNames{m},'Projection (arb. units)'},'FontSize',12);end
        if m==2,xlabel(ax,'Time (s)');end
        if m==1 && d==2,key(ax,hh,cfg);end
    end
end
drawnow;
end
function f=newfig(name,sz,cfg)
f=figure('Name',name,'Color','w','Units','centimeters','Position',[2 2 sz], ...
    'Renderer','painters','Visible',cfg.visible,'NumberTitle','off');
set(f,'PaperUnits','centimeters','PaperSize',sz,'PaperPosition',[0 0 sz],'PaperPositionMode','manual');
end
function h=band(ax,t,y,s,c,cfg)
y=y(:);s=s(:);
fill(ax,[t;flipud(t)],[y-s;flipud(y+s)],c,'FaceAlpha',.22,'EdgeColor','none','HandleVisibility','off');
h=plot(ax,t,y,'Color',c,'LineWidth',cfg.lineWidth);
end
function yl=limits(r,l,sr,sl)
v=[r(:)-sr(:);r(:)+sr(:);l(:)-sl(:);l(:)+sl(:)];
lo=min(v);hi=max(v);span=max(hi-lo,1);yl=[lo-.05*span hi+.18*span];
end
function format(ax,yl,cfg)
set(ax,'FontName','Arial','FontSize',cfg.fontSize,'LineWidth',cfg.axisWidth, ...
    'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1,'TickDir','out', ...
    'TickLength',[.015 .015],'Box','off','XColor','k','YColor','k','Layer','top');
xlim(ax,cfg.xLimits);ylim(ax,yl);xticks(ax,0:4);
end
function marks(ax,d,yl,label,cfg)
go=cfg.goTimes(d);
for event=[0 cfg.sampleOffset go]
    xline(ax,event,'--','Color',[.5 .5 .5],'LineWidth',.6,'HandleVisibility','off');
end
barY=yl(2)-.1*diff(yl);
line(ax,[cfg.sampleOffset go],[barY barY],'Color',[.75 .86 .66],'LineWidth',3,'HandleVisibility','off');
if label,text(ax,mean([cfg.sampleOffset go]),barY+.025*diff(yl),sprintf('%g s',go-cfg.sampleOffset), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','Arial','FontSize',cfg.fontSize);end
end
function key(ax,h,cfg)
legend(ax,h,{'R','L'},'Location','northeast','Box','off','Color','white','FontSize',cfg.fontSize-2);
end

function export_results(A,cfg,F)
out=cfg.outputDir;
if cfg.exportFigures
    save_figures(F,cfg);
end
rows=cell(A.nSessions*4,7);z=0;delays=[2 1];choices={'R','L'};
for s=1:A.nSessions
    S=A.sessions(s);
    for d=1:2,for q=1:2
        z=z+1;rows(z,:)={S.name,char(S.mouse),S.lastNeuron-S.firstNeuron+1,delays(d),choices{q},S.nTrials(d,q),S.nTrials(d,q+2)};
    end,end
end
T=cell2table(rows,'VariableNames',{'session','mouse','ALM_units','delay_s','choice_label','correct_trials_plotted','error_trials_mode_fit_only'});
writetable(T,fullfile(out,'session_trial_counts.csv'));
rows=cell(12,6);z=0;
for e=1:3
    E=A.panelB.examples(e);
    for d=1:2,for q=1:2
        z=z+1;rows(z,:)={sprintf('N%d',e),E.globalNeuronID,E.session,delays(d),choices{q},E.nTrials(d,q)};
    end,end
end
writetable(cell2table(rows,'VariableNames',{'example','global_neuron_ID','session','delay_s','choice','n_trials'}),fullfile(out,'Fig3b_trial_counts.csv'));
tb=table(A.time(:),'VariableNames',{'time_s'});tc=tb;td=tb;
for d=1:2,for q=1:2
    tag=sprintf('%ds_%s',delays(d),choices{q});
    for e=1:3
        E=A.panelB.examples(e);base=sprintf('N%d_%s',e,tag);
        tb.([base '_mean_Hz'])=E.mean{d,q}(:);tb.([base '_SD_Hz'])=E.sd{d,q}(:);
    end
    tc.(['delay_' tag '_mean_Hz'])=A.panelC.mean{d,q}(:);
    for m=1:2
        modes={'ramp','choice'};base=[modes{m} '_' tag];
        td.([base '_mean'])=A.panelD.mean(:,m,d,q);td.([base '_SD'])=A.panelD.sd(:,m,d,q);
    end
end,end
writetable(tb,fullfile(out,'Fig3b_curves.csv'));writetable(tc,fullfile(out,'Fig3c_curves.csv'));writetable(td,fullfile(out,'Fig3d_curves.csv'));
writematrix(A.rampNeuronIDs,fullfile(out,'ramp_neuron_IDs.csv'));
if cfg.saveProcessedData
    fprintf('Saving processed trial rates and projections (no raw spikes)...\n');
    save(fullfile(out,'Fig3_bcd_processed.mat'),'A','-v7.3');
end
fid=fopen(fullfile(out,'analysis_report.txt'),'w');clean=onCleanup(@()fclose(fid));
fprintf(fid,'Fig. 3b/c/d all-trial analysis\nGenerated: %s\n',datestr(now,31));
fprintf(fid,'%d ALM neurons; %d sessions; %d mouse IDs; %d ramp-classified neurons.\n',A.nNeurons,A.nSessions,A.nMice,numel(A.rampNeuronIDs));
fprintf(fid,'R/L correspond to the original R_hit/L_hit flags; only correct trials are plotted.\n');
fprintf(fid,'All correct and error trials contribute to the mode fits. Ignore trials are excluded.\n');
fprintf(fid,'b: mean +/- sample SD across all correct trials for each example and condition.\n');
fprintf(fid,'c: equal-weight mean of the %d selected neuron condition means; no error band.\n',numel(A.rampNeuronIDs));
fprintf(fid,'c: reference line is %g Hz.\n',cfg.referenceRate);
fprintf(fid,'d: %s\n',A.panelD.definition);
fprintf(fid,'Trial counts differ between sessions and conditions; see session_trial_counts.csv.\n');
fprintf(fid,'No arbitrary across-session trial pairing, resampling, trial truncation, or random seed.\n');
fprintf(fid,'Modes and neuron selection are descriptive fits to the same data, not held-out estimates.\n');
fprintf(fid,'The SD bands are not estimates of between-animal uncertainty.\n');
fprintf('Exports complete: %s\n',out);
end

function save_figures(F,cfg)
for panel={'b','c','d'}
    p=panel{1};f=F.(p);name=fullfile(cfg.outputDir,['Fig3' p '_all_trials']);
    savefig(f,[name '.fig']);
    print(f,[name '.svg'],'-dsvg','-painters');
    exportgraphics(f,[name '.pdf'],'ContentType','vector','BackgroundColor','white');
    exportgraphics(f,[name '.png'],'Resolution',180,'BackgroundColor','white');
end
end

function check_results(A)
cfg=A.cfg; nt=numel(A.time);
for e=1:numel(A.panelB.examples)
    E=A.panelB.examples(e);
    for d=1:2
        for q=1:2
            Y=double(E.rates{d,q});
            assert(size(Y,1)==E.nTrials(d,q) && size(Y,2)==nt);
            assert(numel(E.trialIDs{d,q})==E.nTrials(d,q));
            assert(max(abs(mean(Y,1)-E.mean{d,q}))<1e-4);
            assert(max(abs(std(Y,0,1)-E.sd{d,q}))<1e-4);
            ref=smoothdata(A.meanRate{d}(E.globalNeuronID,:,q),2,'movmean',cfg.displayWindow);
            assert(max(abs(ref-E.mean{d,q}))<1e-8,'Example mean does not use all trials');
        end
    end
end
for d=1:2
    for q=1:2
        totalMean=zeros(nt,2);totalVar=zeros(nt,2);
        for s=1:A.nSessions
            Z=double(A.panelD.sessionTrials{s,d,q});
            assert(size(Z,1)==A.sessions(s).nTrials(d,q));
            totalMean=totalMean+reshape(mean(Z,1),nt,2);
            totalVar=totalVar+reshape(var(Z,0,1),nt,2);
            assert(isempty(intersect(A.sessions(s).trialIDs{1,q},A.sessions(s).trialIDs{2,q})), ...
                'Long and short trials overlap');
        end
        assert(max(abs(totalMean-A.panelD.mean(:,:,d,q)),[],'all')<1e-4);
        assert(max(abs(sqrt(totalVar)-A.panelD.sd(:,:,d,q)),[],'all')<1e-4);
    end
end
assert(all(isfinite(A.panelD.sd),'all') && all(A.panelD.sd>=0,'all'));
fprintf('Checks passed: trial counts, disjoint delays, all-trial means, and session-based SD.\n');
end
