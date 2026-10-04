function figureS5(mode)
% Reproduce the experimental panels of Figure S5 from single-trial PSTHs.
if nargin==0, mode='run'; end
mode=validatestring(mode,{'run','plot'});
root=fileparts(mfilename('fullpath'));
out=fullfile(root,'output_S5_from_psth');
if ~isfolder(out), mkdir(out); end
cfg=default_settings(root);
cfg=load_data_settings(cfg);

longToShort={}; shortToLong={};
for k=1:numel(cfg.sessions)
    delay=double(h5read(cfg.psthFile,['/sessions/' cfg.sessions{k} '/delay_duration_s']));
    isLong=delay>1.5;
    assert(all(delay~=1.5) && nnz(diff(isLong))==1, ...
        'Expected exactly one switch in %s; inspect before analyzing.',cfg.sessions{k});
    if isLong(1) && ~isLong(end)
        longToShort{end+1}=cfg.sessions{k};
    else
        shortToLong{end+1}=cfg.sessions{k};
    end
end
sessions={longToShort,shortToLong};
types=[2 1];
bins={ [1 10;11 20;21 40;41 60;61 100], ...
       [1 10;11 30;31 60;61 90;91 150] };
names={'Long_to_short_exp','Short_to_long_exp'};
figureSize=[15 10]; fontSize=16; lineWidth=1.5;

oldPath=path; cleanup=onCleanup(@() path(oldPath));

for direction=1:2
    file=fullfile(out,[names{direction} '_data.mat']);
    if strcmp(mode,'run')
        requested=sessions{direction};
        present=ismember(requested,cfg.sessions);
        missing=requested(~present);
        if ~isempty(missing)
            warning('S5:MissingSessions','Missing sessions: %s',strjoin(missing,', '));
        end
        assert(any(present),'No requested S5 sessions present');
        result=analyze(cfg,requested(present),types(direction),bins{direction});
        result.requested_sessions=requested;
        result.missing_sessions=missing;
        result.exact_session_coverage=isempty(missing);
        result.input_file=cfg.psthFile;
        exportTables(result,out,names{direction});
        save(file,'result','-v7.3');
    else
        loaded=load(file,'result'); result=loaded.result;
    end
    fig=drawFigure(result,figureSize,fontSize,lineWidth);
    savefig(fig,fullfile(out,[names{direction} '.fig']));
    print(fig,fullfile(out,[names{direction} '.svg']),'-dsvg','-painters');
    exportgraphics(fig,fullfile(out,[names{direction} '.png']),'Resolution',150);
    fprintf('Saved %s.fig, .svg and .png.\n',names{direction});
end
end

function result=analyze(cfg,sessionNames,type_to_plot,trial_bins)
T_cue_aligned_sel=cfg.time;
window_size_plot=100;
ramp_mode_idx=7; top_frac=.10;
nSessions = numel(sessionNames);

session_bin_means = cell(nSessions, 1);     % each session: cell array of traces
session_other_means = cell(nSessions, 1);   % each session: opposite-type last 10 trace
session_bin_labels = cell(nSessions, 1);
session_n_ramping = nan(nSessions,1);
session_neurons=cell(nSessions,1);

plot_tag = '';
other_tag = '';
reward_t = [];
other_reward_t = [];

for isess = 1:nSessions
    sessionName = sessionNames{isess};
    fprintf('\n============================================================\n');
    fprintf('Processing session: %s\n', sessionName);

    k=find(strcmp(cfg.sessions,sessionName));
    S=load_session(cfg,k);
    isLong=S.delay>1.5;
    changes=find(diff(isLong)~=0);
    assert(numel(changes)==1 && all(S.delay~=1.5),'Expected one delay switch');
    assert(isLong(end)==(type_to_plot==1),'Wrong switch direction');
    bases=cell(1,2);
    for d=1:2
        M=zeros(S.nUnits,S.nTime,4);
        for q=1:4
            ids=S.trialIDs{d,q};
            for tr=ids(:)'
                M(:,:,q)=M(:,:,q)+read_trial(S,tr);
            end
            M(:,:,q)=M(:,:,q)/numel(ids);
        end
        bases{d}=fit_modes(M,T_cue_aligned_sel,cfg.goTimes(d));
    end
    basis_long=bases{1}; basis_short=bases{2};
    clear M

    long_ramp_loading  = abs(basis_long(:,  ramp_mode_idx));
    short_ramp_loading = abs(basis_short(:, ramp_mode_idx));

    nNeurons = numel(long_ramp_loading);
    nTop = max(1, ceil(top_frac * nNeurons));

    [~, idx_sort_long]  = sort(long_ramp_loading,  'descend');
    [~, idx_sort_short] = sort(short_ramp_loading, 'descend');

    top_long_idx  = idx_sort_long(1:nTop);
    top_short_idx = idx_sort_short(1:nTop);

    ramping_neurons = intersect(top_long_idx, top_short_idx);

    fprintf('Total neurons: %d\n', nNeurons);
    fprintf('Top 10%% threshold count per basis: %d\n', nTop);
    fprintf('Ramping neurons in BOTH long and short basis: %d\n', numel(ramping_neurons));

    if isempty(ramping_neurons)
        warning('No overlapping ramping neurons found in session %s. Skipping.', sessionName);
        continue;
    end

    session_n_ramping(isess) = numel(ramping_neurons);
    session_neurons{isess}=ramping_neurons;

    X=zeros(S.nTotalTrials,numel(ramping_neurons),S.nTime);
    for tr=1:S.nTotalTrials
        Y=read_trial(S,tr);
        X(tr,:,:)=reshape(Y(ramping_neurons,:),[1 numel(ramping_neurons) S.nTime]);
    end
    allTrials.long.psth=X(isLong,:,:);
    allTrials.short.psth=X(~isLong,:,:);
    clear X Y

    switch type_to_plot
        case 1
            trials_to_use       = allTrials.long.psth;
            nTrialsToUse        = size(allTrials.long.psth, 1);
            reward_t            = 3.3;
            plot_tag            = 'long';

            other_trials_to_use = allTrials.short.psth;
            other_reward_t      = 2.3;
            other_tag           = 'short';

        case 2
            trials_to_use       = allTrials.short.psth;
            nTrialsToUse        = size(allTrials.short.psth, 1);
            reward_t            = 2.3;
            plot_tag            = 'short';

            other_trials_to_use = allTrials.long.psth;
            other_reward_t      = 3.3;
            other_tag           = 'long';

        otherwise
            error('type_to_plot must be 1 (long) or 2 (short).');
    end

    this_session_bin_mean = {};
    this_session_bin_labels = {};

    for ibin = 1:size(trial_bins,1)
        sidx = trial_bins(ibin,1);
        eidx = trial_bins(ibin,2);

        if sidx > nTrialsToUse
            continue;
        end

        eidx = min(eidx, nTrialsToUse);

        this_block = trials_to_use(sidx:eidx, :, :);

        this_mean = squeeze(mean(this_block, 1));   % [nRampNeurons x nTime]
        if isvector(this_mean)
            this_mean = this_mean(:)'; 
        end
        this_pop_mean = mean(this_mean, 1);         % [1 x nTime]

        this_session_bin_mean{end+1} = this_pop_mean; %#ok<AGROW>
        this_session_bin_labels{end+1} = sprintf('%d-%d', sidx, eidx); %#ok<AGROW>
    end

    if isempty(this_session_bin_labels)
        warning('No valid trial bins in session %s. Skipping.', sessionName);
        continue;
    end

    fprintf('Using %d %s-trial bins in %s:\n', numel(this_session_bin_labels), plot_tag, sessionName);
    disp(this_session_bin_labels');

    nOtherTrials = size(other_trials_to_use, 1);
    if nOtherTrials < 1
        warning('No opposite-type trials found in session %s. Skipping overlay.', sessionName);
        this_other_pop_mean = nan(1, numel(T_cue_aligned_sel));
    else
        other_sidx = max(1, nOtherTrials - 9);
        other_eidx = nOtherTrials;

        other_block = other_trials_to_use(other_sidx:other_eidx, :, :);
        other_mean = squeeze(mean(other_block, 1));   % [nRampNeurons x nTime]
        if isvector(other_mean)
            other_mean = other_mean(:)';
        end
        this_other_pop_mean = mean(other_mean, 1);

        fprintf('Overlaying last %d trials from %s type: %d-%d\n', ...
            other_eidx - other_sidx + 1, other_tag, other_sidx, other_eidx);
    end

    session_bin_means{isess}   = this_session_bin_mean;
    session_other_means{isess} = this_other_pop_mean;
    session_bin_labels{isess}  = this_session_bin_labels;

    fprintf('Ramping neuron indices for %s:\n', sessionName);
    disp(ramping_neurons(:)');
end

valid_sessions = find(~cellfun(@isempty, session_bin_means));
if isempty(valid_sessions)
    error('No valid sessions available for averaging.');
end

fprintf('\n============================================================\n');
fprintf('Valid sessions used for cross-session averaging:\n');
disp(sessionNames(valid_sessions)');

fprintf('Number of ramping neurons per valid session:\n');
for i = 1:numel(valid_sessions)
    isess = valid_sessions(i);
    fprintf('%s: %d\n', sessionNames{isess}, session_n_ramping(isess));
end

nBinsPerSession = cellfun(@numel, session_bin_means(valid_sessions));
nBinsToUse = min(nBinsPerSession);

if nBinsToUse < 1
    error('No common trial bins across sessions.');
end

bin_labels = session_bin_labels{valid_sessions(1)}(1:nBinsToUse);

fprintf('\nUsing %d common %s-trial bins across sessions:\n', nBinsToUse, plot_tag);
disp(bin_labels');

bin_mean_ramping = cell(1, nBinsToUse);

valid_weights = session_n_ramping(valid_sessions);
valid_weights = valid_weights(:);   % column vector

for ibin = 1:nBinsToUse
    traces_this_bin = nan(numel(valid_sessions), numel(T_cue_aligned_sel));

    for i = 1:numel(valid_sessions)
        isess = valid_sessions(i);
        traces_this_bin(i,:) = session_bin_means{isess}{ibin};
    end

    weighted_sum = sum(traces_this_bin .* valid_weights, 1, 'omitnan');
    weight_denom = sum((~isnan(traces_this_bin)) .* valid_weights, 1);

    bin_mean_ramping{ibin} = weighted_sum ./ weight_denom;
end

other_traces = nan(numel(valid_sessions), numel(T_cue_aligned_sel));
for i = 1:numel(valid_sessions)
    isess = valid_sessions(i);
    other_traces(i,:) = session_other_means{isess};
end

weighted_sum_other = sum(other_traces .* valid_weights, 1, 'omitnan');
weight_denom_other = sum((~isnan(other_traces)) .* valid_weights, 1);
other_pop_mean = weighted_sum_other ./ weight_denom_other;

result.time_s=T_cue_aligned_sel;
result.mean_Hz=vertcat(bin_mean_ramping{:});
result.plotted_Hz=smoothdata(result.mean_Hz,2,'movmean',window_size_plot);
result.pre_switch_Hz=other_pop_mean;
result.plotted_pre_switch_Hz=smoothdata(other_pop_mean,'movmean',window_size_plot);
result.bin_labels=bin_labels;
result.trial_bins=trial_bins(1:nBinsToUse,:);
result.session_names=sessionNames;
result.n_mice=numel(unique(extractBefore(string(sessionNames),'_')));
result.valid_sessions=valid_sessions;
result.session_n_ramping=session_n_ramping;
result.session_neurons=session_neurons;
result.session_bin_means=session_bin_means;
result.session_bin_labels=session_bin_labels;
result.session_pre_switch_means=session_other_means;
result.post_switch_go_s=reward_t;
result.pre_switch_go_s=other_reward_t;
result.type_to_plot=type_to_plot;
result.smoothing_ms=window_size_plot;
end

function fig=drawFigure(d,sz,fontSize,lineWidth)
fig=figure('Color','w','Units','centimeters','Position',[2 2 sz]);
a=axes(fig,'Position',[.16 .18 .80 .69]); hold(a,'on');
t=d.time_s; n=size(d.plotted_Hz,1);
startColor=[.82 .78 .80]; endColor=[.55 .08 .18];
colors=[linspace(startColor(1),endColor(1),n)' ...
    linspace(startColor(2),endColor(2),n)' ...
    linspace(startColor(3),endColor(3),n)'];
handles=gobjects(1,n+1);
for k=1:n
    pre=t<=d.post_switch_go_s; post=t>=d.post_switch_go_s;
    handles(k)=plot(a,t(pre),d.plotted_Hz(k,pre),'-','Color',colors(k,:), ...
        'LineWidth',lineWidth);
    plot(a,t(post),d.plotted_Hz(k,post),'-','Color',.3*colors(k,:)+.7, ...
        'LineWidth',lineWidth,'HandleVisibility','off');
end
pre=t<=d.pre_switch_go_s; post=t>=d.pre_switch_go_s;
handles(end)=plot(a,t(pre),d.plotted_pre_switch_Hz(pre),'k--','LineWidth',lineWidth);
plot(a,t(post),d.plotted_pre_switch_Hz(post),'--','Color',[.7 .7 .7], ...
    'LineWidth',lineWidth,'HandleVisibility','off');
events=[0 1.3 d.post_switch_go_s d.pre_switch_go_s];
for k=1:4
    ls='--'; if k==3, ls='-'; end
    xline(a,events(k),ls,'Color',[.5 .5 .5],'LineWidth',1,'HandleVisibility','off');
end
peak=max([d.plotted_Hz(:);d.plotted_pre_switch_Hz(:)]);
xlim(a,[-.8 4.8]); ylim(a,[0 1.10*peak]);
xlabel(a,'Time (s)'); ylabel(a,'Average firing rate (Hz)');
set(a,'FontName','Arial','FontSize',fontSize,'LineWidth',1.2, ...
    'Box','off','TickDir','out','LabelFontSizeMultiplier',1);
legend(a,handles,[d.bin_labels(:)' {'-10-0'}], ...
    'Location','northeast','Box','off','FontSize',fontSize-4);
end

function exportTables(d,out,name)
T=array2table([d.time_s(:),d.plotted_pre_switch_Hz(:),d.plotted_Hz'], ...
    'VariableNames',[{'time_s','pre_switch_Hz'}, ...
    arrayfun(@(i)sprintf('post_bin_%d_Hz',i),1:size(d.plotted_Hz,1),'UniformOutput',false)]);
writetable(T,fullfile(out,[name '_curves.csv']));
rows={};
for k=d.valid_sessions(:)'
    for j=1:numel(d.bin_labels)
        bounds=sscanf(d.session_bin_labels{k}{j},'%d-%d');
        rows(end+1,:)={d.session_names{k},d.session_n_ramping(k),j,bounds(1),bounds(2),diff(bounds)+1}; %#ok<AGROW>
    end
end
writetable(cell2table(rows,'VariableNames',{'session','n_neurons','bin','first_trial','last_trial','n_trials'}), ...
    fullfile(out,[name '_sample_sizes.csv']));
fid=fopen(fullfile(out,[name '_coverage.txt']),'w'); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'Sessions: %s\nMissing sessions: %s\nSelected neurons: %d\n', ...
    strjoin(d.session_names,', '),strjoin(d.missing_sessions,', '),sum(d.session_n_ramping(d.valid_sessions)));
fprintf(fid,'Mice: %d\n',d.n_mice);
fprintf(fid,'Displayed bins include all outcomes; neuron selection uses correct/error condition means.\n');
fprintf(fid,'Session curves weighted by selected neuron count; 100-ms moving average.\n');
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
