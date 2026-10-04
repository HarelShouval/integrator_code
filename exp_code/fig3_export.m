function fig3_export(A,cfg,F)
out=cfg.outputDir;
if cfg.exportFigures
    fig3_save_figures(F,cfg);
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
report=fig3_compare_previous(A);
if cfg.saveProcessedData
    fprintf('Saving processed trial rates and projections (no raw spikes)...\n');
    save(fullfile(out,'Fig3_bcd_processed.mat'),'A','report','-v7.3');
end
fid=fopen(fullfile(out,'analysis_report.txt'),'w');clean=onCleanup(@()fclose(fid));
fprintf(fid,'Fig. 3b/c/d all-trial analysis\nGenerated: %s\n',datestr(now,31));
fprintf(fid,'%d ALM neurons; %d sessions; %d mouse IDs; %d ramp-classified neurons.\n',A.nNeurons,A.nSessions,A.nMice,numel(A.rampNeuronIDs));
fprintf(fid,'R/L correspond to the original R_hit/L_hit flags; only correct trials are plotted.\n');
fprintf(fid,'All correct and error trials contribute to the mode fits. Ignore trials are excluded.\n');
fprintf(fid,'b: mean +/- sample SD across all correct trials for each example and condition.\n');
fprintf(fid,'c: equal-weight mean of the %d selected neuron condition means; no error band.\n',numel(A.rampNeuronIDs));
fprintf(fid,'c: reference line is %g Hz (previous plot script used 34 Hz).\n',cfg.referenceRate);
fprintf(fid,'d: %s\n',A.panelD.definition);
fprintf(fid,'Trial counts differ between sessions and conditions; see session_trial_counts.csv.\n');
fprintf(fid,'No arbitrary across-session trial pairing, resampling, trial truncation, or random seed.\n');
fprintf(fid,'Modes and neuron selection are descriptive fits to the same data, not held-out estimates.\n');
fprintf(fid,'The SD bands are not estimates of between-animal uncertainty.\n');
for p={'b','c','d'}
    r=report.(p{1});
    if isstruct(r),fprintf(fid,'Old/new Fig. 3%s: %d/%d mean curves agree within 1e-8; maximum absolute difference %.15g.\n', ...
            p{1},r.matchedCurves,r.totalCurves,r.maxAbsoluteMeanDifference);end
end
fprintf('Exports complete: %s\n',out);
end
