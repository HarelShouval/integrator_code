function Y=fig3_read_trial_psth(S,trial)
Y=double(h5read(S.file,[S.group '/psth_Hz'],[1 1 trial],[S.nTime S.nUnits 1]))';
assert(all(isfinite(Y),'all') && all(Y>=0,'all'),'Invalid PSTH values');
end
