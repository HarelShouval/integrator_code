function fig3_selfcheck(A)
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
