function rates = fig3_rate(counts,cfg,applyDisplaySmoothing)
rates=conv2(counts,ones(1,cfg.rateWindow)/(cfg.rateWindow*cfg.dt),'same');
rates=rates(:,cfg.crop);
if applyDisplaySmoothing
    rates=smoothdata(rates,2,'movmean',cfg.displayWindow);
end
end
