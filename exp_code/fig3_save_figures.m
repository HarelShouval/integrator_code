function fig3_save_figures(F,cfg)
for panel={'b','c','d'}
    p=panel{1};f=F.(p);name=fullfile(cfg.outputDir,['Fig3' p '_all_trials']);
    savefig(f,[name '.fig']);
    print(f,[name '.svg'],'-dsvg','-painters');
    exportgraphics(f,[name '.pdf'],'ContentType','vector','BackgroundColor','white');
    exportgraphics(f,[name '.png'],'Resolution',180,'BackgroundColor','white');
end
end
