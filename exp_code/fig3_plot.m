function F=fig3_plot(A,cfg)
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
