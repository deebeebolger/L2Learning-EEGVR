function plotsingle_tf(hdl,~)

    D=get(hdl,'UserDat');
    time = D{1,1};
    timeindx = D{1,2};
    freq = D{1,3};
    powdata = D{1,4};
    powlims = D{1,5};
    currchan = D{1,6};
    blstring = D{1,7};
    cone_oi = D{1,8};
    
    f1=figure; set(f1,'PaperUnits','normalized','PaperPosition',[0.0235308 0.0272775 0.894169 0.909249],'Color',[1 1 1]);
    
    T = time(timeindx);
    
    ax = axes('parent',f1);
    imagesc('Parent',ax,'XData',T,'YData',freq,'CData',powdata,'CDataMapping','scaled')
    ax.XLim = [T(1) T(end)];
    ax.YLim = [4 40];
    ax.CLim = powlims;
    ax.Layer = 'top';
    ax.XLabel.String = 'Time (seconds)';
    ax.YLabel.String = 'Frequency (Hz)';
    
    hold on
    plot(ax,time,cone_oi,'w--','linewidth',2)
    title([currchan,' : ',blstring]);
    colormap(ax,jet);
    cb = colorbar;
    cb.Label.String = ['ERSP (',blstring(1:2),' )'];

end