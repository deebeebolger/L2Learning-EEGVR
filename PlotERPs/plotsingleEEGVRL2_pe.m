function plotsingleEEGVRL2_pe (hdl,~)

D = get(hdl,'UserDat');
ecnt = D{1,1};
Cond1 = D{1,2}';
Cond2 = D{1,3}';
colrz = D{1,4};
time = D{1,5};
legendnoms = D{1,6};
eois = D{1,7}; 
topermut = D{1,8};

if size(eois,1)==1
    eois =eois';
end

assignin('base','Cond1',Cond1)
assignin('base','Cond2',Cond2)

f1 = figure;
set(f1,'Color',[1 1 1]);

axs = zeros(1,size(topermut,2));
linz = zeros(1,size(topermut,2));
           
[axs1,ploth1] = plot_ConfInt2(Cond1,time,nanmean(Cond1,1)',colrz(1,:));   %Call of plotConfInt2() function to calculate and plot the 95% CI (sem)
hold on
axs(1,1) = axs1;
linz(1,1) = ploth1;

[axs2,ploth2]= plot_ConfInt2(Cond2,time,nanmean(Cond2,1)',colrz(2,:));
hold on
axs(1,2) = axs2;
linz(1,2) = ploth2;


set(gca,'YDir','reverse','XAxisLocation','origin','YAxisLocation','origin','Box','off',...
        'YGrid','off','XGrid','off') % Set the current axis properties
Titre = title(eois{ecnt,1});
Titre.FontSize = 16;
xlabel(gca,'Time (ms)','FontSize',12)
ylabel(gca,'Potential (\muV)','FontSize', 12);

[tvals_curr,pvalue_corr,permout] = plot_perm(topermut,time,8,[time(1) time(2)],'yes');    
time_sig = time(~isnan(permout));
if ~isempty(time_sig)
    lowlims = time_sig(1);
    hilims = time_sig(end);
    ylimits = get(gca,'YLim');
    hold(gca,'on');
    shadcols={[1 1 0.2],[0.7 1 1],[1 0.8 0.8],[0.8 0.8 1],[0.8 1 0.8],[10.8 1],[1 1 0.7],[0.7 1 1],[1 0.8 0.8],[0.8 0.8 1],[0.8 1 0.8],[10.8 1]};
    for icnt=1:length(lowlims)
        yintval =ones(1,length(lowlims(icnt):hilims(icnt)))*ylimits(1);
        ha =fill([lowlims(icnt) lowlims(icnt):hilims(icnt) hilims(icnt)],[ylimits(2) yintval ylimits(2)],shadcols{1,icnt},'FaceAlpha',0.3,'EdgeColor','none');
        assignin('base','ha',ha);
        hold on
    end
end


legend(linz,legendnoms,'FontSize', 14);


end 
