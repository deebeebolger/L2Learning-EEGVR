%% Date: 03-03-2022            Programmed by: Deirdre BOLGER
%  Function to plot ERPs of the selected data.
%
%*****************************************************************


%% LOAD IN EEG DATA 
conds  = {'Filler', 'Test'};
numsuj = 2;
condnames = conds;

data_meansuj = cell(1,size(conds,2));
data_GA      = cell(1,size(conds,2));
trialcnt     = zeros(2,size(conds,2));

for condcnt = 1:length(conds)

    toplotERPDir = fullfile(filesep,'Volumes','deepassport','Projects','Projet-L2-VREEG','Processed_Segmented_Data','VBP',conds{1,condcnt},filesep);
    filenum      = dir(strcat(toplotERPDir,'*.set'));
    filenom      = {filenum.name};
    
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;                %open eeglab session
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    
    EEG = pop_loadset('filename',filenom,'filepath',toplotERPDir);
    
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(EEG.setname),'gui','off'); % current set = xx;
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',char(EEG.setname),'filepath',EEG.filepath);  % Saves a copy of the current resampled dataset to the current directory
    eeglab redraw
    
    
    %% DEFINE THE DATA TO PLOT
    time = EEG.times;
    dataAll = {ALLEEG(:).data}';
    D = cellfun(@(x) mean(x,3), dataAll, 'UniformOutput',false);
    datasuj = zeros(length(EEG.chanlocs), length(time),numsuj);

    for scnt = 1:numsuj

        datasuj(:,:,scnt) = D{scnt,1};

    end

    data_meansuj{1,condcnt} = datasuj;
    GA = squeeze(mean(datasuj,3));
    data_GA{1,condcnt} = GA;
    trialcnt(:, condcnt)      = cell2mat(cellfun(@(y) size(y,3), dataAll, 'UniformOutput',false))';

end


%% CREATE A LIST DIALOGUE TO SELECT THE CHANNELS THAT YOU WISH TO PLOT

chanscurr = EEG.chanlocs;
[Eindx,Etf] = listdlg('PromptString',{'Select electrodes to plot.',...
    'You can select several electrodes at once.',''},...
    'SelectionMode','multiple','ListString',{chanscurr.labels});  %Eindx gives the indices of the selected channels

%% Define the configuration of the ERP plot.

prompt = 'How do you want me to configure the ERP data: "electrode" or "grid" or "grid_intconf"?';
dlgtitre = 'Specify ERP Layout';
dims = [1 80];
answr = inputdlg(prompt,dlgtitre,dims);

%% PLOTTING

allchans = length(Eindx);  % The number of channels to plot;
notmtchans = Eindx;         % The indices of the channels to plot.
chanoi = {chanscurr(Eindx).labels};
t = time;
colrs=[ones(1,3).*[0.4 0.15 0.15];ones(1,3).*[0.4 0.4 0.4];ones(1,3).*[0 0.6 0.5];ones(1,3).*[0.2 0.1 0.8]];

if strcmp(answr{1,1},'electrode')

    % SET UP ELECTRODE CONFIGURATION

    hndl=figure; set(hndl,'PaperUnits','normalized','Position',[680 417 727 681],'Color',[1 1 1]);
    orient portrait; axis ('normal')

    xvals=zeros(allchans,1);
    yvals=zeros(allchans,1);
    pwidth    = 0.8;     % 0.75, width and height of plot array on figure
    pheight   = 0.8;
    axwidth=0.07;
    axheight=0.08;

    %Read in channel locations file
    [elocs,titres,theta,rads,inds]=readlocs(chanscurr(notmtchans));
    channoms= strvcat(chanscurr(1:allchans).labels);
    Th=pi/180*theta;           %convert degrees to radians

    %Convert from polar to cartesian
    [ycart,xcart]=pol2cart(Th,rads);
    xvals(notmtchans) = ycart;
    yvals(notmtchans) = xcart;

    %Find the positions of all the channels
    mtchans = setdiff(1:allchans,notmtchans);        %find the channels indices common to both
    allchans_sqrt = floor(sqrt(allchans))+1;

    for i=1:length(mtchans)

        xvals(mtchans(i))=0.5+0.2*floor((i-1)/74);  %allchans - x axes
        yvals(mtchans(i))=-0.2+mod(i-1,allchans)/74;   %allchans - y axes

    end

    channoms2=channoms(1:allchans,:);
    xvals = xvals(1:allchans);
    yvals = yvals(1:allchans);

    if length(xvals) > 1
        if length(unique(xvals)) > 1
            gcapos = get(gca,'Position'); axis off;
            xvals = (xvals-mean([max(xvals) min(xvals)]))/(max(xvals)-min(xvals)); % this recenters
            xvals = gcapos(1)+gcapos(3)/2+pwidth*xvals;  %this controls width of plot
        end
    end

    gcapos = get(gca,'Position'); axis off;
    yvals = gcapos(2)+gcapos(4)/2+pheight*yvals;  % controls height of plot


    ho=zeros(length(chanoi),1);
    sig_elecs = zeros(length(chanoi),length(t));
    sig_times = cell(length(chanoi),1);
    nosigs=0;
    pvalue_corr = cell(length(chanoi),1);
    tvals = cell(length(chanoi),1);

    %% PLOT THE DATA IN 10-20 CONFIGURATION

    axs     = zeros(2,1);   %initialise the axes
    linz    = zeros(2,1);  %initialise the line object
    toperm  = cell(1,2);
    legnoms = cell(1,2); %initialise the legend names cell string array

    Axes = [];

    for chancnt =1:length(chanoi)

        xcenter(chancnt) = xvals(notmtchans(chancnt));
        ycenter(chancnt) = yvals(notmtchans(chancnt));
        Axes = [Axes axes('Units','normalized','Position', [ycenter(chancnt)-axheight/2 xcenter(chancnt)-axwidth/2 axheight axwidth])];
        hold on;

        for dcnt = 1:size(data_GA,2)     % for each condition
            
            Dcurr           = [];
            Dcurr           = data_GA{1,dcnt}(chancnt,:);
            ploth1          = plot(time',Dcurr','Color',colrs(dcnt,:));   %Call of plotConfInt2() function to calculate and plot the 95% CI (sem)
            axs1            = gca;
            hold on
            axs(dcnt)       = axs1;
            linz(dcnt)      = ploth1;
            toperm{1,dcnt}  = squeeze(data_meansuj{1,dcnt}(chancnt,:,:));
            legnoms{1,dcnt} = condnames{1,dcnt};

        end

        set(axs1,'YDir','reverse','XAxisLocation','origin','YAxisLocation','origin','Box','off',...
            'YGrid','off','XGrid','off') % Set the current axis properties
        title(chanoi{1,chancnt})
        [tvals_curr,pvalue_corr,permout] = plot_PermT(toperm,time,3,'no');  %Call of function to carry out permutation t-test with fdr correction.

        set(hndl,'CurrentAxes',Axes(chancnt));
        set(Axes(chancnt),'HitTest','on','SelectionHighlight','on','UserData',{chancnt, squeeze(data_meansuj{1,1}(chancnt,:,:)),squeeze(data_meansuj{1,2}(chancnt,:,:)),...
            colrs, time, legnoms, chanoi, toperm},'NextPlot','add');
        set(Axes(chancnt),'ButtonDownFcn',@plotsingleEEGVRL2_pe)

        if chancnt == 1
            legend(linz,legnoms,'Position',[0.1 0.75 0.2 0.1],'FontSize',14,'Box','off');
        end

    end  % end of chancnt loop

elseif strcmp(answr{1,1},'grid')

    hndl=figure; set(hndl,'PaperUnits','normalized','Position',[680 417 727 681],'Color',[1 1 1]);
    orient portrait; axis ('normal')

    axs     = zeros(2,1);   %initialise the axes
    linz    = zeros(2,1);  %initialise the line object
    toperm  = cell(1,2);
    legnoms = cell(1,2); %initialise the legend names cell string array
    subcols = 4;
    subrows = ceil(length(chanoi)/subcols);
    

    Axes = [];

    for chancnt =1:length(chanoi)

        Axes(chancnt) = subplot(subcols,subrows,chancnt);


        for dcnt = 1:size(data_GA,2)     % for each condition
            Dcurr = [];
            Dcurr = data_GA{1,dcnt}(chancnt,:)';
            ploth1 = plot(time,Dcurr,'Color',colrs(dcnt,:));   %Call of plotConfInt2() function to calculate and plot the 95% CI (sem)
            axs1 = gca;
            hold on
            axs(dcnt) = axs1;
            linz(dcnt) = ploth1;
            toperm{1,dcnt} = squeeze(data_meansuj{1,dcnt}(chancnt,:,:));
            legnoms{1,dcnt} = condnames{1,dcnt};

        end

        set(axs1,'YDir','reverse','XAxisLocation','origin','YAxisLocation','origin','Box','off',...
            'YGrid','off','XGrid','off') % Set the current axis properties
        title(chanoi{1,chancnt})
        [tvals_curr,pvalue_corr,permout]= plot_PermT(toperm,time,3,'no');   %Call of function to carry out permutation t-test with fdr correction.

        set(hndl,'CurrentAxes',Axes(chancnt));
        set(Axes(chancnt),'HitTest','on','SelectionHighlight','on','UserData',{chancnt, squeeze(data_meansuj{1,1}(chancnt,:,:)),squeeze(data_meansuj{1,2}(chancnt,:,:)),...
            colrs, time, legnoms, chanoi, toperm},'NextPlot','add');
        set(Axes(chancnt),'ButtonDownFcn',@plotsingleEEGVRL2_pe)

        if chancnt == 1
            legend(linz,legnoms,'Position',[0.1 0.75 0.2 0.1],'FontSize',14,'Box','off');
        end

    end  % end of chancnt loop%
elseif strcmp(answr{1,1},'grid_intconf')

    chaninfo = EEG.chanlocs;

    chanoi = [];
    chanoi = {'F3' 'Fz' 'F4' 'FC3' 'FCz' 'FC4' 'C3' 'Cz' 'C4' 'P3' 'Pz' 'P4'};

    uiwait(msgbox({'The selected channels will be overwritten by channels defined in script:'; char(chanoi)},'Attention','warn','modal'));

    f1 = figure;
    set(f1,'Color',[1 1 1]);

    

    allchans = {chaninfo.labels};
    for i = 1:length(chanoi)
        chanindx(i) = find(strcmp(allchans,chanoi{1,i}));
    end

    for chancnt = 1:length(chanoi)

        toperm = cell(1,2);
        legnoms = cell(1,2);
        for dcnt = 1:size(data_GA,2)     % for each condition
            toperm{1,dcnt} = squeeze(data_meansuj{1,dcnt}(chanindx(chancnt),:,:));
            legnoms{1,dcnt} = condnames{1,dcnt};

        end
        
        
        subcols = 4;
        subrows = ceil(length(chanoi)/subcols);
        Axes(chancnt) = subplot(subrows,subcols,chancnt);

        axs  = zeros(1,size(toperm,2));
        linz = zeros(1,size(toperm,2));

        Cond1 = squeeze(data_meansuj{1,1}(chanindx(chancnt),:,:))';
        Cond2 = squeeze(data_meansuj{1,2}(chanindx(chancnt),:,:))';
       

        [axs1,ploth1] = plot_ConfInt2(Cond1,time,nanmean(Cond1,1)',colrs(1,:));   %Call of plotConfInt2() function to calculate and plot the 95% CI (sem)
        hold on
        axs(1,1) = axs1;
        linz(1,1) = ploth1;

        [axs2,ploth2]= plot_ConfInt2(Cond2,time,nanmean(Cond2,1)',colrs(2,:));
        hold on
        axs(1,2) = axs2;
        linz(1,2) = ploth2;

        set(Axes(chancnt),'YDir','reverse','XAxisLocation','origin','YAxisLocation','origin','Box','off',...
            'YGrid','off','XGrid','off') % Set the current axis properties 
        Titre = title(allchans{1,chanindx(chancnt)});
        Titre.FontSize = 16;
        xlabel(gca,'Time (ms)','FontSize',12)
        ylabel(gca,'Potential (\muV)','FontSize', 12);

        if chancnt ==1
            legend(linz,legnoms,'FontSize', 14);
        end
    end
end


T = sum(trialcnt,2);
sprintf('********Trial total for condition 1, %s: %d ***************',conds{1,1},T(1));
sprintf('********Trial total for condition 2,  %s: %d **************',conds{1,2}, T(2));


%%