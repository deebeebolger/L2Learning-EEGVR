%% Date: 03-03-2022     Programmed by: Deirdre BOLGER
%  Script to isolate choosen conditions in already segmented data.
%  This script should run for several datasets.
%
%**********************************************************************

%% LOAD IN ALL THE DATASETS THAT YOU WISH TO SEPARATE.
%  You can select those segmented datasets manually.

segDir = fullfile(filesep,'Volumes','deepassport','Projects','Projet-L2-VREEG','Processed_Segmented_Data',filesep);
current_cond = 'Test';

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;                %open eeglab session
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

EEG = pop_loadset();

[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(EEG.setname),'gui','off'); % current set = xx;
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',char(EEG.setname),'filepath',EEG.filepath);  % Saves a copy of the current resampled dataset to the current directory
eeglab redraw

%% EXTRAIRE VERB-TYPES, BLOCK-NUMBERS AND ITEMS FROM EEG STRUCTURE OF FIRST DATASET (should be same for all datasets).

% Extraire verb-types

type_all = unique({EEG.event.eventlabels});
ix_vtype = contains(type_all,'verb');
vtype = cell(sum(ix_vtype),1);
[vtype{:,1}] = deal(type_all{1,ix_vtype});

% Extraire block numbers

block_all = unique({EEG.event.blocknum});
ix_btype  = contains(block_all,'block');
btype     = cell(sum(ix_btype),1);
[btype{:,1}] = deal(block_all{1,ix_btype});

% Extraire individual items

X = {EEG.event.items};
ix = cell2mat(cellfun(@ischar, X, 'UniformOutput',false));
itypes_all = cell(sum(ix),1);
[itypes_all{:,1}] = X{1,ix};
ittype = unique(itypes_all);
ix_item   = contains(ittype,'32');
itemtype  = cell(sum(ix_item),1);
[itemtype{:,1}] = deal(ittype{ix_item,1});


%% CHOOSE VERB TYPE, BLOCK NUMBER AND ITEMS (this will apply to all the selected datasets)

[indx_verb,tf] = listdlg('SelectionMode','multiple','ListString',vtype);

[indx_block,tf1] = listdlg('SelectionMode','multiple','ListString',btype);

[indx_item,tf2] = listdlg('SelectionMode','multiple','ListString',itemtype);


for count = 1:length(ALLEEG)
    
    %% RERIEVE THE CURRENT DATASET 

    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',count,'study',0); 
    eeglab redraw;
    
    disp(['The current dataset is ',EEG.filename])
    
    %% Determine the indices of the epochs to segment.

    % First check which block to segment.

    blockindx = find(contains({EEG.event.blocknum},btype(indx_block)));

    % For these selected blocks detect the correct verb types
    verbindx = find(contains({EEG.event(blockindx).eventlabels},vtype(indx_verb)));
    VIndex = blockindx(verbindx);

    % For those blocks with the selected verb type, isolate the selected items.


    itemIndx = find(contains({EEG.event(VIndex).items},itemtype(indx_item)));
    SegIndex = VIndex(itemIndx);

    % Verb title
    if length(indx_verb)==length(vtype)
        toseg_verbs = "allverbs-";
    else
        toseg_verbs = vtype{indx_verb} ;
    end

    % Block title
    if length(indx_block) == length(btype)
        toseg_block = "allblocks-";
    elseif length(indx_block) == 2
        toseg_block = strcat(btype{indx_block(1)},'_',btype{indx_block(2)});
    else
        toseg_block = btype{indx_block};
    end

    %% CARRY OUT CONDITION SEPARATION

    EpochIndex = zeros(length([EEG.epoch]),1);
    for i = 1:length([EEG.epoch])

        EpochIndex(i) = sum(ismember(EEG.epoch(i).event, SegIndex));

    end

    EEG = eeg_checkset( EEG );
    EEG = pop_select( EEG, 'trial',find(EpochIndex));
    EEG = eeg_checkset( EEG );

    %% DEFINE THE PATH TO WHICH TO SAVE THE NEW DATASET AND SAVE

    pathsuj = fullfile(segDir,'VBP',filesep);
    D = dir(pathsuj);
    parent_folder = D.folder;

    if sum(strcmp({D.name},current_cond)) == 0
        

        mkdir(fullfile(parent_folder,current_cond));
        savedir = fullfile(parent_folder,current_cond);

    else

        disp('Conditions folder already exists');
        savedir = fullfile(parent_folder,current_cond);
    end
    
    % Define the title of the new dataset. 
    epoched_name = strcat(EEG.filename(1:end-4),'-',toseg_block,toseg_verbs,'-epoched');

    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(epoched_name),'gui','off');
    EEG = eeg_checkset( EEG );
    

    EEG = pop_saveset( EEG, 'filename',char(epoched_name),'filepath',savedir);
    EEG = eeg_checkset( EEG );
    eeglab redraw


end % end of ALLEEG count
