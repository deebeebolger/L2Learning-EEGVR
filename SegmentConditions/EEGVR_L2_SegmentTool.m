%% Most recent version of the Segmentation Tool for the EEGVR_l2Learning project.



list = {'SegmentAll','SeparateConds'};
[indx,tf] = listdlg('PromptString','Select the segmentation type:','SelectionMode','single','ListString',list,...
    'ListSize',[160 150]);

answr = list{1,indx};

switch answr
    case 'SegmentAll'


        %% LOAD IN THE INDIVIDUAL VERB DATA AND ADD TO THE EVENTS FIELD OF THE CURRENT EEG STRUCTURE


        [TrigInfile, TrigInpath] = uigetfile({'*.xlsx'},'Load in the VERB_TRIGS.xlsx file');
        [ParamInfile, ParamInpath] = uigetfile({'*.mat'}, 'Load in the EEGVR_L2_Parameters.mat file');

        trigdataIn = readtable(fullfile(TrigInpath, TrigInfile),'FileType','spreadsheet'); % HERE YOU NEED TO PUT IN THE PATH TO THE TRIG DATA EXCEL FILE
        load(fullfile(ParamInpath, ParamInfile));


        vprompt = 'List the verbs to exclude:';
        vdlgtitle = 'User Input';
        vdims = [1 35];
        vdefineinput = {'deplacer soulever secouer'};
        verbs2excl = inputdlg(vprompt, vdlgtitle, vdims, vdefineinput);

        % Need to parse the verbs defined.
        if length(verbs2excl)>1

            ifind =  find(isspace(verbs2excl{1,1}))+1;
            ifind = cat(2,ifind, length(verbs2excl{1,1})+2);
            ifind = cat(2,1, ifind);

            vIn = verbs2excl{1,1};
            VerbExcl = cell(length(ifind)-1,1);
            fcounter = 0;

            for fcnt = 1:length(VerbExcl)

                VerbExcl{fcnt,1} = vIn(ifind(fcnt):ifind(fcnt+1)-2);

            end

        else
            VerbExcl = [];

        end


        %% Load IN THE CURRENT EEG *.set FILE.

        [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;                %open eeglab session
        [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

        [ALLEEG, EEG] = pop_loadset();
        dataset_len = length(ALLEEG);

        for Ecounter = 1:dataset_len

            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',Ecounter,'gui','off'); % current set = xx;
            EEG = eeg_checkset( EEG );
            EEG = pop_saveset( EEG, 'filename',char(EEG.setname),'filepath',EEG.filepath);  % Saves a copy of the current resampled dataset to the current directory
            eeglab redraw

            %% ADD INDIVIDUAL ITEM NAMES TO THE EEG STRUCTURE

            Trigcodes = str2double({EEG.event.type});

            trigread = trigdataIn{:,2};
            itemread = trigdataIn{:,1};

            items_all = cell(length(Trigcodes),1);
            items_all(:) = {''};

            for tcount = 1:length(trigread)

                ix = [];
                ix = find(Trigcodes == trigread(tcount));
                [items_all{ix}] = deal(itemread{tcount,1}(1:end-4));

            end

            % Add to the events structure.

            for ecount = 1:length(items_all)

                b = EEG.event(ecount).blocknum;

                if ~isempty(b)
                    EEG.event(ecount).items = items_all{ecount};
                end

            end

            [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
            EEG = eeg_checkset( EEG );
            EEG = pop_saveset( EEG, 'filename',EEG.setname,'filepath',EEG.filepath);
            eeglab redraw

            %% DATA SEGMENTATION : EXTRAIRE THE VERB TYPES

            % Extraire verb-types
            E = {EEG.event.eventlabels};
            Emt = find(cell2mat(cellfun(@isempty, E, 'UniformOutput',false)));
            for iemt = 1:length(Emt)

                EEG.event(Emt(iemt)).eventlabels = '';
            end
            type_all = unique({EEG.event.eventlabels});

            ix_vtype = contains(type_all,'verb');
            vtype = cell(sum(ix_vtype),1);
            [vtype{:,1}] = deal(type_all{1,ix_vtype});

            % Extraire block numbers
            Bk = {EEG.event.blocknum};
            Bmt = find(cell2mat(cellfun(@isempty, Bk, 'UniformOutput',false)));
            for ibmt = 1:length(Bmt)

                EEG.event(Bmt(ibmt)).blocknum = '';
            end
            block_all = unique({EEG.event.blocknum});
            ix_btype  = contains(block_all,'block');
            btype     = cell(sum(ix_btype),1);
            [btype{:,1}] = deal(block_all{1,ix_btype});

            % Extraire individual items

            X = {EEG.event.items};
            ix = cell2mat(cellfun(@isempty, X, 'UniformOutput',false));
            ix1 = ~ix;
            itypes_all = cell(sum(ix1),1);
            [itypes_all{:,1}] = X{1,ix1};
            item_type = unique(itypes_all);


            if ~isempty(VerbExcl)
                keepindx = ~contains(item_type,VerbExcl);
                item_type = item_type(keepindx);
            end

            respcorr = {EEG.event.response}
            resptype = unique(respcorr);
            rindx = contains(resptype,'correct')
            resp_type = resptype(1, rindx);


            %% CHOOSE VERB TYPE, BLOCK NUMBER AND ITEMS

            if Ecounter == 1

                [indx_verb,tf] = listdlg('SelectionMode','multiple','ListString',vtype);
    
                [indx_block,tf1] = listdlg('SelectionMode','multiple','ListString',btype);
    
                [indx_item,tf2] = listdlg('SelectionMode','multiple','ListString',item_type);
    
                [indx_resp, tf3] = listdlg('SelectionMode','multiple','ListString',resp_type);
            end

            %% Find the intersection between index_verb, indx_block and indx_time
            %  Determine the indices of the epochs to segment.

            % First check which block to segment.
            blockindx = find(contains({EEG.event.blocknum},btype(indx_block)));

            % For these selected blocks detect the correct verb types
            verbindx = find(contains({EEG.event(blockindx).eventlabels},vtype(indx_verb)));
            VIndex = blockindx(verbindx);

            % For those blocks with the selected verb type, isolate the selected items.
            itemIndx = find(contains({EEG.event(VIndex).items},item_type(indx_item)));
            SegIndex = VIndex(itemIndx);

            % For those blocks with the selected verb and item type,
            % isolate the selected response.
            respIndx = find(strcmp({EEG.event(SegIndex).response}, resp_type(indx_resp)));
            RIndx = SegIndex(respIndx);

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


            %Response title
             if length(indx_resp) == length(resp_type)
                toseg_block = "allresponse-";
           
            else
                toseg_resp = resp_type{indx_resp};
            end


            %% CARRY OUT THE SEGMENTATION

            disp('--------------------Segmentation-----------------------------------');

            % Load in the parameters mat file with segmentation information.

            pstim = mystruct.segmentation.poststim_ms./1000;
            bline = mystruct.segmentation.baseline_ms;

            % Prepare epoched data title.
            s = strfind(EEG.setname,'-');
            epoch_name = strcat(EEG.setname(1:s(1)),toseg_block,toseg_verbs,toseg_resp,'-epoched');

            % Carry out the epoching
            timelims = [bline(1)/1000 pstim(2)];
            [outEEG, indxs ] = pop_epoch(EEG,{},timelims, 'eventindices',RIndx);


            pathsuj = EEG.filepath;
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(epoch_name),'gui','off'); % current set = xx;
            EEG = eeg_checkset( EEG );
            EEG = pop_saveset( outEEG, 'filename',char(epoch_name),'filepath',pathsuj);  % Saves a copy of the current resampled dataset to the current directory
            eeglab redraw

            %% CARRY OUT BASELINE CORRECTION

            disp('--------------------Baseline correction-----------------------------------');
            baseline_low = mystruct.segmentation.baseline_correction_ms(1);  % Change here to change baseline to use for correction.
            baseline_hi  = mystruct.segmentation.baseline_correction_ms(2);

            epochBL_name =strcat(epoch_name,'-bl');
            EEG = pop_rmbase( EEG, [baseline_low baseline_hi]);
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(epochBL_name),'gui','off');
            EEG = eeg_checkset( EEG );

            EEG = pop_saveset( EEG, 'filename',char(epochBL_name),'filepath',pathsuj);
            EEG = eeg_checkset( EEG );
            eeglab redraw

        end %End of Ecounter loop

    case 'SeparateConds'

        segDir = uigetdir(cd,'Choose a folder in which to save data');
        %segDir = fullfile(filesep,'Volumes','deepassport','Projects','Projet-L2-VREEG','Processed_Segmented_Data',filesep);

        [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;                %open eeglab session
        [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

        EEG = pop_loadset();

        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(EEG.setname),'gui','off'); % current set = xx;
        EEG = eeg_checkset( EEG );
        EEG = pop_saveset( EEG, 'filename',char(EEG.setname),'filepath',EEG.filepath);  % Saves a copy of the current resampled dataset to the current directory
        eeglab redraw

        %% EXTRAIRE VERB-TYPES, BLOCK-NUMBERS AND ITEMS FROM EEG STRUCTURE

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
        
        
        %% CHOOSE VERB TYPE, BLOCK NUMBER AND ITEMS

        [indx_verb,tf] = listdlg('SelectionMode','multiple','ListString',vtype);

        [indx_block,tf1] = listdlg('SelectionMode','multiple','ListString',btype);

        [indx_item,tf2] = listdlg('SelectionMode','multiple','ListString',itemtype);

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

        %% Carry out the separation of the conditions 
        EpochIndex = zeros(length([EEG.epoch]),1);
        for i = 1:length([EEG.epoch])

           EpochIndex(i) = sum(ismember(EEG.epoch(i).event, SegIndex));

        end

        EEG = eeg_checkset( EEG );
        EEG = pop_select( EEG, 'trial',find(EpochIndex));
        EEG = eeg_checkset( EEG );

        %%
        
        pathsuj = fullfile(segDir,filesep);
        D = dir(pathsuj);
        parent_folder = D.folder;
        
        if sum(strcmp({D.name},'Conditions')) == 0
        %% 

            mkdir(fullfile(parent_folder,'Conditions'));
            savedir = fullfile(parent_folder,'Conditions');

        else

            disp('Conditions folder already exists');
            savedir = fullfile(parent_folder,'Conditions');
        end

        epoched_name = strcat(EEG.filename(1:end-4),'-',toseg_block,toseg_verbs,'-epoched');

        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(epoched_name),'gui','off');
        EEG = eeg_checkset( EEG );
%% 
         
        EEG = pop_saveset( EEG, 'filename',char(epoched_name),'filepath',savedir);
        EEG = eeg_checkset( EEG );
        eeglab redraw


    otherwise
        disp('Invalid Selection!')
end
