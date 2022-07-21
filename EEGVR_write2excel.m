% EEGVR_write2excel()

savedir = cd;

%% Open EEGLAB. User can manually select the conditions

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

EEG = pop_loadset();
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(EEG.setname),'gui','off'); % current set = xx;
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',char(EEG.setname),'filepath',EEG.filepath);  % Saves a copy of the current resampled dataset to the current directory
eeglab redraw

%% Maybe a dialogue box opens in which to mark "Group" and "Session"
test = 1;

while test == 1

    prompt_study = {'Task (Mismatch/Passive Listening)'}
    dlg_title ='';
    deflts = {''};
    num_lignes=[1];
    studyin = inputdlg(prompt_study,dlg_title,num_lignes,deflts);
    options.resize='on';

    if strcmp(studyin{1,1}, 'Mismatch')
        prompt={'Group (Test/Control): ', 'Session (Pre/Post):'};
        dlg_title ='';
        deflts = {'',''};
        num_lignes=[1;1];
        infoIn = inputdlg(prompt,dlg_title,num_lignes,deflts);
        options.resize='on';
        test = 0;
    elseif strcmp(studyin{1,1}, 'Passive Listening')
        prompt={'Group (Test/Control): ', 'Session (Pre/Post):'};
        dlg_title ='';
        deflts = {'',''};
        num_lignes=[1;1];
        infoIn = inputdlg(prompt,dlg_title,num_lignes,deflts);
        options.resize='on';
        test = 0;
    else
        fprintf('Invalid entry, try again')
        test = 1;

    end
end


%% If it is Match-Mismatch Task data

if strcmp(studyin{1,1}, 'Mismatch')

    % Create a table for the current subject
    condtable = table;

    Group = infoIn{1,1};
    Session = infoIn{2,1};

    %% Determine if the current subject has match or mismatch conditions
    %  Need to loop through the selected subjects

    % Initialize average over ROI variables as there will be one value per
    % participant.
    AvgROI_midline         = zeros(length(ALLEEG),1);
    AvgROI_LH              = zeros(length(ALLEEG),1);
    AvgROI_RH              = zeros(length(ALLEEG),1);
    curr_participant       = cell(length(ALLEEG),1);
    group                  = cell(length(ALLEEG),1);
    session                = cell(length(ALLEEG),1);
    condition              = cell(length(ALLEEG),1);


    for sujcnt = 1:length(ALLEEG)

        % Retrieve the current dataset.
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve', sujcnt,'gui','off'); % current set = xx;
        EEG = eeg_checkset( EEG );
        EEG = pop_saveset( EEG, 'filename',char(EEG.setname),'filepath',EEG.filepath);  % Saves a copy of the current resampled dataset to the current directory
        eeglab redraw




        % Extract the trials corresponding to each condition
        mismatch_indx = zeros(length(EEG.epoch),1);
        match_indx = zeros(length(EEG.epoch),1);
        for epcnt = 1:length(EEG.epoch)
            currepoch = (EEG.epoch(epcnt).eventeventlabels);
            if sum(strcmp(currepoch, 'mismatch_verb'))>0
                mismatch_indx(epcnt) = epcnt;
            elseif sum(strcmp(currepoch, 'match_verb'))>0
                match_indx(epcnt) = epcnt;
            end
        end

        MMindx = mismatch_indx([mismatch_indx >0],1);  % Trial index for mismatch condition
        Mindx  = match_indx([match_indx >0],1);        % Trial index for match condition

        %% Find the average over Match and Mismatch trials

        if sujcnt == 1
            ROI = [];
            ROI.midline = {'Fz', 'FCz', 'Cz', 'CPz', 'Pz'};
            ROI.leftH   = {'F1', 'F3', 'F5', 'FC1', 'FC3', 'FC5', 'C1','C3', 'C5', 'CP1', 'CP3', 'CP5', 'P1', 'P3', 'P5'};
            ROI.rightH  = {'F2', 'F4', 'F6', 'FC2', 'FC4', 'FC6', 'C2', 'C4', 'C6', 'CP2', 'CP4', 'CP6', 'P2', 'P4', 'P6'};

            AllChans = {EEG.chanlocs.labels};
            midline_indx = find(ismember(AllChans, ROI.midline));
            leftH_indx = find(ismember(AllChans, ROI.leftH));
            rightH_indx = find(ismember(AllChans,  ROI.rightH));

            time_wind = [300 600];

        end % end of if sujcnt

        time = EEG.times;
        twind_indx = find([time>= time_wind(1) & time<= time_wind(2)]);

        AllData = EEG.data;

        if ~isempty(Mindx)

            condition{sujcnt,1} = 'Match';           % Assign a condition name
            matchdata_trial = AllData(:,:,Mindx);
            matchAvg = squeeze(mean(matchdata_trial, 3));

            % Find average over each spatial-temporal region of interest
            % Midline ROI
            matchAvg_midline         = matchAvg(midline_indx, twind_indx);
            matchAvg_midchan         = squeeze(mean(matchAvg_midline,1));  % Average over the channels
            matchAvg_ROImidline      = mean(matchAvg_midchan);
            AvgROI_midline(sujcnt,1) = matchAvg_ROImidline;

            % Left-hemisphere ROI
            matchAvg_leftH   = matchAvg(leftH_indx, twind_indx);
            matchAvg_LHchan  = squeeze(mean(matchAvg_leftH,1));  % Average over the channels
            matchAvg_ROILH  = mean(matchAvg_LHchan);
            AvgROI_LH(sujcnt,1) = matchAvg_ROILH;


            % Right-hemisphere ROI
            matchAvg_rightH   = matchAvg(rightH_indx, twind_indx);
            matchAvg_RHchan  = squeeze(mean(matchAvg_rightH,1));  % Average over the channels
            matchAvg_ROIRH  = mean(matchAvg_RHchan);
            AvgROI_RH(sujcnt,1) = matchAvg_ROIRH;

        end

        if ~isempty(MMindx)

            condition{sujcnt,1} = 'Mismatch';         % Assign a condition name
            mismatch_trial = AllData(:,:, MMindx);
            mismatchAvg = squeeze(mean(mismatch_trial, 3));

            % Find average over each spatial-temporal region of interest
            % Midline ROI
            mismatchAvg_midline    = mismatchAvg(midline_indx, twind_indx);
            mismatchAvg_midchan    = squeeze(mean(mismatchAvg_midline,1));  % Average over the channels
            mismatchAvg_ROImidline = mean(mismatchAvg_midchan);
            AvgROI_midline(sujcnt,1) = mismatchAvg_ROImidline;


            % Left-hemisphere ROI
            mismatchAvg_leftH   = mismatchAvg(leftH_indx, twind_indx);
            mismatchAvg_LHchan  = squeeze(mean(mismatchAvg_leftH,1));  % Average over the channels
            mismatchAvg_ROILH   = mean(mismatchAvg_LHchan);
            AvgROI_LH(sujcnt,1) = mismatchAvg_ROILH;

            % Right-hemisphere ROI
            mismatchAvg_rightH   = mismatchAvg(rightH_indx, twind_indx);
            mismatchAvg_RHchan  = squeeze(mean(mismatchAvg_rightH,1));  % Average over the channels
            mismatchAvg_ROIRH   = mean(mismatchAvg_RHchan);
            AvgROI_RH(sujcnt,1) = mismatchAvg_ROIRH;

        end

        curr_participant{sujcnt,1} = EEG.setname;
        group{sujcnt,1}   = Group;
        session{sujcnt,1} = Session;

    end

% Assign data to the table
condtable.Participant = curr_participant;
condtable.Group       = group;
condtable.Session     = session;
condtable.Condition   = condition
condtable.ROImidline  = AvgROI_midline;
condtable.ROILH       = AvgROI_LH;
condtable.ROIRH       = AvgROI_LH;

% Save the table contents in an excel file.
mmfilename = ['MisMatch_',Group,'_',Session,'.xlsx'];
writetable(condtable, fullfile(savedir,mmfilename), 'Sheet', 1);




elseif strcmp(studyin{1,1}, 'Passive Listening')

end



%%




