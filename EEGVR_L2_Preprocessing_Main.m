%% EEGVR_L2_Preprocessing_Main        Programmed: D. Bolger
% Specific pre-processing script for L2LearnEEGVR project. 
% Script to carry out the different steps in the pre-processing of EEG
% data.
% It applies functions from EEGLAB, the PREP pipeline.
% Functions from the ADJUST toolbox are called in the "CREx_ICA_calc()".
% The following functions are called:
% Needs to be filled in...

%% CALL OF FUNCTION TO GENERATE THE PARAMETERS MAT-FILE

EEGVR_L2_create_matfile()

%% OPEN THE PARAMETERS MAT-FILE 
% The parameter mat-file contains the necessary pre-processing parameters.

Path_Params = fullfile(filesep,'Volumes','deepassport','Projects','Projet-L2-VREEG');
load(fullfile(Path_Params,'EEGVR_L2_Parameters.mat'));

prefx = mystruct.participant_data.prefixes;

%% Open EEGLAB 

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;               
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

%% LOAD IN THE RAW FILES (*.BDF) INTO CURRENT EEGLAB SESSION
% It will automatically load all the *.bdf files from the "Raw_Data"
% folder.

pathbase = mystruct.paths.base;
pathraw = mystruct.paths.raw_data;
pathproc = mystruct.paths.processed_data;
pathchan = mystruct.paths.chanconfig;
pathchan_info = mystruct.paths.chaninfo;

allfiles= dir(pathraw);
filenum=dir(strcat(pathraw,'*.bdf')); %find all the *.bdf files in the current folder
filenom= {filenum.name};  % the *.bdf file names

%% LOOP THROUGH THE PREPROCESSING STAGES FOR EACH PARTICIPANT.
% Create a folder for each session-type and participant, if necessary.
% Within the script, the following steps will be carried out.
% 1- Reject any electrodes already defined for rejection.
% 2- Resampling (if necessary)
% 3- Re-referencing
% 4- Filtering


for counter = 1:1  %length(filenom)
    
    %% LOAD IN THE CURRENT SUBJECT DATA
    curr_suj = filenom{1, counter};
    EEG = pop_biosig(fullfile(pathraw,curr_suj), 'channels',[1:72], 'ref', [] ,'refoptions',{'keepref' 'off'} );
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(filenom{1,counter}),'gui','off'); 
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    
    % CREATE NEW FOLDER FOR SESSION CONDITIONS IN PROCESSED DATA FOLDER
    ispresent = cellfun(@(s) contains(curr_suj, s),prefx);
    if ~exist(fullfile(pathproc,prefx{ispresent}),'dir')
        [~,~] = mkdir(pathproc,prefx{ispresent}); 
    end
    
    % CREATE NEW FOLDER FOR CURRENT SUBJECT IN PROCESSED DATA FOLDER
    pathcurr = fullfile(pathproc,prefx{ispresent},filesep);
    if ~exist(fullfile(pathcurr,curr_suj(1:end-4)),'dir')
        [~,~] = mkdir(pathcurr,curr_suj(1:end-4)); 
    end
    
    pathsuj = fullfile(pathcurr,curr_suj(1:end-4),filesep);  % The path for the current subject. 
    fnom_raw = curr_suj(1:end-4); 
    
    % SAVE THE CURRENT SUBJECT DATA AS *.SET FILE
    EEG = pop_saveset( EEG, 'filename',char(fnom_raw),'filepath',pathsuj);  % Saves a copy of the current resampled dataset to the current directory
    eeglab redraw
    
    %% ADD CHANNEL INFORMATION TO THE CURRENT DATASET.
    % Channel coordinates and labels are added to the current dataset and
    % the dataset is saved to the current subject-level directory.
    % The Chaninfo.mat file is loaded as it contains the electrode labels.
    % From EEGLAB plugins, the file, "standard-10-5-cap385.elp" is loaded
    % as this contains the correct coordinates for the 10-20 system used here.
    
    fnom_chans = strcat(fnom_raw,'-chan');
    EEG = pop_chanedit(EEG, 'lookup',pathchan_info);                % Load channel path information
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',char(fnom_chans),'filepath',pathsuj);
    eeglab redraw
    
    %% Correct the channels labels if required.
    fnom_chans = EEG.setname;
    X =char({EEG.chanlocs.labels});
    if strcmp(X(1,1:2),'1-')
        labels_corr = X(:,3:end);

        for cntr = 1:length(EEG.chanlocs)

            EEG.chanlocs(cntr).labels = labels_corr(cntr,:);

        end
        
        [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        EEG = eeg_checkset( EEG );
        EEG = pop_saveset( EEG, 'filename',char(fnom_chans),'filepath',pathsuj);
        eeglab redraw

    else
        
       disp('***************Cool! Channels labels are all good********************$')
        
    end
    
    %% PREPARE INFORMATION TEXT-FILE FOR THE CURRENT SUBJECT.
    % This textfile should be saved in each subject folder.

    fname = strcat(curr_suj,'-info.txt');
    fdir = strcat(pathsuj,fname);
    fid = fopen(fdir,'w');
    fprintf(fid,['---------',curr_suj,'----------\n']);
    
    %% RESAMPLE DATA TO THE SAMPLING RATE DEFINED IN THE PARAMETER FILE (IF REQUIRED)
    % Resamples using the EEG resampling function.
    % If the user has the Matlab signal processing toolbox, it uses the
    % Matlab resample() function.
    % Write information regarding the resampling to the subject-level text
    % file.

    SR_orig = EEG.srate;
    SR_new = mystruct.sampling_rate;

    if SR_orig>SR_new

        fprintf(fid,'\nDownsampled from %dHz to %dHz\n',SR_orig,SR_new);
        disp('***********************************Resampling to 512Hz*******************************************')
        fnom_rs = strcat(fnom_chans,'-rs');

        EEG = pop_resample(EEG, SR);   %resample the data at sampling rate defined, sr.
        EEG =eeg_checkset(EEG);
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(fnom_rs),'gui','off'); % current set = xx;
        EEG = eeg_checkset( EEG );
        EEG = pop_saveset( EEG, 'filename',char(fnom_rs),'filepath',pathsuj);  % Saves a copy of the current resampled dataset to the current directory
        eeglab redraw

    elseif SR_orig==SR_new
        fprintf(fid,'\nSampling rate already at %dHz. No need to resample.\n',SR_new);
        disp('***********************************Resampling can be skipped*******************************************');
        fnom_rs = fnom_chans;
    end
    
    %% APPLY BAND-PASS FILTER BETWEEN THE LOWER AND UPPER LIMITS SPECIFIED IN PARAMETERS FILE.
    % It applies a FIR windowed sinc filter using a blackman filter.
    % The filter frequency response is plotted.
    % The details of the filtering are written to subject information txt file.

    f_low = mystruct.filter_params.hp_lim;
    f_hi =  mystruct.filter_params.lp_lim;
    filter_wind = mystruct.filter_params.filter_window;
    fnom_filt = strcat(fnom_rs, '-filt');

    disp('*********************************Bandpass filtering using a FIR windowed sinc filter***********************************\n')
    disp('This could be a bit slow as filter order is quite high!');
    
    fsamp = EEG.srate;
    fcuts = mystruct.filter_params.band_edges;
    mags  =  mystruct.filter_params.band_amplitude ;
    devs  =  mystruct.filter_params.dev_min;

    [M,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,fsamp);   % Calculate filter order and beta (shape factor)
     M = M + rem(M,2);                                      % Filter order must be even. 
    
    [EEG,~,b] = pop_firws(EEG,'fcutoff',[f_low f_hi],'forder',M,'ftype','bandpass','wtype',filter_wind, 'warg',beta);
    fvtool(b);                                                                                               % Visualise the filter characteristics
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(fnom_filt),'gui','off');   % Save the resampled data as a newdata set.
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',char(fnom_filt),'filepath',pathsuj);
    eeglab redraw

    fprintf(fid,'Band-pass filtered %f2 - %dHz with %d-order fir windowed sinc filter (%s).\n',f_low,f_hi,M,filter_wind);
    
    %% VISUALISE THE REFERENCE AND THEIR SPECTRA CALCULATED USING MULTI-TAPERS.
    % Saves a figure of the spectra of the references to the current
    % folder (in the CREx_SpectCalc_multitap() function.

    refs = mystruct.reference.electrodes_index;
    eegplot(EEG.data(refs,:),'srate',EEG.srate,'eloc_file',EEG.chanlocs([refs(1) refs(2)]),'events',EEG.event,'color',{'g' 'b'},'dispchans',2,...
        'winlength',20,'title','Visualise reference electrodes (10sec per window)');

    specnom_ref = fullfile(pathsuj,strcat(fnom_filt,'-spectref'));

    CREx_SpectCalc_multitap(EEG,[refs],[1 60],specnom_ref,.1);
    
    %% VISUALISE THE SPECTRA OF ALL SCALP ELECTRODES BEFORE RE-REFERENCING.
    % The spectrum is saved as a *.fig file. 

    specnom_scalp = fullfile(pathsuj,strcat(fnom_filt,'-spectscalp'));

    CREx_SpectCalc_multitap(EEG,1:64,[1 60],specnom_scalp,.05);
    
    %% RE-REFERENCE THE DATA TO THE ELECTRODES SPECIFIED IN THE PARAMETERS FILE.
    % The channels used for referencing are generally EXG1 and EXG2,
    % channels 65 and 66, respectively.
    % The details of the re-referencing are written to the information text file.

    disp('***********************Rereference to Defined Channel:  does zero potential exist?*****************************')
    EEG = pop_reref(EEG, refs, 'method','standard','keepref','on');
    fnom_ref = strcat(fnom_filt,'-rref');
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(fnom_ref),'gui','off');   %save the resampled data as a newdata set.
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',char(fnom_ref),'filepath',pathsuj);
    EEG = eeg_checkset( EEG );
    eeglab redraw

    refset = CURRENTSET;   % Mark the current set.

    fprintf(fid,'Rereferenced using channels %s and %s.\n\n',EEG.chanlocs(refs(1)).labels,EEG.chanlocs(refs(2)).labels);
    
    %% REJECT ANY ELECTRODES THAT WERE MARKED FOR GENERAL REJECTION DURING DATA ACQUISITION
    
    if ~isempty(mystruct.quality_check.bad_electrodes_indx)
    
        chans2rej = mystruct.quality_check.bad_electrodes_indx;
        fnom_rej1 = strcat(fnom_ref,'-chanrej1');
        EEG=pop_select(EEG,'nochannel',chans2rej);
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(fnom_rej1),'gui','off'); % current set = xx;
        EEG = eeg_checkset( EEG );
        EEG = pop_saveset( EEG, 'filename',char(fnom_rej1),'filepath',pathsuj);  % Saves a copy of the current resampled dataset to the current directory
        eeglab redraw
    else 
        disp('---------------No electrodes marked for rejection at this point------------');
    end
    
    
end







