%% L2learnEEGVR_EEGPreprocessing_MAIN         Programmed: D. Bolger
% Specific pre-processing script for L2LearnEEGVR project. 
% Script to carry out the different steps in the pre-processing of EEG
% data.
% It applies functions from EEGLAB, the PREP pipeline.
% Functions from the ADJUST toolbox are called in the "CREx_ICA_calc()".
% The following functions are called:
% Needs to be filled in...

%% OPEN THE PARAMETERS TEXTFILE. 
% The parameter textfile contains the necessary pre-processing parameters as well as the paths.

paramfile_nom = 'parameters-L2EEGVR.txt';  % The title of the parameters file. 
% The path to the parameters path. This should be the only path that needs
% to be changed.
paramfile_path = fullfile(filesep,'Users','bolger','Documents','work','Projects','Projet-L2-VREEG','Data',paramfile_nom); % Put the path to your parameters file here. 

fid2 = fopen(paramfile_path);  % Define a file identifier.
mydata = textscan(fid2,'%s %s');  % Scan textfile...

for i = 1:length(mydata{1,1})                     % generate a parameters structure from the parameters text file
    Params.(genvarname(mydata{1,1}{i})) = mydata{1,2}(i);
end

%% OPEN EEGLAB SESSION AND LOAD BDF FILE OF CURRENT SUBJECT.
% A folder is created for the current subject in which all pre-processed
% *.set files will be saved.

pathnom = Params.Basedir{1,1};
data_pathnom = Params.Datadir{1,1};
sujcurr = 'S6MMP';

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;               
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

allfiles= dir(data_pathnom);
filenum=dir(strcat(data_pathnom,'*.bdf'));                      %find all the *.bdf files in the current folder
indx = ismember({filenum(:).name},strcat(sujcurr,'.bdf'));
filenom= filenum(indx).name;
datapath_full = fullfile(data_pathnom,filenom);

[status,error] = mkdir(data_pathnom,sujcurr);     % Create a folder for the processed files of the current subject.
DIRsave_curr = fullfile(data_pathnom,sujcurr,filesep);  % Path to where the processed data will be saved. 

EEG = pop_biosig(datapath_full, 'channels',[1:72], 'ref', [] ,'refoptions',{'keepref' 'off'} );
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(sujcurr),'gui','off'); 
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',char(filenom),'filepath',DIRsave_curr);  % Saves a copy of the current resampled dataset to the current directory
eeglab redraw

%% CALL OF FUNCTION TO ADD TRIGGER INFORMATION THE EEGLAB EEG EVENTS STRUCTURE. 
% This function is specific to L2LearnVR project.
% It completes block number information missing from the EEGLAB data
% structure by consulting the logfiles for each subjects.

fileIn = 'trigs_learning_paradigm.xlsx';  % Excel file defining the individual items and their corresponding triggers.
fileIn_full = fullfile(pathnom,fileIn);
logfilein = 'Subject-logfiles.xlsx';   % Excel file containing log files (from *.dat files) for each subject. There is a subject par sheet.
logfileIn_full = fullfile(pathnom,logfilein);

% Call of function 
EEG = L2LearnVR_trig(EEG, fileIn_full, logfileIn_full, sujcurr);

% Save new event information to a new dataset, with "trigcorr" added. 
fnom_new = strcat(sujcurr,'-trigcorr');
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(fnom_new),'gui','off');   %save the resampled data as a newdata set.
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',char(fnom_new),'filepath',DIRsave_curr);
EEG = eeg_checkset( EEG );
eeglab redraw

%% ADD CHANNEL INFORMATION TO THE CURRENT DATASET.
% Channel coordinates and labels are added to the current dataset and
% the dataset is saved to the current subject-level directory.
% The Chaninfo.mat file is loaded as it contains the electrode labels.
% From EEGLAB plugins, the file, "standard-10-5-cap385.elp" is loaded
% as this contains the correct coordinates for the 10-20 system used here.
chanloc_path = Params.Chanlocdir{1,1};

fnom_chans = strcat(fnom_new,'-chan');
EEG = pop_chanedit(EEG, 'lookup',chanloc_path);                % Load channel path information
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',char(fnom_chans),'filepath',DIRsave_curr);
eeglab redraw

%% PREPARE INFORMATION TEXT-FILE FOR THE CURRENT SUBJECT.
% This textfile should be saved in each subject folder.
    
fname = strcat(sujcurr,'-info.txt');
fdir = strcat(DIRsave_curr,fname);
fid = fopen(fdir,'w');
fprintf(fid,['---------',sujcurr,'----------\n']);

%% RESAMPLE DATA TO THE SAMPLING RATE DEFINED IN THE PARAMETER FILE (IF REQUIRED)
% Resamples using the EEG resampling function.
% If the user has the Matlab signal processing toolbox, it uses the
% Matlab resample() function.
% Write information regarding the resampling to the subject-level text
% file.

SR_orig = EEG.srate;
SR_new = str2double(Params.srate{1,1});

if SR_orig>SR_new

    fprintf(fid,'\nDownsampled from %dHz to %dHz\n',SR_orig,SR_new);
    disp('***********************************Resampling to 512Hz*******************************************')
    fnom_rs = strcat(fnom_chans,'-rs');

    EEG = pop_resample(EEG, SR);   %resample the data at sampling rate defined, sr.
    EEG =eeg_checkset(EEG);
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(fnom_rs),'gui','off'); % current set = xx;
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',char(fnom_rs),'filepath',DIRsave_curr);  % Saves a copy of the current resampled dataset to the current directory
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

f_low = str2double(Params.fc_low{1,1});
f_hi = str2double(Params.fc_hi{1,1});
fnom_filt = strcat(fnom_rs, '-filt');

disp('*********************************Bandpass filtering using a FIR windowed sinc filter***********************************')
[M, dev] = pop_firwsord('blackman',EEG.srate, 2);
[EEG,com,b] = pop_firws(EEG,'fcutoff',[f_low f_hi],'forder',M,'ftype','bandpass','wtype','blackman');
fvtool(b);                                      % Visualise the filter characteristics
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(fnom_filt),'gui','off');   %save the resampled data as a newdata set.
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',char(fnom_filt),'filepath',DIRsave_curr);
eeglab redraw

fprintf(fid,'Band-pass filtered %f2 - %dHz with %d-order fir windowed sinc filter (blackman).\n',f_low,f_hi,M);

%% VISUALISE THE REFERENCE AND THEIR SPECTRA CALCULATED USING MULTI-TAPERS.
% Saves a figure of the spectra of the references to the current
% folder (in the CREx_SpectCalc_multitap() function.

refs = [str2double(Params.references{1,1}(1:2)) str2double(Params.references{1,1}(3:4))];
eegplot(EEG.data(refs,:),'srate',EEG.srate,'eloc_file',EEG.chanlocs(refs(1):refs(2)),'events',EEG.event,'color',{'g' 'b'},'dispchans',2,...
    'winlength',20,'title','Visualise reference electrodes (10sec per window)');

specnom_ref = fullfile(DIRsave_curr,strcat(fnom_filt,'-spectref'));

CREx_SpectCalc_multitap(EEG,65:66,[1 60],specnom_ref,.1);

%% VISUALISE THE SPECTRA OF ALL SCALP ELECTRODES BEFORE RE-REFERENCING.
% The spectrum is saved as a *.fig file. 

specnom_scalp = fullfile(DIRsave_curr,strcat(fnom_filt,'-spectscalp'));

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
EEG = pop_saveset( EEG, 'filename',char(fnom_ref),'filepath',DIRsave_curr);
EEG = eeg_checkset( EEG );
eeglab redraw

here = CURRENTSET;   % Mark the current set.

fprintf(fid,'Rereferenced using channels %s and %s.\n\n',EEG.chanlocs(refs(1)).labels,EEG.chanlocs(refs(2)).labels);

%%  VISUALISE THE SPECTRA OF ALL 64 SCALP SIGNALS CALCULATED USING MULTI-TAPERS.
% Visualise the spectra of all 64 scalp signals after re-referencing.

specnom_scalp = fullfile(DIRsave_curr,strcat(fnom_ref,'-spectscalp'));

CREx_SpectCalc_multitap(EEG,1:64,[1 60],specnom_scalp,.05);

%% RUN FUNCTION TO LOCATE BAD ELECTRODES AND, IF DATA IS SEGMENTED, EPOCHS VISUALLY.
% This works best on segmented data but can also work on continuous also,
% but can be very slow so chose only 10000ms of data.

%EpochChan_dlg(EEG); 

%% CALL OF FUNCTION TO REJECT BAD ELECTRODES
% Remember that it is best to reject the bad electrodes before carrying out
% ICA. 

% CREx_RejBadChans();

%% CALL OF FUNCTION TO CARRY OUT ICA ON CONTINUOUS DATA
% It uses functions from the ADJUST toolbox to detect thos ICs that
% correspond to artifacts.
% It can also apply a PCA before carrying out ICA to reduce the number of components and speed up ICA computation for continuous data.
% The number of PCA components is calculated automatically based on
% explained variance (99% explained variance).

% dopca = 1; % or 0 if not doing PCA before ICA.
% CREx_ICA_calc(dopca);

