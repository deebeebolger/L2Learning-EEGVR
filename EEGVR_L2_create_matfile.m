%% Date: July 26th 2021       Programmed by: Deirdre Bolger (ILCB)
%  Script to prepare the *.mat file with the pre-processing and analysing
%  parameters for the EEG-VR_LearnL2 project.
%*********************************************************************************

function EEGVR_L2_create_matfile()

mystruct = [];
mat_path = fullfile(filesep,'Volumes','deepassport','Projects','Projet-L2-VREEG'); %Path to save mat-file
mat_title = 'EEGVR_L2_Parameters.mat';

%% Generate the mystructure.

mystruct.sampling_rate = 512;

% Define kaiser window filter parameters. 
% The filter order is calculated on the basis of the parameters. 
mystruct.filter_params.hp_lim = 0.1;
mystruct.filter_params.lp_lim = 40;
mystruct.filter_params.filter_type = 'FIR';
mystruct.filter_params.filter_window = 'kaiser';
mystruct.filter_params.band_edges = [0.05 0.1 40 41];
mystruct.filter_params.band_amplitude = [0 1 0];
mystruct.filter_params.dev_min = [0.01 0.05 0.01];

mystruct.reference.info = 'average of mastoids';
mystruct.reference.electrodes_labels = {'mastoidL' 'mastoidR'};
mystruct.reference.electrodes_index = [67 69];   % Need to verify this.

mystruct.segmentation.baseline_ms = [- 150 0];
mystruct.segmentation.poststim_ms = [0 1100];
mystruct.segmentation.baseline_correction_ms = [-150 0];

mystruct.paths.base = fullfile(filesep,'Volumes','deepassport','Projects','Projet-L2-VREEG',filesep);
mystruct.paths.chanconfig = fullfile(filesep,'Volumes','deepassport','Projects','Projet-L2-VREEG','Chaninfo.mat');
mystruct.paths.chaninfo = fullfile(filesep,'Users','bolger','Documents','MATLAB','eeglab2020_0','plugins','dipfit','standard_BESA','standard-10-5-cap385.elp');
mystruct.paths.raw_data = fullfile(filesep,'Volumes','deepassport','Projects','Projet-L2-VREEG','Raw_Data',filesep);
mystruct.paths.processed_data = fullfile(filesep,'Volumes','deepassport','Projects','Projet-L2-VREEG','Processed_Data',filesep);

mystruct.participant_data.number = [];  % the number of participants
mystruct.participant_data.groups = {};  % the participant groups
mystruct.participant_data.ages = [];
mystruct.participant_data.sex = [];
mystruct.participant_data.handedness = {};

%Subject naming system.
mystruct.participant_data.prefixes = {'MM' 'VB' 'MMP' 'VBP'};
mystruct.participant_data.prefixes_long = {'Mismatch_Pretest' 'Verb_Pretest' 'Mismatch_Posttest' 'Verb_Posttest'};

mystruct.acquisition_system.name = 'Biosemi Actiview2';
mystruct.acquisition_system.electrode_num = 64;
mystruct.acquisition_system.electexterne_labels = {'vertical_eye' 'horizontal_eye' 'left_mastoid' 'right_mastoid'};
mystruct.acquisition_system.electexterne_indx = [65 66 67 69];

% Electrodes and subjects marked as bad during acquisition
mystruct.quality_check.bad_electrodes_label = {};
mystruct.quality_check.bad_electrodes_indx = [];
mystruct.quality_check.bad_participants_indx = {};

% Information on experimental setup.
mystruct.experiment.tasks = {'Passive Listening' 'Match-Mismatch'};
mystruct.experiment.sessions = {'Pre-learning' 'Post-learning'};
mystruct.experiment.trial_types = {'match' 'mismatch'};
mystruct.experiment.stimuli_type = 'Verbs';
mystruct.experiment.stimuli_modality = 'Auditory';


mystruct.experiment.Passive_Listening.test_verbs = [101:112];
mystruct.experiment.Passive_Listening.test_blocks = [31:33];
mystruct.experiment.Passive_Listening.filler_verbs = [125:136];
mystruct.experiment.Passive_Listening.filler_blocks = [41:43];
mystruct.experiment.Passive_Listening.pre_verb = 96;
mystruct.experiment.Passive_Listening.end_trial = 97;
mystruct.experiment.Passive_Listening.verb_trialord = {'pre_verb' 'test_verbs' 'test_blocks' 'end_trial'};
mystruct.experiment.Passive_Listening.filler_trialord = {'pre_verb' 'filler_verbs' 'filler_blocks' 'end_trial'};

mystruct.experiment.Match_Mismatch.match_verbs = [101:112];
mystruct.experiment.Match_Mismatch.match_blocks = [31:33];
mystruct.experiment.Match_Mismatch.mismatch_verbs = [201:212];
mystruct.experiment.Match_Mismatch.mismatch_blocks = [41:43];
mystruct.experiment.Match_Mismatch.pre_verb = 95;
mystruct.experiment.Match_Mismatch.end_trial = 98;
mystruct.experiment.Match_Mismatch.verbmatch_trialord = {'pre_verb' 'match_verbs' 'match_blocks' 'end_trial'};
mystruct.experiment.Match_Mismatch.verbmismatch_trialord = {'pre_verb' 'mismatch_verbs' 'mismatch_blocks' 'end_trial'};
mystruct.experiment.Match_Mismatch.match_verb_rep = 3;      % Every verb is visualised 3 times as match
mystruct.experiment.Match_Mismatch.mismatch_verb_rep = 3;   % Every verb is visualised 3 times as mismatch.

mystruct.analysis_info.electrodes_of_interest.midline = {'Fz' 'FCz' 'Cz' 'CPz' 'Pz'};
mystruct.analysis_info.electrodes_of_interest.RH = {'F2' 'F4' 'F6' 'FC2' 'FC4' 'FC6' 'C2' 'C4' 'C6' 'CP2' 'CP4' 'CP6' 'P2' 'P4' 'P6'};
mystruct.analysis_info.electrodes_of_interest.LH = {'F1' 'F3' 'F5' 'FC1' 'FC3' 'FC5' 'C1' 'C3' 'C5' 'CP1' 'CP3' 'CP5' 'P1' 'P3' 'P5'};



%% SAVE THE MAT-File using the path defined at the start.

save(fullfile(mat_path,mat_title),'mystruct');

end
