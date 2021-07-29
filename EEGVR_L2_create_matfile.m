%% Date: July 26th 2021       Programmed by: Deirdre Bolger (ILCB)
%  Script to prepare the *.mat file with the pre-processing and analysing
%  parameters for the EEG-VR_LearnL2 project.
%*********************************************************************************

struct = [];

struct.sampling_rate = 512;

struct.filter_params.hp_lim = 0.1;
struct.filter_params.lp_lim = 40;
struct.filter_params.filter_type = 'FIR';
struct.filter_params.filter_window = 'Kaiser';
struct.filter_params.order = 1408;

struct.reference.info = 'average of mastoids';
struct.reference.electrodes_labels = {'mastoidL' 'mastoidR'};
struct.reference.electrodes_index = [67 69];   % Need to verify this.

struct.segmentation.baseline_ms = [- 150 0];
struct.segmentation.poststim_ms = [0 1100];
struct.baseline_correction_ms = [-150 0];

struct.paths.chanconfig = fullfile('');
struct.paths.raw_data = fullfile('');
struct.paths.processed_data = fullfile('');

struct.participant_data.number = [];  % the number of participants
struct.participant_data.groups = {};  % the participant groups
struct.participant_data.ages = [];
struct.participant_data.sex = [];
struct.participant_data.handedness = {};

struct.acquisition_system.name = 'Biosemi Actiview2';
struct.acquisition_system.electrode_num = 64;
struct.acquisition_system.electexterne_labels = {'vertical_eye' 'horizontal_eye' 'left_mastoid' 'right_mastoid'};
struct.acquisition_system.electexterne_indx = [65 66 67 69];

% Electrodes and subjects marked as bad during acquisition
struct.quality_check.bad_electrodes_label = {};
struct.quality_check.bad_electrodes_indx = [];
struct.quality_check.bad_participants_indx = {};

% Information on experimental setup.
struct.experiment.tasks = {'Passive Listening' 'Match-Mismatch'};
struct.experiment.sessions = {'Pre-learning' 'Post-learning'};
struct.experiment.trial_types = {'match' 'mismatch'};
struct.experiment.stimuli_type = 'Verbs';
struct.experiment.stimuli_modality = 'Auditory';

struct.experiment.tasks.Passive_Listening.test_verbs = 101:112;
struct.experiment.tasks.Passive_Listening.test_blocks = 31:33;
struct.experiment.tasks.Passive_Listening.filler_verbs = 125:136;
struct.experiment.tasks.Passive_Listening.filler_blocks = 41:43;
struct.experiment.tasks.Passive_Listening.pre_verb = 96;
struct.experiment.tasks.Passive_Listening.end_trial = 97;
struct.experiment.tasks.Passive_Listening.verb_trialord = {'pre_verb' 'test_verbs' 'test_blocks' 'end_trial'};
struct.experiment.tasks.Passive_Listening.filler_trialord = {'pre_verb' 'filler_verbs' 'filler_blocks' 'end_trial'};

struct.experiment.tasks.Match_Mismatch.match_verbs = 101:112;
struct.experiment.tasks.Match_Mismatch.match_blocks = 31:33;
struct.experiment.tasks.Match_Mismatch.mismatch_verbs = 201:212;
struct.experiment.tasks.Match_Mismatch.mismatch_blocks = 41:43;
struct.experiment.tasks.Match_Mismatch.pre_verb = 95;
struct.experiment.tasks.Match_Mismatch.end_trial = 98;
struct.experiment.tasks.Match_Mismatch.verbmatch_trialord = {'pre_verb' 'match_verbs' 'match_blocks' 'end_trial'};
struct.experiment.tasks.Match_Mismatch.verbmismatch_trialord = {'pre_verb' 'mismatch_verbs' 'mismatch_blocks' 'end_trial'};
struct.experiment.tasks.Match_Mismatch.match_verb_rep = 3;      % Every verb is visualised 3 times as match
struct.experiment.tasks.Match_Mismatch.mismatch_verb_rep = 3;   % Every verb is visualised 3 times as mismatch. 






