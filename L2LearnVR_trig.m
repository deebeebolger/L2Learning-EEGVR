function EEG = L2LearnVR_trig(EEG, full_fileIn, full_logfileIn, curr_suj)


XIn = readtable(full_fileIn,'sheet',2, 'ReadVariableNames',1);  %Read in trigger information from excel file as a table.
LogsIn = readtable(full_logfileIn, 'sheet',curr_suj, 'ReadVariableNames',1);

allverbs = [XIn.Match_Verb', XIn.Mismatch_Verb']';
allverb_trigs = [XIn.Match_Verbtrig',XIn.Mismatch_Verbtrig']';

%% Isolate the event list for comparison with the logfile.

data_events = [EEG.event.type]';
codes = LogsIn{:,1};
code_comm = intersect(unique(codes)',unique([EEG.event.type])); % Find the codes that are common to both the log file and the data in EEGLAB

iblockm = find(sum(codes == [XIn.Match_Block1(1), XIn.Match_Block2(1), XIn.Match_Block3(1)],2));
iblockmm = find(sum(codes == [XIn.Mismatch_Block1(1), XIn.Mismatch_Block2(1), XIn.Mismatch_Block3(1)],2));
iall_blocks = sort(cat(1,iblockm,iblockmm));

imVerb_eeg = find(sum(data_events == [XIn.Match_Verbtrig',XIn.Mismatch_Verbtrig'],2)); %Indices of triggers marking matching verbs in EEGLAB events file. 
imVerb_codes = find(sum(codes == [XIn.Match_Verbtrig', XIn.Mismatch_Verbtrig'],2));

[mlen, i] = min([length(imVerb_eeg),length(imVerb_codes)]);
difflen = abs(diff([length(imVerb_eeg),length(imVerb_codes)]));

if i ==2 %Implies that there should be two extra verb trigger codes in the eeg data recorded.
    
    ivmatch = data_events(imVerb_eeg-1) == 95;
    ibadm = find(ivmatch==0);
    imv_corr = imVerb_eeg(find(ivmatch));  % The indices of the correct verb trigger codes.  
    
end


for counter = 1:length(imv_corr)
    
    EEG.event(imv_corr(counter)).verb = allverbs{allverb_trigs == EEG.event(imv_corr(counter)).type};
    EEG.event(imv_corr(counter)).blocknum = codes(iall_blocks(counter));
    
end 

%% Find the correct trigger codes for match and mismatch trials.
% For match trials, correct ==> 100 after a trigger code between 101-112
% For mismatch trials, correct ==> 200 after a trigger code between
% 201-212.
% Incorrect is 100 before 202-212 verb trigger for mismatch and 200 before
% 101-112 for match trials.
% Maybe locate the 98 after the verb and the in/correct trigger follows
% that. 
% Also identifies match and mismatch trials.

icorr_test = [EEG.event(imv_corr+2).type];
% Look for 98s...
ix = find(icorr_test==98);
icorr_test(ix) = [EEG.event(imv_corr(ix)+3).type];

%Compare each in/correct trigger to the verb before
%[EEG.event(imv_corr).type]

T = [[EEG.event(imv_corr).type]', icorr_test'];
vtrig_str = num2str(T(:,1));
rtrig_str = num2str(T(:,2));

for counter2 = 1:length(imv_corr)
    
    vid_trig = EEG.event(imv_corr(counter2)-2).type;
    EEG.event(imv_corr(counter2)-2).verb = strcat('vid-',XIn.Match_Verb{vid_trig});
     
    if strcmp(vtrig_str(counter2,1),'2')
        EEG.event(imv_corr(counter2)).verb = strcat(EEG.event(imv_corr(counter2)).verb,'-mismatch');
    elseif strcmp(vtrig_str(counter2,1),'1')
        EEG.event(imv_corr(counter2)).verb = strcat(EEG.event(imv_corr(counter2)).verb,'-match');
    end
    
    if strcmp(vtrig_str(counter2,1),rtrig_str(counter2,1))
        EEG.event(imv_corr(counter2)).response = 'correct';
        
    elseif ~strcmp(vtrig_str(counter2,1),rtrig_str(counter2,1))
        EEG.event(imv_corr(counter2)).response = 'incorrect';
    end
end


end


