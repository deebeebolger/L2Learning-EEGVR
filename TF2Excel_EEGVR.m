%%*************************************************************************************************************
% Date: 19-09-2022            Programmed by: d. Bolger
% Script to take load in time-frequency data from the *timefreq.mat files
% and to average TF data over defined frequency bands and time-windows.
% The data needs to be saved in an excel file.
% This is carried out for each participant.

function TF2Excel_EEGVR(timefreq_structIn, time_wind, scount, projpath, datapath)

%% Load in the verb and trigger data from VERB_TRIGS_July2022.xlsx file.


filename = 'VERB_TRIGS_July2022.xlsx';

ds = spreadsheetDatastore(fullfile(projpath, filename));   % this will load in all the data from all the sheets of the excel file.
dataIn = readall(ds);
ttype  = {dataIn.Type};
match1 = find(cell2mat(cellfun(@strcmp, ttype,{'Match'},'Uniformoutput',false)));        % First element of match1 is start of match
mmatch1 = find(cell2mat(cellfun(@strcmp, ttype, {'Mismatch'}, 'Uniformoutput', false))); % First element of mmatch1 is start of mismatchg
vb1 = find(cell2mat(cellfun(@strcmp, ttype, {'Verbs'}, 'UniformOutput',false)));         % First element of vb1 is start of Passive_listening
fill1 = find(cell2mat(cellfun(@strcmp, ttype, {'Fillers'}, 'UniformOutput',false)));     % First element of fill1 is start of Fillers
verb_test = dataIn.Items(vb1);
filler_all = dataIn.Items(fill1);


%%

Groups        = {'Control', 'Test'};
Sessions      = {'Pre_test', 'Post_test'};
Verbtype      = {'Learned', 'Filler'};
freqbands_nom = {'mu', 'beta1', 'beta2'};
blocks        = {'block1', 'block2', 'block3'};
freqbands_hz  = {[8 12], [13 20], [21 32]};

ROInames = {'Midline', 'Lateral_left', 'Lateral_right'};
chanois = {'FCz', 'Cz', 'CPz'; 'FC3', 'C3', 'CP3'; 'FC4','C4', 'CP4'};  % Each row corresponds to each defined ROI.

% Need to find the indices of the channels in the 64-channel context.
chaninfo  = load(fullfile(projpath,'Chaninfo.mat'));
chans_all = {chaninfo.Chaninfo};
chans64    = chans_all{1,1};
ROIs64_indx = cell(size(chanois,1),1);
for roi_cnt = 1:size(chanois,1)
    ROIs64_indx{roi_cnt,1} = find(ismember(chans64, chanois(roi_cnt,:)));
end


%% Set up the structure of the output Excel file
%  Participant; Group; Session; Verb; Block; item; ROI; Electrode; mV

% Calculate the total data length for the current participant.
channum = size(chanois,1)*size(chanois,2);
totallength = length(Verbtype)*length(Sessions)*length(blocks)*length(ROInames)*(length(verb_test)+length(filler_all))*channum*length(freqbands_nom);

datastruct = [];
currsuj    = timefreq_structIn.subject;     % Current participant name to find out group and session
isctrl     = regexp(currsuj, 'C\w*', 'match');    % Test for control group
istest     = regexp(currsuj, 'S\w*', 'match');    % Test for test group.
ispretest  = regexp(currsuj, '\w*B\>','match');  % Parse current participant name to test for pre_test
isposttest = regexp(currsuj, '\w*P\>', 'match'); % Parse current participant name to test for post_test

groupidx = ~[isempty(isctrl), isempty(istest)];
sessidx  = ~[isempty(ispretest), isempty(isposttest)];
currGroup = Groups{find(groupidx)};       % Determine the current group (Control or Test)
currSession = Sessions{sessidx};    % Determine the current session (Pre or Post)

%% Define the group, subject and session columns for the current dataset
% The columns will be of length totallength.

Group = cellstr(repmat(string(currGroup), [totallength, 1]));
Participant = cellstr(repmat(string(currsuj), [totallength, 1]));
Session = cellstr(repmat(string(currSession), [totallength, 1]));

%% Extract the current event data: verb-type, block and items

currevents = timefreq_structIn.eventinfo;
getvb_types = {currevents.eventlabels};
getblocks   = {currevents.blocknum};
getitems    = {currevents.items};

findmt_not   = find(~ismember(getvb_types, ''));  % Find the empty eventlabel entries.
vbtypes_curr = getvb_types(findmt_not);
blocks_curr  = getblocks(findmt_not);
items_curr   = getitems(findmt_not);
datacurr     = timefreq_structIn.TFepochs;
fillerIdx    = find(ismember(vbtypes_curr, 'filler_verb'));   % Index of filler verbs
testvIdx     = find(ismember(vbtypes_curr, 'test_verb'));     % Index of test verbs

%% Create the frequency band column
Freqband = repmat(freqbands_nom', [totallength/length(freqbands_nom'),1]);
ECol = [];
for rcount = 1:length(ROInames)

    elec_curr = reshape(repmat(chanois(rcount,:),[length(freqbands_nom), 1]), [length(freqbands_nom)*size(chanois,2), 1]);
    ECol = cat(1,ECol,elec_curr);

end
Rcurr = reshape(repmat(ROInames, [length(freqbands_nom)*size(chanois,2), 1]), [(length(freqbands_nom)*size(chanois,2))*length(ROInames), 1]);
Rcurr2 = repmat(Rcurr, [length(verb_test),1]);
Rcurr3 = repmat(Rcurr2, [length(filler_all),1]);
RCurr4 = repmat(Rcurr3, [length(blocks),1]);
ROI    = repmat(RCurr4, [length(Verbtype),1]);

ECol2  = repmat(ECol, [length(verb_test),1]);
ECol3  = repmat(ECol2, [length(filler_all),1]);
ECol4  = repmat(ECol3, [length(blocks),1]);
Electrodes = repmat(ECol4, [length(Verbtype),1]);

ICol1 = reshape(repmat(verb_test', [length(Rcurr2),1]),[length(Rcurr2)*length(verb_test),1]);
ICol2 = repmat(ICol1,[length(blocks),1]);
ICol1_fill = reshape(repmat(filler_all', [length(Rcurr2),1]), [length(Rcurr2)*length(filler_all),1]);
ICol2_fill = repmat(ICol1_fill, [length(blocks),1]);
Items = cat(1,ICol2, ICol2_fill);

BCol1 = reshape(repmat(blocks, [length(ICol1_fill), 1]), [length(ICol1_fill)*length(blocks),1]);
Block = repmat(BCol1, [length(Verbtype),1]);

Verbtype = reshape(repmat(Verbtype, [length(BCol1),1]), [length(BCol1)*length(Verbtype),1]);

%% Calculate the power values to assign to each row in the table.
%  Need to find the correct trials to average ==> find verbtype, blocknum, item,
% Note that "Learned" implies "test_verb" in the data.
EventInfo = timefreq_structIn.eventinfo;
Curr_PowData = timefreq_structIn.TFepochs;
PowData   = zeros(length(Verbtype),1);
Time = timefreq_structIn.time;

for bigcounter =1:length(Verbtype)

    fprintf('Find current verbtype, block and item: %s,\t%s,\t%s\n',Block{bigcounter,1}, Verbtype{bigcounter,1}, Items{bigcounter,1});
    if strcmp(Verbtype{bigcounter,1}, 'Learned')
        curr_voi = 'test_verb';
    else
        curr_voi = 'filler_verb';
    end

    vtype_epidx = find(ismember({EventInfo.eventlabels}, curr_voi));
    epochIndx = [EventInfo(vtype_epidx).epoch];
    epochsOI  = ismember({EventInfo(vtype_epidx).items}, Items{bigcounter,1}(1:end-4));
    if sum(epochsOI) ==0
        fprintf('The current item, %s\t, is not in the current dataset.', Items{bigcounter,1})
        PowData(bigcounter,1) = nan;

    else
        %fprintf('Calculating power data for the current item, %s\n', Items{bigcounter,1});
        % Use the data TFepochs field of the timefreq_structIn structure.
        % Check the current electrode to determine which cell to extract.
        % Extract the time-interval defined.
        % Extract the power data for the current frequency-band of
        % interest.
        epochsIndx_curr = find(epochsOI);
        Xdata = find(ismember(timefreq_structIn.chansoi{1,1}, Electrodes(bigcounter,1)));
        assignin('base', "Xdata", Xdata)
        D = Curr_PowData{1, Xdata};   % Extract the data for the current electrode of interest (freq X time X epochs)
        assignin('base','D', D)

        % Find the indices of the time interval of interest.
        twois = str2double(time_wind{1,1})./1000;
        tidx = find(Time >=twois(1) & Time<= twois(2));

        % Find the indices of the current frequency band of interest.
        allfreqs = timefreq_structIn.freqs;
        currfreq = Freqband{bigcounter,1};
        curr_fb =  freqbands_hz{1,strcmp(freqbands_nom, currfreq)};
        fidx    = find(allfreqs>=curr_fb(1) & allfreqs<=curr_fb(2));

        Xcurr1 = D(fidx,tidx,epochsIndx_curr);
        Xmean_ep = squeeze(mean(Xcurr1,3));
        Xmean_t  = squeeze(mean(Xmean_ep,2));
        Xmean_f  = mean(Xmean_t);
        PowData(bigcounter,1) = Xmean_f;
    end
end


%% Assign the data to a table.

powdata_short = arrayfun(@(x)  sprintf('%.2f', x), PowData, 'UniformOutput', false);
currtable = table(Group, Participant, Session, Verbtype, Block, Items, ROI, Electrodes, Freqband, str2double(powdata_short);

%% Need to search for the excel file and if it does not exist, create one.
%  The script will search in the current directory.
fprintf('Checking if the excel file with time-frequency data is in %s\n', datapath);
filelist = dir(fullfile(datapath,'*Timefreq4Stats.xlsx'));
assignin('base','filelist', filelist);
if isempty(filelist)

    % Create the file.
    xlname = inputdlg('Enter the title of the excel file to save time-frequency data (Timefreq4Stats.xlsx)', 'Enter excel name',[1 80], ...
        {'_Timefreq4Stats.xlsx'});
    save_xlname = xlname{1,1};
    writetable(currtable, fullfile(datapath,save_xlname), 'Sheet','AllSubjects', 'WriteVariableNames',true);


else
    % If the file already exists append the currtable to the existing excel
    % file on the sheet defined by "currGroup".
    writetable(currtable, fullfile(datapath,filelist.name), 'Sheet','AllSubjects', 'WriteMode','Append', 'WriteVariableNames', false);

end


end



