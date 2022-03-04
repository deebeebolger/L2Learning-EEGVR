

data2seg = EEG.data;
% The events to use as the segmentation reference are expressed as samples
% (latency)

toseg = 'match_verb';
toseg_indx = find(strcmp({EEG.event.eventlabels}, toseg));
eventsamps = [EEG.event(toseg_indx).latency];
timelims = [-1.1 1.9];

[epocheddata, newtime, indices, rerefevent, rereflatencies ] = epoch(data2seg, toseg_indx,timelims, 'srate',512);

[outEEG, indxs ] = pop_epoch(EEG,{},timelims, 'eventindices',toseg_indx);

s = strfind(EEG.setname,'-');
epoch_name = strcat(EEG.setname(1:s(1)),toseg,'-','epoch');
pathsuj = EEG.filepath;
EEG = pop_saveset( outEEG, 'filename',char(epoch_name),'filepath',pathsuj);  % Saves a copy of the current resampled dataset to the current directory
eeglab redraw

%%

x = size(EEG.data,2);
Fs = EEG.srate;

fb = cwtfilterbank('SignalLength',x,'SamplingFrequency',Fs,...
    'FrequencyLimits',[4 80]);
figure
freqz(fb)

[wt, f, coi] = cwt(EEG.data(48,:,1),'FilterBank',fb);

