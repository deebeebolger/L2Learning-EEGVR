%function CREx_ERSP_MorletComplex()

%% Script to carry out time-frequency analysis of segmented EEG data.
% The time-frequency analysis is carried out using complex Morlet wavelets.
% The spectral bandwidth is defined by defining the fwhm and the empirical
% fwhm is estimated for each frequency of interest and expressed in time
% and frequency domain.
% The resulting time-frequency map is converted to ERSP by baseline
% correction.
% Baseline correction is carried out by calling the function
% CREX_TF_baseline(). This function gives the possibility of carry
% the following baseline corrections:
% 1. decibel conversion
% 2. Percentage change
% 3. Express post-stimulus activity as z-score in relation to baseline
% interval.
% This script carries out the time-frequency decomposition for a sin
% subject.
% This script carries out baseline decomposition for a single subject. The
% time-frequency matrix (ERSP) is saved as a matfile using the current
% filepath defined for the current dataset.
% The script expects a segmented EEGLAB *.set file.
% The ERSP is plotted for all electrodes of interest.
% Programmed by: D. Bolger                 Date: 28-09-2021
%***************************************************************************
%% Look up the current participant data to load.
Sujindx = 1;
Groups = {'Control', 'Test'};

% Define the paths.
base_fpath = '/Users/bolger/Documents/work/Projects/L2Learn_EEG';
project_path = '/Users/bolger/matlab/Projects/L2Learning_EEGVR_project';

%% Present a dialogue box so that the user can entre the temporal window
% limits (in seconds) over which the data will be averaged.

prompt = {'Enter the min and max in milliseconds of each time window','Name the time windows'};
dlgtitle = 'Define time windows';
dims     = [10 25];
tans     = inputdlg(prompt, dlgtitle, dims);
time_winds = string(tans{1,1});
twind_name = cellstr(tans{2,1});

Twinds = cellfun(@split, time_winds, 'UniformOutput',false);   % stored as a cell array of strings.


%% Open EEGLAB

for gcnt = 1:numel(Groups)

    group_fpath = fullfile(base_fpath, Groups{1,gcnt});
    DGroup = dir(fullfile(group_fpath, '*.set'));
    sujnames = {DGroup.name};


    %% Define the parameters.
    for scounter = 1:numel(Sujindx)

        findsujs = cell2mat(strfind(sujnames, strcat('S',num2str(Sujindx(scounter)),'VB')));
        SIndx = [findsujs>0];

        if isempty(findsujs)

            fprintf('Datasets for subject %d do not exist\n', Sujindx(scounter))

        else

            suj2use =  sujnames(find(SIndx));
            % Want to have the order: VB, VBP
            Odre = [];
            Odre(1) = find(~contains(sujnames, 'VBP'));
            Odre(2) = find(contains(sujnames, 'VBP'));
            sujname_sort = cell(2,1);
            sujname_sort{1,1} = suj2use{1,Odre(1)};
            sujname_sort{2,1} = suj2use{1,Odre(2)};


            [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
            [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

            EEG = pop_loadset('filename',sujname_sort, 'filepath',group_fpath);
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(EEG.setname),'gui','off'); % current set = xx;
            EEG = eeg_checkset( EEG );
            EEG = pop_saveset( EEG, 'filename',char(EEG.setname),'filepath',EEG.filepath, 'savemode', 'resave');  % Saves a copy of the current resampled dataset to the current directory
            eeglab redraw

            for ds_counter = 1:length(ALLEEG)

                DIn = ALLEEG(ds_counter).data;
                T = ALLEEG(ds_counter).times./1000;                        % Change times from ms to seconds.
                chansoi = {'FCz', 'Cz', 'CPz', 'FC3', 'C3', 'CP3', 'FC4', 'C4', 'CP4'};
                allchans = {ALLEEG(ds_counter).chanlocs.labels};
                % chansoi  =  {allchans{1:64}};
                chanindx = find(ismember(allchans, chansoi));
                eindx = chanindx;                               % Electrode indices.
                chans = chansoi;

                %% Set wavelet parameters.

                wavet = -2:1/ALLEEG(ds_counter).srate:2;                          % The wavelet time vector.
                midp = dsearchn(wavet',0);                                                   % Find the midpoint of wavelet length.
                datalen = size(DIn,2);                                                            % data length (length of trial).
                wavlen = size(wavet,2);                                                         % wavelet length.
                nConv = wavlen+datalen-1;                                                    % convolution length.
                wavlen_half = floor(wavlen/2)+1;                                            % half wavelet length.
                freqs = ALLEEG(ds_counter).srate*(0:(datalen/2))/datalen;   % Create the frequency vector.
                foi = freqs(freqs>=4 & freqs<=40);                                           % Define the frequency band of interest.

                % Define the fwhm as a function of frequency of interest.
                fwhm = linspace(0.4, 0.1, numel(foi));

                %% Calculate the complex Morlet wavelet for each center frequency of interest.

                % initialize time intervals of interest in seconds.
                tindx = [T>=-0.15 & T<=1.3];
                blindx = [T>=-0.15 & T<=0];
                Tnew = T(tindx);
                tpostindx = [Tnew>=0 & Tnew<=1.3];

                % Initialize variables
                tf = zeros(length(foi),length(T),size(DIn,3));
                bl = zeros(length(foi),length(find(blindx)),size(DIn,3));
                meanTF = cell(1,numel(eindx));
                meanBL = cell(1, numel(eindx));
                meanTFBL = cell(1, numel(eindx));
                meanTFBL_v2 = cell(1, numel(eindx));
                TFepBL   = cell(1, numel(eindx));

                empfwhm = zeros(length(foi),size(DIn,3),2);  % spectral empirical FWHM
                empfwhmt = zeros(length(foi),size(DIn,3),2);  % temporal empirical FWHM
                tmr = zeros(1,size(DIn,3));

                for ecnt = 1:length(eindx)  % For each electrode

                    for trcnt = 1:size(DIn,3)  % For each trial

                        for fcnt = 1:length(foi)  % For each frequency..ugly coding Dee!

                            Dcurr = DIn(eindx(ecnt),:,trcnt);
                            DIn_fft = fft(Dcurr,nConv, 2);                                                                      % fft of current signal (single trial).

                            % Define the gaussian in the time domain.
                            gwin = exp( (-4*log(2)*wavet.^2) ./ fwhm(fcnt)^2 );                                        % Calculate the gaussian using the current fwhm value.
                            emp_fwhmt(fcnt,trcnt,1) = wavet(midp-1+dsearchn(gwin(midp:end)',.5));      % Calculate the empirical fwhm (time domain)
                            emp_fwhm(fcnt,trcnt,2) = 1/emp_fwhmt(fcnt,1);                                               % Calculate the empirical fwhm (spectral domain)

                            wlet_foi = fft( exp(2*1i*pi*foi(fcnt)*wavet).*gwin,nConv );                              % Calculate the complex Morlet wavelet for frequencies of interest (foi), given FWHM
                            wletfoi_norm = wlet_foi./max(wlet_foi);                                                          % Normalize

                            DInConv = ifft(wletfoi_norm.*DIn_fft);                                                         % Convolve (find ifft of product signal and wavelet).
                            tf_pow = abs(DInConv).^2;                                                                          % Calculate the power
                            tf(fcnt,:,trcnt) = tf_pow(wavlen_half:end-wavlen_half+1);                            % Trim and reshape
                            bl(fcnt,:,trcnt) = tf(fcnt,blindx,trcnt);                                                              % Extract the baseline.

                        end

                    end
                    mean_tfpow = squeeze(mean(tf,3));
                    mean_tfpow = mean_tfpow(:,tindx);
                    meanTF{1,ecnt} = mean_tfpow;
                    
                    % Calculate the mean baseline for the current subject averaged over all trials and conditions (for current subject).
                    mean_blpow = squeeze(mean(bl,3));
                    meanBL{1,ecnt} = mean_blpow;

                    % Call of function to carry out baseline correction. Using the mean of baseline.
                    [tfpow_blc, bltype] = CREX_TF_baseline(mean_tfpow,mean_blpow,'dbels');
                    meanTFBL{1,ecnt} = tfpow_blc;

                    [tfpow_ep_blc, bltype] = CREX_TF_baseline(tf,mean_blpow,'dbels');
                    TFepBL{1,ecnt} = tfpow_ep_blc;
                    meanTFBL_v2{1,ecnt} = squeeze(mean(tfpow_ep_blc,3));

                end

                %% Save the current time-frequency matrix in a structure as a mat-file.

                s = strfind(ALLEEG(ds_counter).setname,'-');
                matname = [ALLEEG(ds_counter).setname(1:s(1)-1),'-chanois-timefreq.mat'];
                eventdata = EEG.event;

                timefreq_results = [];
                timefreq_results.subject = ALLEEG(ds_counter).setname(1:s(1)-1);
                timefreq_results.time  = Tnew;
                timefreq_results.freqs = foi;
                timefreq_results.meanTF = meanTF;
                timefreq_results.ersp = meanTFBL;
                timefreq_results.erspv2 =meanTFBL_v2;
                timefreq_results.baseline = meanBL;
                timefreq_results.TFepochs = TFepBL;
                timefreq_results.bltype = bltype;
                timefreq_results.chansoi = {chansoi};
                timefreq_results.eventinfo = eventdata;

                save(fullfile(ALLEEG(ds_counter).filepath,matname),'timefreq_results'); % Save the time-frequency data to a mat-file in current subject folder.

                %% Carry out continuous wavelet transform (to get the COI (Cone of Influence))

                [wt,f,coi] = cwt(squeeze(DIn(48,:,8)),'amor',ALLEEG(ds_counter).srate);    % Based on single trial for Cz
                hf = figure;
                AX = axes('parent',hf);
                imagesc('Parent',AX,'XData',T,'YData',f,'CData',abs(wt),'CDataMapping','scaled');
                AX.Layer = 'top';
                AX.YDir = 'normal';
                AX.YScale = 'log';
                AX.XLim = [T(1) T(end)];
                hold on
                plot(AX,T,coi,'w--','linewidth',2);

                %% Call of function to write the current time-frequency power values to an excel file.

                %TF2Excel_EEGVR(timefreq_results, Twinds, ds_counter, project_path, base_fpath);

                %% Plot the result of applying the complex Morlet wavelet transform.

                hf1 = figure;
                %     set(hf1,'NumberTitle', 'off', ...
                %         'Name', sprintf('Participant: %s',ALLEEG(scounter).setname(1:s(1)-1)));
                %     scr_siz = get(0,'ScreenSize') ;
                %     pos = floor([scr_siz(3) scr_siz(4) scr_siz(3) scr_siz(4)]);
                %set(hf1,'Position',pos,'Color',[1 1 1]);
                row = 3;
                clim = [-10 10];
                bline_t = T(blindx);

                for ecounter = 1:length(eindx)

                    ax1(ecounter) = subplot(row,ceil(length(chans)/row),ecounter);
                    imagesc('Parent',ax1(ecounter),'XData',T(tindx),'YData',foi,'CData',meanTFBL_v2{1,ecounter},'CDataMapping','scaled')
                    xlim([bline_t(1) Tnew(end)]);
                    ax1(ecounter).YLim = [4 40];                    %Set Y-axis limites
                    ax1(ecounter).CLim = clim;                      %Set time-frequency power limits.
                    ax1(ecounter).Layer = 'top';
                    ax1(ecounter).XLabel.String = 'Time (seconds)'; %Set x-axis label
                    ax1(ecounter).YLabel.String = 'Frequency (Hz)'; %Set y-axis label

                    hold on
                    plot(ax1(ecounter),T,coi,'w--','linewidth',2)
                    title([chans{ecounter},' : ',bltype]);
                    colormap(ax1(ecounter),jet);
                    cb(ecounter) = colorbar;                                 %Colorbar
                    cb(ecounter).Label.String = ['ERSP (',bltype(1:2),' )']; %Colorbar title defined
                    set(ax1(ecounter),'HitTest','on','SelectionHighlight','on','UserData',{T,tindx,foi,meanTFBL_v2{1,ecounter},clim,chans{ecounter},bltype,coi},'Nextplot','replace');
                    set(ax1(ecounter),'ButtonDownFcn',@plotsingle_tf)

                end

            end
        end

    end % end of subject index counter

end % end of Group counter
