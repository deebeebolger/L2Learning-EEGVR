

%% Open EEGLAB 

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;               
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

EEG = pop_loadset();
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(EEG.setname),'gui','off'); % current set = xx;
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',char(EEG.setname),'filepath',EEG.filepath);  % Saves a copy of the current resampled dataset to the current directory
eeglab redraw

%% Define the parameters.

DIn = EEG.data;
T = EEG.times./1000;                        % Change times from ms to seconds.
eindx = 1:62;                               % Electrode indices.
chans = {EEG.chanlocs.labels};

%% Set wavelet parameters.

wavet = -2:1/EEG.srate:2;                   % The wavelet time vector.
midp = dsearchn(wavet',0);                  % Find the midpoint of wavelet length.
datalen = size(DIn,2);                      % data length (length of trial).
wavlen = size(wavet,2);                     % wavelet length.
nConv = wavlen+datalen-1;                   % convolution length.
wavlen_half = floor(wavlen/2)+1;            % half wavelet length.
freqs = EEG.srate*(0:(datalen/2))/datalen;  % Create the frequency vector.
foi = freqs(freqs>=4 & freqs<=80);          % Define the frequency band of interest.

% Define the fwhm as a function of frequency of interest.
fwhm = linspace(0.4, 0.1, numel(foi));

%% Calculate the complex Morlet wavelet for each center frequency of interest.

% initialize time intervals of interest.
tindx = [T>=-0.5 & T<=1.3];
blindx = [T>=-0.5 & T<=0];
Tnew = T(tindx);
tpostindx = [Tnew>=0 & Tnew<=1.3];

% Initialize variables.
tf = zeros(length(foi),length(T),size(DIn,3));
bl = zeros(length(foi),length(find(blindx)),size(DIn,3));
meanTF = cell(1,numel(eindx));
meanBL = cell(1, numel(eindx));
meanTFBL = cell(1, numel(eindx));

empfwhm = zeros(length(foi),size(DIn,3),2);
tmr = zeros(1,size(DIn,3));

wb = waitbar(0,'Ready...','Name','Calculating time-frequency...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(wb,'canceling',0);

tStart = tic; 
for ecnt = 1:length(eindx)
    tic;
    for trcnt = 1:size(DIn,3)
        
        for fcnt = 1:length(foi)
            Dcurr = DIn(eindx(ecnt),:,trcnt);
            DIn_fft = fft(Dcurr,nConv, 2);
            gwin = exp( (-4*log(2)*wavet.^2) ./ fwhm(fcnt)^2 );               % Calculate the gaussian using the current fwhm value.
            empfwhm(fcnt,trcnt,1) = wavet(midp-1+dsearchn(gwin(midp:end)',.5));     % Calculate the empirical fwhm (time domain)
            empfwhm(fcnt,trcnt,2) = 1/empfwhm(fcnt,1);                              % Calculate the empirical fwhm (spectral domain)
            
            waveX = fft( exp(2*1i*pi*foi(fcnt)*wavet).*gwin,nConv );
            waveX_norm = waveX./max(waveX);                                   % Normalize
            DInConv = ifft(waveX_norm.*DIn_fft);                              % Convolve
            tf_pow = abs(DInConv).^2;                                         % Calculate the power
            tf(fcnt,:,trcnt) = tf_pow(wavlen_half:end-wavlen_half+1);         % Trim and reshape
            bl(fcnt,:,trcnt) = tf(fcnt,blindx,trcnt);                         % Extract the baseline.
            
        end
        
    end
    mean_tfpow = squeeze(mean(tf,3));
    mean_tfpow = mean_tfpow(:,tindx);
    meanTF{1,ecnt} = mean_tfpow;
    
    mean_blpow = squeeze(mean(bl,3));
    meanBL{1,ecnt} = mean_blpow;
    
    % Call of function to carry out baseline correction. Using the mean of
    % baseline.
    [tfpow_blc, bltype] = CREX_TF_baseline(mean_tfpow,mean_blpow,'dbels');
    meanTFBL{1,ecnt} = tfpow_blc;
    
    % Update waitbar and message
    waitbar(ecnt/length(eindx),wb,sprintf('%s',chans{1,ecnt}))
    tmr(trcnt) = toc;
end

ttotal = sum(tmr);
disp(ttotal)
delete(wb)

%% Save the current time-frequency matrix in a structure as a mat-file.

s = strfind(EEG.setname,'-');
matname = [EEG.setname(1:s(1)-1),'-timefreq.mat'];

timefreq_results = [];
timefreq_results.subject = EEG.setname(1:s(1)-1);
timefreq_results.freqs = foi;
timefreq_results.ersp = meanTFBL;
timefreq_results.baseline = meanBL;
timefreq_results.bltype = bltype;
timefreq_results.chansoi = {chans(eindx)};

save(fullfile(EEG.filepath,matname),'timefreq_results'); % Save the time-frequency data to a mat-file in current subject folder.

%% Carry out continuous wavelet transform (to get the COI)

[wt,f,coi] = cwt(squeeze(DIn(48,:,8)),'amor',EEG.srate);
hf = figure;
AX = axes('parent',hf);
imagesc('Parent',AX,'XData',T,'YData',f,'CData',abs(wt),'CDataMapping','scaled');
AX.Layer = 'top';
AX.YDir = 'normal';
AX.YScale = 'log';
AX.XLim = [T(1) T(end)];
hold on
plot(AX,T,coi,'w--','linewidth',2);

%% Plot the result of applying the complex Morlet wavelet transform.

hf1 = figure;
set(hf1,'NumberTitle', 'off', ...
    'Name', sprintf('Participant: %s',EEG.setname(1:s(1)-1)));
scr_siz = get(0,'ScreenSize') ;
pos = floor([scr_siz(3) scr_siz(4) scr_siz(3) scr_siz(4)]);
set(hf1,'Position',pos,'Color',[1 1 1]);
row = 8;
clim = [-8 8];

for ecounter = 1:length(eindx)

    ax1(ecounter) = subplot(row,ceil(64/row),ecounter);
    imagesc('Parent',ax1(ecounter),'XData',T(tindx),'YData',foi,'CData',meanTFBL{1,ecounter},'CDataMapping','scaled')
    xlim([Tnew(1) Tnew(end)]);
    ax1(ecounter).YLim = [4 40];
    ax1(ecounter).CLim = clim;
    ax1(ecounter).Layer = 'top';
    
    hold on
    plot(ax1(ecounter),T,coi,'w--','linewidth',2)
    title([chans{ecounter},' : ',bltype]);
    colormap(ax1(ecounter),jet);
    colorbar;
    set(ax1(ecounter),'HitTest','on','SelectionHighlight','on','UserData',{T,tindx,foi,meanTFBL{1,ecounter},clim,chans{ecounter},bltype,coi},'Nextplot','replace');
    set(ax1(ecounter),'ButtonDownFcn',@plotsingle_tf)

end

