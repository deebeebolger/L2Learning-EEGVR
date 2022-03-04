function [blcorr, bltitre] = CREX_TF_baseline(activity, bline, choice)
% Function to carry out baseline correction method of choice on the
% time-frequency data for a single electrode.
% activity => post-stimulus time-frequency matrix (frequencyXtime). This
% could be the mean over all trials or a single trial.
% bline => base-line time-frequency matrix (frequencyXtime). This may be
% the mean over all trials or a single trial.
% choice => the baseline correction approach to apply (string). 
%**************************************************************************
    meanbl = squeeze(mean(bline,2));

    switch choice
        case 'dbels'    % Decibel conversion
            
            blcorr = 10*log10(activity./meanbl);
            bltitre = 'dB conversion';
            
        case 'pcent'    % Percentage change
            
            
            blcorr = 100*((activity - meanbl)./meanbl);
            bltitre = '% change';
            
        case 'zscore'   % Z-transform
            
           numer = activity - meanbl;
           denom1 = (bline - meanbl).^2;
           denom2 = sum(denom1,2)./size(denom1,2);
           denom = sqrt(denom2);
           blcorr = numer./denom;
           bltitre = 'zscore';
        
        otherwise 
            disp('****Not a value baseline correction choice*****')
            blcorr = [];
    end

end