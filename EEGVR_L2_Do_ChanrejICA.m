% Date: 17-08-2021            Programmed by: D. Bolger
% Script to call the channel rejection and ICA functions.
% Each section can be run separately by selecting the individual section
% and running "Run Section". 
% For each operation (channel rejection and ICA), the functions allow you
% to manually choose the dataset that you want to process. 


%% CALL OF FUNCTION TO REJECT BAD ELECTRODES
% Remember that it is best to reject the bad electrodes before carrying out
% ICA. 

CREx_RejBadChans();

%% CALL OF FUNCTION TO CARRY OUT ICA ON CONTINUOUS DATA
% It uses functions from the ADJUST toolbox to detect thos ICs that
% correspond to artifacts.
% It can also apply a PCA before carrying out ICA to reduce the number of components and speed up ICA computation for continuous data.
% The number of PCA components is calculated automatically based on
% explained variance (99% explained variance).

ICAcalc = 1;   % 1 = to calculate the ICA components, 0 = ICA components are already calculated and want to reject components.

if ICAcalc == 1
    dopca = 1; % or 0 if not doing PCA before ICA.
    doICA = 1; % or 0 if not doing ICA calculation (and only rejection of components).
    doRej = 0; % or 1 if only rejecting components and carrying out back-projection of retained components
    doVis = 1; % or 0 if you do not want to visualise and data before and after IC rejection. 
elseif ICAcalc == 0
    dopca = 0; % or 0 if not doing PCA before ICA.
    doICA = 0; % or 0 if not doing ICA calculation (and only rejection of components).
    doRej = 1; % or 1 if only rejecting components and carrying out back-projection of retained components
    doVis = 1; % or 0 if you do not want to visualise and data before and after IC rejection. 
end
        
EEGVR_L2_CREx_ICA_calc(dopca, doICA, doRej,doVis);

