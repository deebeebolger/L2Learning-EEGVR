"""Script to calculate the ERSP for segmented data files from EEGLAB
    Also carries out statistical analysis.
    Date: April 2023
    Programmed by: Deirdre Bolger
"""
import mne
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from mne.stats import permutation_cluster_test, permutation_cluster_1samp_test
import glob
import os
import pandas as pd
import seaborn as sns
import math
import pickle5 as pkl5

def define_triggers(trigData, epoch_curr, titlecurr):


    print(f"The current file name is {titlecurr}")
    eventID = epoch_curr.event_id
    trigcodes = list(eventID.keys())
    valueID   = list(eventID.values())
    allevents = epoch_curr.events
    Trigs = list(trigData['Triggers'])
    Items = list(trigData['Items'])
    Verbs = list(trigData['Type'])

    # Check title of current file and determine of pre- or post-test
    if 'VB_' in titlecurr:
        test_type = '/pretest'
    elif 'VBP_' in titlecurr:
        test_type = '/posttest'
    elif 'VB _' in titlecurr:
        test_type = '/pretest'
    elif 'VBP _' in titlecurr:
        test_type = '/posttest'

    VItem = [vindx+'/'+iindx for vindx, iindx in zip(Verbs, Items)]
    # test-type to the key of each item in the dictionary.
    VItem2 = [vi +test_type for vi in VItem]
    EventDict = dict(zip(VItem2, Trigs))

    return EventDict, test_type


def createEpochs(fullfiles_In, currentcond_path, savepath_cond, AllEpochs, trigInfo_In):

    CIndx = []
    posttest_indx = []
    pretest_indx = []
    for count, Fcurr in enumerate(fullfiles_In):
        Epochcurr = []
        evindx = []
        Epochcurr = mne.read_epochs_eeglab(Fcurr)
        curr_events = Epochcurr.event_id
        # Extract name of current dataset and check if pre- or post-test
        fname_split = Fcurr.split("/")
        fnom = fname_split[-1:]
        fnom = fnom[0]  #convert to string


        # Load in the triggers textfile for the current subject.
        txtname_split = Fcurr.split("/")
        txtname_curr = 'Trigs_' + txtname_split[-1][0:-4] + '.txt'
        currpath_txt = os.path.join(currentcond_path, txtname_curr)
        trigdataIn = pd.read_csv(currpath_txt, sep=",", header=0)
        Items = list(trigdataIn['items'].values)
        Nindices = [Nindx for Nindx, N in enumerate(Items) if type(N) == str]

        # Recreate the events array
        col2 = np.squeeze(np.zeros((len(trigdataIn), 1)))
        col1 = trigdataIn['latency'].values
        col3 = trigdataIn['type'].values
        EventsCurr = np.column_stack([col1[Nindices].astype(int), col2[Nindices].astype(int), col3[Nindices].astype(int)])
        Epochcurr.events = EventsCurr

        new_eventID, ttype = define_triggers(trigInfo_In, Epochcurr, fnom) # Call of function.
        ev = list(new_eventID.values())
        kev = list(new_eventID.keys())
        missing_events = set(np.unique(ev))-set(np.unique(EventsCurr[:,2]))
        print(f'Those events missing from data are {missing_events}\n')
        # Find indices of the missing events
        I = [idx for idx, Evnts in enumerate(ev) if Evnts in missing_events]
        for ev_counter in I:
            print(f'Removing key item {kev[ev_counter]}\n')
            new_eventID.pop(kev[ev_counter], None)

        Epochcurr.event_id = new_eventID
        Epochcurr.selection = np.arange(0, len(EventsCurr))
        Epochcurr.drop_log = tuple([] for _ in range(len(EventsCurr)))
        #Evokedcurr = Epochcurr.copy().average()
        AllEpochs.append(Epochcurr)
        #AllEvoked.append(Evokedcurr)
        CIndx.append(evindx)

        if ttype == '/pretest':
            pretest_indx.append(count)
        elif ttype == '/posttest':
            posttest_indx.append(count)

    return AllEpochs, EventsCurr, pretest_indx, posttest_indx

def doPadding(epochdataIn, padur, srate, pret, postt):
    """
    Function to carry out padding at the beginning and the end of each trial of the input epoched data.
    The input epoched data is a list of EEGLAB Epoched objects for all participants for a given condition.
    :param epochdataIn: List of EEGLAB epoch objects for all participants for a given condition.
    :param padur: The padding duration in seconds.
    :param srate: sampling rate
    :return:
    """
    pad_nsamps = int(srate*padur)
    # recreate the time vector.
    T = epochdataIn[0].times
    bline = np.arange(abs(T[0])+padur, 0, (1/srate)*-1)
    bline = bline*-1
    postim = np.arange(0, T[-1]+padur, 1/srate)
    Tnew = np.concatenate((bline, postim), axis=0)
    NewEpochs = []

    for count, currepoch in enumerate(epochdataIn):
        # Get the epochs object data to create new epoch array
        print(count)
        if count in pret:
            print('The current dataset is a pre-test')
            getdata_cond = 'pretest'
        elif count in postt:
            print('The current dataset is a post-test.')
            getdata_cond = 'posttest'

        new_info = mne.create_info(ch_names=currepoch.info['ch_names'], sfreq=currepoch.info['sfreq'])
        new_events = epochdataIn[count].events
        new_event_id = epochdataIn[count].event_id
        chansoi = epochdataIn[count].ch_names

        # Get the 3D data array.
        currdata = epochdataIn[count].get_data(chansoi, getdata_cond)   # Get the data matrix for the current epoch object.
        data_shape = np.shape(currdata)   # Get all epochs of the current epochs object as a 3D array - trials X channels X time
        print(data_shape)

        # Create channels X pad_nsamps array to use for padding
        new_epoch = []
        Enew = []
        for trialcnt in range(0,data_shape[0]):
            print(trialcnt)
            data_start = currdata[trialcnt, :, 0]
            data_end = currdata[trialcnt, :, data_shape[2]-1]
            startpad = np.tile(data_start,pad_nsamps).reshape(data_shape[1], pad_nsamps)   # To concatenate to start of each trial
            endpad = np.tile(data_end, pad_nsamps).reshape(data_shape[1], pad_nsamps)      # To concatenate to end of each trial
            trialcurr = currdata[trialcnt, :,:]
            trialcurr_ex1 = np.concatenate((startpad, trialcurr), axis=1)
            trialcurr_ex2 = np.concatenate((trialcurr_ex1, endpad), axis=1)
            trialcurr_ex2 = trialcurr_ex2[:,:-1]
            new_epoch.append(trialcurr_ex2)

        # Create a new epochs array
        Epochnew = np.array(new_epoch)  # Convert the list of epochs to a 3D array.
        Enew = mne.EpochsArray(Epochnew, info=new_info, events=new_events, tmin=Tnew[0], event_id=new_event_id)
        NewEpochs.append(Enew)

    return NewEpochs


def doTFAnalysis(freqs, cyclesN, chan_picks, EpochsIn, tlims):
    """

    :param freqs:
    :param cyclesN:
    :param chan_picks:
    :param EpochsIn:
    :param Tshort:
    :param tlims:
    :return:
    """

    ERSPfillers_all = []
    ERSPverbs_all = []
    BLfill = []
    for counter, currep in enumerate(EpochsIn):
        TFRcurr = mne.time_frequency.tfr_morlet(EpochsIn[counter], freqs=freqs, n_cycles=cyclesN, use_fft=False,
                                                return_itc=False, average=False,
                                                decim=1, picks=chan_picks)
        # Calculate the average baseline
        currTF_data = TFRcurr.data  # Extract data from the TFRcurr object: trials X chans X freqs X time
        currrTF_short = currTF_data[:,: ,: ,tlims[0]:tlims[1]]  # Get the real data from the original time vector. (-150ms to 1000ms)

        TFRcurr_avg = TFRcurr.average()  # Calculate the average over all epochs for the current subject
        currTFavg_data = TFRcurr_avg.data  # Extract data from the average TFR object: chans X freqs X time
        currTFavg_short  = currTFavg_data[:,:,tlims[0]:tlims[1]]

        print(f'The size the trial-level TF data is {np.shape(currrTF_short)}')
        print(f'The size the trial-average TF data is {np.shape(currTFavg_short)}')

        #Now extract the baseline interval.
        T = TFRcurr.times
        baseline_lims = [T[tlims[0]], 0]
        baseline_indx = TFRcurr.time_as_index(baseline_lims)  # Get the indices of the baseline time limits.
        avg_baseline = currTFavg_data[:, :, baseline_indx[0]:baseline_indx[1]] # This is average over fillers and test verbs.

        print(f'The baseline limits are {baseline_lims} in seconds')
        print(f'The real limits of the baseline are {T[tlims[0]]}ms and 0')
        print(f'The baseline min and max of the baseline indices are {baseline_indx[0]} and {baseline_indx[-1]}')
        print(f'The shape of the average baseline is {np.shape(avg_baseline)}')

        # Extract the trials corresponding to conditions (Fillers, Verbs)
        TFdata_fillers = TFRcurr['Fillers'].data
        TFdata_verbs = TFRcurr['Verbs'].data
        # Get the data over the real trial range as defined by tlims.
        TFdata_fillshort = TFdata_fillers[:,:,:,tlims[0]:tlims[1]]
        TFdata_verbshort = TFdata_verbs[:,:,:,tlims[0]:tlims[1]]

        dshape_f = np.shape(TFdata_fillshort)
        dshape_v = np.shape(TFdata_verbshort)
        print(f'The shape of the filler TF data is {dshape_f}')
        print(f'The shape of the verb TF data is {dshape_v}')

        # Set up for loop to loop through the channels
        ERSP_fillers = np.empty([dshape_f[1], dshape_f[0], dshape_f[2], dshape_f[3]])
        ERSP_verbs = np.empty([dshape_v[1], dshape_v[0], dshape_v[2], dshape_v[3]])

        print(f'The shape of the ERSP_fillers variable is {np.shape(ERSP_fillers)}')
        print(f'The shape of the ERSP verbs variable is {np.shape(ERSP_verbs)}')

        # Call of a fuction to carry out baseline correction at individual trial level.
        bline_type = 'dB'
        for cntr, chanidx in enumerate(chan_picks):
            print(cntr)
            ERSP_fillers[cntr, :, :, :], BLcurr_filler = doBaselineCorrect(TFdata_fillshort, avg_baseline, cntr, bline_type, dshape_f)
            ERSP_verbs[cntr, :, :, :], BLcurr_verbs = doBaselineCorrect(TFdata_verbshort, avg_baseline, cntr, bline_type, dshape_v)

            ERSPfillers_mean = np.mean(ERSP_fillers, 1)  # Average over epochs for current channel
            ERSPverbs_mean = np.mean(ERSP_verbs, 1) # Average over epochs for current channel

        ERSPfillers_all.append(ERSPfillers_mean)
        ERSPverbs_all.append(ERSPverbs_mean)
        BLfill.append(BLcurr_filler)

    return ERSPfillers_all, ERSPverbs_all, TFRcurr

def doBaselineCorrect(TFdataIn, BLavg, counter, BLtype, dshapeIn):
    """
    Function to carry out trial-level baseline correction using the average
    baseline. Average over all conditions for the participant.
    :param TFdataIn: data for current channel: trials X freqs X time
    :param BLavg: baseline for the current channel: freqs X time
    :param counter: counter for current channel index
    :param BLtype: defines the baseline correction type (pcent= percentage, dB=10*log10)
    :param dshapeIn: the shape of the current TF data.
    :return: erspArray: the numpy array with baseline corrected data for current condition/subject.
    """
    erspArray = np.empty([dshapeIn[0], dshapeIn[2], dshapeIn[3]])
    X = np.squeeze(TFdataIn[:, counter, :, :])  # data for current channel: trials X freqs X time
    bl = np.squeeze(BLavg[counter, :, :])  # baseline for the current channel: freqs X time
    meanbl = np.squeeze(np.mean(bl, 1))
    stdBL = np.std(bl,1)

    print(f'The size of the current data is {np.shape(X)}')
    print(f'The size the baseline data is {np.shape(bl)}')
    print(f'The size of the mean baseline is {np.shape(meanbl)}')
    print(f'The size of the standard deviation of the baseline frequencies is {np.shape(stdBL)}')

    for epcntr in range(0, dshapeIn[0]):

        if BLtype == 'zscore':
            ersp_curr = np.squeeze(X[epcntr, :, :].T) - meanbl
            ersp_bl = ersp_curr/stdBL
        elif BLtype == 'pcent':  # percentage (based on the gain model)
            ersp_curr = np.squeeze(X[epcntr, :, :].T)/meanbl
            ersp_bl = ersp_curr * 100
        elif BLtype == "dB":  #log ratio (based on gain model)
            ersp_curr = np.squeeze(X[epcntr, :, :].T)/meanbl
            #ersp_curr = ersp_curr*100
            ersp_log = np.log10(ersp_curr)*10
            ersp_bl = ersp_log

        erspArray[epcntr, :, :] = ersp_bl.T
        return erspArray, meanbl

def findSigERSP(chanois, freqs, time, CondIn, alpha_value, Vmin, Vmax, clevels, n_perms, tail, threshT):
    """
    Function to find significant ERSP at the individual channel level using a cluster-based permutation test.
    The function applied for the significance testing is a non-parametric cluster-level paired t-test and FDR correction is applied.
    The input to the statistical test is the ERSP expressed either as percentage or dBs.



    :return:
    """

    fig, axes = plt.subplots(3, 3)
    Cntr = 0
    clusterpvals = []
    Tvals = []
    ttstats_all = []
    for ecount, ax in enumerate(axes):
        for ecount2, ax2 in enumerate(ax):
            print(Cntr)
            mask = []
            CondIn_erspv = np.squeeze(CondIn[0])
            CondIn_erspf = np.squeeze(CondIn[1])
            condInERSP = [np.squeeze(CondIn_erspv[:,Cntr,:,:]), np.squeeze(CondIn_erspf[:,Cntr,:,:])]
            CondIn_GAv = np.squeeze(CondIn_erspv[:,Cntr,:,:]).mean(0)
            CondIn_GAf = np.squeeze(CondIn_erspf[:,Cntr,:,:]).mean(0)
            CondIn1 = np.squeeze(CondIn_erspv[:, Cntr, :, :])
            CondIn2 = np.squeeze(CondIn_erspf[:, Cntr, :, :])
            CondIn2test = CondIn2-CondIn1
            CondIn_diff = CondIn_GAf-CondIn_GAv   # Condition 2 - Condition 1

            print(f'The shape of CondIn1 is {np.shape(CondIn1)}')
            CondIn1_dstack = np.dstack(CondIn1)
            CondIn2_dstack = np.dstack(CondIn2)
            print(f'The shape of CondIn1 dstacked is {np.shape(CondIn1_dstack)}')

            # Plot the time-frequency results
            T_obs, clusters, cluster_p_values, H0 = permutation_cluster_1samp_test(CondIn2test, threshold=threshT, n_permutations=n_perms,
                                               tail=0, out_type='mask', verbose=True)
            clusterpvals.append(cluster_p_values)
            Tvals.append(T_obs)
            print(np.shape(T_obs))

            ttest_results = sp.stats.ttest_rel(CondIn1_dstack, CondIn2_dstack, axis=2, nan_policy='propagate', alternative='two-sided')
            tt_stats = ttest_results.statistic
            pvals = ttest_results.pvalue
            dof   = ttest_results.df
            print(f'The shape of the ttest-statistic is {np.shape(tt_stats)}')
            ttstats_all.append(tt_stats)

            T_obs_plot = np.nan * np.ones_like(T_obs)
            for c, p_val in zip(clusters, cluster_p_values):
                if p_val <= alpha_value:
                    print('Significant data to plot')
                    T_obs_plot[c] = T_obs[c]

            # f_idx, t_idx = np.unravel_index(np.nanargmax(np.abs(T_obs_plot)), CondIn_GA.shape[1:])
            vmax2 = np.max(np.abs(T_obs))
            vmin2 = -vmax2
            vmax1 = np.max(np.abs(CondIn_diff))
            vmin1 = -vmax1
            ax2.imshow(CondIn_diff, cmap=plt.cm.RdBu_r,
                       extent=[time[0], time[-1], freqs[0], freqs[-1]],
                       aspect='auto', origin='lower', vmin=vmin1, vmax=vmax1, alpha=0.5)
            ax2.imshow(T_obs_plot, cmap=plt.cm.coolwarm,
                       extent=[time[0], time[-1], freqs[0], freqs[-1]],
                       aspect='auto', origin='lower', vmin=vmin2, vmax=vmax2, alpha=1)
            if ecount ==len(axes)-1:
                ax2.set_xlabel('Time (ms)')
            else:
                ax2.set_xlabel(' ')

            if ecount2==0:
                ax2.set_ylabel('Frequency (Hz)')
            else:
                ax2.set_ylabel(' ')
            ax2.axvline(0, linewidth=1, color="black", linestyle=":")
            ax2.set_title(chanois[Cntr], fontsize=10)
            Cntr+=1

    return clusterpvals, Tvals, ttstats_all

def Tvals_save(path, filename, tdata, time, freq, Ttestvals, filenom_ttest):
    """
    Function to save the Tvalue array to textfile.

    :param path: path in which to save current textfile.
    :param filename: Title of text file (string)
    :param tdata: Arrays of t-values (List of numpy arrays)
    :param Ttestvals: Array of t-values from paired t-test.
    :return:
    """
    fullpath = os.path.join(path, filename)
    np.savetxt(fullpath, tdata)

    fullpath_ttest = os.path.join(path, filenom_ttest)
    np.savetxt(fullpath_ttest, Ttestvals)

    # save the time and frequency vectors in a matfile
    tfinfo_dict = {'Time' : time, 'Frequency' : freq}
    mat_nom = 'TF_time_freq.mat'
    matpath_full = os.path.join(path, mat_nom)
    sp.io.savemat(matpath_full, {'tfinfo_dict' : tfinfo_dict})


path_triginfo = '/Users/bolger/Documents/work/Projects/L2learn_EEG/VERB_TRIGs_July2022.xlsx'
path_to_data = '/Users/bolger/Documents/work/Projects/L2learn_EEG/Data_for_TF/HighPerformers'
groups = ['Control', 'Test']
path_control = os.path.join(path_to_data, groups[0])
path_test = os.path.join(path_to_data, groups[1])
savepath_control = os.path.join(path_control, 'mne_epoched')
savepath_test = os.path.join(path_test, 'mne_epoched')

df = pd.read_excel(path_triginfo, sheet_name='Passive_Listening')
trigData = pd.DataFrame(df, columns=['Items', 'Triggers', 'Type'])

controlF = os.listdir(path_control)
testF = os.listdir(path_test)

control_files = [x for x in controlF if x.endswith('.set')]
test_files = [x for x in testF if x.endswith('.set')]

controlF_full = [os.path.join(path_control, f) for f in control_files]
testF_full = [os.path.join(path_test, f) for f in test_files]
TFAll = []
EpochAll_control = []
EpochAll_test = []
EvokedAll_control = []
EvokedAll_test = []

kwargs = dict(n_permutations=1000, step_down_p=0.05, seed=1,
              buffer_size=None, out_type='mask')  # for cluster test

EpochAll_control, control_events, pretestIC, posttestIC  = createEpochs(controlF_full, path_control, savepath_control, EpochAll_control, trigData)
EpochAll_test, test_events, pretestIT, posttestIT = createEpochs(testF_full, path_test, savepath_test, EpochAll_test, trigData)

paddur=1.0
sfreq = EpochAll_control[0].info['sfreq']
EpochAllC_ext = doPadding(EpochAll_control, paddur, sfreq, pretestIC, posttestIC)
EpochAllT_ext = doPadding(EpochAll_test, paddur, sfreq, pretestIT, posttestIT)

EpochControl_pre = [EpochAllC_ext[preI] for preI in pretestIC]
EpochControl_post = [EpochAllC_ext[postI] for postI in posttestIC]
EpochTest_pre = [EpochAllT_ext[preIt] for preIt in pretestIT]
EpochTest_post = [EpochAllT_ext[postIt] for postIt in posttestIT]

## Carry out time frequency analysis using wavelets.
Freqs = np.arange(5, 40, 0.5)
cycles_n = 7   #np.round(np.linspace(5,10,len(Freqs)))
fwhm = mne.time_frequency.fwhm(Freqs, cycles_n)
AllChans = EpochAll_control[0].ch_names
ChanOI = ('FC3', 'FCz', 'FC4', 'C3', 'Cz', 'C4', 'CP3', 'CPz', 'CP4')
ChanoiIndx = [chidx for chidx, chanc in enumerate(AllChans) if chanc in ChanOI]

newtimes = EpochAllC_ext[0].times
lowlim_indx = np.argmin(np.abs(newtimes--0.15))
hilim_indx = np.argmin(np.abs(newtimes-1.0))
lims_indx = [lowlim_indx, hilim_indx]
Tin = newtimes[lims_indx[0]:lims_indx[1]]

# Call of function to carry out time-frequency analysis and baseline correction.
ERSPcontrol_filler, ERSPcontrol_verb, TFR_ob = doTFAnalysis(Freqs, cycles_n, ChanoiIndx, EpochControl_post, lims_indx)
ERSPtest_filler, ERSPtest_verb, TFRob = doTFAnalysis(Freqs, cycles_n, ChanoiIndx, EpochTest_post, lims_indx)
ERSPcontrol_filler_pre, ERSPcontrol_verb_pre, TFR_ob = doTFAnalysis(Freqs, cycles_n, ChanoiIndx, EpochControl_pre, lims_indx)
ERSPtest_filler_pre, ERSPtest_verb_pre, TFRob = doTFAnalysis(Freqs, cycles_n, ChanoiIndx, EpochTest_pre, lims_indx)

# Need to create an array for each condition and verb-type (verb, filler) and channel
ConFiller_arraypre = np.asarray(ERSPcontrol_filler_pre)
ConVerb_arraypre = np.asarray(ERSPcontrol_verb_pre)
TestFiller_arraypre = np.asarray(ERSPtest_filler_pre)
TestVerb_arraypre = np.asarray(ERSPtest_verb_pre)

ConFiller_array = np.asarray(ERSPcontrol_filler)
ConVerb_array = np.asarray(ERSPcontrol_verb)
TestFiller_array = np.asarray(ERSPtest_filler)
TestVerb_array = np.asarray(ERSPtest_verb)

Time = EpochAll_control[0].times
shape_control = np.shape(ConFiller_array)
shape_test = np.shape(TestVerb_array)

alpha = 0.025
vmin, vmax = -1.5, 1.5  # set min and max ERDS values in plot
permN = 2000
Tail = 0
levels = np.linspace(vmin,vmax, 50)
cnorm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)  # min, center & max ERDS
dof_ctrl = []
dof_test = []
dof_ctrl = shape_control[0] - 1
dof_test = shape_test[0]-1

t_thresh_ctrl = sp.stats.t.ppf(1 - alpha / 2, df=dof_ctrl)
t_thresh_test = sp.stats.t.ppf(1-alpha/2, df=dof_test)

Ctrl_comptype = [ConVerb_arraypre, ConFiller_arraypre]  # Pretest
Test_comptype = [TestVerb_arraypre, TestFiller_arraypre] # Pretest
PCtrl_comptype = [ConVerb_array, ConFiller_array]  # Post-test
PTest_comptype = [TestVerb_array, TestFiller_array]  # Post-test

Ctrlfiller_prepost = [ConFiller_arraypre, ConFiller_array]   # For post-test - pre-test difference
Testfiller_prepost = [TestFiller_arraypre, TestFiller_array]
Ctrlverb_prepost = [ConVerb_arraypre, ConVerb_array]
Testverb_prepost = [TestVerb_arraypre, TestVerb_array]

pvals_clusterPreC, Tmap, ttestvals1 = findSigERSP(ChanOI, Freqs, Tin, Ctrl_comptype, alpha, vmin, vmax, levels, permN, Tail, t_thresh_ctrl) # Call of function to carry out significance testing on ERSP for individual conditions.
pvals_clusterPreT, Tmap, ttestvals2 = findSigERSP(ChanOI, Freqs, Tin, Test_comptype, alpha, vmin, vmax, levels, permN, Tail, t_thresh_test)
pvals_clusterPC, Tmap, ttestvals3 = findSigERSP(ChanOI, Freqs, Time, PCtrl_comptype, alpha, vmin, vmax, levels, permN, Tail, t_thresh_ctrl)
pvals_clusterPT, Tmap, ttestvals4 = findSigERSP(ChanOI, Freqs, Time, PTest_comptype, alpha, vmin, vmax, levels, permN, Tail, t_thresh_test)

pval_testfill, Tmap, ttestvals5 = findSigERSP(ChanOI, Freqs, Tin, Testfiller_prepost, alpha, vmin, vmax, levels, permN, Tail, t_thresh_test)
pval_ctrlverb, Tmap, ttestvals6 = findSigERSP(ChanOI, Freqs, Tin, Ctrlverb_prepost, alpha, vmin, vmax, levels, permN, Tail, t_thresh_ctrl)
pval_testverb, Tmap, ttestvals7 = findSigERSP(ChanOI, Freqs, Tin, Testverb_prepost, alpha, vmin, vmax, levels, permN, Tail, t_thresh_test)
pval_ctrlfiller, Tmap, ttestvals8 = findSigERSP(ChanOI, Freqs, Tin, Ctrlfiller_prepost, alpha, vmin, vmax, levels, permN, Tail, t_thresh_ctrl)


# Save the array with t-values to a text file
fpath = path_test
filenom = 'HP_testverb_Tmap_Cz.txt'
filenom_ttest = 'HP_testverb_Ttest_Cz_vals.txt'
Tvals_save(fpath, filenom, Tmap[4], Tin, Freqs, ttestvals7[4], filenom_ttest)  # Call of function to save the array of tvalues to file.

T = ttestvals7[4]
kde = sp.stats.gaussian_kde(ttestvals7[4].T)
#visualize KDE
fig = plt.figure()
ax = fig.add_subplot(111)
x_eval = np.linspace(-2, 2, num=200)
ax.plot(x_eval, kde(x_eval), 'k-')

isnorm = sp.stats.shapiro(ttestvals7[5].T)
plt.hist(ttestvals7[5].T, edgecolor='black',density=True, bins=50)
plt.title('Histogram of t-test values (paired t-test):\n test statistic for normality (Shapiro-Wilkes) = '+ str(isnorm.statistic))
plt.show()

plt.hist(Tmap[4], edgecolor='black', bins=20)
plt.title('Histogram of test-test values (paired t-test from non-parametric cluster test)')
plt.show()