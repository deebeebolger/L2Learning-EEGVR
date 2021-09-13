function EEG = EEGVR_L2_addconds(EEG, mys)

sessiontypes = {'Match_Mismatch','Passive_Listening','Match_Mismatch','Passive_Listening'};
sujtitle = EEG.setname;
titletypes = mys.participant_data.prefixes;
titletypes = cellfun(@(t) t(1:end-1),titletypes, 'UniformOutput',false);
ispresent = cellfun(@(s) contains(sujtitle, s),titletypes);
Stypes = sessiontypes(ispresent);
curr_struct = mys.experiment.(matlab.lang.makeValidName(Stypes{1,1}));
allevents = {EEG.event.type};

 %% Erase the word "condition" from the event trigger labels. Leave only the number following the word "condition".
        
fx = find(contains({EEG.event.type},'condition'));
isindx = cellfun(@(E) isspace(E),{EEG.event(fx).type},'UniformOutput',false);

for ecnt = 1:length(fx)

    spacem = find(isindx{1,ecnt});
    EEG.event(fx(ecnt)).type = EEG.event(fx(ecnt)).type(spacem+1:end);

end

%% Now add "tstart" (start of trial) and "tend" (end of trial) labels to the "eventlabels column.
        
tstart_trig = curr_struct.pre_verb;
tend_trig = curr_struct.end_trial;

tstart_indx = find(contains({EEG.event.type},string(tstart_trig)));
tend_indx = find(contains({EEG.event.type},string(tend_trig)));

for tscnt = 1:length(tstart_indx)
    EEG.event(tstart_indx(tscnt)).eventlabels = 'trial_start';
end

for tecnt = 1:length(tend_indx)
    EEG.event(tend_indx(tecnt)).eventlabels = 'trial_end';
end

%% Different actions to be taken depending on whether the data is Passive_Listening or Match_Mismatch

switch Stypes{1,1}
    case 'Passive_Listening'
        %% Now add the labels "verb" and "filler" to the events structure.
        % Create a new field called "eventlabels" with the labels corresponding
        % to the trigger codes in the "type" column.
        
        for testcnt = 1:length(curr_struct.test_verbs)
            
            x = cellstr(string(curr_struct.test_verbs(testcnt)));
            testindx = cell2mat(cellfun(@(teststr) contains(allevents,teststr), x, 'UniformOutput',false));
            testIndx = find(testindx);
            for icnt = 1:length(testIndx)
                EEG.event(testIndx(icnt)).eventlabels = 'test_verb';
            end
        end
        
        for fillercnt = 1:length(curr_struct.filler_verbs)
            fx = cellstr(string(curr_struct.filler_verbs(fillercnt)));
            fillindx = cell2mat(cellfun(@(fillstr) contains(allevents,fillstr), fx, 'UniformOutput',false));
            fillIndx = find(fillindx);
            for icnt2 = 1:length(fillIndx)
                EEG.event(fillIndx(icnt2)).eventlabels = 'filler_verb';
            end
        end
        
        %% Now tidy up the triggers:
        % - take out the second test_verb or "filler_verb" after the
        % start-trial trigger.
        % - take out those triggers that correspond to trial-start but have
        % neither a "test_verb" nor a "filler_verb" following.
        % - assign the trials to block 1, 2 or 3
        
        
        for ievent = 2:length({EEG.event.type})
            
            if strcmp(EEG.event(ievent).eventlabels, 'test_verb') && strcmp(EEG.event(ievent-1).eventlabels, 'trial_start')
                fprintf('%s \n','Test_verb correctly positioned.');
                
            elseif strcmp(EEG.event(ievent).eventlabels, 'test_verb') && ~strcmp(EEG.event(ievent-1).eventlabels, 'trial_start')
                fSpec = 'Test_verb is incorrectly positioned after %s rather than  or trial_start\n';
                fprintf(fSpec,EEG.event(ievent-1).eventlabels);
                
                EEG.event(ievent).eventlabels = '';
                
            elseif strcmp(EEG.event(ievent).eventlabels, 'filler_verb') && strcmp(EEG.event(ievent-1).eventlabels, 'trial_start')
                fprintf('%s \n','Filler_verb correctly positioned.');
                
            elseif strcmp(EEG.event(ievent).eventlabels, 'filler_verb') && ~strcmp(EEG.event(ievent-1).eventlabels, 'trial_start')
                fSpec = 'filler_verb is incorrectly positioned after %s rather than  or trial_start\n';
                fprintf(fSpec,EEG.event(ievent-1).eventlabels);
                
                EEG.event(ievent).eventlabels = '';
                
            elseif strcmp(EEG.event(ievent).eventlabels, 'trial_start') && ievent == length({EEG.event.type})
                
                EEG.event(ievent).eventlabels = '';
                
            elseif strcmp(EEG.event(ievent).eventlabels, 'trial_start') && ~ischar(EEG.event(ievent+1).eventlabels)
                
                EEG.event(ievent).eventlabels = '';
            end
        end
        
        %% Assign the test and filler verbs to blocks 1, 2 or 3
        
        % Start with the test verbs
        alltv_indx = find(strcmp({EEG.event.eventlabels},'test_verb'));    % Find the index of the test_verbs
        
        alltv = unique({EEG.event(alltv_indx).type});
        
        for tvcnt = 1:length(alltv)
            
            currindx_tv = find(strcmp({EEG.event.type},alltv{1,tvcnt}));
            
            if length(currindx_tv)>3
                
                itv =  arrayfun(@(x) ismember(alltv_indx,x), currindx_tv,'UniformOutput',false);
                corrtv_indx = find(cell2mat(cellfun(@sum, itv,'UniformOutput',false)));
                Ctv = currindx_tv(corrtv_indx);
                
                EEG.event(Ctv(1)).blocknum = 'block1';
                EEG.event(Ctv(2)).blocknum = 'block2';
                EEG.event(Ctv(3)).blocknum = 'block3';
                
            else
                
                EEG.event(currindx_tv(1)).blocknum = 'block1';
                EEG.event(currindx_tv(2)).blocknum = 'block2';
                EEG.event(currindx_tv(3)).blocknum = 'block3';
            end
            
        end
        
        %Now for the filler verbs
        allfv_indx = find(strcmp({EEG.event.eventlabels},'filler_verb'));   % Find the index of the filler_verbs
        allfv = unique({EEG.event(allfv_indx).type});
        
        for fvcnt = 1:length(allfv)
            
            currindx = find(strcmp({EEG.event.type},allfv{1,fvcnt}));
            
            if length(currindx)>3
                
                i =  arrayfun(@(x) ismember(allfv_indx,x), currindx,'UniformOutput',false);
                corr_indx = find(cell2mat(cellfun(@sum, i,'UniformOutput',false)))
                C = currindx(corr_indx);
                
                EEG.event(C(1)).blocknum = 'block1';
                EEG.event(C(2)).blocknum = 'block2';
                EEG.event(C(3)).blocknum = 'block3';
                
            else
                
                EEG.event(currindx(1)).blocknum = 'block1';
                EEG.event(currindx(2)).blocknum = 'block2';
                EEG.event(currindx(3)).blocknum = 'block3';
            end
            
        end
        
    case 'Match_Mismatch'
        
          %% Now add the labels "Match" and "Mismatch" to the events structure.
        % Create a new field called "eventlabels" with the labels corresponding
        % to the trigger codes in the "type" column.
        
        for matchcnt = 1:length(curr_struct.match_verbs)
            x = cellstr(string(curr_struct.match_verbs(matchcnt)));
            matchindx = cell2mat(cellfun(@(matchstr) contains(allevents,matchstr), x, 'UniformOutput',false));
            matchIndx = find(matchindx);
            for icnt = 1:length(matchIndx)
                EEG.event(matchIndx(icnt)).eventlabels = 'match_verb';
            end
        end
        
        for mismcnt = 1:length(curr_struct.mismatch_verbs)
            fx = cellstr(string(curr_struct.mismatch_verbs(mismcnt)));
            mismindx = cell2mat(cellfun(@(mismstr) contains(allevents,mismstr), fx, 'UniformOutput',false));
            mismIndx = find(mismindx);
            for icnt2 = 1:length(mismIndx)
                EEG.event(mismIndx(icnt2)).eventlabels = 'mismatch_verb';
            end
        end
        %% Now tidy up the triggers:
        % - take out the second match_verb or "mismatch_verb" after the
        % start-trial trigger.
        % - take out those triggers that correspond to trial-start but have
        % neither a "match_verb" nor a "mismatch_verb" following.
        % - assign the trials to block 1, 2 or 3
        
        
        for ievent = 2:length({EEG.event.type})
            
            if strcmp(EEG.event(ievent).eventlabels, 'match_verb') && strcmp(EEG.event(ievent-1).eventlabels, 'trial_start')
                fprintf('%s \n','match_verb correctly positioned.');
                
            elseif strcmp(EEG.event(ievent).eventlabels, 'match_verb') && ~strcmp(EEG.event(ievent-1).eventlabels, 'trial_start')
                fSpec = 'Test_verb is incorrectly positioned after %s rather than  or trial_start\n';
                fprintf(fSpec,EEG.event(ievent-1).eventlabels);
                
                EEG.event(ievent).eventlabels = '';
                
            elseif strcmp(EEG.event(ievent).eventlabels, 'mismatch_verb') && strcmp(EEG.event(ievent-1).eventlabels, 'trial_start')
                fprintf('%s \n','Mismatch_verb correctly positioned.');
                
            elseif strcmp(EEG.event(ievent).eventlabels, 'mismatch_verb') && ~strcmp(EEG.event(ievent-1).eventlabels, 'trial_start')
                fSpec = 'mismatch_verb is incorrectly positioned after %s rather than  or trial_start\n';
                fprintf(fSpec,EEG.event(ievent-1).eventlabels);
                
                EEG.event(ievent).eventlabels = '';
                
            elseif strcmp(EEG.event(ievent).eventlabels, 'trial_start') && ievent == length({EEG.event.type})
                
                EEG.event(ievent).eventlabels = '';
                
            elseif strcmp(EEG.event(ievent).eventlabels, 'trial_start') && ~ischar(EEG.event(ievent+1).eventlabels)
                
                EEG.event(ievent).eventlabels = '';
            end
        end      
        
         %% Assign the match and mismatch verbs to blocks 1, 2 or 3
        
        % Start with the test verbs
        alltv_indx = find(strcmp({EEG.event.eventlabels},'match_verb'));    % Find the index of the test_verbs
        
        alltv = unique({EEG.event(alltv_indx).type});
        
        for tvcnt = 1:length(alltv)
            
            currindx_tv = find(strcmp({EEG.event.type},alltv{1,tvcnt}));
            
            if length(currindx_tv)>3
                
                itv =  arrayfun(@(x) ismember(alltv_indx,x), currindx_tv,'UniformOutput',false);
                corrtv_indx = find(cell2mat(cellfun(@sum, itv,'UniformOutput',false)));
                Ctv = currindx_tv(corrtv_indx);
                
                EEG.event(Ctv(1)).blocknum = 'block1';
                EEG.event(Ctv(2)).blocknum = 'block2';
                EEG.event(Ctv(3)).blocknum = 'block3';
                
            else
                
                EEG.event(currindx_tv(1)).blocknum = 'block1';
                EEG.event(currindx_tv(2)).blocknum = 'block2';
                EEG.event(currindx_tv(3)).blocknum = 'block3';
            end
            
        end
        
        %Now for the filler verbs
        allfv_indx = find(strcmp({EEG.event.eventlabels},'mismatch_verb'));   % Find the index of the filler_verbs
        allfv = unique({EEG.event(allfv_indx).type});
        
        for fvcnt = 1:length(allfv)
            
            currindx = find(strcmp({EEG.event.type},allfv{1,fvcnt}));
            
            if length(currindx)>3
                
                i =  arrayfun(@(x) ismember(allfv_indx,x), currindx,'UniformOutput',false);
                corr_indx = find(cell2mat(cellfun(@sum, i,'UniformOutput',false)))
                C = currindx(corr_indx);
                
                EEG.event(C(1)).blocknum = 'block1';
                EEG.event(C(2)).blocknum = 'block2';
                EEG.event(C(3)).blocknum = 'block3';
                
            else
                
                EEG.event(currindx(1)).blocknum = 'block1';
                EEG.event(currindx(2)).blocknum = 'block2';
                EEG.event(currindx(3)).blocknum = 'block3';
            end
            
        end
        
        %% Assign trials as either correct or incorrect for both match and mismatch.
        % Match trial correct ==> 200
        % Mismatch trial correct ==> 100
        % Match trial incorrect ==> 100
        % Mismatch trial incorrect ==> 200
        
        match_trials = find(strcmp({EEG.event.eventlabels}, 'match_verb'));
        
        for tcounter = 1:length(match_trials)
            
            if strcmp(EEG.event(match_trials(tcounter)+1).eventlabels, 'trial_end')
                corrtrig = str2double(EEG.event(match_trials(tcounter)+2).type);
                if corrtrig == 200
                    EEG.event(match_trials(tcounter)).response = 'correct';
                elseif corrtrig == 100
                    EEG.event(match_trials(tcounter)).response = 'incorrect';
                end
                
            elseif isempty(EEG.event(match_trials(tcounter)+1).eventlabels)
                corrtrig = str2double(EEG.event(match_trials(tcounter)+3).type);
                
                if corrtrig == 200
                    EEG.event(match_trials(tcounter)).response = 'correct';
                elseif corrtrig == 100
                    EEG.event(match_trials(tcounter)).response = 'incorrect';
                end
                
            end
        end
         
        mismatch_trials = find(strcmp({EEG.event.eventlabels}, 'mismatch_verb'));
        
        for tcounter = 1:length(mismatch_trials)
            
            if strcmp(EEG.event(mismatch_trials(tcounter)+1).eventlabels, 'trial_end')
                corrtrig = str2double(EEG.event(mismatch_trials(tcounter)+2).type);
                if corrtrig == 200
                    EEG.event(mismatch_trials(tcounter)).response = 'incorrect';
                elseif corrtrig == 100
                    EEG.event(mismatch_trials(tcounter)).response = 'correct';
                end
            elseif isempty(EEG.event(mismatch_trials(tcounter)+1).eventlabels)
                corrtrig = str2double(EEG.event(mismatch_trials(tcounter)+3).type);
                if corrtrig == 200
                    EEG.event(mismatch_trials(tcounter)).response = 'incorrect';
                elseif corrtrig == 100
                    EEG.event(mismatch_trials(tcounter)).response = 'correct';
                end
            end
        end
        
      
        
end

end