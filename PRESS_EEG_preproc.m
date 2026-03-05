% --------------------------------------------------------------------- %
%   EEG PREPROCESSING
% --------------------------------------------------------------------- %
%
% Description:  Matlab-script to pre-process Muse EEG data
%             - Study name: PRESS (dossiernr.: 25U-0078)
%             - Measure/instrument: Muse S Headband recordings
%             - Data type: EEG time series with trigger events (time stamps)
%             - Design: within-subjects longitudinal, 3 timepoints (T0 ("M1:") 27th 03/'25, T1 ("M2"): 10th and 11th 04/'25, T1 ("M4"): 24th 04/'25)
%
% Required toolbox:
% EEGLAB (used: version 2024.2)
% Required EEGLAB plugins:
% - clean_rawdata (used: version 2.11)
%
% Notes: 
% 1. Cleaning based on paper: Delorme, A., & Martin, J. A. (2021, December). Automated data cleaning for the Muse EEG. In 2021 IEEE International Conference on Bioinformatics and Biomedicine (BIBM) (pp. 1-5). IEEE. https://doi.org/10.1109/BIBM52615.2021.9669415
% 2. No re-referencing possible due to low no. channels; (or yes?)
% 3. No downsampling: not needed due to low sampling rate;
% 4. No ICA possible due to low no. channels
%
% Date:             July 2025, update March 2026
% Matlab version:   2024b
%
% --------------------------------------------------------------------- %

%% Clear workspace

clear
close all


%% Initialize

% initialize eeglab in nogui mode (no GUI opens, all initializes synchronously)
%cd('/Users/fsmits2/Downloads/eeglab2024.2')
addpath('/Users/fsmits2/Downloads/eeglab2024.2') %('/Applications/eeglab2024.0')
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab('nogui');
% NOTE: If EEGlab is started with GUI, do not continue before eeglab is fully started >> code will not run because EEGlab function start ups may still run asynchronously [in the background]
% % with GUI:
% [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;




%% Set paths and subject IDs

% change directory to where script is
cd('/Users/fsmits2/MATLAB/Projects/PRESSsandbox')

% find path name to research folder structure (RFS)
path2RFS = '/Users/fsmits2/Networkshares/Onderzoek/Groep Geuze/25U-0078_PRESS/'; % Enter your path to RFS. End with slash ('/' on Mac, '\' on Windows)

% set other paths 
path2data    = [path2RFS 'E_ResearchData/2_ResearchData/'];
path2EEGsets = [path2data '1. Verwerkte data/Muse/EEGLAB sets with markers/'];
path2save    = [path2data '1. Verwerkte data/Muse/'];
%%%path2scripts = [path2RFS 'F_DataAnalysis/1f_DataAnalysisScripts/';

% enter subject names (extract from key file where subject IDs are linked)
key_filename = [path2data 'SubjectID_koppelbestand_Castor-PVT-Garmin.csv'];
subj_tab     = readtable( key_filename, 'ReadVariableNames', 1);
subj_list    = table2array( subj_tab(:,1) );

% enter session and experimental task names
sessions = {'M1', 'M2', 'M4'};
tasks    = {'rust', 'startle'};


%% Load data, re-reference, filter

%  %% * Example filename:
% filename = 'sub-110027_ses-M1-rust_task-Default_run-001_eeg.xdf.set'; 

% Loop over subjects, sessions and tasks (rest and startle)
for subj_i = 1:length(subj_list)
    for sess_i = 1:length(sessions)
        for task_i = 1:length(tasks)

            fprintf('\n****\nStart processing subject %i session %s %s\n****\n\n', subj_list(subj_i), sessions{sess_i}, tasks{task_i});
            filename = sprintf('sub-%i_ses-%s-%s_task-Default_run-001_eeg.xdf.set', subj_list(subj_i), sessions{sess_i}, tasks{task_i});

            % -- Check if file exists --
            fullPath = fullfile(path2EEGsets, filename);
            if ~exist(fullPath, 'file')
                fprintf('File %s not found. Skipping.\n', fullPath);
                continue;  % skip to next iteration of your subject/session/task loop
            end

            % -- Load raw EEG set --
            EEG      = pop_loadset('filename',filename,'filepath',path2EEGsets);

            %  %% * if needed: View set characteristics or plot for visual inspection
            % disp(EEG);
            % pop_eegplot( EEG, 1, 1, 1); % Inspect data

            % -- Add info to EEG structure --
            EEG.filename  = filename;
            EEG.setname   = filename;
            EEG.subject   = subj_list(subj_i);
            EEG.session   = sess_i;
            EEG.condition = tasks{task_i};

            % -- Select only brain activity channels (AF7,AF8,TP9,TP10) (i.e., remove flat line right AUX channel)
            EEG = pop_select( EEG, 'channel',{'TP9','AF7','AF8','TP10'});

            % -- Rename trigger events --
            % 0) Rename 'habituation_probe' (startle task) to 'habituation_sound' to avoid matching 'probe'
            if strcmp(tasks{task_i}, 'startle') %this only needed for startle task
                hab_idx = find(contains({EEG.event.type}, 'habituation_probe'));
                [EEG.event(hab_idx).type] = deal('habituation_sound');
            end
            % 1) Remove number in event name that represented time indication
            conditions = ["open", "close", "end", "pre-pulse", "probe", "habituation"];
            for cond = conditions
                idx = find(contains({EEG.event.type}, cond));  % partial match now safe
                [EEG.event(idx).type] = deal(char(cond));
            end
            % 2) Adapt events of startle task
            if strcmp(tasks{task_i}, 'startle') %this only needed for startle task
                % a) Recode probe events to probes with or without preceding pre-pulse
                % logical index of probe and pre-pulse events
                probe_idx = strcmpi({EEG.event.type}, 'probe');
                pp_idx    = strcmpi({EEG.event.type}, 'pre-pulse');
                % get corresponding latencies
                probe_lats = [EEG.event(probe_idx).latency];
                pp_lats    = [EEG.event(pp_idx).latency];
                % check where probe and pre-pulse latencies match 
                for ev = 1:length(EEG.event)
                    if strcmpi(EEG.event(ev).type, 'probe')
                        probe_lat = EEG.event(ev).latency;
                        if any(pp_lats == probe_lat) %check if any pre-pulses have the same latency 
                            EEG.event(ev).type = 'probe_pp';
                        else
                            EEG.event(ev).type = 'probe_solo';
                        end
                    end
                end
                % b) Subtract 100 ms of the pre-pulse latencies (latency of pre-pulse was recorded at the same time of probe, while presented 100 ms earlier)
                EEG.event(pp_idx) = arrayfun(@(e) setfield(e, 'latency', e.latency - 0.1*EEG.srate), EEG.event(pp_idx));
                % Clean up event structure
                EEG = eeg_checkset(EEG, 'eventconsistency');
            end

            % % -- Re-reference --
            % % Re-reference to avg mastoids when interested in frontal activity - advise by: https://doi.org/10.1101/2021.11.02.466989 Cannard, C., Wahbeh, H., & Delorme, A. (2021, December). Validating the wearable MUSE headset for EEG spectral analysis and Frontal Alpha Asymmetry. IEEE.
            % mastoid1 = find(strcmpi( {EEG.chanlocs.labels}, 'TP9' ));
            % mastoid2 = find(strcmpi( {EEG.chanlocs.labels}, 'TP10' ));
            % EEG      = pop_reref( EEG, [mastoid1 mastoid2]); %re-references to the average of 2 channels << check if this is right, or put the "real" formula here: Vref = ½(vLM + vRM) where Vref is the reference signal and vLM and vRM are the signals for the left and right mastoid channels respectively

            % -- Filter -- 
            % Filter settings: Delorme 2021 uses lowpass filter at 40 Hz. Here, we likewise apply no highpass filter and we lowpass filter at 34 Hz for this is the maximum that we can apply universally to all datasets because some datasets have a reduced sampling rate of 69.9, and filter requires at least double the sampling rate.
            % Output:          'EEG' contains the filtered EEGLAB structure, 'com' contains the history string (the matlab command), 'b' contains the filter coefficients (plot to see the filter)
            [EEG, com, b] = pop_eegfiltnew(EEG,'hicutoff',34); % with highpass filer at 1 hz: EEG = pop_eegfiltnew(EEG, 'locutoff',1,'hicutoff',34); 

            % -- Save processed EEG set --
            fprintf('\n****\nSaving processed subject %i session %s %s\n****\n\n', subj_list(subj_i), sessions{sess_i}, tasks{task_i});
            SaveName = sprintf( '%i-%s-%s_filtered.set', subj_list(subj_i), sessions{sess_i}, tasks{task_i} );
            EEG      = pop_saveset( EEG,'filename',SaveName,'filepath', path2save );

             % -- Clear this EEG set from workspace before iterating to the next --
            clear EEG
            ALLEEG(1:end) = [];
        end
    end
end



%% Epoch the data

% Loop over subjects, sessions and tasks (rest and startle)
for subj_i = 1:length(subj_list)
    for sess_i = 1:length(sessions)
        for task_i = 1:length(tasks)

            % * if needed: Create epoch markers automatically
            if strcmp(tasks{task_i}, 'rust') %here: only needed for rest-EEG
                % To divide long-duration recordings (e.g., resting-state EEG) into epochs, add trigger events after each (f.e.) 1-second period
                conds = {'open','close'};
                for ci = 1:numel(conds)
                    cond = conds{ci};
                    starts = find(strcmpi({EEG.event.type}, cond));
                    for s = 1:length(starts)
                        period     = EEG.event(starts(s)).latency + [0 60*EEG.srate];
                        event_lats = period(1):EEG.srate:(period(2)-EEG.srate);
                        % append events
                        for segi = 1:length(event_lats)
                            new_event = EEG.event(1); %copy template of event to have all the struct-fields EEGlab requires
                            new_event.type    = cond; %overwrite event type
                            new_event.latency = event_lats(segi);  %overwrite event latency, in samples
                            EEG.event(end+1)  = new_event; %append
                        end
                    end
                end
            end
            EEG = eeg_checkset(EEG, 'eventconsistency'); %make sure events are sorted after inserting

            % for startle task: you want epochs from -1800 to 150 ms around startle probes, in order to quantify and baseline-correct the maximum amplitude between 20 and 120 milliseconds after probe onset
            % Approach startle reflex as an eyeblink ("poor man's ocular-EMG") and try to quantify the amplitude of the blinks 20-120 ms after the probe as
            % detected by a blink algortihm for headband EEG such as Chang et al. 2014: https://doi.org/10.1016/j.cmpb.2015.10.011

            % Epoch the data
            EEG  = pop_epoch( EEG, {'open' 'close'},  [0  1], 'epochinfo', 'yes'); %epoch 0-1 seconds around event

            % Save
            fprintf('\n****\nSave epoched data subject %i session %s task %s\n****\n\n', subj_list(subj_i), sessions{sess_i}, tasks{task_i});
            SaveName = sprintf( '%i-%s-%s_epoched.set', subj_list(subj_i), sessions{sess_i}, tasks{task_i} );
            EEG      = pop_saveset( EEG, 'filename',SaveName,'filepath', path2save );

        end
    end
end

% % % % Add event after each 1-second into eyes-open or eyes-closed condition 
% % % % 1) find event latencies of start
% % % open_begin  = sort(find(strcmpi( {EEG.event.type}, 'open' ))); %Find the open/closed eyes onset triggers
% % % close_begin = sort(find(strcmpi( {EEG.event.type}, 'close' )));
% % % 
% % % % eyes open: find start latencies
% % % for i = 1:length(open_begin)
% % %     open_period(i,1) = EEG.event(open_begin(i)).latency;
% % %     open_period(i,2) = EEG.event(open_begin(i)).latency + 60*EEG.srate;
% % % end
% % % % eyes closed: find start latencies
% % % for i = 1:length(close_begin)
% % %     close_period(i,1) = EEG.event(close_begin(i)).latency;
% % %     close_period(i,2) = EEG.event(close_begin(i)).latency + 60*EEG.srate;
% % % end
% % % 
% % % % eyes open: insert events
% % % for i = 1:length(open_begin)
% % %     trggLats = open_period(i,1):EEG.srate:(open_period(i,2)-1*EEG.srate);
% % %     for segi = 1:length(trggLats)
% % %         EEG = pop_editeventvals(EEG,'insert',{ segi ,[],[],[]},...
% % %             'changefield',{ segi ,'type',    'open' },...
% % %             'changefield',{ segi ,'latency', trggLats(segi)/EEG.srate }); % latency of events (EEG.event.latency) is defined in data points, not in time. But when you want to change it, you need to define in seconds.
% % %     end
% % % end
% % % 
% % % % eyes closed: insert events
% % % for i = 1:length(close_begin)
% % %     trggLats = close_period(i,1):EEG.srate:(close_period(i,2)-1*EEG.srate);
% % %     for segi = 1:length(trggLats)
% % %         EEG = pop_editeventvals(EEG,'insert',{ segi ,[],[],[]},...
% % %             'changefield',{ segi ,'type',    'close' },...
% % %             'changefield',{ segi ,'latency', trggLats(segi)/EEG.srate }); % latency of events (EEG.event.latency) is defined in data points, not in time. But when you want to change it, you need to define in seconds.
% % %     end
% % % end
% % % 
% % % 
% % % % Epoch the data 
% % % EEG  = pop_epoch( EEG, {'open' 'close'},  [0  1], 'epochinfo', 'yes'); % epoch 0-1 seconds around event

%% Artifact rejection
% Evaluate and remove if necessary bad channels and epochs 
% Note: 
% No re-referencing possible due to low no. channels;
% No downsampling not needed due to low sampling rate);
% No ICA possible due to low no. channels

% 0) Import a channel location file 
% (needed? for Artifact Subspace Reconstruction, ASR)
EEG = pop_chanedit(EEG, {'lookup','/Users/fsmits2/Downloads/eeglab2024.2/plugins/dipfit/standard_BEM/elec/standard_1005.elc'});


% 1) initiate or load matrix to save noisy channels and rejected epochs
% Check if file exists:
fullPath_badch = fullfile([ path2save '/Overview_badchannels.txt']);
fullPath_rejep = fullfile([ path2save '/Overview_rejected_epochs.txt']);

if ~exist(fullPath_badch, 'file')
    % When not found, create a new matrix
    fprintf('Bad channels file not found: %s. \nCreating new file.\n', fullPath_badch);
    bdchns      = cell(length(subj_list),length(sessions)+1); 
    bdchns(:,1) = num2cell(subj_list');
else
    % When found, read existing matrix
    fprintf('Bad channels file found! %s. \nReading file.\n', fullPath_badch);
    bdchns = table2cell( readtable( fullPath_badch ) );
end

if ~exist(fullPath_rejep, 'file')
    % When not found, create a new matrix
    fprintf('Rejected epochs file not found: %s. \nCreating new file.\n', fullPath_rejep);
    rej_epocs = [subj_list  nan( length(subj_list),length(sessions)*length(tasks) )];
else
    % When found, read existing matrix
    fprintf('Rejected epochs file found! %s. \nReading file.\n', fullPath_rejep);
    rej_epocs = table2cell( readtable( fullPath_rejep ) );
end


% 2) Loop over files
for subj_i = 1:length(subj_list)
    for sess_i = 1:length(sessions)
        for task_i = 1:length(tasks)

            filename = sprintf( '%i-%s-%s_epoched.set', subj_list(subj_i), sessions{sess_i}, tasks{task_i} );

            % Load EEG set
            EEG = pop_loadset('filename', filename , 'filepath', path2save);

            % Detect noisy channels
            pop_eegplot( EEG, 1, 1, 1);
            for chan = 1:EEG.nbchan
                % - - Compute standard deviation - -                
                chansd = std(EEG.data(chan, :)');
                fprintf('%f = STD channel %s.\n', chansd, EEG.chanlocs(chan).labels);
                % - — Compute PSD - - 
                % Concatenate epochs
                signal = reshape(EEG.data(chan,:,:), [], 1); %squeeze over timepoints × trials
                % Compute mean bandpower per channel in PSD 5-34 Hz (34 Hz is based on filter freq, and 5 Hz is based on Delorme, A., & Martin, J. A. (2021, December). Automated data cleaning for the Muse EEG.)
                % Delorme et al. use mean bandpower >25 log10(µV2)/Hz as threshold.
                [pxx,f] = pwelch(signal, window, 0, [], EEG.srate);
                idx = f >= 5 & f <= 34;
                band_power = mean(log10(pxx(idx)));
                fprintf('%f = mean PSD channel %s.\n', band_power, EEG.chanlocs(chan).labels);
            end
            
            % Visually check noisy channels in plot and decide to keep or reject
            m3a = "no"; m3 = -1;
            while m3a ~= "yes"
                m3a = input('Ready to input channels to leave out? Type [yes] ','s');
            end       
            while m3 == -1
                m3 = str2double( input('How many channels to leave out? ','s') );
            end
            badchannels = {[]}; badchannrs = [];
            if m3 > 0
                for badchani = 1:m3
                    badchannels{badchani} = input(['Which channel to leave out? nr: ' num2str(badchani) ' ' ],'s') ;
                    badchannrs(badchani)  = find( strcmpi( badchannels{badchani}, {EEG.chanlocs.labels} ));
                end
                bdchns{subj_i,sess_i+1}  = string(badchannels);
            else
                bdchns{subj_i,sess_i+1}  = '0';
            end
            EEG.eventdescription = { {'Too much noise in channels: '} badchannels };

            % Leave noisy channels (previously detected and saved) out of consideration for epoch rejection
            badchannels = bdchns{subj_i, sess_i+1};
            chanarray   = 1:length(EEG.chanlocs);
            if sum( strcmpi( badchannels, '0') ) < 1
                badchans = regexp(badchannels, ',', 'split');
                badchannrs = [];
                for badchani = 1:length(badchans)
                    badchannrs(badchani) = find( strcmpi( badchans{badchani}, {EEG.chanlocs.labels} ));
                end
                chanarray(badchannrs) = []; % remove the noisy channels from chanarray
            end

            % - - Semi-automatic artifact rejection - -
            % Artifact Subspace Reconstruction (ASR)
            % Method (ASR) based on Delorme, A., & Martin, J. A. (2021, December). Automated data cleaning for the Muse EEG. In 2021 IEEE International Conference on Bioinformatics and Biomedicine (BIBM) (pp. 1-5). IEEE. https://doi.org/10.1109/BIBM52615.2021.9669415
            % Paramaters based on musemonitor code (https://github.com/sccn/eeglab_musemonitor_plugin/blob/master/pop_musemonitor.m)
            % BurstCriterion value: 11 ("The best method is the ASR method with a parameter of 11.")
            channel_crit = 0.8;
            burst_crit   = 11;

            % from musemonitor code:
            EEG = clean_rawdata(EEG, ...
                'FlatlineCriterion','off',...
                'ChannelCriterion',channel_crit,...
                'LineNoiseCriterion',5,...
                'Highpass',[0.25 0.75],...
                'BurstCriterion',burst_crit,...
                'WindowCriterion',0.25,...
                'BurstRejection','on',...
                'Distance','Euclidian',...
                'WindowCriterionTolerances',[-Inf 7] );

            % Delorme paper says: "We varied the WindowCriterionTolerances argument 
            % from 5 to 15 in increments of 1, set the WindowCriterion parameter 
            % to 0 (to automatically reject bad portions of data instead of trying 
            % to correct them), and disabled all other features including BurstRemoval. 
            % We used default values for filtering parameters."
            % Implementing this text in code below:
            EEG = pop_clean_rawdata(EEG, ...
                'FlatlineCriterion','off', ...
                'ChannelCriterion','off', ...
                'LineNoiseCriterion','off', ...
                'Highpass',[0.25 0.75], ...
                'BurstCriterion',burst_crit, ...
                'WindowCriterion',0, ...
                'BurstRejection','off', ...
                'Distance','Euclidian', ...
                'WindowCriterionTolerances',[-Inf 15] );

           % from code David ten Kate:
            EEG_clean = clean_artifacts(EEG, ...
                'FlatlineCriterion','off', ...
                'ChannelCriterion',25, ...
                'LineNoiseCriterion',5, ...
                'Highpass',[0.25 0.75], ...
                'BurstCriterion',11, ...
                'WindowCriterion',0.25, ...
                'BurstRejection','on', ...
                'Distance','Euclidian', ...
                'WindowCriterionTolerances',[-Inf 7]);


            % - - Artifact rejection with ASR - -
            % Create continuous signal from epoched data (required for ASR) - -
            EEG = eeg_epoch2continuous(EEG);
            % Artifact rejection with ASR in repair mode (no removal of data, just get suggestion of what should be removed):
            EEG_asr = clean_artifacts(EEG, ...
                'FlatlineCriterion','off', ...
                'ChannelCriterion','off', ...  % prevent automatic removal of channels
                'LineNoiseCriterion',5, ...
                'Highpass',[0.25 0.75], ...
                'BurstCriterion',11, ... 
                'BurstCriterionRefMaxBadChns','off', ...               
                'WindowCriterion',0.5, ...
                'BurstRejection','off', ... % prevent automatic removal of segments
                'Distance','Euclidian');

            % Visualize artifacts:
            vis_artifacts(EEG_asr, EEG, 'WindowLength', 50, 'DisplayMode', 'both');

            % find bad samples according to ASR:
            bad_samples = find(~EEG_asr.etc.clean_sample_mask);

              
            % View the bad segments in plot
            %   Scale value to 100 and 29 epochs per window. 
            %rej_epocs(subj_i,1+sess_i) = EEG.trials; % Save total number of epochs
            EEG = eeg_checkset( EEG );
            pop_eegplot( EEG, 1, 1, 0); % Plot data with marked epochs but do not immediately reject, only mark as noisy

% ----------------- 04/03/26
            m0 = -1;
            while m0 == -1
                m0 = input('Ready to reject samples? ','s');
                while isempty(m0)
                    m0 = input('Ready to reject samples? [yes]: ','s');
                end
            end

            noisyepocs = find(EEG.reject.rejmanual > 0) % see & save final series of marked epochs
            length(noisyepocs)

            m1 = [] ;
            while isempty(m1)
                m1 = input('How many epochs rejected? [enter number]: ','s');
            end

            rej_epocs(subj_i,3+sess_i) = str2double(m1);
            EEG.epochdescription       = [m1 '/' num2str(rej_epocs(subj_i,1+sess_i)) ' trials rejected'];
            EEG                        = pop_rejepoch( EEG, noisyepocs , 1);
            [ALLEEG EEG CURRENTSET]    = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off');
            EEG = eeg_checkset( EEG );
            [ALLEEG, EEG, CURRENTSET]  = eeg_store( ALLEEG, EEG );

           




            %     Gradient:  Specifies that the absolute difference between two adjacent sample points of data must not exceed a value (artifact of weird spikes). Starting values from Boost tutorial: Gradient: 75 μV
            %     Amplitude: Specifies that the voltage must not  exceed a certain value (artifacts like eye blinks). Starting values from Boost tutorial: Max-Min: 150 μV/200 ms
            %     Diff max-Min: Sets the threshold for the difference between the minimum and maximum voltages within the entire segment (voltage drifts). Starting values from Boost tutorial: Amplitude: -100 μV, +100 μV"
            winpnts = round(200/(1000/EEG.srate)); % points for window of 200 ms segments in the epoch
            winidx  = 1:winpnts:EEG.pnts;
            windiff = nan(1,length(winidx));
            EEG.reject.rejmanual = zeros(1, EEG.trials); % Initialize the array for marked trials
            EEG.reject.rejmanualE = zeros(length(EEG.chanlocs), EEG.trials);

            % Loop over channels and epochs
            for ichan = chanarray
                for itrial = 1:EEG.trials
                    gradient = max( abs( diff(EEG.data(ichan, :, itrial)) ) );
                    ampliMax = max(EEG.data(ichan, :, itrial));
                    ampliMin = min(EEG.data(ichan, :, itrial));
                    for iwin = 1:length(winidx)-1
                        [winmin, winmax] = bounds(EEG.data( ichan, winidx(iwin):winidx(iwin)+winpnts-1, itrial));
                        windiff(iwin) = diff([winmin, winmax]);
                    end
                    diffV = max(windiff);

                    if gradient > 50 || ampliMax > 100 || ampliMin < -100 || diffV > 150  % gradient > 50 || ampliMax > 75 || ampliMin < -75 || diffV > 150
                        EEG.reject.rejmanual(1,itrial) = 1; % Mark the trial when a criterium is met
                        EEG.reject.rejmanualE(ichan,itrial) = 1;
                    end
                end
            end

            % View the marked trials in plot
            %   Scale value to 100 and 29 epochs per window. Pay attention to VEOG.
            %   Reject the epoch around tACS artifact
            rej_epocs(subj_i,1+sess_i) = EEG.trials; % Save total number of epochs
            find(EEG.reject.rejmanual > 0) % See the marked epoch numbers
            EEG = eeg_checkset( EEG );
            pop_eegplot( EEG, 1, 1, 0); % Plot data with marked epochs but do not immediately reject, only mark as noisy

            m0 = -1;
            while m0 == -1
                m0 = input('Ready to reject epocs? ','s');
                while isempty(m0)
                    m0 = input('Ready to reject epocs? [yes]: ','s');
                end
            end

            noisyepocs = find(EEG.reject.rejmanual > 0) % see & save final series of marked epochs
            length(noisyepocs)

            m1 = [] ;
            while isempty(m1)
                m1 = input('How many epochs rejected? [enter number]: ','s');
            end

            rej_epocs(subj_i,3+sess_i) = str2double(m1);
            EEG.epochdescription       = [m1 '/' num2str(rej_epocs(subj_i,1+sess_i)) ' trials rejected'];
            EEG                        = pop_rejepoch( EEG, noisyepocs , 1);
            [ALLEEG EEG CURRENTSET]    = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off');
            EEG = eeg_checkset( EEG );
            [ALLEEG, EEG, CURRENTSET]  = eeg_store( ALLEEG, EEG );

            % Save
            fprintf('\n****\nSave clean data subject %i session %i\n****\n\n', subj_list(subj_i), sessions(sess_i));
            SaveName = [filename '_CleanEEG.set'];
            EEG      = pop_saveset( EEG, 'filename',SaveName,'filepath', path2EEGsets );
            cd(path2save);
            writecell(  bdchns,     [path2save '/Overview_badchannels_'     'PRESS_restEEG' char(datetime('today')) '.txt'], 'Delimiter',',');
            writematrix(rej_epocs,  [path2save '/Overview_rejected_epochs_' 'PRESS_restEEG' char(datetime('today')) '.txt'], 'Delimiter',',');

            m2 = 0;
            while m2 == 0
                m2 = input('Continue? [Y/N] ','s');
                if m2 == 'Y'
                    continue
                else
                    return
                end
            end

            clear EEG
            close all

            % open EEGlab
            cd('/Users/fsmits2/Downloads/eeglab2024.2')
            [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab('nogui');

        end
    end
end

%% Spectral analysis

% FFT over all channels
dataX   = fft(EEG.data,[],2);
dataPow = mean( abs(dataX) ,3);

% vector of frequencies
hz = linspace(0,EEG.srate,floor(EEG.pnts));

% frequency cutoffs in Hz and indices
frex = [ 1 30 ];
fidx = dsearchn(hz',frex');

addpath('/Volumes/home-3/Documents/Administratie/Congressen, cursussen en praatjes/Cursussen/Analyzing Neural Time Series data - alleen online/NTSA_spectral')

figure(2), clf
subplot(121)
plot(hz(fidx(1):fidx(2)),squeeze(dataPow(1,fidx(1):fidx(2))),'b'), hold on
plot(hz(fidx(1):fidx(2)),squeeze(dataPow(2,fidx(1):fidx(2))),'r'), hold on
plot(hz(fidx(1):fidx(2)),squeeze(dataPow(3,fidx(1):fidx(2))),'g'), hold on
plot(hz(fidx(1):fidx(2)),squeeze(dataPow(4,fidx(1):fidx(2))),'k')

subplot(122)
topoplotIndie(dataPow,EEG.chanlocs,'numcontour',0,'electrodes','off','shading','interp');
title('Spectral power 1-30 Hz' )
