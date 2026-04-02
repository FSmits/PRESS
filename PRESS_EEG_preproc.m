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
% 2. No re-referencing due to low no. channels and noisy re-referencing electrodes that also contain brain data, even though re-referencing to TP9 and TP10 is advised when interested in frontal activity by: https://doi.org/10.1101/2021.11.02.466989 Cannard, C., Wahbeh, H., & Delorme, A. (2021, December). Validating the wearable MUSE headset for EEG spectral analysis and Frontal Alpha Asymmetry. IEEE.
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

% Define epoch lengths (in seconds)
epoch_length = 1; %in seconds


%% Load data, trim, filter

%  %% * Example filename:
% filename = 'sub-110027_ses-M1-rust_task-Default_run-001_eeg.xdf.set'; 

% Initiate cell variable to save which datasets are skipped for pre-processing
skipped_sets = cell(2,1); 
skip_i       = 1;

% Loop over subjects, sessions and tasks (rest and startle)
for subj_i = 1:length(subj_list)
    for sess_i = 1:length(sessions)
        for task_i = 1:length(tasks)

            fprintf('\n****\nStart processing subject %i session %s %s\n****\n\n', subj_list(subj_i), sessions{sess_i}, tasks{task_i});
            filename = sprintf('sub-%i_ses-%s-%s_task-Default_run-001_eeg.xdf.set', subj_list(subj_i), sessions{sess_i}, tasks{task_i});

            fullPath = fullfile(path2EEGsets, filename);

            % -- Check if file exists --
            if ~exist(fullPath, 'file')
                fprintf('Not found: File %s.\n', fullPath);
                % Note: some resting-state EEG files are saved as "M1-rest" instead of "M1-rust"
                if task_i == 1 && ~exist(filename, 'file')
                    % redefine full filename
                    filename = sprintf('sub-%i_ses-%s-%s_task-Default_run-001_eeg.xdf.set', subj_list(subj_i), sessions{sess_i}, 'rest');
                    fullPath = fullfile(path2EEGsets, filename);
                end
            end
            % if also not found with name 'rest', then skip
            if ~exist(fullPath, 'file')
                fprintf('Skipping.\n', fullPath);
                skipped_sets{skip_i} = sprintf('Not found: sub-%i_ses-%s-%s_task-Default_run-001_eeg.xdf.set', subj_list(subj_i), sessions{sess_i}, tasks{task_i}); 
                skip_i = skip_i + 1;
                continue;  % skip to next iteration of your subject/session/task loop
            end

            % -- Load raw EEG set --
            EEG      = pop_loadset('filename',filename,'filepath',path2EEGsets);

            % -- Check if there is at least 20s of recorded data --
            if EEG.pnts/EEG.srate < 20
                fprintf('Less than 20s recorded data in %s.\n Skipping.\n', fullPath);
                skipped_sets{skip_i} = sprintf('<20s of data: sub-%i_ses-%s-%s_task-Default_run-001_eeg.xdf.set', subj_list(subj_i), sessions{sess_i}, tasks{task_i}); 
                skip_i = skip_i + 1;
                continue;  % skip to next iteration of your subject/session/task loop
            end

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
                if any(pp_idx > 0 ) % only if pre-pulse events are found in the data
                    EEG.event(pp_idx) = arrayfun(@(e) setfield(e, 'latency', e.latency - 0.1*EEG.srate), EEG.event(pp_idx));
                end
                % Clean up event structure
                EEG = eeg_checkset(EEG, 'eventconsistency');
            end

            % -- Trim the data to time of experimental task --
            % find first and last event latency
            EEG = eeg_checkset(EEG, 'eventconsistency');
            clear start_event; clear first_event; clear last_event;
            if task_i == 1
                start_event = find(strcmp({EEG.event.type}, 'open'),1); %first eyes open condition in resting-state EEG
            else
                start_event = find(strcmp({EEG.event.type}, 'habituation'),1); %first habituation probe in startle paradigm
            end
            first_event = EEG.event(start_event).latency;
            start_time  = first_event / EEG.srate; % convert to seconds
            last_event = EEG.event(end).latency;
            end_time   = last_event / EEG.srate;
            % Define how much time of data to keep after last experimental event (absence of 'end' event)
            if task_i == 1
                extra_time = 60; % resting-state EEG: maintain 60 seconds of data after last event
                if EEG.event(end).type == "end"
                    end_time   = EEG.event(end).latency / EEG.srate;
                    extra_time = 0;
                end
            else
                extra_time = 5; % startle paradigm: maintain 5 seconds of data after last event
            end
            % trim dataset
            EEG = pop_select(EEG, 'time', [start_time-2 end_time+extra_time]); %trim from 2 seconds before first event to 'extra time' seconds after last event

            % -- Filter -- 
            % Filter settings: Delorme 2021 uses lowpass filter at 40 Hz. Here, we likewise apply no highpass filter and we lowpass filter at 34 Hz for this is the maximum that we can apply universally to all datasets because some datasets have a reduced sampling rate of 69.9, and filter requires at least double the sampling rate.
            % Output:          'EEG' contains the filtered EEGLAB structure, 'com' contains the history string (the matlab command), 'b' contains the filter coefficients (plot to see the filter)
            [EEG, com, b] = pop_eegfiltnew(EEG,'hicutoff',34); % with highpass filer at 1 hz: EEG = pop_eegfiltnew(EEG, 'locutoff',1,'hicutoff',34); 

            % -- Save processed EEG set --
            fprintf('\n****\nSaving processed subject %i session %s %s\n****\n\n', subj_list(subj_i), sessions{sess_i}, tasks{task_i});
            SaveName = sprintf( '%i-%s-%s_trimmed_filtered.set', subj_list(subj_i), sessions{sess_i}, tasks{task_i} );
            EEG      = pop_saveset( EEG,'filename',SaveName,'filepath', path2save );

             % -- Clear this EEG set from workspace before iterating to the next --
            clear EEG
            ALLEEG(1:end) = [];
        end
    end
end

% Save cell array with information about skipped datasets
cd(path2save);
writecell(skipped_sets, [path2save '/Overview_missing_EEGdata_PRESS_' char(datetime('today')) '.txt'], 'Delimiter',',');

% change directory to where script is
cd('/Users/fsmits2/MATLAB/Projects/PRESSsandbox')
          

%% Insert extra 'epoch' markers before cleaning

% Loop over subjects, sessions and tasks (rest and startle)
for subj_i = 1:length(subj_list)
    for sess_i = 1:length(sessions)
        for task_i = 1:length(tasks)

            filename = sprintf( '%i-%s-%s_trimmed_filtered.set', subj_list(subj_i), sessions{sess_i}, tasks{task_i} );

             % -- Check if file exists --
            fullPath = fullfile(path2save, filename);
            if ~exist(fullPath, 'file')
                fprintf('Not found: File %s.\n Skipping.\n', fullPath);
                continue;  % skip to next iteration of your subject/session/task loop
            end

            % Load EEG set
            EEG = pop_loadset('filename', filename , 'filepath', path2save);

            % Define condition event names per task
            if strcmp(tasks{task_i}, 'rust')                 
                conds  = {'open','close'};
            elseif strcmp(tasks{task_i}, 'startle') 
                conds = {'habituation','probe_pp','probe_solo'};
            end

            % Make sure events are sorted before inserting new events
            EEG = eeg_checkset(EEG, 'eventconsistency');

            % -- Append events --
            for ci = 1:numel(conds)
                % find events belonging to conditions of interest
                cond = conds{ci};
                starts = find(strcmpi({EEG.event.type}, cond));
                for s = 1:length(starts)
                    current_event = starts(s);
                    % skip if it is the last event in dataset
                    if current_event == length(EEG.event)
                        continue
                    end
                    % find latencies of this event and next event
                    start_lat  = EEG.event(current_event).latency;
                    next_lat   = EEG.event(current_event + 1).latency;
                    % when this condition event is the last event, i.e., no 'end event' was present in original data, compute latency of the end of condition
                    if start_lat == max([EEG.event.latency])
                        end_cond_lat = start_lat+60*EEG.srate-EEG.srate;
                        if end_cond_lat < EEG.pnts
                            event_lats = start_lat:(EEG.srate*epoch_length):end_cond_lat;
                        else
                            % when end of condition latency is beyond latency of last data sample, only add markers until end of data
                            event_lats = start_lat:(EEG.srate*epoch_length):EEG.pnts;
                        end
                    else
                        event_lats = start_lat:(EEG.srate*epoch_length):(next_lat-EEG.srate);
                    end
                    % append events
                    for segi = 2:length(event_lats) %start from 2, because first event already exists (original trigger event)
                        new_event         = EEG.event(1); %copy template of event to have all the struct-fields EEGlab requires
                        new_event.type    = [cond '+' num2str((segi-1)*epoch_length) 'sec']; %overwrite event type and add time lag to original event
                        new_event.latency = event_lats(segi);  %overwrite event latency, in samples
                        EEG.event(end+1)  = new_event; %append
                    end
                end
            end

            % Make sure events are sorted after inserting
            EEG = eeg_checkset(EEG, 'eventconsistency');


            % -- Save processed EEG set --
            fprintf('\n****\nSaving processed subject %i session %s %s\n****\n\n', subj_list(subj_i), sessions{sess_i}, tasks{task_i});
            SaveName = sprintf( '%i-%s-%s_epochsmarked.set', subj_list(subj_i), sessions{sess_i}, tasks{task_i} );
            EEG      = pop_saveset( EEG,'filename',SaveName,'filepath', path2save );

             % -- Clear this EEG set from workspace before iterating to the next --
            clear EEG
            ALLEEG(1:end) = [];


        end
    end
end


%% Clean data
% Evaluate and remove if necessary bad channels and epochs 

% 1) initiate or load matrix to save noisy channels and rejected epochs
% Check if file exists:
fullPath_badch  = fullfile([ path2save '/Overview_badchannels.txt']);
fullPath_rejdat = fullfile([ path2save '/Overview_rejecteddata.txt']);

if ~exist(fullPath_badch, 'file')
    % When not found, create a new matrix
    fprintf('Bad channels file not found: %s. \nCreating new file.\n', fullPath_badch);
    bdchns      = cell(length(subj_list),length(sessions)*length(tasks)+1);
    bdchns(:,1) = num2cell(subj_list');
else
    % When found, read existing matrix
    fprintf('Bad channels file found. %s. \nReading file.\n', fullPath_badch);
    bdchns      = readcell( fullPath_badch);
    idx_missing = cellfun(@(x) isa(x,'missing'), bdchns);
    bdchns(idx_missing) = {[]};
end

if ~exist(fullPath_rejdat, 'file')
    % When not found, create a new matrix
    fprintf('Rejected epochs file not found: %s. \nCreating new file.\n', fullPath_rejdat);
    rej_data           = [subj_list  nan( length(subj_list),length(sessions)*length(tasks) )];
    % Initiate cell variable to save visual evaluation of rejections by ASR
    rej_data_eval      = cell(length(subj_list),length(sessions)*length(tasks)+1);
    rej_data_eval(:,1) = num2cell(subj_list');
    
else
    % When found, read existing matrix
    fprintf('Rejected epochs file found. %s. \nReading file.\n', fullPath_rejdat);
    rej_data = table2array( readtable( fullPath_rejdat ) );
    % Load evaluation cell variable too
    rej_data_eval = readcell( fullfile([ path2save '/Evaluation_rejecteddata.txt']) );
end


% 2) Loop over files
for subj_i = 1:length(subj_list)
    for sess_i = 1:length(sessions)
        for task_i = 1 %1:length(tasks)

            filename = sprintf( '%i-%s-%s_epochsmarked.set', subj_list(subj_i), sessions{sess_i}, tasks{task_i} );

             % -- Check if file exists --
            fullPath = fullfile(path2save, filename);
            if ~exist(fullPath, 'file')
                fprintf('Not found: File %s.\n Skipping.\n', fullPath);
                continue;  % skip to next iteration of your subject/session/task loop
            end

            % Load EEG set
            EEG = pop_loadset('filename', filename , 'filepath', path2save);

            % % Import a channel location file
            % % (needed? for Artifact Subspace Reconstruction, ASR)
            % EEG = pop_chanedit(EEG, {'lookup','/Users/fsmits2/Downloads/eeglab2024.2/plugins/dipfit/standard_BEM/elec/standard_1005.elc'});

             % -- Check if there is at least 20s of recorded data --
            if EEG.pnts/EEG.srate < 20
                fprintf('Less than 20s recorded data in %s.\n Skipping.\n', fullPath);
                continue;  % skip to next iteration of your subject/session/task loop
            end

            % -- Detect noisy channels --
            badchannels = {[]}; badchannrs = [];
            badchan_idx = sess_i + task_i + (sess_i-1); % defines column in 'bdchns' where bad channel should be written
            if ~ismissing(bdchns{subj_i,sess_i+1})
                % Load bad channels if channel rejection already done
                if bdchns{subj_i,sess_i+1} ~= 0
                    badchannels = {bdchns{subj_i,sess_i+1}};
                end
                % Remove bad channels
                EEG = pop_select(EEG,'nochannel', badchannels);
            else
                % Else: check channel signals and reject bad channels
                pop_eegplot( EEG, 1, 1, 1);
                for chan = 1:EEG.nbchan
                    % - Compute standard deviation
                    chansd = std(EEG.data(chan, :)');
                    fprintf('%f = STD channel %s.\n', chansd, EEG.chanlocs(chan).labels);
                    % - Compute PSD
                    % Compute mean bandpower per channel in PSD 5-34 Hz (34 Hz is based on filter freq, and 5 Hz is based on Delorme, A., & Martin, J. A. (2021, December). Automated data cleaning for the Muse EEG.)
                    % Delorme et al. use mean bandpower >25 log10(µV2)/Hz as threshold.
                    % to suppress pop-up window, pass the window numbers to pwelch
                    window      = min(EEG.srate, size(EEG.data,2));
                    winvec      = hann(window);
                    EEG_datauV  = EEG.data * 1e6;   % convert to µV
                    [pxx,f]     = pwelch(EEG_datauV(chan,:)', winvec, 0, floor(EEG.srate), EEG.srate); %remember to transpose EEG.data, because pwelch function expects an array of timepoints x channels, while EEG.data is structured as channels x timepoints.
                    idx = f >= 5 & f <= 34;
                    band_power = mean(log10(pxx(idx)));
                    fprintf('%f log10(µV2)/Hz = mean PSD channel %s.\n Delorme 2021 threshold: 25 log10(µV2)/Hz \n', band_power, EEG.chanlocs(chan).labels);
                end

                % Visually check noisy channels in plot and decide to keep or reject
                m3a = "no"; m3 = -1;
                while m3a ~= "yes"
                    m3a = input('Ready to input channels to leave out? Type [yes] ','s');
                end
                while m3 == -1
                    m3 = str2double( input('How many channels to leave out? ','s') );
                end                
                if m3 > 0
                    for badchani = 1:m3
                        badchannels{badchani} = input(['Which channel to leave out? nr: ' num2str(badchani) ' ' ],'s') ;
                        badchannrs(badchani)  = find( strcmpi( badchannels{badchani}, {EEG.chanlocs.labels} ));
                    end
                    bdchns{subj_i,sess_i+1}  = string(badchannels);
                    EEG.eventdescription = { {'Too much noise in channels: '} badchannels };
                    % Remove bad channels
                    EEG = pop_select(EEG,'nochannel', badchannels);
                else
                    bdchns{subj_i,sess_i+1}  = '0';
                end
            end
            

            % Save original data (before artifact rejection)
            EEG_orig = EEG;

            % -- Semi-automatic artifact rejection --
            % Artifact Subspace Reconstruction (ASR)
            % Choice for ASR method is based on Delorme, A., & Martin, J. A. (2021, December). Automated data cleaning for the Muse EEG. In 2021 IEEE International Conference on Bioinformatics and Biomedicine (BIBM) (pp. 1-5). IEEE. https://doi.org/10.1109/BIBM52615.2021.9669415
            % Paramaters based on musemonitor code (https://github.com/sccn/eeglab_musemonitor_plugin/blob/master/pop_musemonitor.m)
            % BurstCriterion value: 11 (from above paper: "The best method is the ASR method with a parameter of 11.")
            
            % Note: ASR requires continuous signal (if epoched, use: EEG = eeg_epoch2continuous(EEG); )
            % Check ASR - step 1: Run once with reconstruction: Artifact rejection with ASR in repair mode (see what segments are detected):
            % Note: this will not produce the exact same result as final clean_artifacts code below to segments, because of difference in some parameters here that prevent data rejection (needed for visualisation).
            EEG_rec = clean_artifacts(EEG, ...
                'FlatlineCriterion','off', ...
                'ChannelCriterion','off', ...  % prevent automatic removal of channels
                'LineNoiseCriterion','off', ... % alternative (more stringent): 'LineNoiseCriterion',5, ...
                'Highpass',[0.25 0.75], ...
                'BurstCriterion',11, ... % alternative (less stringent): 'BurstCriterion',25, ...
                'BurstCriterionRefMaxBadChns','off', ...               
                'WindowCriterion','off', ...
                'BurstRejection','off', ...  
                'Distance','Euclidian');

            % % Check ASR - step 2 (optional): Visually check detected bad segments in plot (just verification of the ASR cleaning method)
            % % Visualize artifacts and input your evaluation:
            % m4 = 'N/A';
            % if EEG_rec.pnts == EEG_orig.pnts
            %     vis_artifacts(EEG_rec, EEG_orig, 'WindowLength', 50, 'DisplayMode', 'both');
            %     % Evaluate ASR segment rejection (check if or not ASR also rejects clean segments, f.e. high amplitude alpha oscillations can be confused for noise)
            %     while ~strcmpi(m4,'yes') && ~strcmpi(m4,'no')
            %         m4 = input('Marked segments rejection is OK? Type [yes] or [no]: ','s');
            %     end
            %     EEG.epochdescription          = sprintf('%.1f percent of data rejected by ASR. Were rejections evaluated as acceptable by visual inspection? > %s', perc_rejected, m4);
            %     rej_data_eval{subj_i,rej_idx} = EEG.epochdescription;
            % else
            %     sprintf('No visualization possible: reconstructed data differs %i samples from original data', EEG_orig.pnts - EEG_rec.pnts);
            % end            

            % Check ASR - step 3: Run again without reconstruction (segment rejection) 
            EEG = clean_artifacts(EEG, ...
                'FlatlineCriterion','off', ...
                'ChannelCriterion','off', ...  % prevent automatic removal of channels method 1
                'LineNoiseCriterion','off', ... % prevent automatic removal of channels method 2; alternative (more stringent): 'LineNoiseCriterion',5, ...
                'Highpass',[0.25 0.75], ...
                'BurstCriterion',11, ... % alternative (less stringent): 'BurstCriterion',25, ...
                'BurstCriterionRefMaxBadChns','off', ...               
                'WindowCriterion',0.5, ...
                'BurstRejection','on', ... 
                'Distance','Euclidian');
            
            % Save percentage and evaluation of rejected data by ASR:
            bad_samples              = find(~EEG.etc.clean_sample_mask);
            perc_rejected            = length(bad_samples)/length(EEG.etc.clean_sample_mask) * 100;
            perc_rejected            = round(perc_rejected,1); % round to 1 decimal
            rej_idx                  = sess_i + task_i + (sess_i-1); % defines column in 'rej_data' where percentage should be written
            rej_data(subj_i,rej_idx) = perc_rejected;


            % -- Epoch the data --
            % For data with rejected segments:
            % select the 'epoch'-markers ('open+..sec' and 'close+..sec' for resting-state EEG)
            idx     = startsWith({EEG.event.type}, {'open', 'close'});
            markers = {EEG.event(idx).type};
            % create epochs around markers
            EEG     = pop_epoch(EEG, markers, [0 epoch_length]);

            % Remove epochs with boundary (if not already done during epoching)
            bad_epochs = cellfun(@(epcs) any(strcmp(epcs,'boundary')), {EEG.epoch.eventtype});
            EEG        = pop_rejepoch(EEG, bad_epochs, 0);

            % For reconstructed data:
            % select the 'epoch'-markers ('open+..sec' and 'close+..sec' for resting-state EEG)
            idx     = startsWith({EEG_rec.event.type}, {'open', 'close'});
            markers = {EEG_rec.event(idx).type};
            % create epochs around markers
            EEG_rec = pop_epoch(EEG_rec, markers, [0 epoch_length]);

            % Remove epochs with boundary (if not already done during epoching)
            bad_epochs = cellfun(@(epcs) any(strcmp(epcs,'boundary')), {EEG_rec.epoch.eventtype});
            EEG_rec    = pop_rejepoch(EEG_rec, bad_epochs, 0);


            % -- Save clean dataset --
            fprintf('\n****\nSave clean data subject %i session %s task %s\n****\n\n', subj_list(subj_i), sessions{sess_i}, tasks{task_i});
            % For data with rejected segments:
            SaveName = sprintf( '%i-%s-%s_cleanRejectedsegments.set', subj_list(subj_i), sessions{sess_i}, tasks{task_i} );
            EEG      = pop_saveset( EEG, 'filename',SaveName,'filepath', path2EEGsets );
            % For reconstructed data:
            SaveName_rec = sprintf( '%i-%s-%s_cleanReconstructed.set', subj_list(subj_i), sessions{sess_i}, tasks{task_i} );
            EEG_rec      = pop_saveset( EEG_rec, 'filename',SaveName_rec,'filepath', path2EEGsets );
            % Save bad channels and rejected data percentages
            cd(path2save);
            writecell(       bdchns, [path2save '/Overview_badchannels.txt'], 'Delimiter',',');
            writematrix(   rej_data, [path2save '/Overview_rejecteddata.txt'], 'Delimiter',',');
            writecell(rej_data_eval, [path2save '/Evaluation_rejecteddata.txt'], 'Delimiter',',');

            % m2 = 0;
            % while m2 == 0
            %     m2 = input('Continue? [Y/N] ','s');
            %     if m2 == 'Y'
            %         continue
            %     else
            %         return
            %     end
            % end

            clear EEG
            close all

            % open EEGlab
            cd('/Users/fsmits2/Downloads/eeglab2024.2')
            [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab('nogui');

        end
    end
end

% change directory to where script is
cd('/Users/fsmits2/MATLAB/Projects/PRESSsandbox')