clear; close all

addpath '/Users/fsmits2/Downloads/'
addpath '/Users/fsmits2/Downloads/eeglab2024.2'

% SubjID recorded with Muse 2022 or Muse 2025 device at T0 (Meting 1), T1 (meting2) and t3 (meting 4) (list by copy because change directory (cd) function does not work for RFS so cannot read files?): 
field1 = 'muse2022';
value1 = {[110002 110004 110005 110007 110009 110011 110013 110016 110018 ...
    110019 110021 110023 110025 110027 110030]; % subjIDs M1 muse2022
    [110003 110004 110006 110007 110008 110009 110011 110023 110024 ...
    110026 110027 110028]; % subjIDs M2 muse2022
    [110004 110007 110014 110018 110019 110021 110022 110023 110024 ...
    110025 110027 110028 110029 110030]}; % subjIDs M4 muse2022
field2 = 'muse2025';
value2 = {[110001 110003 110006 110008 110010 110012 110014 110015 110017 ...
    110020 110022 110024 110026 110028 110029]; % subjIDs M1 muse2025
    [110005 110010 110013 110014 110015 110017 110018 110020 110021 ...
    110022 110025 110029]; % subjIDs M2 muse2025
    [110003 110005 110006 110008 110009 110010 110011 110012 110013 ...
    110015 110016 110017 110020 110026];}; % subjIDs M4 muse2025
muse_group = struct(field1, value1, field2, value2);

% enter subject names (extract from key file where subject IDs are linked)
path2RFS     = '/Users/fsmits2/Networkshares/Onderzoek/Groep Geuze/25U-0078_PRESS/'; % Enter your path to RFS. End with slash ('/' on Mac, '\' on Windows)
key_filename = [path2RFS 'E_ResearchData/2_ResearchData/' 'SubjectID_koppelbestand_Castor-PVT-Garmin.csv'];
subj_tab     = readtable( key_filename, 'ReadVariableNames', 1);
subj_list    = table2array( subj_tab(:,1) );

% define session names
sessions            = {'M1', 'M2', 'M4'};
session_foldernames = {'T0-Meting1', 'T1-Meting2', 'T3-Meting4'};

% define the recordings (measurements)
recordings = {{'rust' 'rest'}, {'startle'}};

% change filename according to the file you want to load
prefix = '/Users/fsmits2/Networkshares/Onderzoek/Groep Geuze/25U-0078_PRESS/E_ResearchData/2_ResearchData/0. Ruwe data (niet in werken)/Muse/';

% Search all files in all subfolders (depending on connection, this may take a while)
files = dir(fullfile(prefix,'**','*._task-Default_run-001_eeg.xdf')); % example for EEGLAB files

% Loop over subjects and recordings
for subj_i = 1:length(subj_list)
    for sess_i = 1:length(sessions)
        for rec_i = 1:length(recordings)

            % find current subject-id
            subjid = subj_list(subj_i);

            % find in which muse-subfolder dataset is stored
            if find(subjid == muse_group(sess_i).muse2022) > 0
                muse_folder = ['Muse 2022_' session_foldernames{sess_i} '/'];
            elseif find(subjid == muse_group(sess_i).muse2025) > 0
                muse_folder = ['Muse 2025_' session_foldernames{sess_i} '/'];
            else
                disp('Problem: this subj id is not found in any Muse folder')
            end

            % define name of selected recording
            recname = ['ses-' sessions{sess_i} '-' recordings{rec_i}{1}];

            % print current subject-id and recording to command window
            disp([num2str(subjid) ' ' recname])

            % define full filename
            filename = fullfile( prefix, muse_folder, ['sub-' num2str(subjid)], recname, 'eeg', ['sub-' num2str(subjid) '_' recname '_task-Default_run-001_eeg.xdf']);

            % check if file exists
            if ~exist(filename, 'file')
                fprintf('Not found: File %s.\n', filename);
                % Note: some resting-state EEG files are saved as "M1-rest" instead of "M1-rust"
                if rec_i == 1 && ~exist(filename, 'file')
                    recname = ['ses-' sessions{sess_i} '-' recordings{rec_i}{2}];
                    % redefine full filename
                    filename = fullfile( prefix, muse_folder, ['sub-' num2str(subjid)], recname, 'eeg', ['sub-' num2str(subjid) '_' recname '_task-Default_run-001_eeg.xdf']);
                end
            end
            if ~exist(filename, 'file')
                fprintf('Skipping.\n', filename);
                continue;  % skip to next iteration of your subject/session/task loop
            end


            % Load the streams from xdf file
            streams = load_xdf(filename);

            % find which stream is which
            stream_muse = NaN;
            stream_marker = NaN;
            for i = 1:length(streams)
                if streams{i}.info.name == "Muse"
                    stream_muse = i;
                elseif streams{i}.info.name == "MyMarkerStream"
                    stream_marker = i;
                end
            end

            % Muse has timestamps in unix time (absolute time)
            % Markerstream has timestamps in local time ->> don't use those timestamps.
            % To sync with Muse, Marerstream has unix time of each marker in Marker names (stream{2}.timeseries)
            % In this for-loop, we read the marker unix times from the marker names and save those in an extra created struct in streams{1} called time_stamps_unix
            for i = 1:length(streams{stream_marker}.time_series)
                unix_stamps(i,:) = strsplit( streams{stream_marker}.time_series{i} , ':' );
                streams{stream_marker}.time_stamps_unix(i) = str2double(cell2mat( unix_stamps(i,2) ));
            end

            Markers_time = streams{stream_marker}.time_stamps_unix(1);
            Muse_time2 = str2double(streams{stream_muse}.info.first_timestamp);
            delta = Markers_time - Muse_time2; %in seconds

            timestamps_delta = streams{stream_marker}.time_stamps_unix - Muse_time2;

            % start EEGlab
            eeglab('nogui');
            % If using eeglab GUI:
            % - Stream name to import: Muse
            % - datatset name: [free choice]

            % Now first import .xdf file in eeglab with .xdf plugin
            EEG = pop_loadxdf(filename , 'streamname', 'Muse', 'streamtype', 'EEG', 'exclude_markerstreams', {});

            % Save the
            setname = strsplit(filename, '/');
            EEG.setname = setname{end};

            % Read the EFFECTIVE sampling rate from Muse.
            % NOTE: IMPORTANT to not set equal srate for all files but read this from the muse stream file, because the sampling rate differs per recording. Usually srate is ~256Hz, sometimes errors during recording and srate is lower (~70Hz)
            srate = streams{stream_muse}.info.effective_srate;
            % Change the sampling rate to the effective sampling rate in the eegset
            EEG = pop_editset(EEG, 'srate', srate);
            % Add the marker timestamps based on the effective sampling rate
            for i = 1:length(timestamps_delta)
                EEG.event(i).latency = round(timestamps_delta(i) * srate); %EEGLAB wants latency in sampling points so multiply by sampling rate
            end

            % % view:
            % pop_eegplot( EEG, 1, 1, 1);

            % save as eeglab dataset including the correctly times markers:'
            setname2save = [EEG.setname '.set'];
            path2save = '/Users/fsmits2/Networkshares/Onderzoek/Groep Geuze/25U-0078_PRESS/E_ResearchData/2_ResearchData/0. Ruwe data (niet in werken)/Muse/EEGLAB sets with markers/';
            EEG = pop_saveset( EEG, 'filename',setname2save,'filepath',path2save);

            % clear EEG and close EEGLAB
            clear EEG; close all

        end
    end
end