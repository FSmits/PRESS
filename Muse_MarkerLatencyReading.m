clear; close all

addpath '/Users/fsmits2/Downloads/'
addpath '/Users/fsmits2/Downloads/eeglab2024.2'

% SubjID recorded with Muse 2022 device at T0 (Meting 1) (list by copy because change directory (cd) function does not work for RFS):
id_t0_muse2022 = ["sub-110002" "sub-110004" "sub-110005" "sub-110007" "sub-110009" ...
    "sub-110011" "sub-110013" "sub-110016" "sub-110018" "sub-110019" "sub-110021" ...
    "sub-110023" "sub-110025" "sub-110027" "sub-110030"];
id_t0_muse2025 = ["sub-110001" "sub-110003" "sub-110006" "sub-110008" "sub-110010" ...
    "sub-110012" "sub-110014" "sub-110015" "sub-110017" "sub-110020" "sub-110022"...
    "sub-110024" "sub-110026" "sub-110028" "sub-110029"];
% merge
id = [id_t0_muse2022 id_t0_muse2025];

% define the recordings (measurements)
recordings = {{'rust' 'rest'}, {'startle'}};

% change filename according to the file you want to load
prefix = '/Volumes/heronderzoek-1/Groep Geuze/25U-0078_PRESS/E_ResearchData/2_ResearchData/Muse/';

% Loop over subjects and recordings
for subj_i = 16:length(id)

    for rec_i = 1:length(recordings)

        % find current subject-id
        subjid = id(subj_i);

        % Find the Muse folder with subject's data (2022 or 2025)
        if any(strcmpi(subjid, id_t0_muse2022))
            muse_folder = 'Muse 2022_T0-Meting1/';
        elseif any(strcmpi(subjid, id_t0_muse2025))
            muse_folder = 'Muse 2025_T0-Meting1/';
        else
            print('Problem: this subj id is not found in any Muse folder')
        end

        % define name of selected recording
        recname = ['ses-M1-' recordings{rec_i}{1}];

        % print current subject-id and recording to command window
        disp([subjid, recname])

        % define full filename
        filename = fullfile( prefix, muse_folder, subjid, recname, 'eeg', [char(subjid) '_' recname '_task-Default_run-001_eeg.xdf']);

        % check if file exists
        if exist(filename, 'file') > 0
            disp('file exists')
            % Note: some resting-state EEG files are saved as "M1-rest" instead of "M1-rust"
        elseif rec_i == 1 && exist(filename, 'file') == 0
            recname = ['ses-M1-' recordings{rec_i}{2}];
            % define full filename
            filename = fullfile( prefix, muse_folder, subjid, recname, 'eeg', [char(subjid) '_' recname '_task-Default_run-001_eeg.xdf']);
            if exist(filename, 'file') > 0
                disp('file exists')
            else
                disp('file does not exist')
            end
        else
            disp('file does not exist')
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
        eeglab
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
        for i = 1:size(EEG.event,2)-1
            EEG.event(i).latency = round(timestamps_delta(i) * srate); %EEGLAB wants latency in sampling points so multiply by sampling rate
        end

        % % view:
        % pop_eegplot( EEG, 1, 1, 1);

        % save as eeglab dataset including the correctly times markers:'
        setname2save = [EEG.setname '.set'];
        path2save = '/Volumes/heronderzoek-1/Groep Geuze/25U-0078_PRESS/E_ResearchData/2_ResearchData/Muse/EEGLAB sets with markers/';
        EEG = pop_saveset( EEG, 'filename',setname2save,'filepath',path2save);

        % clear EEG and close EEGLAB
        clear EEG; close all

    end
end