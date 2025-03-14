addpath('/Users/kcmggz/Downloads/load_xdf.m')

% change filename according to the file you want to load
filename = '/Users/kcmggz/Documents/CurrentStudy/sub-P20/ses-4/eeg/sub-P20_ses-4_task-Default_run-001_eeg.xdf';
streams = load_xdf(filename);

% Muse has timestamps in unix time (absolute time)
% Markerstream has timestamps in local time ->> don't use those timestamps.
% To sync with Muse, Marerstream has unix time of each marker in Marker names (stream{1}.timeseries)
% In this for-loop, we read the marker unix times from the marker names and save those in an extra created struct in streams{1} called time_stamps_unix
for i = 1:length(streams{1}.time_series)
    unix_stamps(i,:) = strsplit( streams{1}.time_series{i} , ':' );
    streams{1}.time_stamps_unix(i) = str2double(cell2mat( unix_stamps(i,2) ));
end

Markers_time = streams{1}.time_stamps_unix(1);
%Muse_time = streams{2}.time_stamps(1); %DON'T USE, THIS IS THE WRONG 1stSAMPLE TIME
Muse_time2 = str2double(streams{2}.info.first_timestamp);
delta = Markers_time - Muse_time2; %in seconds

timestamps_delta = streams{1}.time_stamps_unix - Muse_time2;

% Now first import .xdf file in eeglab with .xdf plugin
EEG = pop_loadxdf(filename , 'streamname', 'Muse', 'streamtype', 'EEG', 'exclude_markerstreams', {});

% Save the 
setname = strsplit(filename, '/');
EEG.setname = setname{end};

% - Stream name to import: Muse
% - datatset name: [free choice]
srate = 256; % sampling rate Muse
for i = 1:size(EEG.event,2)-1
    EEG.event(i).latency = round(timestamps_delta(i) * srate); %EEGLAB wants latency in sampling points so multiply by sampling rate
end

