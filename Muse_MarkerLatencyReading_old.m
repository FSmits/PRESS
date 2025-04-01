addpath('/Users/kcmggz/Downloads/')

% change filename according to the file you want to load
filename = '/Users/kcmggz/Documents/CurrentStudy/sub-P20/ses-4/eeg/sub-P20_ses-4_task-Default_run-001_eeg.xdf';
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
% To sync with Muse, Marerstream has unix time of each marker in Marker names (stream{1}.timeseries)
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

% view:
pop_eegplot( EEG, 1, 1, 1);
