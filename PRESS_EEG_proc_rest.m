% --------------------------------------------------------------------- %
%   EEG RESTING-STATE EEG PROCESSING
% --------------------------------------------------------------------- %
%
% Description:  Matlab-script to process Muse resting-state EEG data
%             - Study name: PRESS (dossiernr.: 25U-0078)
%             - Measure/instrument: Muse S Headband recordings
%             - Data type: EEG time series with trigger events (time stamps)
%             - Design: within-subjects longitudinal, 3 timepoints (T0 ("M1:") 27th 03/'25, T1 ("M2"): 10th and 11th 04/'25, T1 ("M4"): 24th 04/'25)
%
% Required toolbox:
% EEGLAB (used: version 2024.2)
%
% Notes: 
% Cleaning based on script 'PRESS_EEG_preproc.m'
%
% Date:             March 2026
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

% find path name to research folder structure (RFS)
path2RFS = '/Users/fsmits2/Networkshares/Onderzoek/Groep Geuze/25U-0078_PRESS/'; % Enter your path to RFS. End with slash ('/' on Mac, '\' on Windows)

% set other paths 
path2data    = [path2RFS 'E_ResearchData/2_ResearchData/'];
path2save    = [path2data '1. Verwerkte data/Muse/'];

% enter subject names (extract from key file where subject IDs are linked)
key_filename = [path2data 'SubjectID_koppelbestand_Castor-PVT-Garmin.csv'];
subj_tab     = readtable( key_filename, 'ReadVariableNames', 1);
subj_list    = table2array( subj_tab(:,1) );

% enter session and experimental task names
sessions = {'M1', 'M2', 'M4'};
tasks    = {'rust', 'startle'};


%% Epoch the data

% Loop over subjects and sessions
for subj_i = 1:length(subj_list)
    for sess_i = 1:length(sessions)

        filename = sprintf( '%i-%s-%s_clean.set', subj_list(subj_i), sessions{sess_i}, 'rust' );

        % Check if file exists
        fullPath = fullfile(path2save, filename);
        if ~exist(fullPath, 'file')
            fprintf('Not found: File %s.\n Skipping.\n', fullPath);
            continue;  % skip to next iteration of your subject/session loop
        end

        % Load EEG set
        EEG = pop_loadset('filename', filename , 'filepath', path2save);

        % Epoch the data
        EEG  = pop_epoch( EEG, {'open' 'close'},  [0  2], 'epochinfo', 'yes'); %epoch 0-2 seconds around event

        % Compute spectral power with Welch's method
        [pxx,f] = pwelch(EEG.data, window, 0, [], EEG.srate);
        log_power = log10(pxx);

                % % % Save
                % % fprintf('\n****\nSave epoched data subject %i session %s task %s\n****\n\n', subj_list(subj_i), sessions{sess_i}, tasks{task_i});
                % % SaveName = sprintf( '%i-%s-%s_epoched.set', subj_list(subj_i), sessions{sess_i}, tasks{task_i} );
                % % EEG      = pop_saveset( EEG, 'filename',SaveName,'filepath', path2save );

    end
end

            % for startle task: you want epochs from -1800 to 150 ms around startle probes, in order to quantify and baseline-correct the maximum amplitude between 20 and 120 milliseconds after probe onset
            % Approach startle reflex as an eyeblink ("poor man's ocular-EMG") and try to quantify the amplitude of the blinks 20-120 ms after the probe as
            % detected by a blink algortihm for headband EEG such as Chang et al. 2014: https://doi.org/10.1016/j.cmpb.2015.10.011

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
