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

            % for startle task: you want epochs from -1800 to 150 ms around startle probes, in order to quantify and baseline-correct the maximum amplitude between 20 and 120 milliseconds after probe onset
            % Approach startle reflex as an eyeblink ("poor man's ocular-EMG") and try to quantify the amplitude of the blinks 20-120 ms after the probe as
            % detected by a blink algortihm for headband EEG such as Chang et al. 2014: https://doi.org/10.1016/j.cmpb.2015.10.011


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

% Define epoch lengths (in seconds)
epoch_length = 1; %in seconds


%% Extract spectral readouts

% Pre-specify a matrix (dataframe) to save outcomes in the size: Subject x session x conditions x channels x frequency points 
no_frqpoints  = 256 * 2 + 1; % max. sampling frequency * 2 + 1
no_channels   = 2; % AF7, AF8
no_conditions = 2; % eyes open, eyes closed
datafr_rso    = nan(length(subj_list), length(sessions), no_conditions, no_channels, no_frqpoints); %for resting-state EEG data eyes-open
datafr_rsc    = nan(length(subj_list), length(sessions), no_conditions, no_channels, no_frqpoints); %eyes-closed


% Loop over subjects and sessions
for subj_i = 1:length(subj_list)
    for sess_i = 1:length(sessions)

        filename = sprintf( '%i-%s-%s_cleanRejectedsegments.set', subj_list(subj_i), sessions{sess_i}, 'rust' );

        % Check if file exists
        fullPath = fullfile(path2save, filename);
        if ~exist(fullPath, 'file')
            fprintf('Not found: File %s.\n Skipping.\n', fullPath);
            continue;  % skip to next iteration of your subject/session loop
        end

        % -- Load EEG set -- 
        EEG = pop_loadset('filename', filename , 'filepath', path2save);


        % -- Select channels: AF7, AF8 --
        EEG = pop_select(EEG, 'channel', {'AF7', 'AF8'});
        

        % -- Prepare Fourier Transform --
        % Make sure only one event per epoch is 'read' in case of more events than epochs.
        events_to_keep = [];
        for epoci = 1:length({EEG.epoch.eventtype})
            ev_idx         = EEG.epoch(epoci).event;  % indices into EEG.event
            events_to_keep = [events_to_keep ev_idx(1)];
        end
        % Split dataset into eyes-open and eyes-closed
        openidx  = startsWith({EEG.event(events_to_keep).type}, 'open');
        closeidx = startsWith({EEG.event(events_to_keep).type}, 'close');
        % Put EEG data in double precision for good computation performance
        EEG.datax_o = double( EEG.data(:,:,openidx) ); 
        EEG.datax_c = double( EEG.data(:,:,closeidx) ); 
        nfft        = round(EEG.srate) * 4; % For zero-padding (upsampling) and overlapping in Welch's method
        nOverlap    = size(EEG.datax_c,2)/2; % For no overlap: 0;  for 50% overlap: size(EEG.datax,2)/2
        % Create Hann window to taper the data with.
        hannw       = hann(EEG.pnts); % = .5 * (1 - cos(2*pi*linspace( 0, 1, size(EEG.datax,2) ) ));
        % Pre-specify power spectrum variable to save FFT results
        powspec_o   = nan( EEG.nbchan, round(EEG.srate)*2+1, size(EEG.datax_o,3) ); 
        powspec_c   = nan( EEG.nbchan, round(EEG.srate)*2+1, size(EEG.datax_c,3) ); 


        % -- Do the FFT --
        fprintf('\n****\nCompute power spectrum - subject %i session %s resting-state\n****\n\n', subj_list(subj_i), sessions{sess_i});
        % Loop over channels
        for chani = 1:EEG.nbchan
            % Do the FFT for each frequency and epoch using Welch's method (MATLAB's function pwelch)
            % Only if enough data: >25 1-sec epochs
            if sum(openidx) > 25
                [powspec_o(chani,:,:), hz] = pwelch( squeeze( EEG.datax_o(chani,:,:) ), hannw, nOverlap, nfft, round(EEG.srate) );
            end
            if sum(closeidx) > 25
                [powspec_c(chani,:,:), hz] = pwelch( squeeze( EEG.datax_c(chani,:,:) ), hannw, nOverlap, nfft, round(EEG.srate) );
            end
        end

        % average PSD over trials
        powspec_o = mean(powspec_o, 3);
        powspec_c = mean(powspec_c, 3);

        % log transform the power spectra
        logpowspec_o = log10(powspec_o);
        logpowspec_c = log10(powspec_c);

        % Save powespectra in the dense matrix:
        fprintf('\n*** Save dense power spectrum - subject %i session %s\n', subj_list(subj_i), sessions{sess_i});
        datafr_rso(subj_i, sess_i, 1, 1:EEG.nbchan , 1:size(logpowspec_o,2)) = logpowspec_o;
        datafr_rsc(subj_i, sess_i, 2, 1:EEG.nbchan , 1:size(logpowspec_c,2)) = logpowspec_c;

    end
end

% -- Write to file --
cd(path2save)

filename_hz = 'hz_rust.mat';
save(filename_hz,'hz');

filename_o = 'logpowerspectra_rust_open.mat';
save(filename_o,'datafr_rso', '-v7.3');

filename_c = 'logpowerspectra_rust_closed.mat';
save(filename_c,'datafr_rsc', '-v7.3');



% change directory to where script is
cd('/Users/fsmits2/MATLAB/Projects/PRESSsandbox')