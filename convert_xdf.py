import mne
import pyxdf  # Import pyxdf directly

# Load the XDF file
xdf_file = "sub-P009_ses-S004_task-Default_run-001_eeg.xdf"  # Replace with your actual file
streams, header = pyxdf.load_xdf(xdf_file)

# Identify EEG and Marker streams
eeg_stream = None
marker_stream = None

for stream in streams:
    name = stream['info']['name'][0]
    if 'Muse' in name:
        eeg_stream = stream
    elif 'Marker' in name or 'Markers' in name:
        marker_stream = stream

if eeg_stream:
    # Extract EEG data
    sfreq = float(eeg_stream['info']['nominal_srate'][0])  # Sampling rate
    ch_names = [ch['label'][0] for ch in eeg_stream['info']['desc'][0]['channels'][0]['channel']]
    ch_types = ['eeg'] * len(ch_names)  # Assuming all channels are EEG

    # Convert to MNE Raw object
    raw = mne.io.RawArray(
        data=eeg_stream['time_series'].T, 
        info=mne.create_info(ch_names, sfreq, ch_types)
    )

    # Save EEG data as .set file for EEGLAB
    raw.save("sub-P009_ses-S004_task-Default_run-001_eeg.set", fmt="eeglab", overwrite=True)
    print("Saved EEG data as .set file!")

if marker_stream:
    # Extract event markers and timestamps
    event_times = marker_stream['time_stamps']
    event_values = [int(m[0]) for m in marker_stream['time_series']]  # Convert to integers

    # Convert events to MNE format
    events = [[int(t * sfreq), 0, v] for t, v in zip(event_times, event_values)]
    events = mne.events_from_annotations(raw)  # Create MNE event structure

    # Save events separately
    with open("events.txt", "w") as f:
        for t, v in zip(event_times, event_values):
            f.write(f"{t},{v}\n")

    print("Saved events separately as events.txt!")


