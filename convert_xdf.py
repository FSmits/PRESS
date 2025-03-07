import mne
from mne.externals.pyxdf import load_xdf  # Correct way to load XDF files in MNE

# Load the XDF file
xdf_file = "sub-P009_ses-S004_task-Default_run-001_eeg.xdf"  # Replace with your actual file
streams, _ = load_xdf(xdf_file)

# Find EEG and Marker streams
eeg_stream = None
marker_stream = None

for stream in streams:
    name = stream['info']['name'][0]
    if 'Muse' in name:
        eeg_stream = stream
    elif 'Marker' in name or 'Markers' in name:
        marker_stream = stream

# Convert EEG data to MNE format
if eeg_stream:
    sfreq = float(eeg_stream['info']['nominal_srate'][0])  # Sampling rate
    ch_names = [ch['label'][0] for ch in eeg_stream['info']['desc'][0]['channels'][0]['channel']]
    ch_types = ['eeg'] * len(ch_names)  # Assume all are EEG channels

    raw = mne.io.RawArray(eeg_stream['time_series'].T, mne.create_info(ch_names, sfreq, ch_types))

    # Save as EEGLAB file (.set)
    raw.save("sub-P009_ses-S004_task-Default_run-001_eeg.set", overwrite=True)
    print("Saved EEG data as .set file!")

# Handle event markers
if marker_stream:
    event_times = marker_stream['time_stamps']
    event_values = [int(m[0]) for m in marker_stream['time_series']]  # Convert to integers
    events = mne.events_from_annotations(raw)  # Create MNE event structure

    # Save events separately if needed
    with open("events.txt", "w") as f:
        for t, v in zip(event_times, event_values):
            f.write(f"{t},{v}\n")
    print("Saved events separately as events.txt!")


