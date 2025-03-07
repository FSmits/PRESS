import mne
import pyxdf  # Import pyxdf directly
from mne.export import export_raw  # Import export_raw for saving .set files

# Load the XDF file
xdf_file = "sub-P009_ses-S004_task-Default_run-001_eeg.xdf"  # Replace with your actual file
streams, header = pyxdf.load_xdf(xdf_file)

# Identify EEG and Marker streams
eeg_stream = None
marker_stream = None

for stream in streams:
    name = stream['info']['name'][0]
    if 'EEG' in name:
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

    # Check if there is a marker stream
    if marker_stream:
        # Extract event markers and timestamps
        event_times = marker_stream['time_stamps']
        event_values = [int(m[0]) for m in marker_stream['time_series']]  # Convert markers to integers

        # Convert event times to sample indices
        event_samples = [int(t * sfreq) for t in event_times]  # Convert time (s) to sample index

        # Create MNE event array (n_events x 3 format: [sample, 0, event_id])
        events = [[sample, 0, value] for sample, value in zip(event_samples, event_values)]
        events = mne.EventsArray(events, raw.info["sfreq"])

        # Add events to raw object
        annotations = mne.Annotations(onset=event_times, duration=[0] * len(event_times), description=[str(v) for v in event_values])
        raw.set_annotations(annotations)

        print("Added event markers to EEG data!")

    # Save EEG data with events as .set file for EEGLAB
    raw.export("sub-P009_ses-S004_task-Default_run-001_eeg.set", fmt="eeglab", overwrite=True)
    print("Saved EEG data with events as .set file!")

    # Save events separately (optional)
    with open("events.txt", "w") as f:
        for t, v in zip(event_times, event_values):
            f.write(f"{t},{v}\n")

    print("Saved events separately as events.txt!")

