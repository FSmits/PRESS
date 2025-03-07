from mne.io import read_raw_xdf

# Load the XDF file
raw = read_raw_xdf("your_recording.xdf")

# Convert to EEGLAB format and save as .set file
raw.export("your_recording.set", fmt="eeglab")
