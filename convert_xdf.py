from mne.io import read_raw_xdf

# Load the XDF file
raw = read_raw_xdf("sub-P009_ses-S004_task-Default_run-001_eeg.xdf")

# Convert to EEGLAB format and save as .set file
raw.export("sub-P009_ses-S004_task-Default_run-001_eeg.set", fmt="eeglab")
