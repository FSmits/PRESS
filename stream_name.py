from pyxdf import load_xdf

# Load XDF file
xdf_file = "/Users/kcmggz/Documents/CurrentStudy/sub-P009/ses-S004/eeg/sub-P009_ses-S004_task-Default_run-001_eeg.xdf"  # Replace with your file
streams, _ = load_xdf(xdf_file)

# Print available streams
for stream in streams:
    print(f"Stream Name: {stream['info']['name'][0]}")
