from pyxdf import load_xdf

# Load XDF file
xdf_file = "your_recording.xdf"  # Replace with your file
streams, _ = load_xdf(xdf_file)

# Print available streams
for stream in streams:
    print(f"Stream Name: {stream['info']['name'][0]}")
