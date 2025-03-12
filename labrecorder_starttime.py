from pylsl import StreamInfo, StreamInlet, StreamOutlet, resolve_stream, local_clock
import time

# Resolve the EEG stream
print("Looking for EEG stream...")
streams = resolve_stream('type', 'EEG')  # Change 'EEG' to whatever your EEG stream type is

# Get the first EEG inlet
inlet = streams[0]

# Wait for LabRecorder to start recording
print("Waiting for the first EEG sample...")
sample, timestamp = inlet.pull_sample(timeout=10)  # Adjust timeout as necessary

# When the first sample is received, LabRecorder has started
print(f"LabRecorder started at {timestamp:.6f} seconds")

# Now send the 'StartRecording' marker with the correct timestamp
send_marker("StartRecording")

# Continue with other markers and EEG data
