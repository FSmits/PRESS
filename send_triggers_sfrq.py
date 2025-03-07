import time
from pylsl import StreamInfo, StreamOutlet, StreamInlet, resolve_stream, local_clock

# Define EEG sampling rate (must match Muse's actual rate, usually 256 Hz)
SFREQ = 256  

# Step 1: Create the marker stream FIRST
info = StreamInfo(name='Markers', type='Markers', channel_count=1, nominal_srate=0,
                  channel_format='int32', source_id='marker_stream')
outlet = StreamOutlet(info)

print("Marker stream created. Waiting for EEG stream...")

# Step 2: Wait for EEG stream to appear
streams = resolve_stream('type', 'EEG')  # Wait for EEG to start
eeg_inlet = StreamInlet(streams[0])  # Connect to EEG stream

# Step 3: Get the first EEG sample timestamp as the start time
first_sample, first_timestamp = eeg_inlet.pull_sample()
eeg_start_time = first_timestamp  # This is the exact start time of EEG recording
print(f"EEG Recording Started at LSL Time: {eeg_start_time}")

print("Sending triggers... Press Ctrl+C to stop.")

try:
    while True:
        # Get current LSL time
        current_time = local_clock()

        # Compute sample index relative to EEG start
        sample_index = int((current_time - eeg_start_time) * SFREQ)
        trigger_value = 1  # Change this based on your experiment
        
        # Send the event with corrected sample index
        outlet.push_sample([sample_index])
        print(f"Trigger Sent: {trigger_value} at Sample {sample_index}")

        time.sleep(5)  # Adjust as needed
except KeyboardInterrupt:
    print("Stopped sending triggers.")
