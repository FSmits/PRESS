import time
from pylsl import StreamInfo, StreamOutlet, StreamInlet, resolve_stream, local_clock

# Define EEG sampling rate (must match Muse's actual rate, usually 256 Hz)
SFREQ = 256  

# Wait for the Muse EEG stream to appear on LSL
print("Looking for EEG stream...")
streams = resolve_stream('type', 'EEG')
eeg_inlet = StreamInlet(streams[0])  # Connect to the EEG stream

# Get the first EEG sample timestamp as the start time
eeg_start_time, _ = eeg_inlet.pull_sample()
print(f"EEG Recording Start Time: {eeg_start_time}")

# Create an LSL marker stream
info = StreamInfo(name='Markers', type='Markers', channel_count=1, nominal_srate=0,
                  channel_format='int32', source_id='marker_stream')
outlet = StreamOutlet(info)

print("Sending triggers... Press Ctrl+C to stop.")

try:
    while True:
        # Get current LSL time
        current_time = local_clock()

        # Convert to sample index relative to EEG stream
        sample_index = int((current_time - eeg_start_time) * SFREQ)
        trigger_value = 1  # Change this based on your experiment
        
        # Send the event with corrected sample index
        outlet.push_sample([sample_index])
        print(f"Trigger Sent: {trigger_value} at Sample {sample_index}")

        time.sleep(5)  # Adjust as needed
except KeyboardInterrupt:
    print("Stopped sending triggers.")
