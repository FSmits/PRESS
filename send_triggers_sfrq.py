import time
from pylsl import StreamInfo, StreamOutlet, local_clock

# Define EEG sampling rate (must match Muse's actual rate, usually 256 Hz)
SFREQ = 256  

# Create an LSL marker stream
info = StreamInfo(name='Markers', type='Markers', channel_count=1, nominal_srate=0,
                  channel_format='int32', source_id='marker_stream')
outlet = StreamOutlet(info)

print("Sending triggers... Press Ctrl+C to stop.")

# Track the experiment start time
start_time = local_clock()

try:
    while True:
        # Compute sample index relative to the start of the recording
        current_time = local_clock() - start_time
        sample_index = int(current_time * SFREQ)  # Convert seconds to samples
        trigger_value = 1  # Change this based on your experiment
        
        # Send sample index instead of time
        outlet.push_sample([sample_index])
        print(f"Trigger Sent: {trigger_value} at Sample {sample_index}")

        time.sleep(5)  # Adjust as needed
except KeyboardInterrupt:
    print("Stopped sending triggers.")
