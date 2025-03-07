import time
from pylsl import StreamInfo, StreamOutlet, StreamInlet, resolve_stream, local_clock

# Define EEG sampling rate (256 Hz is typical for Muse)
SFREQ = 256  

# Step 1: Create the marker stream
info = StreamInfo(name='Markers', type='Markers', channel_count=1, nominal_srate=0,
                  channel_format='int32', source_id='marker_stream')
outlet = StreamOutlet(info)

print("Marker stream created. Waiting for EEG stream...")

# Step 2: Wait for EEG stream to appear
streams = resolve_stream('type', 'EEG')  # Wait for EEG to start
eeg_inlet = StreamInlet(streams[0])  # Connect to EEG stream

# Step 3: Wait until EEG stream is ready
print("Waiting for EEG data to start...")
first_sample, first_timestamp = eeg_inlet.pull_sample()
eeg_start_time = first_timestamp  # Use this as the reference time for EEG
print(f"EEG Recording Started at LSL Time: {eeg_start_time}")

print("Sending triggers... Press Ctrl+C to stop.")

try:
    # Start sending triggers after ensuring synchronization
    while True:
        # Get the current LSL time
        current_time = local_clock()

        # Calculate sample index relative to the EEG start time
        sample_index = int((current_time - eeg_start_time) * SFREQ)

        # Skip sending events if sample_index is negative (to avoid out-of-bounds)
        if sample_index < 0:
            continue

        # Define trigger value (change based on experiment)
        trigger_value = 1

        # Send the trigger event with the calculated sample index
        outlet.push_sample([trigger_value])
        print(f"Trigger Sent: {trigger_value} at Sample {sample_index}")

        time.sleep(5)  # Adjust as needed for your experiment's trigger interval
except KeyboardInterrupt:
    print("Stopped sending triggers.")
