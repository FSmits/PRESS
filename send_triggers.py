import time
import pandas as pd
from pylsl import StreamInfo, StreamOutlet, resolve_stream, StreamInlet

# Create an LSL stream for event markers
info = StreamInfo('Markers', 'Markers', 1, 0, 'int32', 'marker_stream')
marker_outlet = StreamOutlet(info)

# Resolve EEG stream
print("Looking for an EEG stream...")
streams = resolve_stream('type', 'EEG')
eeg_inlet = StreamInlet(streams[0])

data_list = []

try:
    while True:
        # Get EEG data sample
        eeg_sample, timestamp = eeg_inlet.pull_sample()

        # Send trigger randomly (Example: Every 5 seconds)
        if int(time.time()) % 5 == 0:
            marker_outlet.push_sample([1])  # Send event marker "1"
            event_marker = 1
        else:
            event_marker = 0

        # Store data in a list
        data_list.append([timestamp] + eeg_sample + [event_marker])

except KeyboardInterrupt:
    print("Stopping recording...")

# Save data to CSV
df = pd.DataFrame(data_list, columns=['Timestamp', 'EEG1', 'EEG2', 'EEG3', 'EEG4', 'Marker'])
df.to_csv('muse_eeg_with_markers.csv', index=False)
print("Data saved as 'muse_eeg_with_markers.csv'")
