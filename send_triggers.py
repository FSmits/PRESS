import time
from pylsl import StreamInfo, StreamOutlet

# Create a marker stream for event triggers
info = StreamInfo(name='Markers', type='Markers', channel_count=1, nominal_srate=0, 
                  channel_format='int32', source_id='marker_stream')
outlet = StreamOutlet(info)

# Send triggers at specific intervals
print("Sending triggers... Press Ctrl+C to stop.")

try:
    while True:
        trigger_value = 1  # Change this based on your experiment
        outlet.push_sample([trigger_value])
        print(f"Trigger Sent: {trigger_value}")
        time.sleep(5)  # Adjust time interval as needed
except KeyboardInterrupt:
    print("Stopped sending triggers.")
