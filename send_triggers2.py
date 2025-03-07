import time
import csv
from pylsl import StreamInfo, StreamOutlet, local_clock

# Create a marker stream for event triggers
info = StreamInfo(name='Markers', type='Markers', channel_count=1, nominal_srate=0, 
                  channel_format='int32', source_id='marker_stream')
outlet = StreamOutlet(info)

# Open a CSV file to save triggers
with open('triggers.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["Timestamp", "Marker"])  # Header

    print("Sending triggers... Press Ctrl+C to stop.")

    try:
        while True:
            trigger_value = 1  # Change this based on your experiment

            # Get the accurate LSL timestamp
            timestamp = local_clock()

            # Send trigger with timestamp
            outlet.push_sample([trigger_value], timestamp)

            # Save timestamp and trigger to CSV
            writer.writerow([timestamp, trigger_value])
            print(f"Trigger Sent: {trigger_value} at {timestamp}")

            time.sleep(5)  # Adjust interval as needed

    except KeyboardInterrupt:
        print("Stopped sending triggers.")
