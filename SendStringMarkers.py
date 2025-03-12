"""Example program to demonstrate how to send string-valued markers into LSL."""

import random
import time
import keyboard
import datetime

from pylsl import StreamInfo, StreamOutlet, StreamInlet, resolve_stream,local_clock

def wait_for_eeg_stream():
    """Wait for an EEG stream to appear and return its inlet."""
    print("Waiting for EEG stream...")
    streams = resolve_stream('type', 'EEG')  # Search for EEG streams
    inlet = StreamInlet(streams[0])  # Connect to the first found EEG stream
    print("EEG stream found!")
    return inlet

def main():
    # first create a new stream info (here we set the name to MyMarkerStream,
    # the content-type to Markers, 1 channel, irregular sampling rate,
    # and string-valued data) The last value would be the locally unique
    # identifier for the stream as far as available, e.g.
    # program-scriptname-subjectnumber (you could also omit it but interrupted
    # connections wouldn't auto-recover). The important part is that the
    # content-type is set to 'Markers', because then other programs will know how
    #  to interpret the content
    info = StreamInfo('MyMarkerStream', 'Markers', 1, 0, 'string', 'myuidw43536')

    # next make an outlet
    outlet = StreamOutlet(info)

    print("Marker stream created. Waiting for EEG stream...")

    # Wait for EEG stream before sending triggers
    eeg_inlet = wait_for_eeg_stream()

    # Get EEG stream timestamp to synchronize markers
    _, first_eeg_timestamp = eeg_inlet.pull_sample(timeout=5)
    if first_eeg_timestamp is None:
        print("Error: Could not retrieve EEG timestamp.")
        return

    eeg_start_time, _ = eeg_inlet.pull_sample()
    
    # Get EEG stream timestamp to synchronize markers (previously by local_clock: start_time = local_clock() )
    start_time = time.time()
    # Calculate the difference (delta) between the sample's timestamp and the current LSL time
    delta = first_eeg_timestamp - start_time
    # Convert the LSL timestamp to the actual time by adjusting with the delta
    real_time = time.time() - delta  # Real-world time in seconds

    markernames = ['Test', 'Blah', 'Marker', 'XXX', 'Testtest', 'Test-1-2-3']
    srate = 256

    # Wait for user to start the experiment
    input("Press Enter to start sending markers...")

    print("Sending triggers... Press Ctrl+C to stop.")
   
    while True:

        # get a new sample (you can also omit the timestamp part if you're not interested in it)
        # Get EEG stream timestamp to synchronize markers (previously by local_clock: timestamp = local_clock() )
        timestamp = time.time()
        elapsed_time = timestamp - first_eeg_timestamp
        latency = srate * elapsed_time
        markername = [random.choice(markernames)]
        current_time = local_clock()
        # sample_index = (current_time - eeg_start_time) * srate
        print(eeg_start_time, markername, timestamp, elapsed_time, first_eeg_timestamp, latency)
      
        # Combine marker name and latency into a single string
        marker_data = f"{markername}:{timestamp:.3f}"

        print(f"Sent Marker: {marker_data}, Timestamp: {timestamp}")

        # Send as a single-element list
        outlet.push_sample([latency])
        
        time.sleep(1.7)


if __name__ == '__main__':
    main()
