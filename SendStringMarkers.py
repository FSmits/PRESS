"""Example program to demonstrate how to send string-valued markers into LSL."""

import random
import time
import keyboard

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

    # Wait for user to start the experiment
    input("Press Enter to start sending markers...")

    # Get EEG stream timestamp to synchronize markers
    _, first_eeg_timestamp = eeg_inlet.pull_sample(timeout=5)
    if first_eeg_timestamp is None:
        print("Error: Could not retrieve EEG timestamp.")
        return

    print("Sending triggers... Press Ctrl+C to stop.")

    markernames = ['Test', 'Blah', 'Marker', 'XXX', 'Testtest', 'Test-1-2-3']
    start_time = local_clock()
    srate = 256
    while True:

        # get a new sample (you can also omit the timestamp part if you're not interested in it)
        timestamp = local_clock()
        elapsed_time = timestamp - start_time
        since_eeg_time = first_eeg_timestamp - start_time
        latency = srate * (elapsed_time + since_eeg_time)
        markername = [random.choice(markernames)]
        print(markername, timestamp, elapsed_time, since_eeg_time, latency)
      
        # Combine marker name and latency into a single string
        marker_data = f"{markername}:{latency:.3f}"

        print(f"Sent Marker: {marker_data}, Timestamp: {timestamp}")

        # Send as a single-element list
        outlet.push_sample([marker_data], timestamp)
        
        time.sleep(3.1)


if __name__ == '__main__':
    main()
