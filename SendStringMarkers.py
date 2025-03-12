"""Example program to demonstrate how to send string-valued markers into LSL."""

import random
import time

from pylsl import StreamInfo, StreamOutlet, resolve_streamt


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

    # Step 2: Wait for EEG stream to appear
    streams = resolve_stream('type', 'EEG')  # Wait for EEG to start
    eeg_inlet = StreamInlet(streams[0])  # Connect to EEG stream

    # Step 3: Wait until EEG stream is ready
    print("Waiting for EEG data to start...")
    first_sample, first_timestamp = eeg_inlet.pull_sample()
    start_time = first_timestamp  # Use this as the reference time for EEG
    print(f"EEG Recording Started at LSL Time: {eeg_start_time}")

    print("Sending triggers... Press Ctrl+C to stop.")

    markernames = ['Test', 'Blah', 'Marker', 'XXX', 'Testtest', 'Test-1-2-3']
    srate = 256
    while True:

        # get a new sample (you can also omit the timestamp part if you're not interested in it)
        timestamp = local_clock()
        elapsed_time = local_clock() - start_time
        latency = srate * elapsed_time
        print(timestamp, elapsed_time, latency)
      
        # pick a sample to send an wait for a bit
        outlet.push_sample([random.choice(markernames)], latency)
        time.sleep(3.1)


if __name__ == '__main__':
    main()
