"""Example program to demonstrate how to send string-valued markers into LSL."""

import random
import time

from pylsl import StreamInfo, StreamOutlet, resolve_stream,local_clock


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

    print("Sending triggers... Press Ctrl+C to stop.")

    markernames = ['Test', 'Blah', 'Marker', 'XXX', 'Testtest', 'Test-1-2-3']
    start_time = local_clock()
    srate = 256
    while True:

        # get a new sample (you can also omit the timestamp part if you're not interested in it)
        timestamp = local_clock()
        elapsed_time = local_clock() - start_time
        latency = srate * elapsed_time
        print(timestamp, elapsed_time, latency)
      
        # pick a sample to send an wait for a bit
        outlet.push_sample([random.choice(markernames)], timestamp)
        time.sleep(3.1)


if __name__ == '__main__':
    main()
