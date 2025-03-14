"""Example program to demonstrate how to send string-valued markers into LSL."""

import random
import time
import keyboard
import datetime
import simpleaudio as sa

from pylsl import StreamInfo, StreamOutlet, StreamInlet, resolve_stream,local_clock

def play_beep():  
    """Plays a beep sound from beep.wav"""
    try:
        wave_obj = sa.WaveObject.from_wave_file("beep.wav")
        play_obj = wave_obj.play()
        play_obj.wait_done()  # Wait for the sound to finish
    except Exception as e:
        print(f"Error playing beep: {e}")

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
    _, first_eeg_timestamp = eeg_inlet.pull_sample(timeout=10)
    if first_eeg_timestamp is None:
        print("Error: Could not retrieve EEG timestamp.")
        return

    markernames = ['open', 'close', 'open', 'close', 'end']
    markersleep = [60, 60, 60, 60, 0.01]
    num_markers = 5
    srate = 256

    # Wait for user to start the experiment
    input("Press Enter to start sending markers...")

    print("Sending triggers... Press Ctrl+C to stop.")
   
    try:
        for i in range(num_markers):

            # get a new sample (you can also omit the timestamp part if you're not interested in it)
            timestamp = time.time()
            elapsed_time = timestamp - first_eeg_timestamp
            latency = srate * elapsed_time

            # Select marker based on count
            markername = markernames[i % len(markernames)]
      
            # Combine marker name and latency into a single string
            marker_data = f"{markername}:{timestamp:.3f}"

            print(f"Sent Marker: {marker_data}, Timestamp: {timestamp}")

            # Send as a single-element list
            outlet.push_sample([marker_data])

            play_beep()  # Play beep sound

            # Select sleep time based on predefined list
            time2sleep = markersleep[i % len(markersleep)]
            print(f"Sleeping for {time2sleep} seconds before next beep...")
            time.sleep(time2sleep)

        print("Experiment complete. 20 startle probes sent.")

    except KeyboardInterrupt:
        print("\nExperiment stopped.")


if __name__ == '__main__':
    main()
