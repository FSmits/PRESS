import time
import random
import numpy as np
import simpleaudio as sa
from pylsl import StreamInfo, StreamOutlet, StreamInlet, resolve_stream, local_clock

def generate_startle_probe(duration_ms=50, sample_rate=44100, amplitude=1.0):
    """
    Generates a 50 ms white noise burst (100 dB) as an audio buffer.
    """
    num_samples = int(sample_rate * (duration_ms / 1000.0))
    noise = np.random.uniform(-amplitude, amplitude, num_samples).astype(np.float32)
    
    # Convert to 16-bit PCM format
    audio_data = (noise * 32767).astype(np.int16)
    return sa.play_buffer(audio_data, 1, 2, sample_rate)

def wait_for_eeg_stream():
    """Wait for an EEG stream to appear and return its inlet."""
    print("Waiting for EEG stream...")
    streams = resolve_stream('type', 'EEG')  # Search for EEG streams
    if not streams:
        print("Error: No EEG stream found!")
        return None
    inlet = StreamInlet(streams[0])  # Connect to the first found EEG stream
    print("EEG stream found!")
    return inlet

def main():
    # Create LSL Marker Stream
    info = StreamInfo('MyMarkerStream', 'Markers', 1, 0, 'string', 'myuid_startle')
    outlet = StreamOutlet(info)

    print("Marker stream created. Waiting for EEG stream...")

    # Wait for EEG stream
    eeg_inlet = wait_for_eeg_stream()
    if eeg_inlet is None:
        return  # Exit if no EEG stream is found

    # Get EEG stream timestamp to synchronize markers
    sample, first_eeg_timestamp = eeg_inlet.pull_sample(timeout=10)
    if first_eeg_timestamp is None:
        print("Error: Could not retrieve EEG timestamp.")
        return

    print("Press Enter to start the 5-minute startle probe experiment...")
    input()

    num_probes = 20  # Total number of startle probes
    markername = "probe"

    print("Experiment started. Sending 20 probes over 5 minutes...")

    try:
        for i in range(num_probes):
            timestamp = time.time()
            elapsed_time = timestamp - first_eeg_timestamp
            
            print(f"Probe {i+1}/{num_probes}: {markername}, Timestamp: {timestamp:.3f}, Clock: {clocktime:.3f}")

            # Send LSL marker
            outlet.push_sample([f"{markername}:{timestamp:.3f}"])

            # Play the startle probe (100 dB white noise for 50 ms)
            generate_startle_probe()

            # Wait for next probe (random ITI between 10-20 seconds)
            iti = random.uniform(10, 20)
            time.sleep(iti)

        print("End of Startle paradigm.")

    except KeyboardInterrupt:
        print("\nExperiment stopped.")

if __name__ == '__main__':
    main()
