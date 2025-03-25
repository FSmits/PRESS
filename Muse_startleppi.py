import time
import random
import numpy as np
import simpleaudio as sa
from pylsl import StreamInfo, StreamOutlet, StreamInlet, resolve_stream, local_clock

def generate_startle_probe(duration_ms=50, sample_rate=44100, amplitude=1.0):
    """
    Generates a white noise burst at the specified amplitude (default: 100 dB).
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

def habituation_phase(outlet):
    """
    Runs 6 habituation startle probes before the experiment begins.
    """
    num_habituation_probes = 6
    print("Starting habituation phase: 6 startle probes before main experiment.")

    for i in range(num_habituation_probes):
 
        # Intertrial interval (ITI) (10-20s)
        iti = random.uniform(10, 20)

        if i > 5:  
            iti = 0.5  

        print(f"Waiting {iti:.2f} seconds for next habituation probe...\n")
        time.sleep(iti)

        timestamp = time.time()
        print(f"Habituation Probe {i+1}/{num_habituation_probes}, Timestamp: {timestamp:.3f}")

        # Send LSL marker
        outlet.push_sample([f"habituation_probe:{timestamp:.3f}"])

        # Play startle probe (100 dB)
        generate_startle_probe(amplitude=2.99)

    print("Habituation phase complete.\n")

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

    print("Press Enter to start habituation startles")
    input()

    # Run habituation phase
    habituation_phase(outlet)

    print("Press Enter to start the 5-minute startle probe experiment...")
    input()

    num_probes = 20  # Total number of startle probes
    ppi_trials = random.sample(range(num_probes), 10)  # Select 10 trials for PPI
    print(f"PPI trials will occur at indices: {ppi_trials}")

    print("Experiment started. Sending 20 probes over 5 minutes...")

    # Wait for next probe (random ITI between 10-20 seconds)
    iti = random.uniform(10, 20)
    print(f"Waiting {iti:.2f} seconds for next trial...\n")
    time.sleep(iti)

    try:
        for i in range(num_probes):
            timestamp = time.time()
            elapsed_time = timestamp - first_eeg_timestamp

            if i in ppi_trials:
                # Pre-Pulse (70 dB, 50 ms) occurs 100 ms before the main probe
                print(f"Trial {i+1}/{num_probes}: Pre-Pulse (70 dB) + Startle Probe (100 dB), Timestamp: {timestamp:.3f}")
                outlet.push_sample([f"pre-pulse:{timestamp:.3f}"])
                generate_startle_probe(amplitude=0.025)  # 70 dB pre-pulse
                time.sleep(0.12)  # 120 ms gap before main probe

            # Main Startle Probe (100 dB, 50 ms)
            print(f"Trial {i+1}/{num_probes}: Startle Probe (100 dB), Timestamp: {timestamp:.3f}")
            outlet.push_sample([f"probe:{timestamp:.3f}"])
            generate_startle_probe(amplitude=2.99)  # 100 dB startle probe

            # Wait for next probe (random ITI between 10-20 seconds)
            iti = random.uniform(10, 20)
            print(f"Waiting {iti:.2f} seconds for next trial...\n")
            time.sleep(iti)

        print("End of Startle paradigm.")

    except KeyboardInterrupt:
        print("\nExperiment stopped.")

if __name__ == '__main__':
    main()
