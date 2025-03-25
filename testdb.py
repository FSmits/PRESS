import numpy as np
import simpleaudio as sa
import time

def generate_white_noise(duration_s, sample_rate=44100, amplitude=1.0):
    """
    Generates and plays white noise for a given duration.
    
    Parameters:
    - duration_s: Duration of the noise in seconds.
    - sample_rate: Sampling rate in Hz.
    - amplitude: Amplitude of the noise (0.316 for ~70 dB, 1.0 for ~100 dB).
    """
    num_samples = int(sample_rate * duration_s)
    noise = np.random.uniform(-amplitude, amplitude, num_samples).astype(np.float32)
    
    # Convert to 16-bit PCM format
    audio_data = (noise * 32767).astype(np.int16)
    
    # Play sound
    play_obj = sa.play_buffer(audio_data, 1, 2, sample_rate)
    play_obj.wait_done()

# Play 70 dB noise for 5 seconds
print("Playing 70 dB white noise...")
generate_white_noise(duration_s=5, amplitude=0.316)

# Pause before next sound
time.sleep(2)

# Play 100 dB noise for 5 seconds
print("Playing 100 dB white noise...")
generate_white_noise(duration_s=5, amplitude=1.0)

print("Test complete.")
