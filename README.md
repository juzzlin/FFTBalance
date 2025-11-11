# FFTBalance

FFTBalance is an **experimental** command-line tool for balancing the frequency spectrum of audio files using FFT (Fast Fourier Transform). It takes an input audio file and a reference audio file, analyzes their frequency content, and applies a gain adjustment to the input file to match the spectral characteristics of the reference file.

## Prerequisites

To build and run FFTBalance, you need the following libraries installed on your system:

*   **libsndfile**: A C library for reading and writing files containing sampled sound (like WAV, AIFF, FLAC, Ogg/Vorbis).
*   **FFTW3**: The Fastest Fourier Transform in the West, a C subroutine library for computing the discrete Fourier transform (DFT).
*   **CMake**: A cross-platform free software tool designed to simplify the process of compilation of source code.
*   **A C++ compiler**: (e.g., GCC, Clang)

## Build Instructions

1.  **Clone the repository:**
    ```bash
    git clone <repository_url>
    cd FFTBalance
    ```

2.  **Create a build directory and navigate into it:**
    ```bash
    mkdir build && cd build
    ```

3.  **Run CMake to configure the project:**
    ```bash
    cmake ..
    ```

4.  **Build the project:**
    ```bash
    make
    ```
    This will compile the `fft_balance` executable in the `build` directory.

## Run Instructions

After building, you can run the `fft_balance` executable from the `build` directory.

```bash
./fft_balance <input_audio_file> <reference_audio_file> <output_audio_file> [num_bands] [max_boost_db] [max_cut_db]
```

*   `<input_audio_file>`: The audio file to be processed.
*   `<reference_audio_file>`: The audio file whose spectral balance will be used as a reference.
*   `<output_audio_file>`: The path where the balanced audio file will be saved.
*   `[num_bands]`: Optional. The number of frequency bands to use for analysis (default: 8, range: 2-30).
*   `[max_boost_db]`: Optional. The maximum gain in dB that can be applied to a band (default: 3.0, must be non-negative).
*   `[max_cut_db]`: Optional. The maximum cut in dB that can be applied to a band (default: 3.0, must be non-negative).

**Examples:**

```bash
# Basic usage with default settings
./fft_balance Input.wav Reference.wav Balanced.wav

# Using 16 frequency bands
./fft_balance Input.wav Reference.wav Balanced_16_bands.wav 16

# Using 10 frequency bands, with a max boost of 6 dB and max cut of 3 dB
./fft_balance Input.wav Reference.wav Balanced_10_bands_6_3.wav 10 6.0 3.0
```

```
Loading input track: Input.wav...
Loading reference track: Reference.wav...
Spectral analysis parameters: Sample Rate=48000 Hz, Channels=2, Frames=10563584, Bands=8
Spectral analysis complete. Calculating gains...
Input Overall Avg Band RMS: 122306
Reference Overall Avg Band RMS: 886528
Applying level match factor (Ref * 0.137961) to normalize spectra before EQ.
Using Dynamic Amplitude Floor: 0.0886528
Using Gain Limits: Boost=3 dB (factor 1.41254), Cut=3 dB (factor 0.707946).
Band 1 (20 Hz - 47.4275 Hz): Gain = -3 dB
Band 2 (47.4275 Hz - 112.468 Hz): Gain = 2.14909 dB
Band 3 (112.468 Hz - 266.704 Hz): Gain = 2.52302 dB
Band 4 (266.704 Hz - 632.456 Hz): Gain = -0.771232 dB
Band 5 (632.456 Hz - 1499.79 Hz): Gain = -3 dB
Band 6 (1499.79 Hz - 3556.56 Hz): Gain = 0.02975 dB
Band 7 (3556.56 Hz - 8433.93 Hz): Gain = 2.14274 dB
Band 8 (8433.93 Hz - 20000 Hz): Gain = 0.525492 dB
Applying EQ to 2 channels...
Applying final output scaling factor: 1
Successfully matched spectral balance and wrote to: Balanced.wav
```

