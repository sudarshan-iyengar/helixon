## Code Description for the SDM method and the Binaural Rendering

The "Code SDM and Auralization" folder is created starting from the SDM toolbox in Matlab.
This repository contains MATLAB scripts for analyzing Spatial Room Impulse Responses (SRIRs) using the Spatial Decomposition Method (SDM), specifically tailored for the GRASVI50 dataset. The scripts are an adaptation and extension of the standard SDM toolbox, developed for research on non-cooperative binaural auralization and SRIR interpolation techniques.

### Prerequisites

- MATLAB (R2021a or later recommended)
- SOFA MATLAB/Octave API
- GRASVI50 impulse response dataset

### Installation

- Clone or download this repository.
- Install the SOFA API for MATLAB/Octave from [SOFA API](https://github.com/sofacoustics/SOFAtoolbox).
- Place the GRASVI50 impulse response dataset in a known directory.
### Repository Structure
Code SDM and Auralization: Contains modified MATLAB files from the SDM toolbox.
- **spai_srir.m**: Calculates and visualizes the SRIR from our dataset.
- **spai_srir_b_auralization.m**: Extends spai_srir.m for binaural rendering.
- **combining_auralization.m**: Combines all auralized SRIRs into one WAV file.
- **equalize_norm.m**: Equalizes the samples.
- **extract_json.m**: Extracts interpolated SRIRs from JSON files.
- **filterData.m**: Filters data for DOA and pressure of interpolated SRIRs.
- **TestingSingleAuralization.m**: Tests the auralization with one SRIR.
### Folders:
- **Measurements**: Contains 6 measured RIRs at every discrete position.
- **Interpolated_points**: Stores interpolated SRIRs in JSON format.
- **Figures**: Includes spatio-temporal visualizations for each section plane.
- **AurSRIR**: Contains auralized SRIR WAV files and combined file (AurSRIR.wav).
- **AurInterpolated30cm**: Auralized SRIRs with 30cm resolution and combined WAV.
- **AurInterpolated210cm**: Auralized SRIRs with 210cm resolution and combined WAV.

### WAV Files:
- **120 BPM - ROCK.wav**: Music signal used during measurements.
- **inv_chirp_signal.wav**: Inverse chirp signal.
- **measured_signal.wav**: Measured music signal.

### Usage
- Set your MATLAB folder to the repository's location.
- Modify currentFolder variable in scripts to point to your dataset location.
- Run the desired scripts to perform SRIR calculations, visualizations, and auralizations.
- Results and visualizations are stored in respective output directories.

### Contributing

We welcome contributions to enhance and extend this toolbox. Please maintain the existing structure for consistency.
