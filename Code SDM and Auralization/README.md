## Code Description for the SDM method and the Binaural Rendering

The "Code SDM and Auralization" folder is created starting from the SDM toolbox in Matlab.
It contains all files from the Matlab toolbox in which changes were made to adapt it to our needs. 
Furthermore, it also contains new Matlab files, folders with the data of the measurements, figures etc. and wav files. 
### New Matlab files
- spai_srir.m: Calculates the SRIR from our measured dataset and visualizes it
- spai_srir_b_auralization.m: Extension of the spai_srir.m file, calculates the SRIR and does the binaural rendering for each SRIR
- combining_auralization.m: Combines all auralized SRIRs to one big wav file
- equalize_norm.m: Equalizes the the samples 
- extract_json.m: Extracts the json files where the interpolated SRIRs are stored in
- filterData.m: Filters the data to get the DOA and the pressure of the interpolated SRIRs
- TestingSingleAuralization.m: Test file to test the auralization with one SRIR
### Folders
- Measurements: Contains the 6 measured RIRs at every discrete position 
- Interpolated_points: Contains the interpolated SRIRs for a resolution of respectively 30cm and 210cm in json format
- Figures: Contains the spatio temporal visualization for each section plane (Transverse, Median and Lateral)
- AurSRIR: Contains the auralized SRIRs wav files without interpolation and the combined wav file (AurSRIR.wav)
- AurInterpolated30cm: Contains the auralized SRIRs wav files with interpolation (resolution 30cm) and the combined wav file (AurInterpolated30cm.wav)
- AurInterpolated210cm: Contains the auralized SRIRs wav files with interpolation (resolution 30cm) and the combined wav file (AurInterpolated210cm.wav)
### Wav files
- 120 BPM - ROCK.wav: Music signal used during the measurements
- inv_chirp_signal.wav: Inverse chirp signal 
