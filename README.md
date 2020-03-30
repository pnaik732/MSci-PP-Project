# MSci-PP-Project
This reposition contains the code, installation and project lab notes for the final year project.

## Get started

### Getting the source
Clone with git
```
git clone http://github.com/tz16559/MSci-PP-Project/ --recursive
```
### Pre-requirements

Get AmpGen (version 1.2) installed

Get python3

Get ROOT, uproot, scipy, numpy, pandas, matplotlib, subprocess installed

## How to use the code

### MC data
To generate MC data, go to the AmpGen/build/bin/
```
cd B0_event_generator
```
You can specify the range of the factor to multiply to the amplitude and phase of the resonance. Change the maximum seed number if necessary (value of b). More amplitudes and phases used would require a larger number of seeds. Start generating events by typing:
```
python run_generator.py
```

### Analyse MC data
Copy the generated event directories to B0_event_analysis, open run_program.py and make sure the part in between ='s is the same as in run_generator.py. start analysing by typing:
```
python run_program.py
```

### Analyse LHCb data
In the MSci-PP-Project directory,
```
mkdir data_new
```
put the data in data_new, by default, the names should be
```
Data_sig_tos_weights-Run1.pkl
Data_sig_tos_weights-Run2.pkl
Data_sig_tis_weights-Run1.pkl
Data_sig_tis_weights-Run2.pkl
```
Go to B0_event_analysis_LHCb, start the analysis by typing
```
python run_program.py
```
Note the root file generated in the results directory is used for amplitude analysis by copy them to B0_amplitude_plot

### Amplitude analysis  - LHCb data
Copy B0_amplitude_plot to AmpGen/build/bin/

Go to B0_amplitude_plot, run the SignalOnltFitter and plot the graph by typing
```
python run_compare.py
```
You can change the data type by changing the names and save it to a differnet name (line 23 and 37) inside the run_compare.py



