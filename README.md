# Telescope_PyAnalysis

Prerequisits:
- TelescopeEvent lib (to read the EUDAQ data files directly)
- Corryvreckan libs (to read the EUDAQ-->Corryvreckan data files or AllPix2-->Corryvreckan simulation files)

Setup:
- Setup ROOT
- `export LD_LIBRARY_PATH=/path/to/TelescopeEvent/libs:$LD_LIBRARY_PATH`
- `export LD_LIBRARY_PATH=/path/to/Corryvreckan/lib:$LD_LIBRARY_PATH`
- put data files somewhere with enough space...
- change config file as needed

Run noise scan:
- change `doNoiseScan` in the config to 1
- `python3 serial_analyzer.py -conf conf/config_file_name.txt`
- change `doNoiseScan` in the config back to 0

Run analysis:
- run noise scan (see above)
- `python3 serial_analyzer.py -conf conf/config_file_name.txt`
- `python3 multiproc_analyzer.py -conf conf/config_file_name.txt`
- to see event displays (fits...):
  - change `doplot` in the config to 1
  - run with `serial_analyzer.py` as above
  - kill the process after as many fits as desired
  - change `doplot` in the config back to 0

Run alignment with cosmics:
- run noise scan (see above)
- step 1: `python3 multiproc_analyzer.py -conf conf/config_file_name.txt`
- step 2: `python3 alignment_analyzer.py -conf conf/config_file_name.txt -det ALPIDE_1`
- step 3: look at the resulting histograms to see if the minimum is contained in all 3 ranges and change as needed
- step 4: put the resulting misalinment values in the corresponding config file for the specific detector
- setp 5: repeat step 1 and then compare the residuals and the chi2 histograms
- step 6: repeat steps 1-5 for another detector, until all "middle" detectors are aligned
