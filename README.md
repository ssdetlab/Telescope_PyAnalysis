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

Run analysis:
- `python3 multiproc_analyzer.py -conf conf/config_cosmics_sim_thr120e.txt`
- `python3 serial_analyzer.py -conf conf/config_source_data_dv9.txt`

Run alignment with cosmics:
- step 1: `python3 multiproc_analyzer.py -conf conf/config_cosmics_sim_thr120e.txt`
- step 2: `python3 alignment_analyzer.py -conf conf/config_cosmics_sim_thr120e.txt`
- step 3: put the resulting misalinment values in the corresponding config file
- setp 4: repeat step 1 and then compare the residuals and the chi2 histograms
