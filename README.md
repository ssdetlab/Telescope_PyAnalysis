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

Run:
- `python3 multiproc_analyzer.py -conf conf/config_cosmics_sim_thr120e.txt`
- `python3 serial_analyzer.py -conf conf/config_source_data_dv9.txt`
