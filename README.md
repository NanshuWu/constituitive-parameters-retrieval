# Constituitive Parameters Retrieval

This is a set of scripts for retrieving constituitive parameters (e.g. relative permeability and permittivity), in the rectangular waveguide system, from measured S parameters. Three are three algorithms named NRW, NNI, and NIST, whose pros and cons are listed in comments to help choose a suitable algorithm for different samples. 

# How to use

1. Connect coax-to-waveguide convertors to corresponding coax ports.
2. Implement the TRL calibration on VNA. 
3. Prepare sample to a uniform rectangular plate featuring the size of rectangular waveguide WR-90 (0.9inches X 0.4inches, or 22.86mm X 10.16mm). 
4. Measure and extract S parameters with .S2P format. 
5. Choose a suitable script to get results.

# Tips

## Data preparation

If you want to test scripts with simulation data, say from CST, be sure to show original data with 'real/imag' format. Then, export S11 and S21 to two individual files. After choosing the script, replace the `filename_1` with your S11 file, and `filename_2` with S21 file. Finally set `is_simu=1`. Or just export files with touchstone format, and the process is as the same as the measurement scenario, as showing below.

If data is from measurement, be sure to save them with `.S2P` format to make sure that four S parameters are included. For experiment file, you only have to replace `filename_1` with the file path, and set `is_simu=0`.

## Choose algorithm according to material properties

| Material length/magnetic properties  | Methods | Speed | Accuracy |
|--------------------------------------|---------|-------|----------|
| Lossy solids + short + whatever magnetics | NRW     | Fast  | Medium   |
| Lossy solids + short + whatever magnetics | ITER     | Slow  | Medium   |
| low loss solids + whatever length + non-magnetics   | NIST    | Slow  | Good     |
| low loss solids + whatever length + non-magnetics   | NNI     | Fast  | Good     |
