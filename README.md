# Constituitive Parameters Retrieval

This is a set of scripts for retrieving constituitive parameters (e.g. relative permeability and permittivity), in the rectangular waveguide system, from measured S parameters. Three are three algorithms named NRW, NNI, and NIST, whose pros and cons are listed in comments to help choose a suitable algorithm for different samples. 

# How to use

1. Connect coax-to-waveguide convertors to corresponding coax ports.
2. Implement the TRL calibration on VNA. 
3. Prepare sample to a uniform rectangular plate featuring the size of rectangular waveguide WR-90 (0.9inches X 0.4inches, or 22.86mm X 10.16mm). 
4. Measure and extract S parameters with .S2P format. 
5. Choose a suitable script to get results.
