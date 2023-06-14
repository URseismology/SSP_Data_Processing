# SSP_Data_Processing

Users of this package should cite:
> [![DOI](https://zenodo.org/badge/653740988.svg)](https://zenodo.org/badge/latestdoi/653740988)
> 
> Olugboji, T. M., Park, J., Karato, S.-I., & Shinohara, M. (2016). Nature of the seismic lithosphere-asthenosphere boundary within normal oceanic mantle from high-resolution receiver functions. Geochemistry, Geophysics, Geosystems, 17(4), 1265–1282.

BBOBS and borehole station data can be retrieved from the University of Rochester Seismology Lab server:
http://repovibranium.earth.rochester.edu/JapaneseProject.zip

### Notes on the usage of the data and processing scripts:
All relevant data and scripts associated with Olugboji et al. (2016) are included in the provided zip file `JanpaneseProject.zip`. The content of this package can be found in the `DataScripts/` directory. All processed SAC files can be found in the `ProcData/` directory. Major processing steps of the SAC files are performed using `DataQCwithOBSpy_c.ipynb`.

<img width="1045" alt="image" src="https://github.com/URseismology/SSP_Data_Processing/assets/66632382/1056e3c5-9257-407e-a798-dcc5d4908b55">
Figure. Location of the stations provided. 2 borehole stations (circles) and 17 ocean-bottom stations (triangles). The stations are color-coded based on the number of years deployed (red, 1 year deployments; white, 2 year deployments; black, 3 year deployments). Adapted from Figure 2 of Olugboji et al. (2016).
