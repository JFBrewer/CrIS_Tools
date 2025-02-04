# CrIS_Tools
This is a repository where I store the custom code I've written for handling the CrIS outputs. It contains several different subpackages

# CHEEREIO Pre-processing Tools
This contains the custom CrIS Filtering tool I built originally for use pre-processing the CrIS-Ethane files for CHEEREIO. It is written to be generally applicable to the various different ROCRv2 species (ethane, HCN, methanol, ethene, ethyne, etc.).

The code written contains a custom implementation of the xESMF regridder based upon that present in the GEOS-Chem GCpy python package. Where GCpy's regridder requires the use of a 'conservative' regridding algorithm, mine allows for customization, something necessary when dealing with the large amount of missing data.

It also contains a 'MapFilter' class, which implements several methods important for the preprocessing of CrIS data for use in CHEEREIO. The basic concept of this class is that it can maintain the complex 4-D structure of the raw CrIS data (lat, lon, ANN, P90) through the process of landmasking and time-averaging the CrIS retrievals.

Finally, it contains a python notebook giving my implementation of this filtering for my own work. I include this to clarify how I intend the tools to be used, for the reference of those for whom its use might be relevant.
