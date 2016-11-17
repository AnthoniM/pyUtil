# pyUtil
A bunch of tools to analyse data in python

`corefunctions.py` contains usefull functions

`coreclasses.py` contains : 
  - `SweepData` : a class to load data file (uses `pyHegel` convension to read the shape of files).
  - `Plot` : a class that implement some basic plot configurations.
  - `ExploreParameters` : a class to manipulate the arguments of a list of function and see the outcome. 
  
`dynamicplot.py` : a class to visualise data along different dimension by interactivelly selecting curves associated
                   with a given set of experimental parameters (very usefull when data is more than 3D). 

`nlfit.py` : a class to fit make multidimentional fit with a set of functions that can share some parameters. 
             (taken from Jean Oliver Simoneau : https://github.com/JeanOlivier)
