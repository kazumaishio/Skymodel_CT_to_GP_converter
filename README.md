# SkyModelConverter
This repository provides the conversion of a skymodel file, from a ctool format (xml) file to a gammapy format (yaml) file.

## Contents of this repository
- ```skymodelconverter```
  contains the conversion library.
- ```scripts```
  contains the example scripts to use the conversion library.

# Usage
After downloading this repository, edit and run one of the scripts on your ctool skymodel file(s). 

## Necessary packages
- astropy
- gammapy
- xmltodict
  ```
  conda create --name ctadatachallenge -c conda-forge gammapy xmltodict astropy notebook jupyter nb_conda_kernels 
  ```

## Necessary input files
- xml skymodel file (for ctool) 
- external files: the set of files referred from the skymodel file (mostly fits format files)
### NOTE: phasecurve files in gps need correction. Use phasecurvecorrector.py!

## Run a script
1. Determine your preferrable input directory and output directory.
2. Place the input files in the input directory. 
3. Open one of the example scripts in the "scripts" directory of this repository, and edit input and output paths, and run the script.

## Output files
The script writes the gammapy skymodel file(s) in the specified path.

# Change log (of significant revision)
- 20240314: bugfix for (1) wrong pivot energy in some spectral models, (2) problematic conversion from DiffuseMapCube
- 20240314: phasecurvecorrector
- 20230507: introducing README
- 20230507: unit change from power of day to power of sec in reading the phaseogram parameters
- 20221122: fist commit: test commit 
