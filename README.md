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
- For gammalib 2.0
  - python 3.11 (>=3.12 does not work for gammalib2.0 because of the use of the obsolete 'imp' package)
  - gammapy (latest = 2.0)
  - gammalib 2.0 ( follow the compilation installation, including cfitsio)
  - xmltodict (installation possible via pip, not via conda)

### you can create the dedicated conda environment for avoiding package confliction as follows;
  ```
  conda create --name ctadatachallenge 
  conda activate ctadatachallenge 
  conda install -c conda-forge gammapy
  (installation process for cfitsio)
  (installation process for gammalib)
  pip install xmltodict
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
- 20260105: adjustment to the gammapy 2.0, including the update of NodeFunction conversion.
  - gammapy version: 2.0.1
  - gammalib version: 2.0.0
  - python version 3.11.14
  - cfitsio version 4.5.0
- 20240314: bugfix for (1) wrong pivot energy in some spectral models, (2) problematic conversion from DiffuseMapCube
- 20240314: phasecurvecorrector
- 20230507: introducing README
- 20230507: unit change from power of day to power of sec in reading the phaseogram parameters
- 20221122: fist commit: test commit 
