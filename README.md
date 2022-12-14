<!-- PROJECT SHIELDS -->
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![Apache-2.0 License][license-shield]][license-url]

# PICO_Fortran
This is a Fortran implementation of the Potsdam Ice-shelf Cavity mOdel (PICO) from <a href="https://tc.copernicus.org/articles/12/1969/2018/">Reese et al. 2018</a>. 

## Dependencies
The implementation is dependent on a NetCDF library to be installed with a Fortran compiler, with the Fortran bindings for NetCDF. PICO_Fortran also depends on the ncio library that is resolved automatically with FPM (see below). 

## Test case
Currently implemented with a test case from the ISOMIP data. 

Original data found in: https://gmd.copernicus.org/articles/9/2471/2016/gmd-9-2471-2016.html and included from https://doi.org/10.5880/PIK.2016.002 into the data directory of PICO_Fortran.
Please refer to the original source when using. License of the dataset: CC BY-SA 4.0

## FPM
PICO_Fortran is packaged for the <a href="https://github.com/fortran-lang/fpm">Fortran Package Manager</a>, using a toml description file. Currently the build and app/run are functional.
