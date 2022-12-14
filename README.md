<!-- README updated loosely based on: https://github.com/othneildrew/Best-README-Template -->

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

## Test case
Currently implemented with a test case from the ISOMIP data. 

Original data found in: https://gmd.copernicus.org/articles/9/2471/2016/gmd-9-2471-2016.html and included from https://doi.org/10.5880/PIK.2016.002 into the data directory of PICO_Fortran.
Please refer to the original source when using. License of the dataset: CC BY-SA 4.0

## FPM
PICO_Fortran is packaged for the <a href="https://github.com/fortran-lang/fpm">Fortran Package Manager</a>, using a toml description file. Currently the build and app/run are functional.

Since there are external dependencies, you need to specify where they lies if they are not in your standard path. For example, setting the following environnement variables for the NetCDF libraries:

```
export NETCDF_CFLAGS=""
export NETCDF_FFLAGS="-I/usr/include"
export FPM_FFLAGS="${NETCDF_CFLAGS} ${NETCDF_FFLAGS}"

export NETCDF_CLIBS="" 
export NETCDF_FLIBS="-L/usr/lib/x86_64-linux-gnu"
export FPM_LDFLAGS="${NETCDF_CLIBS} ${NETCDF_FLIBS}"
```

Then you can run the build with fpm:
```
fpm build
```

With successful built, you can run the standard example:
```
fpm run
```


<!-- Other stuff taken to get the shields correctly -->

[contributors-shield]: https://img.shields.io/github/contributors/dmr-dj/PICO_Fortran
[contributors-url]: https://github.com/dmr-dj/PICO_Fortran/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/dmr-dj/PICO_Fortran
[forks-url]: https://github.com/dmr-dj/PICO_Fortran/network/members
[stars-shield]: https://img.shields.io/github/stars/dmr-dj/PICO_Fortran
[stars-url]: https://github.com/dmr-dj/PICO_Fortran/stargazers
[issues-shield]: https://img.shields.io/github/issues/dmr-dj/PICO_Fortran
[issues-url]: https://github.com/dmr-dj/PICO_Fortran/issues
[license-shield]: https://img.shields.io/github/license/dmr-dj/PICO_Fortran
[license-url]: https://github.com/dmr-dj/PICO_Fortran/blob/master/LICENSE.txt
