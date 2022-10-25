<h1 align="center">
<img src="docs/source/_static/logo.svg" width="300">
</h1><br>


TUV-x: A photolysis rate constant calculator.

[![License](https://img.shields.io/github/license/NCAR/tuv-x.svg)](https://github.com/NCAR/tuv-x/blob/main/LICENSE)
[![CI Status](https://github.com/NCAR/tuv-x/actions/workflows/test.yml/badge.svg)](https://github.com/NCAR/tuv-x/actions/workflows/test.yml)
[![codecov](https://codecov.io/gh/NCAR/tuv-x/branch/main/graph/badge.svg?token=H46AAEAQF9)](https://codecov.io/gh/NCAR/tuv-x)
[![DOI](https://zenodo.org/badge/396946468.svg)](https://zenodo.org/badge/latestdoi/396946468)

Copyright (C) 2020 National Center for Atmospheric Research

Please see the [TUV-x documentation](https://ncar.github.io/tuv-x/) for detailed
installation and usage instructions.

# Build and run (Docker version)

To build and run the stand-alone version of TUV-x, you must have [Docker Desktop](https://www.docker.com/get-started) installed and running. With Docker Desktop running, open a terminal window and run the following command to start the TUV-x container:

```
docker run -it ghcr.io/ncar/tuv-x:release bash
```

Inside the container, you can run the TUV-x tests from the `/build/` folder:

```
cd build/
make test
```

# Build and run (local build version)

To build and run TUV-x locally, you must have the following libraries available:

- [json-fortran](https://github.com/jacobwilliams/json-fortran)
- [NetCDF](https://www.unidata.ucar.edu/software/netcdf/) (both C and Fortran libraries)
- [nc4fortran](https://github.com/geospace-code/nc4fortran)

You must also have CMake installed on your machine. To build and run TUV-x locally,
open a terminal window, navigate to a folder where you would like the TUV-x files to exist,
and run the following commands:

```
git clone --recurse-submodules https://github.com/NCAR/tuv-x.git
cd tuv-x
mkdir build
cd build
ccmake ..
make
./tuv-x
```

# Citation

The following bibtex can be used to cite the work that originally developed
this tool.

A recommended citation is 

> Madronich, Sasha, and Siri Flocke (1999), The role of solar radiation in atmospheric chemistry, in Handbook of Environmental Chemistry, edited by P. Boule, pp. 1-26, Springer-Verlag, Heidelberg.

However, you are encouraged to use the format that best matches the style
you prefer.

```
@incollection{madronich_role_1999,
	address = {Berlin, Heidelberg},
	series = {The {Handbook} of {Environmental} {Chemistry}},
	title = {The {Role} of {Solar} {Radiation} in {Atmospheric} {Chemistry}},
	isbn = {978-3-540-69044-3},
	url = {https://doi.org/10.1007/978-3-540-69044-3_1},
	language = {en},
	booktitle = {Environmental {Photochemistry}},
	publisher = {Springer},
	author = {Madronich, Sasha and Flocke, Siri},
	editor = {Boule, Pierre},
	year = {1999},
	doi = {10.1007/978-3-540-69044-3_1},
	keywords = {Earth-Sun geometry., photolysis rate coefficients, radiative transfer, solar radiation, spectral actinic flux},
	pages = {1--26},
}
```

The TUV-x software can be cited with

```
@software{acom.software.tuvx,
  author       = {Matt Dawson and
                  Kyle Shores and
                  Stacy Walters},
  title        = {NCAR/tuv-x: Version 0.2.0},
  month        = sep,
  year         = 2022,
  publisher    = {Zenodo},
  version      = {v0.2.0},
  doi          = {10.5281/zenodo.7126040},
  url          = {https://doi.org/10.5281/zenodo.7126040}
}
```

and the specific verstion of TUV-x that you are using can be found by
clicking on the zenodo banner above. Choose the appropraite version there
and use the citation provided by Zenodo.