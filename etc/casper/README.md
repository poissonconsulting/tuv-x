# Building TUV-x on CASPER

## Get the source code

- Copy the build script (`build_tuvx_casper_gnu.sh`) to GLADE.

- Log in to CASPER and start an interactive session (BASH shell)

- Create a directory to build TUV-x in:

```
mkdir my-tuvx-build
```

- Create an environment variable named `TUVX_HOME` pointing to the absolute path of your build directory:

```
export TUVX_HOME=/path/to/my-tuvx-build
```

## Build TUV-x

Replace `/path/to/` with the path to the directory you copied the build script to, in the following:

```
cd $TUVX_HOME
. /path/to/build_tuvx_casper_gnu.sh
```

## Run TUV-x
- Whenever you go to run TUV-x after it has been built, make sure you have the correct environment modules loaded:

```
module purge
module load gnu/11.2.0
module load ncarenv/1.3
module load ncarcompilers/0.5.0
module load cmake/3.22.0
module load netcdf/4.8.1
```

- Run a test. The tests use the python numpy package, and we access this through the conda environment module on CASPER:

```
module load conda/latest
conda activate npl
cd $TUVX_HOME/tuv-x/build
make test
```

- Run a configuration (this does not require any Python packages)

```
cd $TUVX_HOME/tuv-x/build
./tuv-x examples/full_config.json
```

