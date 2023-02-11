# Building TUV-x on modeling2

## Set up your library paths

- Add the following lines to your .bash_profile

...
PATH=$PATH:$HOME/.local/bin:$HOME/bin
PATH="/opt/local/bin:$PATH"
export PATH

export PKG_CONFIG_PATH="/opt/local/lib/pkgconfig:$PKG_CONFIG_PATH"
export LD_LIBRARY_PATH="/opt/local/lib64:/opt/local/lib:$LD_LIBRARY_PATH"
...

## Get the source code

- Copy the build script (`build_tuvx_modeling2_gnu.sh`) to modeling2

- Log in to modeling2 (BASH shell)

- Create a directory to build TUV-x in:

```
mkdir my-tuvx-build
```

- Create an environment variable named `TUVX_HOME` pointing the the absolute path of your build directory:

```
export TUVX_HOME=/path/to/my-tuvx-build
```

## Build TUV-x

Replace `/path/to/` with the path to the directory you copied the build script to, in the following:

```
cd $TUVX_HOME
. /path/to/build_tuvx_modeling2_gnu.sh
```

## Run TUV-x
- Run the tests

```
cd $TUVX_HOME/tuv-x/build
make test
```

- Run a configuration

```
cd $TUVX_HOME/tuv-x/build
./tuv-x examples/full_config.json
```

