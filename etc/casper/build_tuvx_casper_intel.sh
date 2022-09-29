# Downloads and builds TUV-x and its dependencies on CASPER using GNU compilers
#
# The TUVX_HOME environment variable must be set to the directory to build TUV-x
# in prior to calling this script


module purge
module load intel/2022.1
module load openmpi/3.1.6
module load ncarenv/1.3
module load ncarcompilers/0.5.0
module load cmake/3.22.0
module load netcdf/4.8.1

if [[ -z "${TUVX_HOME}" ]]; then
  echo "You must set the TUVX_HOME environment variable to the directory where TUV-x should be build."
  return
fi

if [[ ! -d "${TUVX_HOME}" ]]; then
  echo "TUVX_HOME must point to an existing directory"
  return
fi

echo "Building TUV-x"

# get the source code
cd ${TUVX_HOME}
curl -LO https://github.com/jacobwilliams/json-fortran/archive/8.2.1.tar.gz
git clone --recurse-submodules https://github.com/NCAR/tuv-x.git

# extract
cd ${TUVX_HOME}
tar -zxf 8.2.1.tar.gz

INSTALL_ROOT=$TUVX_HOME/install
mkdir -p $INSTALL_ROOT

# json-fortran
JSON_FORTRAN_ROOT=$TUVX_HOME/json-fortran-8.2.1
export JSON_FORTRAN_HOME=$INSTALL_ROOT/jsonfortran-intel-8.2.1
cd $JSON_FORTRAN_ROOT
sed -i 's/\-C $<CONFIG>//' CMakeLists.txt
mkdir -p build
cd build
cmake -D CMAKE_Fortran_COMPILER=ifort \
      -D SKIP_DOC_GEN:BOOL=TRUE \
      -D CMAKE_INSTALL_PREFIX=$INSTALL_ROOT \
      ..
make install
mkdir -p $JSON_FORTRAN_HOME/lib/shared
mv $JSON_FORTRAN_HOME/lib/*.so* $JSON_FORTRAN_HOME/lib/shared

# TUV-x
TUVX_ROOT=$TUVX_HOME/tuv-x
cd $TUVX_ROOT
git checkout release
git submodule update
mkdir -p build
cd build
cmake -D CMAKE_Fortran_COMPILER=ifort \
      -D CMAKE_BUILD_TYPE=release \
      -D NETCDF_INCLUDE_DIR=$NCAR_INC_NETCDF \
      -D NETCDF_C_LIB=$NCAR_LDFLAGS_NETCDF/libnetcdf.so \
      -D NETCDF_FORTRAN_LIB=$NCAR_LDFLAGS_NETCDF/libnetcdff.so \
      -D ENABLE_MPI=OFF \
      -D ENABLE_COVERAGE=OFF \
      ..
make