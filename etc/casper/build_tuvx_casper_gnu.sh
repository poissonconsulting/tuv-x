# Downloads and builds TUV-x and its dependencies on CASPER using GNU compilers
#
# The TUVX_HOME environment variable must be set to the directory to build TUV-x
# in prior to calling this script


module purge
module load gnu/11.2.0
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
curl -LO https://github.com/geospace-code/nc4fortran/archive/refs/tags/v1.4.2.tar.gz
git clone --recurse-submodules https://github.com/NCAR/tuv-x.git

# extract
cd ${TUVX_HOME}
tar -zxf 8.2.1.tar.gz
tar -zxvf v1.4.2.tar.gz

INSTALL_ROOT=$TUVX_HOME/install
mkdir -p $INSTALL_ROOT

# json-fortran
JSON_FORTRAN_ROOT=$TUVX_HOME/json-fortran-8.2.1
export JSON_FORTRAN_HOME=$INSTALL_ROOT/jsonfortran-gnu-8.2.1
cd $JSON_FORTRAN_ROOT
sed -i 's/\-C $<CONFIG>//' CMakeLists.txt
mkdir -p build
cd build
cmake -D CMAKE_Fortran_COMPILER=gfortran \
      -D SKIP_DOC_GEN:BOOL=TRUE \
      -D CMAKE_INSTALL_PREFIX=$INSTALL_ROOT \
      ..
make install
mkdir -p $JSON_FORTRAN_HOME/lib/shared
mv $JSON_FORTRAN_HOME/lib/*.so* $JSON_FORTRAN_HOME/lib/shared

# nc4fortran
NC4_ROOT=$TUVX_HOME/nc4fortran-1.4.2
cd $NC4_ROOT
mkdir build
cd build
cmake ..
make

# TUV-x
TUVX_ROOT=$TUVX_HOME/tuv-x
cd $TUVX_ROOT
mkdir -p build
cd build
cmake -D CMAKE_Fortran_COMPILER=gfortran \
      -D CMAKE_BUILD_TYPE=release \
      -D NETCDF_INCLUDE_DIR=$NCAR_INC_NETCDF \
      -D NETCDF_C_LIB=$NCAR_LDFLAGS_NETCDF/libnetcdf.so \
      -D NETCDF_FORTRAN_LIB=$NCAR_LDFLAGS_NETCDF/libnetcdff.so \
      -D NC4F_INCLUDE_DIR=$NC4_ROOT/build/include \
      -D NC4F_LIB=$NC4_ROOT/build/libnc4fortran.a \
      -D ENABLE_COVERAGE=OFF \
      ..
make
