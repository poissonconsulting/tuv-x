# Downloads and builds TUV-x and its dependencies on modeling2 using GNU compilers
#
# The TUVX_HOME environment variable must be set to the directory to build TUV-x
# in prior to calling this script

if [[ -z "${TUVX_HOME}" ]]; then
  echo "You must set the TUVX_HOME environment variable to the directory where TUV-x should be build."
  return
fi

if [[ ! -d "${TUVX_HOME}" ]]; then
  echo "TUVX_HOME must point to an existing directory"
  return
fi

echo "Building TUV-x"

export PATH="/opt/local/bin:${PATH}"
export LD_LIBRARY_PATH="/opt/local/lib64:/opt/local/lib:/usr/bin:/usr/lib:usr/lib64:usr/local/bin:usr/local/lib:usr/local/lib64"

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
export JSON_FORTRAN_HOME=$INSTALL_ROOT/jsonfortran-gnu-8.2.1
cd $JSON_FORTRAN_ROOT
sed -i 's/\-C $<CONFIG>//' CMakeLists.txt
mkdir -p build
cd build
cmake3 -D CMAKE_Fortran_COMPILER=/opt/local/bin/gfortran \
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
cmake3 -D CMAKE_Fortran_COMPILER=/opt/local/bin/gfortran \
       -D CMAKE_BUILD_TYPE=release \
       -D NETCDF_INCLUDE_DIR=/opt/local/include \
       -D NETCDF_C_LIB=/opt/local/lib/libnetcdf.so \
       -D NETCDF_FORTRAN_LIB=/opt/local/lib/libnetcdff.so \
       -D ENABLE_COVERAGE=OFF \
       ..
make
