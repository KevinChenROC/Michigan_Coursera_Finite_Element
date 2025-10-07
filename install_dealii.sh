
# --- prerequisites (Ubuntu 22.04) ---
sudo apt update
sudo apt install -y \
  build-essential cmake gfortran pkg-config \
  wget curl tar xz-utils \
  libblas-dev liblapack-dev libhdf5-dev \
  libsuitesparse-dev libboost-all-dev \
  zlib1g-dev libbz2-dev
#  ================================ download dealii without parallelization ================================


wget https://github.com/dealii/dealii/releases/download/v8.2.1/dealii-8.2.1.tar.gz
tar -xf dealii-8.2.1.tar.gz
cd dealii-8.2.1

mkdir build ; cd build

# one header to rule them all
cat > dealii_compat.hpp <<'EOF'
#include <cmath>
#include <algorithm>
#include <limits>
using std::isnan;
using std::isinf;
using std::isfinite;
EOF

cmake .. \
  -DCMAKE_INSTALL_PREFIX=$HOME/local/dealii-8.2.1 \
  -DCMAKE_BUILD_TYPE=Release \
  -DBoost_NO_BOOST_CMAKE=ON \
  -DDEAL_II_WITH_MPI=OFF \
  -DDEAL_II_WITH_THREADS=OFF \
  -DDEAL_II_WITH_TBB=OFF \
  -DDEAL_II_WITH_PETSC=OFF \
  -DDEAL_II_WITH_TRILINOS=OFF \
  -DDEAL_II_WITH_HDF5=ON \
  -DDEAL_II_WITH_NETCDF=OFF \
  -DDEAL_II_WITH_UMFPACK=ON \
  -DDEAL_II_WITH_LAPACK=ON \
  -DDEAL_II_WITH_METIS=OFF \
  -DDEAL_II_WITH_ARPACK=OFF \
  -DCMAKE_CXX_STANDARD=11 \
  -DCMAKE_CXX_EXTENSIONS=ON \
  -DCMAKE_CXX_FLAGS="-include $PWD/dealii_compat.hpp" \
  -DHDF5_INCLUDE_DIRS=/usr/include/hdf5/serial \
  -DHDF5_LIBRARIES="/usr/lib/x86_64-linux-gnu/libhdf5_serial.so;/usr/lib/x86_64-linux-gnu/libhdf5_serial_hl.so" \
  -DHDF5_IS_PARALLEL=OFF

make -j8
make install & make test
echo -e '\nexport DEAL_II_DIR=$HOME/local/dealii-8.2.1\nexport CMAKE_PREFIX_PATH=$DEAL_II_DIR:$CMAKE_PREFIX_PATH' >> ~/.zshrc