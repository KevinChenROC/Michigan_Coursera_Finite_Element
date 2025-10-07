
# --- prerequisites (Ubuntu 22.04) ---
sudo apt update
sudo apt install -y \
  build-essential cmake gfortran pkg-config \
  wget curl tar xz-utils \
  libblas-dev liblapack-dev \
  libsuitesparse-dev libboost-all-dev \        # UMFPACK for SparseDirectUMFPACK
  zlib1g-dev libbz2-dev
#  ================================ download dealii without parallelization ================================


wget https://github.com/dealii/dealii/releases/download/v8.2.1/dealii-8.2.1.tar.gz
tar -xf dealii-8.2.1.tar.gz
cd dealii-8.2.1

mkdir build ; cd build

# Create the math compatibility shim directly from the shell
echo '#include <cmath>'                            >  dealii_math_compat.hpp
echo 'using std::isnan;'                          >> dealii_math_compat.hpp
echo 'using std::isinf;'                          >> dealii_math_compat.hpp
echo 'using std::isfinite;'                       >> dealii_math_compat.hpp

cmake .. \
  -DCMAKE_INSTALL_PREFIX=$HOME/local/dealii-8.2.1 \
  -DCMAKE_BUILD_TYPE=Release \
  -DBoost_NO_BOOST_CMAKE=ON \
  -DDEAL_II_WITH_MPI=OFF \
  -DDEAL_II_WITH_THREADS=OFF \
  -DDEAL_II_WITH_TBB=OFF \
  -DDEAL_II_WITH_PETSC=OFF \
  -DDEAL_II_WITH_TRILINOS=OFF \
  -DDEAL_II_WITH_HDF5=OFF \
  -DDEAL_II_WITH_NETCDF=OFF \
  -DDEAL_II_WITH_UMFPACK=ON \
  -DDEAL_II_WITH_LAPACK=ON \
  -DDEAL_II_WITH_METIS=OFF \
  -DDEAL_II_WITH_ARPACK=OFF \
  -DCMAKE_CXX_STANDARD=11 \
  -DCMAKE_CXX_EXTENSIONS=ON \
  -DCMAKE_CXX_FLAGS="-include $PWD/dealii_math_compat.hpp -include algorithm -include limits"

  # -DCMAKE_CXX_FLAGS="-std=gnu++11 -D_GLIBCXX_USE_C99_MATH_TR1 -include math.h"

make -j8 
make install & make test