DEPENDENCY_STORAGE_DIR=${HOME}/opt-rcwa7/depstore
DEPENDENCY_INSTALL_DIR=${HOME}/opt-rcwa7

rm -rf ${DEPENDENCY_STORAGE_DIR} ./lib ./include ./bin
mkdir ${DEPENDENCY_STORAGE_DIR}

cd ${DEPENDENCY_STORAGE_DIR}
wget http://github.com/xianyi/OpenBLAS/archive/v0.2.20.tar.gz

tar xzvf v0.2.20.tar.gz
cd OpenBLAS-0.2.20
make
make PREFIX=${DEPENDENCY_INSTALL_DIR} install

cd ${DEPENDENCY_STORAGE_DIR}
wget http://fftw.org/fftw-3.3.8.tar.gz

tar xzvf fftw-3.3.8.tar.gz
cd fftw-3.3.8
./configure --prefix=${DEPENDENCY_INSTALL_DIR}
make
make install

suitesparse=0
if [ $suitesparse -eq 1 ]
then
    cd ${DEPENDENCY_STORAGE_DIR}
    wget http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-5.2.0.tar.gz

    tar xzvf SuiteSparse-5.2.0.tar.gz
    cd SuiteSparse
    make library
    make install INSTALL=${DEPENDENCY_INSTALL_DIR}
fi
