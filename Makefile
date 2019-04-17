export PETSC_DIR=${HOME}/petsc-3.6.4
export PETSC_ARCH=arch-mumps-opt
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
CLEANFILES = *.o *_exec

CC=mpicc
CXX=mpicxx
PKGDIR=${HOME}/opt-rcwa7

NLOPT_INC = -I /usr/local/include
NLOPT_LIB = /usr/local/lib/libnlopt.a

RCWA_INC = -I ${PKGDIR}/rcwa
RCWA_LIB = ${PKGDIR}/rcwa/rcwa.a

PARRT_INC = -I ${PKGDIR}/parrt
PARRT_LIB = ${PKGDIR}/parrt/parrt.a

BLASLPK_INC = -I ${PKGDIR}/include
BLASLPK_LIB = -L ${PKGDIR}/lib -lopenblas

FFTW3_INC = -I ${PKGDIR}/include
FFTW3_LIB = -L ${PKGDIR}/lib -lfftw3

PTHREAD_INC = -DHAVE_UNISTD_H
PTHREAD_LIB = -lpthread

CFLAGS   += -O3 -Wall -march=native -fcx-limited-range -fno-exceptions
CXXFLAGS += -std=c++11 -O3 -Wall -march=native -fcx-limited-range -fno-exceptions
CPPFLAGS = -I. ${RCWA_INC} ${RCWA_INC}/RNP ${PARRT_INC} ${NLOPT_INC} ${PETSC_CC_INCLUDES}
INCFLAGS = -DHAVE_BLAS -DHAVE_LAPACK ${BLASLPK_INC} -DHAVE_FFTW3 $(FFTW3_INC) -DHAVE_LIBPTHREAD $(PTHREAD_INC)

LIBS=$(RCWA_LIB) $(PARRT_LIB) $(NLOPT_LIB) $(BLASLPK_LIB) $(FFTW3_LIB) $(PTHREAD_LIB) $(PETSC_LIB)

all: lens_exec

lens.o: lens.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCFLAGS) $< -o $@
lens_exec: lens.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@ $(LIBS)
