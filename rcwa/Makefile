CC=gcc
CXX=g++
PKGDIR=${HOME}/opt-rcwa7

BLASLPK_INC = -I ${PKGDIR}/include
BLASLPK_LIB = -L ${PKGDIR}/lib -lopenblas

FFTW3_INC = -I ${PKGDIR}/include
FFTW3_LIB = -L ${PKGDIR}/lib -lfftw3

PTHREAD_INC = -DHAVE_UNISTD_H
PTHREAD_LIB = -lpthread

CFLAGS   += -O3 -Wall -march=native -fcx-limited-range -fno-exceptions
CXXFLAGS += -std=c++11 -O3 -Wall -march=native -fcx-limited-range -fno-exceptions
CPPFLAGS = -I. -IRNP
INCFLAGS = -DHAVE_BLAS -DHAVE_LAPACK ${BLASLPK_INC} -DHAVE_FFTW3 $(FFTW3_INC) -DHAVE_LIBPTHREAD $(PTHREAD_INC)

LIBS=$(BLASLPK_LIB) $(FFTW3_LIB) $(PTHREAD_LIB)

LIBS4OBJS=numalloc.o fft_iface.o sort.o gsel.o rcwa.o wrap.o objfuncs.o Eigensystems.o 

lib: rcwa.a $(LIBS4OBJS)

rcwa.a: $(LIBS4OBJS)
	$(AR) crvs $@ $(LIBS4OBJS) 

numalloc.o: numalloc.c
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
sort.o: sort.c
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
gsel.o: gsel.c
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
fft_iface.o: fft_iface.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
rcwa.o: rcwa.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
wrap.o: wrap.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
objfuncs.o: objfuncs.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
Eigensystems.o: RNP/Eigensystems.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@

clean:
	rm -rf *.a *.o
