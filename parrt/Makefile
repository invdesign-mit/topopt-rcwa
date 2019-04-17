export PETSC_DIR=${HOME}/petsc-3.6.4
export PETSC_ARCH=arch-mumps-opt
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
CLEANFILES = *.o

NLOPT_INC = -I /usr/local/include
RCWA_INC = -I ../rcwa/

CC=mpicc
CXX=mpicxx
CFLAGS   += -O3 -Wall -march=native -fcx-limited-range -fno-exceptions
CXXFLAGS += -std=c++11 -O3 -Wall -march=native -fcx-limited-range -fno-exceptions
CPPFLAGS = -I. ${RCWA_INC} ${NLOPT_INC} ${PETSC_CC_INCLUDES}

LIBS4OBJS=filters.o input.o output.o optitemp.o optvfield.o scatter.o specs.o near2far.o multiobj.o

lib: parrt.a $(LIBS4OBJS)

parrt.a: $(LIBS4OBJS)
	cp ../rcwa/*.o ./ && $(AR) crvs $@ *.o && rm -rf *.o

filters.o: filters.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
input.o: input.c
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
output.o: output.c
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
optitemp.o: optitemp.c
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
optvfield.o: optvfield.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
scatter.o: scatter.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
specs.o: specs.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
near2far.o: near2far.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
multiobj.o: multiobj.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $< -o $@
