ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
ROOTGLIBS  = $(shell root-config --glibs) 
CXXFLAGS  += $(ROOTCFLAGS)
GLIBS      = $(ROOTGLIBS)
GXX	   = /usr/bin/g++ -Wall -O3

ROOTFLAGS   = $(ROOTCFLAGS) $(ROOTLIBS) $(ROOTGLIBS) 
P5640FLAGS  = -L${P5640LIB}/lib -lP5640  -I${P5640LIB}
GSLFLAGS    = -I${EBROOTGSL}/include/gsl  -I/usr/include/gsl -lgsl -lgslcblas

all: Scattering.cpp

Scattering: Scattering.cpp ${P5640LIB}/lib/libP5640.a
	$(GXX) -o Scattering Scattering.cpp $(ROOTFLAGS) $(P5640FLAGS) $(GSLFLAGS)

clean: 
	rm -f Scattering *~
