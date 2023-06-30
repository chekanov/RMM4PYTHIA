# S.Chekanov

# define here PYTHIA and HEPMC directories
# PYTHIADIR=/share/sl6/pythia8

ifndef PYTHIADIR
$(error PYTHIADIR env variable is not set. Run setup.sh first)
endif

ifndef ROOTSYS
$(error ROOTSYS env variable is not set. Run setup.sh first)
endif

ifndef FASTJET
$(error FASTJET env variable is not set. Run setup.sh first)
endif

include ${ROOTSYS}/etc/Makefile.arch

# Root variables
ROOTCFLAGS    = $(shell root-config --nonew --cflags)
ROOTLIBS      = $(shell root-config --nonew --libs)
ROOTGTTLIBS   = $(shell root-config --nonew --glibs)
CXXFLAGS     += $(ROOTCFLAGS)

LIBDIRARCH=lib/
OutPutOpt     = -o  
LIBS         += -L$(PYTHIADIR)/$(LIBDIRARCH) -lpythia8
LIBS         += -L$(FASTJET)/lib -lfastjet

SOURCE_FILES1 := $(shell ls -1 main.cxx)
SOURCE_FILES1 += $(shell ls -1 src/*.cxx)


INCLUDE1=-I./src
INCLUDE2=-I. -I./inc 
INCLUDE3=-I$(FASTJET)/include
INCLUDE4=-I$(HEPMC)/include
INCLUDE5=-I$(PYTHIADIR)/include

# Rui's version
OPT=-O0 -Wall -Wextra -fsanitize=address -lasan

# Sergei version
OPT=-O0 -Wall -Wextra


# build object files 
objects1       = $(patsubst %.cxx,%.o,$(SOURCE_FILES1))


%.o: %.cxx
	$(CXX) $(OPT) $(CXXFLAGS) $(INCLUDE1) $(INCLUDE2) $(INCLUDE3) $(INCLUDE5) -o $@ -c $<

Tasks:     clean main.exe

mydict: inc/LParticle.h
	rm -f src/CParticle_dict*
	@rm -f inc/CParticle_dict*
	@echo "Generating dictionary for CParticle"
	@rootcint src/CParticle_dict.cxx -c inc/CParticle.h
	@rm -f src/LParticle_dict*
	@rm -f inc/LParticle_dict*
	@echo "Generating dictionary for LParticle"
	@rootcint src/LParticle_dict.cxx -c inc/LParticle.h


LIBOBJS = $(patsubst %.cxx,%.o,$(SOURCE_FILES))

main.exe: $(objects1)
	$(LD) $(OPT) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@

clean:
	        @rm -f *.o *~ main.exe src/*.o ;  echo "Clear.." 
