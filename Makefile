os = $(shell uname -s)

INCFLAGS      = -I$(EICDIRECTORY)/include
INCFLAGS      += -I$(ROOTSYS)/include
#INCFLAGS      +=  -I$(FASTJETDIR)/include
INCFLAGS      += -I./src
INCFLAGS      += -I/direct/eic+u/eickolja/software/createJetTrees/StEpSimuJetMaker

LIBPATH       = -L$(EICDIRECTORY)/lib
ROOTLIBS      = $(shell root-config --libs)
LIBPATH       += $(ROOTLIBS) -L$(FASTJETDIR)/lib
LIBS          = -lfastjet -lfastjettools  -lConstituentSubtractor  -lRecursiveTools

# LIBPATH       += -L$(PYTHIA8DIR)/lib -L$(STARPICOPATH)
# LIBS          += -lpythia8  -llhapdfdummy -lTStarJetPico

LIBPATH       += -L./lib
# LIBS	      += -lMyJetlib

LIBPATH       += -L/afs/rhic.bnl.gov/eic/lib -L/eic/u/eickolja/software/createJetTrees/StEpSimuJetMaker
LIBS          += -leicsmear -lStEpSimuJetEvent


#INCFLAGS      += -I$(PYTHIA8DIR)/include -I$(PYTHIA8DIR)/include/Pythia8/ -I$(STARPICOPATH)

#INCFLAGS      += -I/afs/rhic.bnl.gov/eic/PACKAGES/EicRoot/eic-smear/include/
#INCFLAGS      += -I/afs/rhic.bnl.gov/eic/include/
#INCFLAGS      += -I/afs/rhic.bnl.gov/eic/env/new/include/


ifeq ($(os),Linux)
CXXFLAGS      = -std=c++11 -O -fPIC $(shell root-config --libs)
else
CXXFLAGS      = -O -fPIC -pipe -Wall -Wno-deprecated-writable-strings -Wno-unused-variable -Wno-unused-private-field -Wno-gnu-static-float-init
endif

CXXFLAGS      += -g

ifeq ($(os),Linux)
LDFLAGS       =
LDFLAGSS      = --shared
else
LDFLAGS       = -O -Xlinker -bind_at_load -flat_namespace
LDFLAGSS      = -flat_namespace -undefined suppress
LDFLAGSSS     = -bundle
endif


ifeq ($(os),Linux)
CXX          = g++ 
else
CXX          = clang
endif

# for cleanup
SDIR          = src
ODIR          = src/obj
BDIR          = bin


###############################################################################
################### Remake when these headers are touched #####################
###############################################################################
INCS = $(SDIR)/JetAnalyzer.hh $(SDIR)/AnalysisParameters.hh $(SDIR)/EicEnums.hh
INCS += $(SDIR)/EicAnalysis.hh

###############################################################################
# standard rules
$(ODIR)/%.o : $(SDIR)/%.cxx $(INCS)
	@echo 
	@echo COMPILING
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -c $< -o $@

$(BDIR)/%  : $(ODIR)/%.o 
	@echo 
	@echo LINKING
	$(CXX) $^ -o $@ $(LDFLAGS) $(LIBPATH) $(LIBS) 

###############################################################################

###############################################################################
############################# Main Targets ####################################
###############################################################################
all    : $(BDIR)/RunTest \
	 lib/libMyJetlib.so \
	 $(ODIR)/EicAnalysis.o
#	 doxy

## Individual rules
$(ODIR)/EicAnalysis.o 	 	: $(SDIR)/EicAnalysis.cxx $(INCS) $(SDIR)/EicAnalysis.hh

$(BDIR)/RunTest	:		$(ODIR)/RunTest.o 	$(ODIR)/EicAnalysis.o 	lib/libMyJetlib.so

lib/libMyJetlib.so	: $(ODIR)/JetAnalyzer.o
	@echo 
	@echo MAKING LIBRARY
	$(CXX) -shared $(LDFLAGS) $(LIBPATH) $(LIBS) $^ -o $@

lib/libMyJetlib.a	: $(ODIR)/JetAnalyzer.o
	@echo 
	@echo MAKING LIBRARY
	ar -rcs $@ $^

###############################################################################
##################################### MISC ####################################
###############################################################################


doxy: html/index.html

html/index.html : $(INCS) src/* Doxyfile
#	doxygen
	@echo 
	@echo Updating documentation
	( cat Doxyfile ; echo "QUIET=YES" ) | doxygen -

clean :
	@echo 
	@echo CLEANING
	rm -vf $(ODIR)/*.o
	rm -vf $(BDIR)/*
	rm -vf lib/*


.PHONY : clean doxy
