#--------------------------------------------------------------------------

# Compile extra debugging code (slight performance impact)
export WITH_DEBUG = 1

# Compile debug version
export DEBUG = 1

ARCH         := $(shell root-config --arch)
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLDFLAGS  := $(shell root-config --ldflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

ifdef DEBUG
  OPT         = -g
  OPT2        = -g
else
  OPT         = -O
  OPT2        = -O
endif

ifeq ($(ARCH),linux)
# Linux with egcs, gcc 2.9x, gcc 3.x
CXX           = g++
CXXFLAGS      = $(OPT2) -Wall -fPIC -Woverloaded-virtual
LD            = g++
LDFLAGS       = $(OPT2)
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxx8664gcc)
# AMD Opteron and Intel EM64T (64 bit mode) Linux with gcc 3.x
CXX           = g++
CXXFLAGS      = $(OPT2) -Wall -fPIC -Woverloaded-virtual
LD            = g++
LDFLAGS       = $(OPT2)
SOFLAGS       = -shared
endif

ifeq ($(ARCH),macosx64)
# MacOS X >= 10.4 with gcc 64 bit mode (GNU gcc 4.*)
# Only specific option (-m64) comes from root-config
MACOSX_MINOR := $(shell sw_vers | sed -n 's/ProductVersion://p' | cut -d . -f 2)
MACOSXTARGET := MACOSX_DEPLOYMENT_TARGET=10.$(MACOSX_MINOR)
CXX           = g++
CXXFLAGS      = $(OPT2) -pipe -Wall -W -Woverloaded-virtual -Wno-deprecated
LD            = $(MACOSXTARGET) g++
LDFLAGS       = $(OPT2)
# The SOFLAGS will be used to create the .dylib,
# the .so will be created separately
ifeq ($(subst $(MACOSX_MINOR),,1234),1234)
DllSuf        = so
else
DllSuf        = dylib
endif
SOFLAGS       = -dynamiclib -single_module -undefined dynamic_lookup
endif

ifeq ($(CXX),)
$(error $(ARCH) invalid architecture)
endif

INCLUDES      = $(ROOTCFLAGS)

CXXFLAGS     += $(INCLUDES)
LDFLAGS      += $(ROOTLDFLAGS)
LIBS         += $(ROOTLIBS) $(SYSLIBS)  
GLIBS        += $(ROOTGLIBS) $(SYSLIBS)

MAKEDEPEND    = gcc

ifdef WITH_DEBUG
CXXFLAGS     += -DWITH_DEBUG
endif

#-------------------------------------------------------------------------

SRC           = TLab.C TSim.C TTheory.C TRunInfo.C Messages.C Compton.C

OBJ           = $(SRC:.C=.o)
HDR           = $(SRC:.C=.h)
DEP           = $(SRC:.C=.d)
OBJS          = $(OBJ) simLabDict.o

LIBSIMLAB     = libSIMLAB.so
PROGRAMS      = simLab

all:            $(PROGRAMS)

$(LIBSIMLAB):	$(OBJS)
		$(LD) $(LDFLAGS) $(SOFLAGS) -o $@ $^
		@echo "$@ done"

simLab:		simLab.o $(LIBSIMLAB)
		$(LD) $(LDFLAGS) simLab.o -L$(CURDIR) -lSIMLAB $(GLIBS) -lMinuit -o $@

clean:
		rm -f *.o *.d *.so $(PROGRAMS)

realclean:	clean
		rm -f *.d *~ core
		rm -f simLabDict.*

simLabDict.C: 	$(HDR) SIMLAB_LinkDef.h
		@echo "Generating dictionary simLabDict..."
		$(ROOTSYS)/bin/rootcint -f $@ -c $(INCLUDES) $(HDR) \
		 SIMLAB_LinkDef.h

#-------------------------------------------------------------------------

.SUFFIXES:
.SUFFIXES: .c .cc .cpp .C .o .d

%.o:	%.C
	$(CXX) $(CXXFLAGS) -o $@ -c $<

%.d:	%.C
	@echo Creating dependencies for $<
	@$(SHELL) -ec '$(MAKEDEPEND) -MM $(INCLUDES) -c $< \
		| sed '\''s%^.*\.o%$*\.o%g'\'' \
		| sed '\''s%\($*\)\.o[ :]*%\1.o $@ : %g'\'' > $@; \
		[ -s $@ ] || rm -f $@'

###

-include $(DEP)