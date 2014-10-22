
#LOCATION=rmbp;

###############################################################
#
#  First, choose a C and C++ compilers, and set compiler flags.
#  This is done by setting the variables CXX, CC, CFLAGS, and
#  CXXFLAGS.
#
###############################################################
CC=gcc
# A C or C++ compiler, e. g., gcc, cc, xlc.
# There are a few components written in C which may be compiled
# under C++, but C compilers tend to generate slightly better code.
# The C++ and C compilers (if different) must generate
# compatible code.

CXX=g++
# A C++ compiler, e.g., g++, CC, xlC


CFLAGS=-O2
# Flags for the C compiler
# Some useful flags:
#   -O2   --  recommended level of optimization
#   -m64  --  needed to get 64-bit longs on some platforms
#   -g    --  debugging
#   -mcpu=v8  --  needed with gcc/g++ to get the full instruction set
#                 of a SparcStation
#   -qarch=ppc  -- needed with xlc/C to get the full instruction set
#                  of a PowerPC


CXXFLAGS=$(CFLAGS)
# Flags for the C++ compiler (usually the same as CFLAGS)
# Some useful flags:
#   -+  -- needed with xlC to make it treat .c files as C++ files


# Here are a few more variables for standard utilities.
# The defaults should almost always be OK.


AR=ar
# command to make a library

ARFLAGS=ruv
# arguments for AR

RANLIB=ranlib
# set to echo if you want to disable it completely

LDFLAGS=
# arguments for linker for C programs

LDFLAGS_CXX=$(LDFLAGS)
# libraries for linking C++ programs

LDLIBS=-lm
# libraries for linking C programs

LDLIBS_CXX=$(LDLIBS)
# libraries for linking C++ programs

CPPFLAGS=
# arguments for the C preprocessor

LIBTOOL=libtool
# libtool command

###############################################################
#
#  Second, specify your NTL (the Number Theory Library),
#  define the variables NTL_OPT_INCDIR, NTL_OPT_LIBDIR, 
#  NTL_OPT_LIB below.
#
###############################################################
ifdef LOCATION
	NTL_PREFIX=/Users/David/software/ntl-6.1.0
else	
	NTL_PREFIX=/home/david/ntl-6.1.0
endif

NTL_INCDIR=$(NTL_PREFIX)/include

NTL_LIBDIR=$(NTL_PREFIX)/lib

NTL_OPT_INCDIR=-I$(NTL_INCDIR) # NTLI
NTL_OPT_LIBDIR=-L$(NTL_LIBDIR) # NTLL
NTL_OPT_LIB=-lntl # NTL


###############################################################
#
#  Third, if you want to use GMP (the GNU Multi-Precision library),
#  define the variables GMP_OPT_INCDIR, GMP_OPT_LIBDIR, GMP_OPT_LIB below.
#  You also will have to set either NTL_GMP_LIP or NTL_GMP_HACK
#  in the config.h file.
#
#  Using GMP can lead to significant performance gains on some
#  platforms.  You can obtain GMP from http://www.swox.com/gmp.
#  Once you unpack it into a directory, just execute
#     ./configure; make
#  in that directory.
#
###############################################################
ifdef LOCATION
	GMP_PREFIX=/Users/David/software/gmp-6.0.0
else
	GMP_PREFIX=/home/david/gmp-6.0.0
endif

GMP_INCDIR=$(GMP_PREFIX)/include
# directory containing gmp.h if using GMP 

GMP_LIBDIR=$(GMP_PREFIX)/lib
# directory containing libgmp.a if using GMP 

GMP_OPT_INCDIR=-I$(GMP_INCDIR) # GMPI
GMP_OPT_LIBDIR=-L$(GMP_LIBDIR) # GMPL
GMP_OPT_LIB=-lgmp # GMP
# uncomment these if using GMP


###############################################################
#
#  Fourth, if you want to use gf2x (a library for fast
#  multiplication over GF(2)[X]), you need to
#  define the variables GF2X_OPT_INCDIR, GF2X_OPT_LIBDIR, GF2X_OPT_LIB below.
#  You also will have to set NTL_GF2X_LIB 
#  in the config.h file.
#
###############################################################
ifdef LOCATION
	GF2X_PREFIX=/Users/David/software/gf2x-1.1
else
	GF2X_PREFIX=/home/david/gf2x-1.1
endif

GF2X_INCDIR=$(GF2X_PREFIX)/include
# directory containing gf2x.h if using gf2x

GF2X_LIBDIR=$(GF2X_PREFIX)/lib
# directory containing libgf2x.a

GF2X_OPT_INCDIR=-I$(GF2X_INCDIR) # GF2X
GF2X_OPT_LIBDIR=-L$(GF2X_LIBDIR) # GF2X
GF2X_OPT_LIB=-lgf2x # GF2X
# uncomment these if using gf2x

#################################################################
#
#  Finally, you need to indicate where your source files and 
#  header files locate.
#
#################################################################


#################################################################
#
#  That's it!  You can ignore everything else in this file,
#  unless you have added new source files.
# 
#################################################################

# object files
O01=main.o
O02=$(O01) GoppaCode.o
O03=$(O02) McElieceCryptosystem.o 
O04=$(O03) NiederreiterCryptosystem.o
OBJ=$(O04)


# your own source files
S01=main.cpp
SRC=$(S01)

# your own header files
IN01=
INCL=$(IN01)


#################################################################
#
#  Rules for compiling your source files
#
#################################################################
NTL_OPT = $(NTL_OPT_LIB) $(NTL_OPT_INCDIR) $(NTL_OPT_LIBDIR) 
GMP_OPT = $(GMP_OPT_LIB) $(GMP_OPT_INCDIR) $(GMP_OPT_LIBDIR) 
GF2X_OPT = $(GF2X_OPT_LIB) $(GF2X_OPT_INCDIR) $(GF2X_OPT_LIBDIR) 

COMPILE = $(CC) $(CPPFLAGS) $(CFLAGS) $(NTL_OPT) $(GMP_OPT) $(GF2X_OPT) -g -c 
COMPILE_CXX = $(CXX) $(CPPFLAGS) $(CFLAGS) $(NTL_OPT) $(GMP_OPT) $(GF2X_OPT) -g -c 




MAIN.OUT: $(OBJ)
	  $(CXX)  $(OBJ) -g -o hello $(NTL_OPT) $(GMP_OPT) $(GF2X_OPT) 		

main.o: main.cpp
	$(COMPILE_CXX) main.cpp

GoppaCode.o: GoppaCode.cpp
	$(COMPILE_CXX) GoppaCode.cpp

McElieceCryptosystem.o: McElieceCryptosystem.cpp
	$(COMPILE_CXX) McElieceCryptosystem.cpp

NiederreiterCryptosystem.o: NiederreiterCryptosystem.cpp
	$(COMPILE_CXX) NiederreiterCryptosystem.cpp

#################################################################
#
#  Rules for cleaning up
#
#  make clobber removes *everything* created by make,
#    but it does not restore config.h to its default.
#
#  make clean tidies up a bit
#
#################################################################

clean:
	rm -f *.o 
	rm -rf small


