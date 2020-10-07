# number of cores for customary operations
# does not change anything for arpack etc -- change it separately
OMP_NUM_THREADS ?= 1

# !!! THE MOST IMPROTANT THING YOU MAY NEED TO CHANGE !!!
# armadillo library (has to have libarmadillo.so; by default the main armadillo rep)
LIB_ARMADILLO_DIR = /home/kplekhanov/opt/armadillo-9.900.1

# where to build the boject files; automatically created if necessary
BUILD_DIR = build

# reps for source files and headers
# should match the actual directory names
SRC_DIR = src
LIB_DIR = lib

# gcc compiler and options -- usual stuff
CCCOM = g++
CCOPT = -O2 -Wno-write-strings -std=c++11 -Wno-deprecated -march=native -m64
# automatic creation of dependencies for .o files -- pretty cool
CCOPT += -MMD -MP

# imports the constant into cpp -- check headers.hpp
CCDEFS = -DOMP_NUM_THREADS=$(OMP_NUM_THREADS)
CCFLAGS = $(CCOPT) $(CCDEFS)

# doesn't use armadillo vanilla tools by default; uses openblas and lapack instead
LIBDIRS = -L$(LIB_ARMADILLO_DIR)
LIBFLAGS = -larmadillo -DARMA_DONT_USE_WRAPPER -lopenblas -llapack
LDFLAGS = $(LIBDIRS) $(LIBFLAGS)
