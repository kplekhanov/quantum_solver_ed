#--- up to you to change

# number of cores for customary operations
# does not change anything for arpack etc -- change it separately
OMP_NUM_THREADS ?= 16
LIB_ARMADILLO_DIR = /home/kplekhanov/opt/armadillo-9.900.1
BUILD_DIR = build

# should match the actual directory name though
SRC_DIR = src
LIB_DIR = lib

#--- do not change

CCCOM = g++
CCOPT = -O2 -Wno-write-strings -std=c++11 -Wno-deprecated -march=native -m64
# automatic creation of dependencies for .o files -- pretty cool
CCOPT += -MMD -MP
# imports the constant into cpp; check headers.hpp
CCDEFS = -DOMP_NUM_THREADS=$(OMP_NUM_THREADS)
CCFLAGS = $(CCOPT) $(CCDEFS)

LIBDIRS = -L$(LIB_ARMADILLO_DIR)
LIBFLAGS = -larmadillo -DARMA_DONT_USE_WRAPPER -lopenblas -llapack
LDFLAGS = $(LIBDIRS) $(LIBFLAGS)
