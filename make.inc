#MPI:
CC=mpicc
FC=mpif90
 
# with TAU
INST_CC=tau_cc.sh
INST_FC=tau_f90.sh 

  
VCFLAGS = -g -O0

# Intel compiler flags
c5hdn_CFLAGS = -g -O3 -xAVX -fno-alias -openmp

# Intel compiler flags at Stampede
stampede_CFLAGS = -g -O3 -xAVX -xhost -fno-alias -openmp

# Intel compiler flags at Emmy
emmy_CFLAGS = -g -O3 -xAVX -fno-alias -openmp# -mcmodel=large

# Intel compiler flags at Harpertown
harpertown1_CFLAGS = -g -O3 -openmp -funroll-loops -xHOST -no-prec-div

# GCC compiler flags
gcc_CFLAGS = -g -O3


DP=0
USE_SPLIT_STRIDE=0
INST=0
VERBOSE=0
DEBUG=0


# Detect machine
MACHINE=$(shell uname -n)
ifeq ($(findstring local,$(MACHINE)),local)
	CFLAGS=$(gcc_CFLAGS)
	ARCH=mac_laptop
else
ifeq ($(findstring c5hdn,$(MACHINE)),c5hdn)
	CFLAGS=$(c5hdn_CFLAGS)
	ARCH=Hatem_machine
else
ifeq ($(findstring stampede,$(MACHINE)),stampede)
	CFLAGS=$(stampede_CFLAGS)
	ARCH=Stampede
else
ifeq ($(findstring emmy,$(MACHINE)),emmy)
	CFLAGS=$(emmy_CFLAGS) ${LIKWID_INC} -DLIKWID -DLIKWID_PERFMON
	VCFLAGS += ${LIKWID_INC} -DLIKWID -DLIKWID_PERFMON
        LDFLAGS = ${LIKWID_LIB} -llikwid
	ARCH=Emmy
else
ifeq ($(findstring ivyep1,$(MACHINE)),ivyep1)
	CFLAGS=$(emmy_CFLAGS) ${LIKWID_INC} -DLIKWID_PERFMON
	VCFLAGS += ${LIKWID_INC} -DLIKWID_PERFMON
        LDFLAGS = ${LIKWID_LIB} -llikwid
	ARCH=Ivyep1
else
ifeq ($(findstring sandyep4,$(MACHINE)),sandyep4)
	CFLAGS=$(emmy_CFLAGS) ${LIKWID_INC} -DLIKWID_PERFMON
	VCFLAGS += ${LIKWID_INC} -DLIKWID_PERFMON
        LDFLAGS = ${LIKWID_LIB} -llikwid
	ARCH=sandyep4
else
ifeq ($(findstring harpertown1,$(MACHINE)),harpertown1)
	CFLAGS=$(harpertown1_CFLAGS)
	ARCH=Harpertown1
else
	ARCH=Unknown
endif
endif
endif
endif
endif
endif
endif

ifeq ($(INST), 1)
	override CC = $(INST_CC)
	override FC = $(INST_FC)
endif

ifeq ($(VERBOSE), 1)
	CFLAGS+= -openmp_report -vec-report3
	FFLAGS+= -openmp_report -vec-report3
endif

ifeq ($(DEBUG), 1)
	CFLAGS+= -fp-model strict
	FFLAGS+= -fp-model strict
endif

# Selecting single or double precision build
ifeq ($(DP), 1)
	CFLAGS+=-DDP=1
	FFLAGS+=-DDP=1
	VCFLAGS+=-DDP=1
	BUILD_DIR=build_dp
else
	BUILD_DIR=build
endif

# Selecting time or space parallel multi-core wavefront
ifeq ($(USE_SPLIT_STRIDE), 1)
	CFLAGS+=-DUSE_SPLIT_STRIDE=1
	FFLAGS+=-DUSE_SPLIT_STRIDE=1
	VCFLAGS+=-DUSE_SPLIT_STRIDE=1
endif

