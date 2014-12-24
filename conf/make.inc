VCFLAGS += -g -O0

DP=0
USE_SPLIT_STRIDE=0
INST=0
VERBOSE=0
DEBUG=0

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

