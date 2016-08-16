ifdef MIC
  include conf/make.conf.mic
  KERNEL=mwd_kernel.mic
  BUILD_POST=_mic
else
  include conf/make.conf
  KERNEL=mwd_kernel
  BUILD_POST=
endif


VCFLAGS += -g -O0

#DP=0
USE_SPLIT_STRIDE=0
INST=0
VERBOSE=0
DEBUG=0
VTUNE=0

ifeq ($(VERBOSE), 1)
	CFLAGS+= -openmp_report -vec-report3
	FFLAGS+= -openmp_report -vec-report3
endif

ifeq ($(VTUNE), 1)
  CFLAGS+= -DUSE_VUTNE
  FFLAGS+= -DUSE_VTUNE
  LIBS+= $(VTUNE_AMPLIFIER_XE_2016_DIR)/lib64/libittnotify.a
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
	BUILD_DIR=build_dp$(BUILD_POST)
else
	BUILD_DIR=build$(BUILD_POST)
endif

# Selecting time or space parallel multi-core wavefront
ifeq ($(USE_SPLIT_STRIDE), 1)
	CFLAGS+=-DUSE_SPLIT_STRIDE=1
	FFLAGS+=-DUSE_SPLIT_STRIDE=1
	VCFLAGS+=-DUSE_SPLIT_STRIDE=1
endif


