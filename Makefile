include conf/make.conf

KERNEL=mwd_kernel
VCFLAGS += -O3 -qopenmp
BUILD_DIR=build

KERNELS_SRC = $(wildcard src/kernels/*.c)
KERNELS_OBJ = $(patsubst src/kernels/%.c, ${BUILD_DIR}/%.o, $(KERNELS_SRC))


.PHONY: all, dp, split_stride, clean, debug, verbose, vtune

all: ${BUILD_DIR} ${BUILD_DIR}/${KERNEL}

dp: COMFLAGS+=-DDP=1
dp: VCFLAGS+=-DDP=1
dp: all
dp:
	mv ${BUILD_DIR} ${BUILD_DIR}_dp

# Experimental feature
split_stride: COMFLAGS+=-DUSE_SPLIT_STRIDE=1
split_stride: VCFLAGS+=-DUSE_SPLIT_STRIDE=1
split_stride: all

debug: COMFLAGS+= -fp-model strict
debug: all

verbose: COMFLAGS+= -qopt-report=5
verbose: all

vtune: COMFLAGS+= -DUSE_VUTNE
vtune: LIBS+= $(VTUNE_AMPLIFIER_XE_2016_DIR)/lib64/libittnotify.a
vtune: all

CFLAGS+= ${COMFLAGS}
FFLAGS+= ${COMFLAGS}




${BUILD_DIR}:
	mkdir -p ${BUILD_DIR}

${BUILD_DIR}/%.o: src/kernels/%.c
	${CC} ${CFLAGS} -Isrc -o $@ -c $<

${BUILD_DIR}/utils.o: src/utils.c
	${CC} ${CFLAGS} -o $@ -c $<

${BUILD_DIR}/verification.o: src/verification.c 
	${CC} ${VCFLAGS} -o $@ -c $<

${BUILD_DIR}/mpi_utils.o: src/mpi_utils.c 
	${CC} ${CFLAGS} -o $@ -c $<

${BUILD_DIR}/performance.o: src/performance.c 
	${CC} ${CFLAGS} -o $@ -c $< 

${BUILD_DIR}/driver.o: src/driver.c 
	${CC} ${CFLAGS} -o $@ -c $< 

${BUILD_DIR}/${KERNEL}: ${KERNELS_OBJ} ${BUILD_DIR}/driver.o ${BUILD_DIR}/utils.o ${BUILD_DIR}/performance.o  ${BUILD_DIR}/verification.o ${BUILD_DIR}/mpi_utils.o
	${CC} ${CFLAGS} -o $@ $^ ${LDFLAGS} ${LIBS}

clean:
	-rm -rf build build_dp error_snapshot reference_snapshot target_snapshot 


