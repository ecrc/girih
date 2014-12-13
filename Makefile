include conf/make.conf
include conf/make.inc

KERNEL=mwd_kernel

KERNELS_SRC = $(wildcard src/kernels/*.c)
KERNELS_OBJ = $(patsubst src/kernels/%.c, ${BUILD_DIR}/%.o, $(KERNELS_SRC))


.PHONY: all, clean 

all: ${BUILD_DIR} ${BUILD_DIR}/${KERNEL}

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


