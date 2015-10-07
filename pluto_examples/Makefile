SRC=3d7pt
DST=lbpar_$(SRC)
DST2=debug_$(SRC)

CC=icc
PLC=polycc

PLC_FLAGS = -q --tile --parallel --pet --partlbtile
PLC_DEBUG_FLAGS = 

CFLAGS = -O3 -xHost -ansi-alias -ipo -openmp
DEBUG_FLAGS = 
LIKWID_FLAGS = -DLIKWID -DLIKWID_PERFMON  ${LIKWID_INCLUDE} ${LIKWID_LIB} -llikwid
#LDFLAGS = -lm

.PHONY: all, clean 

all: $(DST)

$(SRC).lbpar.c:  $(SRC).c
	$(PLC) $< $(PLC_FLAGS) -o $@

$(DST): $(SRC).lbpar.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS) $(LIKWID_FLAGS)


$(SRC).debug.c:  $(SRC).c
	$(PLC) $< $(PLC_DEBUG_FLAGS) -o $@

$(DST2): $(SRC).c
	$(CC) $(DEBUG_FLAGS) $< -o $@ $(LDFLAGS) $(LIKWID_FLAGS)


clean:
	rm -f *.lbpar.c  *.cloog lbpar_* debug_* .dumpstmttxt

