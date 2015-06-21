SRC=3d7pt
DST=lbpar_$(SRC)

CC=icc
PLC=polycc

PLC_FLAGS = -q --tile --parallel --pet --partlbtile

CFLAGS = -O3 -xHost -ansi-alias -ipo -openmp
LIKWID_FLAGS = -DLIKWID -DLIKWID_PERFMON  ${LIKWID_INCLUDE} ${LIKWID_LIB} -llikwid
#LDFLAGS = -lm

.PHONY: all, clean 

all: $(DST)

$(SRC).lbpar.c:  $(SRC).c
	$(PLC) $(SRC).c $(PLC_FLAGS) -o $@

$(DST): $(SRC).lbpar.c
	$(CC) $(CFLAGS) $< -o $@  $(LDFLAGS) $(LIKWID_FLAGS)

clean:
	rm -f *.lbpar.c  *.cloog lbpar_*

