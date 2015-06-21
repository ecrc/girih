SRC=3d7pt
GEN_SRC=$(SRC)
DST=lbpar_$(GEN_SRC)

CC=icc
PLC=polycc

PLC_FLAGS = -q --tile --parallel --pet --partlbtile

CFLAGS = -O3 -xHost -ansi-alias -ipo -openmp
LIKWID_FLAGS = -DLIKWID -DLIKWID_PERFMON  ${LIKWID_INCLUDE} ${LIKWID_LIB} -llikwid
#LDFLAGS = -lm

.PHONY: all, clean 

all: $(DST)

$(GEN_SRC).lbpar.c:  $(SRC).c
	$(PLC) $< $(PLC_FLAGS) -o $@

$(DST): $(GEN_SRC).lbpar.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS) $(LIKWID_FLAGS)

clean:
	rm -f *.lbpar.c  *.cloog lbpar_*

