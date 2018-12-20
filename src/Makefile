IDIR =../include
CC=mpic++
CFLAGS=-I$(IDIR)

ODIR=obj

_DEPS = problem.h waveform.h WFR.h WFR_GS.h WFR_JAC.h WFR_NEW.h WFR_synch.h problem_toy.h problem_bigger.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS)) 

_OBJ = waveform.o WFR_GS.o WFR_JAC.o WFR_synch.o WFR_NEW.o problem_toy.o problem_bigger.o

_OBJ_TOY_GS = $(_OBJ) main_toy_GS.o
OBJ_TOY_GS = $(patsubst %,$(ODIR)/%,$(_OBJ_TOY_GS))

_OBJ_TOY_JAC = $(_OBJ) main_toy_JAC.o
OBJ_TOY_JAC = $(patsubst %,$(ODIR)/%,$(_OBJ_TOY_JAC))

_OBJ_TOY_NEW = $(_OBJ) main_toy_NEW.o
OBJ_TOY_NEW = $(patsubst %,$(ODIR)/%,$(_OBJ_TOY_NEW))

_OBJ_BIGGER_GS = $(_OBJ) main_bigger_GS.o
OBJ_BIGGER_GS = $(patsubst %,$(ODIR)/%,$(_OBJ_BIGGER_GS))

_OBJ_BIGGER_JAC = $(_OBJ) main_bigger_JAC.o
OBJ_BIGGER_JAC = $(patsubst %,$(ODIR)/%,$(_OBJ_BIGGER_JAC))

_OBJ_BIGGER_NEW = $(_OBJ) main_bigger_NEW.o
OBJ_BIGGER_NEW = $(patsubst %,$(ODIR)/%,$(_OBJ_BIGGER_NEW))

all : TOY_GS TOY_JAC TOY_NEW BIGGER_GS BIGGER_JAC BIGGER_NEW

toy: TOY_GS TOY_JAC TOY_NEW
bigger: BIGGER_GS BIGGER_JAC BIGGER_NEW

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

TOY_GS: $(OBJ_TOY_GS)
	$(CC) -o $@ $^ $(CFLAGS)
TOY_JAC: $(OBJ_TOY_JAC)
	$(CC) -o $@ $^ $(CFLAGS)
TOY_NEW: $(OBJ_TOY_NEW)
	$(CC) -o $@ $^ $(CFLAGS)

BIGGER_GS: $(OBJ_BIGGER_GS)
	$(CC) -o $@ $^ $(CFLAGS)
BIGGER_JAC: $(OBJ_BIGGER_JAC)
	$(CC) -o $@ $^ $(CFLAGS)
BIGGER_NEW: $(OBJ_BIGGER_NEW)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ TOY_* BIGGER_*
