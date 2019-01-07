IDIR =../../include
CC=mpic++
CFLAGS=-I$(IDIR)

ODIR = obj

_BASEDEPS = problem.h waveform.h WFR.h WFR_GS.h WFR_JAC.h WFR_NEW.h WFR_synch.h
BASEDEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

LOCALDEPS = problem_toy.h

_OBJ = waveform.o WFR_GS.o WFR_JAC.o WFR_synch.o WFR_NEW.o problem_toy.o

_OBJ_TOY_GS = $(_OBJ) main_toy_GS.o
OBJ_TOY_GS = $(patsubst %,$(ODIR)/%,$(_OBJ_TOY_GS))

_OBJ_TOY_JAC = $(_OBJ) main_toy_JAC.o
OBJ_TOY_JAC = $(patsubst %,$(ODIR)/%,$(_OBJ_TOY_JAC))

_OBJ_TOY_NEW = $(_OBJ) main_toy_NEW.o
OBJ_TOY_NEW = $(patsubst %,$(ODIR)/%,$(_OBJ_TOY_NEW))

all : GS JAC NEW

$(ODIR)/%.o: ../%.cpp $(BASEDEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/%.o: %.cpp $(LOCALDEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

GS: $(OBJ_TOY_GS)
	$(CC) -o $@ $^ $(CFLAGS)
JAC: $(OBJ_TOY_JAC)
	$(CC) -o $@ $^ $(CFLAGS)
NEW: $(OBJ_TOY_NEW)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ GS JAC NEW
