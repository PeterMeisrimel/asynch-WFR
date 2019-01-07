IDIR =../../include
CC=mpic++
CFLAGS=-I$(IDIR)

ODIR = obj

_BASEDEPS = problem.h waveform.h WFR.h WFR_GS.h WFR_JAC.h WFR_NEW.h WFR_synch.h
BASEDEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

LOCALDEPS = problem_bigger.h

_OBJ = waveform.o WFR_GS.o WFR_JAC.o WFR_synch.o WFR_NEW.o problem_bigger.o

_OBJ_BIGGER_GS = $(_OBJ) main_bigger_GS.o
OBJ_BIGGER_GS = $(patsubst %,$(ODIR)/%,$(_OBJ_BIGGER_GS))

_OBJ_BIGGER_JAC = $(_OBJ) main_bigger_JAC.o
OBJ_BIGGER_JAC = $(patsubst %,$(ODIR)/%,$(_OBJ_BIGGER_JAC))

_OBJ_BIGGER_NEW = $(_OBJ) main_bigger_NEW.o
OBJ_BIGGER_NEW = $(patsubst %,$(ODIR)/%,$(_OBJ_BIGGER_NEW))

all : GS JAC NEW

$(ODIR)/%.o: ../%.cpp $(BASEDEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/%.o: %.cpp $(LOCALDEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

GS: $(OBJ_BIGGER_GS)
	$(CC) -o $@ $^ $(CFLAGS)
JAC: $(OBJ_BIGGER_JAC)
	$(CC) -o $@ $^ $(CFLAGS)
NEW: $(OBJ_BIGGER_NEW)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ GS JAC NEW
