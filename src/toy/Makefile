IDIR =../../include
CC=mpic++
CFLAGS=-I$(IDIR) -Wall

ODIR = obj

_BASEDEPS = problem.h waveform.h waveform_locking.h WFR.h WFR_GS.h WFR_JAC.h WFR_NEW.h utils.h
BASEDEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

LOCALDEPS = problem_toy.h

_OBJ = waveform.o waveform_locking.o WFR.o WFR_GS.o WFR_JAC.o WFR_NEW.o utils.o problem_toy.o main_toy.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

all : TOY

$(ODIR)/%.o: ../%.cpp $(BASEDEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/%.o: %.cpp $(LOCALDEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

TOY: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ TOY
