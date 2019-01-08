IDIR =../../include
CC=mpic++
CFLAGS=-I$(IDIR)

ODIR = obj

_BASEDEPS = problem.h waveform.h WFR.h WFR_GS.h WFR_JAC.h WFR_NEW.h WFR_synch.h input_reader.h
BASEDEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

LOCALDEPS = problem_toy.h

_OBJ = waveform.o WFR_GS.o WFR_JAC.o WFR_synch.o WFR_NEW.o input_reader.o problem_toy.o main_toy.o
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
