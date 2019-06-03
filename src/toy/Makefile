IDIR =../../include
CC=mpic++
CFLAGS=-I$(IDIR)

ODIR = obj

_BASEDEPS = problem.h waveform.h waveform_locking.h WFR.h WFR_GS.h WFR_JAC.h WFR_NEW.h WFR_synch.h input_reader.h
BASEDEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

LOCALDEPS = problem_toy.h

_OBJ = waveform.o waveform_locking.o WFR.o WFR_GS.o WFR_JAC.o WFR_synch.o WFR_NEW.o input_reader.o problem_toy.o main_toy.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

_OBJ_LOG = waveform.o waveform_locking.o WFR.o WFR_GS.o WFR_JAC.o WFR_synch.o WFR_NEW.o input_reader.o problem_toy.o main_toy_log.o
OBJ_LOG = $(patsubst %,$(ODIR)/%,$(_OBJ_LOG))

all : TOY TOY_LOG

$(ODIR)/%.o: ../%.cpp $(BASEDEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/%.o: %.cpp $(LOCALDEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

TOY: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

TOY_LOG: $(OBJ_LOG)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ TOY TOY_LOG
