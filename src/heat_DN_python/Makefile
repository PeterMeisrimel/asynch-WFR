CC=mpic++
CFLAGS=-I../../include/ -Wall
CFLAGS_py1=-I../../include/ -I/usr/include/python3.6m/ -Wall
CFLAGS_py2=-lpython3.6m

ODIR = obj

_BASEDEPS = problem.h waveform.h waveform_locking.h WFR.h WFR_GS.h WFR_JAC.h WFR_NEW.h utils.h
BASEDEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

LOCALDEPS = problem_heat_python.h input_reader_heat.h

_OBJ = waveform.o waveform_locking.o WFR.o WFR_GS.o WFR_JAC.o WFR_NEW.o utils.o problem_heat_python.o main_heat_python.o input_reader_heat.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

all : HEAT_PYTHON

$(ODIR)/%.o: ../%.cpp $(BASEDEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

$(ODIR)/%.o: %.cpp $(LOCALDEPS)
	$(CC) $(CFLAGS) $(CFLAGS_py1) -c -o $@ $< $(CFLAGS_py2)

HEAT_PYTHON: $(OBJ)
	$(CC) $(CFLAGS) $(CFLAGS_py1) -o $@ $^ $(CFLAGS_py2)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ HEAT_PYTHON -r __pycache__
