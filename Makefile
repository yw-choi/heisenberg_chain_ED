.SUFFIXES: .f .F .o .a  .f90 .F90
.PHONY: fdf
COMPILER = ifort
LDFLAGS = -mkl

PROGRAM = main
OBJS = $(PROGRAM).o 

default: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(FC) -o $(PROGRAM).x $(LDFLAGS) $(OBJS) $(LIBS) 

clean:
	rm -f *.x *.o  *.a
	rm -f *.mod

%.o:%.mod

.F.o:
	$(FC) -c $(FFLAGS)  $(DEFS) $<
.f.o:
	$(FC) -c $(FFLAGS)   $<
.F90.o:
	$(FC) -c $(FFLAGS)  $(DEFS) $<
.f90.o:
	$(FC) -c $(FFLAGS)   $<
