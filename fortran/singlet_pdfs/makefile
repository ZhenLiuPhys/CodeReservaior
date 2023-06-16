# Makefile for ppqq.f
# Platform: macOS X
# Compiler: f2c fortran
FC = gfortran

FFLAGD = -g -extend_source -check bounds -check overflow -u -warn unused
FFLAGS = -w -O -g

#LIBFLAGS = -L/sw/lib -lf2c
#LIBFLAGS = -L/sw/lib -lpdflib -lmathlib -lpacklib
#LIBFLAGS = -L/usr/local/cern/2001/lib -lpdflib -lmathlib -lpacklib 

OBJqq = Cteq6Demo.o hanlib.o Cteq6Pdf-2010.o

EXEqq = Cteq6Demo

# ---------------------------------------------------------
# rules:

$(EXEqq): $(OBJqq) 
	$(FC) -o $(EXEqq) $(FFLAGS) $(OBJqq) $(LIBFLAGS)

.f.o:
	$(FC)  $(FFLAGS)  -c  $*.f
