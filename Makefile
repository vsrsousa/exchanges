 PROG = bin/exchanges.x 
 SRCDIR = src
 BUILDDIR = build
 BINDIR = bin

#gfortan Mac Os
# FC = gfortran 
# FFLAGS = -fopenmp -Og -O2 -g
# LIBS = -fopenmp -framework Accelerate #-L/usr/local/opt/lapack/lib/ -llapack -lblas

# Intel fortran linux
 FC = ifort 
 FFLAGS = -qopenmp -O2 -g -mavx -I$(SRCDIR)
 LIBS = -qopenmp -lmkl_intel_lp64  -lmkl_sequential -lmkl_core
#uncomment for debug:
# FFLAGS = -qopenmp -O0 -g -mavx -traceback -check

# gfortran linux
# FC = gfortran 
# FFLAGS = -fopenmp -Og -O2 -g
# LIBS = -L/usr/local/opt/lapack/lib/ -llapack -lblas


LFLAGS =

 OBJ = $(BUILDDIR)/parameters.o $(BUILDDIR)/general.o $(BUILDDIR)/iomodule.o $(BUILDDIR)/find_nnbrs.o $(BUILDDIR)/green_function.o

all: $(PROG)
 
$(PROG):  $(BUILDDIR)/exchanges.o $(OBJ) | $(BINDIR)
	$(FC) $(LFLAGS) -o $@ $(BUILDDIR)/exchanges.o $(OBJ) $(LIBS)

clean:
	rm -rf $(BUILDDIR) $(BINDIR) *.mod


$(BUILDDIR)/%.o: $(SRCDIR)/%.f90 | $(BUILDDIR)
	$(FC) -c $(FFLAGS) -o $@ $<

$(BUILDDIR)/exchanges.o: $(SRCDIR)/exchanges.f90 | $(BUILDDIR)
	$(FC) -c $(FFLAGS) -o $@ $<

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(BINDIR):
	mkdir -p $(BINDIR)
	
	
