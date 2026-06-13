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
 FFLAGS = -qopenmp -O2 -g -mavx -I$(BUILDDIR) -I$(SRCDIR)

# place compiled module files into $(BUILDDIR)
ifeq ($(findstring ifort,$(FC)),ifort)
	FFLAGS += -module $(BUILDDIR)
else
	FFLAGS += -J$(BUILDDIR)
endif
 LIBS = -qopenmp -lmkl_intel_lp64  -lmkl_sequential -lmkl_core
# enable additional runtime checks for debugging (bounds, uninitialized)
DEBUG_FLAGS = -traceback -check bounds -check uninit
# To enable debug checks, append $(DEBUG_FLAGS) to FFLAGS when needed.

# gfortran linux
# FC = gfortran 
# FFLAGS = -fopenmp -Og -O2 -g
# LIBS = -L/usr/local/opt/lapack/lib/ -llapack -lblas


LFLAGS =

 OBJ = $(BUILDDIR)/parameters.o $(BUILDDIR)/general.o $(BUILDDIR)/iomodule.o $(BUILDDIR)/find_nnbrs.o $(BUILDDIR)/green_function.o $(BUILDDIR)/meminfo.o $(BUILDDIR)/mesh_mod.o $(BUILDDIR)/diag_mod.o $(BUILDDIR)/exchange_utils.o $(BUILDDIR)/streaming_mod.o $(BUILDDIR)/io_mod.o $(BUILDDIR)/hamiltonian_mod.o $(BUILDDIR)/setup_mod.o $(BUILDDIR)/input_mod.o $(BUILDDIR)/calc_mod.o

# Choose compile command: for ifort keep normal compile (it uses -module $(BUILDDIR)),
# for other compilers compile from inside $(BUILDDIR) so any default .mod/.o end up there.
ifeq ($(findstring ifort,$(FC)),ifort)
COMPILE = $(FC) -c $(FFLAGS) -o $@ $<
else
COMPILE = cd $(BUILDDIR) && $(FC) -c $(FFLAGS) -o $(notdir $@) ../$<
endif

all: $(PROG)
 
$(PROG):  $(BUILDDIR)/exchanges.o $(OBJ) | $(BINDIR)
	$(FC) $(LFLAGS) -o $@ $(BUILDDIR)/exchanges.o $(OBJ) $(LIBS)

clean:
	rm -rf $(BUILDDIR) $(BINDIR) *.mod


$(BUILDDIR)/%.o: $(SRCDIR)/%.f90 | $(BUILDDIR)
	$(COMPILE)

$(BUILDDIR)/general.o: $(SRCDIR)/general.f90 $(BUILDDIR)/parameters.o | $(BUILDDIR)
	$(COMPILE)

$(BUILDDIR)/iomodule.o: $(SRCDIR)/iomodule.f90 $(BUILDDIR)/general.o $(BUILDDIR)/parameters.o | $(BUILDDIR)
	$(COMPILE)

$(BUILDDIR)/green_function.o: $(SRCDIR)/green_function.f90 $(BUILDDIR)/general.o $(BUILDDIR)/parameters.o | $(BUILDDIR)
	$(COMPILE)

$(BUILDDIR)/diag_mod.o: $(SRCDIR)/diag_mod.f90 $(BUILDDIR)/green_function.o | $(BUILDDIR)
	$(COMPILE)

$(BUILDDIR)/find_nnbrs.o: $(SRCDIR)/find_nnbrs.f90 $(BUILDDIR)/general.o | $(BUILDDIR)
	$(COMPILE)

$(BUILDDIR)/exchange_utils.o: $(SRCDIR)/exchange_utils.f90 $(BUILDDIR)/diag_mod.o $(BUILDDIR)/general.o | $(BUILDDIR)
	$(COMPILE)

$(BUILDDIR)/exchanges.o: $(SRCDIR)/exchanges.f90 $(BUILDDIR)/parameters.o $(BUILDDIR)/general.o $(BUILDDIR)/iomodule.o $(BUILDDIR)/meminfo.o $(BUILDDIR)/mesh_mod.o $(BUILDDIR)/exchange_utils.o $(BUILDDIR)/streaming_mod.o $(BUILDDIR)/io_mod.o $(BUILDDIR)/diag_mod.o $(BUILDDIR)/hamiltonian_mod.o $(BUILDDIR)/setup_mod.o $(BUILDDIR)/input_mod.o $(BUILDDIR)/calc_mod.o | $(BUILDDIR)
	$(COMPILE)

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(BINDIR):
	mkdir -p $(BINDIR)
	
	
