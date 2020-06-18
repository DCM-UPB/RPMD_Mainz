MF = Makefile
FC = gfortran
LIBS = -lfftw3
FC += -Wall -fbounds-check

FFLAGS = -O2
LFLAGS = $(FFLAGS)

EXE = wff_forces.x

SRC = module_f_env.f90 main.f90 setup_options.f90 md_forces.f90 md_tools.f90 random.f90 fourier.f90 io.f90 ljones.f90 buckingham.f90 ewald_cubic.f90 ewald_noncubic.f90 ewald_smeared.f90 epsr.f90 intra_morse_harm.f90 linked_cell.f90 module_thermostat.f90 md_pathi.f90 neigh_list.f90 correct_tem.f90 blas_lapack.f90

.SUFFIXES: .f90 .o

OBJ = $(SRC:.f90=.o)

.f90.o:
	$(FC) $(FFLAGS) -c $<

all: $(EXE)

$(EXE):	$(OBJ)
	$(FC) $(LFLAGS) -o $@ $(OBJ) $(LIBS)

$(OBJ): $(MF)

clean:
	rm -f $(OBJ) $(EXE)
	rm -f *.mod
