MF = Makefile
FC = gfortran
LIBS = -lfftw3
FC += -Wall -fbounds-check -O2 #Flags for the FORTRAN compiler
#FC += -Wall -pedantic -Waliasing -Wcharacter-truncation \
         -Wsurprising -Wintrinsic-shadow -Warray-temporaries \
         -g -fbacktrace -fbounds-check      #Flags for the FORTRAN compiler (DEBUG)
FFLAGS = -cpp -DEPSR_STARTUP_READ
LFLAGS = $(FFLAGS)

EXE = wff_forces.x

SRC = setup_options.f90 module_f_env.f90 linked_cell.f90 neigh_list.f90 fourier.f90 ljones.f90 buckingham.f90 correct_tem.f90 ewald_cubic.f90 ewald_noncubic.f90 ewald_smeared.f90 epsr.f90 e3b.f90 intra_harmonic.f90 intra_morse.f90 intra_morse_harm.f90 io.f90 md_forces.f90 main.f90

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
