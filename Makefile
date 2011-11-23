MF=	Makefile
FC=	gfortran
LIBS =  -L/usr/local/cluster/fftw3//fftw3_gcc-64/lib/ -lfftw3 
#LIBS = -lfftw3 
FFLAGS=	-O3
LFLAGS=	$(FFLAGS)

EXE= qmd.x

SRC= \
	module_thermostat.f90 \
	module_barostat.f90 \
	module_gle.f90 \
        main.f90 \
	setup_positions.f90\
	setup_ice.f90 \
        setup_interface.f90 \
        setup_options.f90 \
	md_eq.f90 \
	md_melt.f90 \
	md_static.f90 \
	md_dynamics.f90 \
	md_pathi.f90 \
	md_forces.f90 \
	md_tools.f90 \
	md_rattle.f90 \
	md_pressure.f90 \
	evolve.f90 \
	evolve_basic.f90 \
	evolve_cl.f90 \
	evolve_pi.f90 \
	evolve_rigid.f90 \
	random.f90 \
	kinetic.f90 \
	fourier.f90 \
	io.f90 \
	ljones.f90 \
	buckingham.f90 \
	ewald_cubic.f90 \
	ewald_noncubic.f90 \
	ewald_smeared.f90 \
	intra_harmonic.f90 \
	intra_morse.f90 \
	intra_morse_harm.f90 \
	intra_spcfw.f90 \
	linked_cell.f90 \
	neigh_list.f90 \
	correct_tem.f90 \
	blas_lapack.f90 \

###############################

.SUFFIXES:
.SUFFIXES: .f90 .o

OBJ=	$(SRC:.f90=.o)

.f90.o:
	$(FC) $(FFLAGS) -c $<

all:	$(EXE)

$(EXE):	$(OBJ)
	$(FC) $(LFLAGS) -o $@ $(OBJ) $(LIBS)

$(OBJ):	$(MF)

tar:
	tar cvf $(EXE).tar $(MF) $(SRC)

clean:
	rm -f $(OBJ) $(EXE) core
