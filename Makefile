MF=	Makefile
FC=	gfortran -cpp -DCP2K_BINDING #-DPARALLEL_BINDING
#LIBS = -lfftw3
LIBS = -L/home/cp2k/trunk/cp2k/lib/Linux-x86-64-gfortran/sopt -lcp2k_lib -lcp2k_base_lib -lcp2k_fft_lib -lcp2k_ma_lib -lcp2k_dbcsr_lib \
				-L/usr/lib -llapack -lblas -lstdc++ -lfftw3\
				/home/grizzly/Programme_fuer_CP2K/libint-1.1.4/lib/libderiv.a \
        /home/grizzly/Programme_fuer_CP2K/libint-1.1.4/lib/libint.a \

ifeq ($(wildcard /home/cp2k/trunk/cp2k/lib/Linux-x86-64-gfortran/sopt),)
FC=	gfortran -cpp
LIBS = -lfftw3
NOCP2K="WARNING: CP2K_BINDINGS NOT COMPILED"
endif
# warnings that could result in wrong code
FC+=	-Wall -pedantic -Waliasing -Wcharacter-truncation -Wconversion -Wsurprising -Wintrinsic-shadow

# speed warnings
FC+=	-Warray-temporaries

FFLAGS=	-O3 
LFLAGS=	$(FFLAGS)

EXE= qmd.x

SRC= \
	module_thermostat.f90 \
	module_barostat.f90 \
	module_gle.f90 \
	module_f_env.f90 \
        main.f90 \
	RPMDDFT.f90\
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
	evolve.f90 \
	evolve_basic.f90 \
	evolve_cl.f90 \
	evolve_pi.f90 \
	evolve_rigid.f90 \
	evolve_cl_RPMDDFT.f90 \
	evolve_pi_RPMDDFT.f90 \
	evolve_pi_rc_RPMDDFT.f90 \
	evolve_cl_pi_RPMDDFT.f90 \
	random.f90 \
	kinetic.f90 \
	fourier.f90 \
	instvacint.f90 \
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
	md_pressure.f90 \
###############################

.SUFFIXES:
.SUFFIXES: .f90 .o

OBJ=	$(SRC:.f90=.o)

.f90.o:
	echo $(NOCP2K)
	$(FC) $(FFLAGS) -c $<

all:	$(EXE)

$(EXE):	$(OBJ)
	echo $(NOCP2K)
	$(FC) $(LFLAGS) -o $@ $(OBJ) $(LIBS)

$(OBJ):	$(MF)

tar:
	tar cvf $(EXE).tar $(MF) $(SRC)

clean:
	rm -f $(OBJ) $(EXE) core
