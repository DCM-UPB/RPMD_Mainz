MF=	Makefile src/globals.inc src/2DPMF.inc
FC=	gfortran 
CC=     gcc
FCFLAGS= -ccp -DCP2K_BINDING #-DPARALLEL_BINDING
LIBS =  -L/root/cp2k/lib/Linux-x86-64-gfortran/sopt \
        -lcp2k_lib -lcp2k_base_lib -lcp2k_fft_lib -lcp2k_ma_lib \
        -lcp2k_dbcsr_lib \
        -L/usr/lib -llapack -lblas -lstdc++ -lfftw3\
        -L/usr/lib64/libint/ -lderiv -lint
INCLUDE= -I/root/obj/Linux-x86-64-gfortran/sopt/ \
	 -I/usr/include/libint
#-----------------------------------------------------
#check whether a cp2k build is actually available
#-----------------------------------------------------
ifeq ($(wildcard /home/cp2k/trunk/cp2k/lib/Linux-x86-64-gfortran/sopt),)
FC=	gfortran 
FCFLAGS= -cpp
INCLUDE= ""
LIBS = -lfftw3
NOCP2K="WARNING: CP2K_BINDINGS NOT COMPILED"
endif

# warnings that could result in wrong code
FCFLAGS+=	-Wall -pedantic -Waliasing -Wcharacter-truncation -Wconversion -Wsurprising -Wintrinsic-shadow

# speed warnings
FCFLAGS+=	-Warray-temporaries

FCFLAGS+= 	-fbounds-check -g 

# Read EPSR only on startup
FCFLAGS+=	-DEPSR_STARTUP_READ

FFLAGS= -O2
LFLAGS=	$(FFLAGS)
FCFLAGS+= $(FFLAGS)

EXE:= qmd.x
WD=$(shell pwd)
SRCDIR:=src
OBJDIR:=build
SRCPATH=$(addprefix $(WD)/,$(SRCDIR))
SRCFILES=\
	module_thermostat.f90 \
	module_barostat.f90 \
	module_gle.f90 \
	module_f_env.f90 \
	intmod.f90 \
	trajavmod.f90 \
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
	io.f90 \
	ljones.f90 \
	epsr.f90 \
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
        e3b.f90 \
        main.f90 
SRC=$(addprefix $(SRCPATH)/,$(SRCFILES) sockets.c)
###############################

.SUFFIXES:
.SUFFIXES: .f90 .o
.PHONY: directories
OBJPATH=$(addprefix $(WD)/,$(OBJDIR))
OBJ=	$(addprefix $(OBJPATH)/,$(SRCFILES:.f90=.o) sockets.o)

all:	directories $(EXE)

directories:
	mkdir -p build
$(EXE):	$(OBJ)
	echo $(NOCP2K)
	$(FC) $(LFLAGS) $^ $(LIBS) -o $@

#$(OBJ):	$(MF) 
$(OBJPATH)/%.o: $(SRCPATH)/%.f90 $(MF)
	$(FC) -c $(FCFLAGS) $(INCDIR) $< -J $(OBJPATH) -o $@	
$(OBJPATH)/%.o: $(SRCPATH)/%.c $(MF)
	$(CC) -O2 -c $< -o $@

zip:
	zip $(EXE).zip $(MF) $(SRC)

tar:
	tar cvf $(EXE).tar $(MF) $(SRC)

clean:
	rm -rf $(EXE) $(OBJPATH)
