FC=      gfortran                           #which FORTRAN compiler to use
CC=      gcc                                #which C comiler to use
CPFLAGS= -DCP2K_BINDING -DEPSR_STARTUP_READ \
         #-DPARALLEL_BINDING                #Flags for the C preprocessor

CCFLAGS= -O2                                #Flags for the C compiler
FCFLAGS= -cpp -O2 -Wall -pedantic -Waliasing -Wcharacter-truncation \
         -Wsurprising -Wintrinsic-shadow -Warray-temporaries \
         -funroll-loops                     #Flags for the FORTRAN compiler
LDFLAGS=                                    #Linker Flags
EXE=     qmd.x                              #Name of the executable

CP2K_DIR=/home/kristof/software/cp2k/branch-4.1#root directrory of cp2k install
LIBINT_LIB=/usr/lib64/libint                #Path to libint library
LIBINT_INC=/usr/include/libint              #Path to libint include files

LIBS =   -L$(CP2K_DIR)/lib/Linux-x86-64-gfortran/sopt \
         -lcp2kbase -lcp2kfft -ldbcsr \
         -lcp2kstart\
         -L/usr/lib -llapack -lblas -lstdc++ -lfftw3 \
         -L$(LIBINT_LIB) -lderiv -lint
#        -lcp2kmain -lcp2kmc \ #-lcp2k -lcp2kma

INCLUDE= -I$(CP2K_DIR)/obj/Linux-x86-64-gfortran/sopt/ \
         -I$(LIBINT_INC)/libint

