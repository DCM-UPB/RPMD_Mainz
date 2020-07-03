FC=      gfortran                           #which FORTRAN compiler to use
CC=      gcc                                #which C comiler to use
CPFLAGS= -DEPSR_STARTUP_READ #-DPARALLEL_BINDING -DENABLE_IPI_DRIVER  #Flags for the C preprocessor

CCFLAGS=                                    #Flags for the C compiler
FCFLAGS= -cpp -Wall -pedantic -Waliasing -Wcharacter-truncation \
         -Wsurprising -Wintrinsic-shadow -Warray-temporaries \
         -ggdb3 -fbacktrace -fcheck=all \
         -ftrap-fpe=invalid,overflow,zero   #Flags for the FORTRAN compiler
LDFLAGS=                                    #Linker Flags
EXE=     qmd-dbg.x                          #Name of the executable

LIBS =   -lfftw3 

INCLUDE= 

