FC=      gfortran                           #which FORTRAN compiler to use
CC=      gcc                                #which C comiler to use
CPFLAGS= -DEPSR_STARTUP_READ #-DPARALLEL_BINDING  #Flags for the C preprocessor

CCFLAGS= -O2                                #Flags for the C compiler
FCFLAGS= -cpp -O2 -Wall -pedantic -Waliasing -Wcharacter-truncation \
         -Wsurprising -Wintrinsic-shadow -Warray-temporaries \
         -funroll-loops                     #Flags for the FORTRAN compiler
LDFLAGS=                                    #Linker Flags
EXE=     qmd.x                              #Name of the executable

LIBS =   -lfftw3 

INCLUDE= 

