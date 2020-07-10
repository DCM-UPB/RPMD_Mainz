Special version of RPMD code, adapted for variable projection force matching. All required source files are identical to the normal code (just symlinks to the main src dir), except main.f90 and md_forces.f90.
Due to the code being cut to the minimum required for force calculation, a lot of input variables are useless and left out of the input list. See the force matching example for details.

