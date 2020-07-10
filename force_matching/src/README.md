Special version of RPMD code, adapted for variable projection force matching. All required source files are identical to the normal code (just symlinks to the main src dir), except main.f90, md_forces.f90, intra_morse.f90 and intra_morse_harm.f90.
Due to the code being cut to the minimum required for force calculation, a lot of input variables are useless and left out of the input list. See the force matching example for details.

NOTES:
- fixed nb=1 setting
- rigid mode not supported
- RPMDDFT force mode not supported
- EPSR correction not supported
- E3B correction not supported 
- purely harmonic intramolecular force split not implemented
