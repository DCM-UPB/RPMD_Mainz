subroutine setup_pi(nb,beta,nbdf1,nbdf2,sig)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Setup of path integral parameters
  ! ------------------------------------------------------------------
  integer nb,nbdf1,nbdf2,nbdf3,rpmddft
  real(8) beta,om,sig
  character*4 type
  common /path_i/ om,type
	common /RPMDDFT/ rpmddft,nbdf3
	
  write(6,*)
  write(6,62) type,nb,om
62 format( ' Path Integral Parameters : ' /1x,&
          '---------------------------' /1x, &
          'type   =      ',a4 /1x, &
          'nb     = ',i9,' beads'/1x, &
          'omega  = ',f12.6,' t**(-1)'/1x)

  ! RP Contraction Options
  ! ------------------------

  if ((nbdf1.gt.0).or.(nbdf2.gt.0)) then
     sig = sig / ToA
     call setup_contraction(nbdf1,nbdf2,nb,sig)
  endif
	if (nbdf3.gt.0) then
		 write(6,'(a,i2,a)') &
     ' * RP contraction of CP2K-Force with ', nbdf3 , ' modes'
	endif
  return
end subroutine setup_pi

subroutine setup_contraction(nbdf1,nbdf2,nb,sig)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Checks the RP contraction settings
  ! ------------------------------------------------------------------
  real(8) alpha,ecut,wm,wh,sig
  integer nb,nbdf1,nbdf2
  common /ew_param/ alpha,ecut,wm,wh  

  write(6,*) 'Ring Polymer Contraction Settings : '
  write(6,*) '------------------------------------'

  if (nbdf1.gt.0) then
     write(6,'(a,i2,a)') &
       ' * RP contraction of Ewald with ', nbdf1 , ' modes'
     if (mod((nbdf1-1),2).ne.0) then
        write(6,*) ' * RP contraction beads must be ODD :',nbdf1
        stop
     endif
     if (nbdf1.gt.nb) then
        write(6,*) ' * RP contraction beads must not exceed nb:', &
                        nbdf1,nb
        stop
     endif
  else
     write(6,*) '* Standard MTS evolution for Ewald'
  endif

  if (nbdf2.gt.0) then
     write(6,'(a,i2,a)') &
          ' * RP contraction of LJ with ', nbdf2 , ' modes'
     if (mod((nbdf2-1),2).ne.0) then
        write(6,*) ' * RP contraction beads must be ODD :',nbdf2
        stop
     endif
     if (nbdf2.gt.nb) then
        write(6,*) ' * RP contraction beads must not exceed nb:', &
             nbdf2,nb
        stop
     endif
  else
     write(6,*) '* Standard MTS evolution for LJ'
  endif

  if (nb.eq.1) then
     if ((nbdf1.gt.0).or.(nbdf2.gt.0)) then
        write(6,*) ' * Beadiabatic cannot be used in classical'
        stop
     endif
  endif
  
  ! Electrostatic contraction scheme

  if (sig.gt.0.d0) then
     write(6,*) '* Electrostatic RP-contraction scheme used'
     write(6,'(a,f8.3,a)') '   with sigma =', sig , ' bohr'
     if ((wm.gt.0).or.(wh.gt.0)) then
        write(6,*)
        write(6,*) &
             'Cannot use RP-contraction with smeared charges'
        write(6,'(2f8.4)') wh, wm
        stop
     endif
  endif
  
  return
end subroutine setup_contraction

subroutine setup_ewald(na,boxlxyz)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Setup Ewald parameters
  !
  ! cvar is a separation parameter which should be determined
  ! to maximize speed.
  !
  ! The parameter ecut controls the accuracy of the calculation.
  ! A sensible value of ecut for MD simulations is exp(-pi*pi).
  ! This value is used by default if ecut <= 0 on entry.
  ! ------------------------------------------------------------------
  integer na,kmax
  real(8) alpha,ecut,wm,wh,rkmax,walpha,wrcut
  real(8) cvar,pi,s,rshort,rlong,boxlxyz(3)
  logical iamcub
  common /ew_param/ alpha,ecut,wm,wh
  common /ew_setup/ wrcut,walpha,rkmax,kmax
  common /symmetry/ iamcub

  write(6,*) 'Ewald Parameters : '
  write(6,*) '-------------------'

  pi = dacos(-1.d0)

  if ((boxlxyz(1).eq.boxlxyz(2)).and.(boxlxyz(1).eq.boxlxyz(3)) &
          .and.(iamcub)) then

     ! Cubic setup - parameters set-up for unit box.

     iamcub = .true.
     cvar = 1.2d0

     wrcut = min(0.5d0,cvar*na**(-1.d0/6.d0))
     if (ecut .le. 0.d0) then
        walpha = pi/wrcut
        kmax = int(walpha)
        rkmax = 2.d0*pi*walpha
     else
        s = dsqrt(-dlog(ecut))
        walpha = s/wrcut
        kmax = int(s*walpha/pi)
        rkmax = 2.d0*s*walpha
     endif
     write(6,*) '* Cubic Ewald with : '
     write(6,'(a,f9.4)')   ' ecut  = ', ecut
     write(6,'(a,f9.4,a)') ' wrcut  = ', wrcut*boxlxyz(1),' bohr '
     write(6,'(a,f9.4,a)') ' walpha = ', walpha/boxlxyz(1), &
                                             ' bohr**(-1)'
     write(6,'(a,f9.4,a)') ' rkmax  = ', rkmax/boxlxyz(1), &
                                             ' bohr**(-1)'
     write(6,'(a,i8,a)') ' kmax   = ', kmax , '  vectors'
     write(6,*) '(Approximation to error function used)'
  else

     ! Non-cubic setup - parameters set-up for non-unit box

     if (iamcub) then
        iamcub = .false.
     endif
     cvar = 1.2d0
     
     ! Shortest side length

     rshort = min(boxlxyz(1),boxlxyz(2),boxlxyz(3))
     wrcut = rshort*min(0.5d0,cvar*na**(-1.d0/6.d0))
     walpha = pi/wrcut

     ! Longest side length

     rlong = max(boxlxyz(1),boxlxyz(2),boxlxyz(3))
     rkmax = 2.d0*pi*walpha
     kmax = int(walpha*rlong)

     write(6,*) '* Non-Cubic Ewald with : '
     write(6,'(a,f9.4)')   ' ecut  = ', ecut
     write(6,'(a,f9.4,a)') ' wrcut  = ', wrcut,' bohr '
     write(6,'(a,f9.4,a)') ' walpha = ', walpha,' bohr**(-1)'
     write(6,'(a,f9.4,a)') ' rkmax  = ', rkmax,' bohr**(-1)'
     write(6,'(a,i8,a)') ' kmax   = ', kmax , '  vectors'
  endif
  write(6,*)

  return
end subroutine setup_ewald

subroutine setup_box_size(lattice,rho,nm,boxlxyz,wmass)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Converts cell numbers and density into box dimensions
  ! ------------------------------------------------------------------
  integer ncellxyz(3),nm,nlat,nflag
  real(8) boxlxyz(3),wmass,rho,vbox,rlat,alat,blat,clat
  real(8) boxl,vcell,totcell
  character*3 lattice
  logical iamcub
  common /lattice/ ncellxyz,nlat
  common /symmetry/ iamcub

  ! Lattice type

  if (lattice.eq.'CUB') then
     nlat = 1             ! Cubic Lattice
  else if (lattice.eq.'FCC') then
     nlat = 2             ! Face Centred Cubic Lattice
  else if (lattice.eq.'ICE') then
     nlat = 3             ! Ice lattice
  else
     write(6,*) 'Main : Invalid lattice choice'
     write(6,*) '*Options - CUB, FCC, ICE'
     write(6,*) lattice
     stop
  endif

  if (nlat.eq.1) then
     nm = ncellxyz(1)*ncellxyz(2)*ncellxyz(3)   ! Cubic Lattice
     nflag = 1
  else if (nlat.eq.2) then
     nm = 4*ncellxyz(1)*ncellxyz(2)*ncellxyz(3) ! FCC Lattice
     nflag = 1
  else if (nlat.eq.3) then
     nm = 8*ncellxyz(1)*ncellxyz(2)*ncellxyz(3) ! Ice lattice
     nflag = 2
  else
     write(6,*) 'Box_size - invalid nlat parameter'
  endif
  
  vbox = (nm*emass*wmass)/(rho*1000.d0) ! Tot. Vol. in m**3
  vbox = vbox / (ToA*1d-10)**3          ! Tot. Vol. in ao**3

  if (nflag.eq.1) then

     ! Cubic box setup
     totcell = dble(product(ncellxyz))
     boxl = (vbox/totcell)**(1.d0/3.d0)
     boxlxyz(1) = boxl*ncellxyz(1)
     boxlxyz(2) = boxl*ncellxyz(2)
     boxlxyz(3) = boxl*ncellxyz(3)
  else

     ! Orthorhombic hexagonal ice box setup
     ! Volume of unit cell :

     vcell = vbox / (ncellxyz(1)*ncellxyz(2)*ncellxyz(3))
     rlat = ((3.d0*dsqrt(3.d0)/64.d0) * vcell)**(1.d0/3.d0)
     
     ! Lattice Parameters
     alat = dsqrt(8.d0/3.d0) * rlat
     blat = dsqrt(8.d0) * rlat
     clat = (8.d0/3.d0) * rlat

     ! Boxlengths
     boxlxyz(1) = alat*ncellxyz(1)
     boxlxyz(2) = blat*ncellxyz(2)
     boxlxyz(3) = clat*ncellxyz(3)

     if (iamcub) then
        write(6,*) '** Cannot do ice in cubic box **'
        iamcub = .false.
     endif

  endif

  return
end subroutine setup_box_size

