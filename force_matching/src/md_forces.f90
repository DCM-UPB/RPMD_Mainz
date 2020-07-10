! Edited FFM version:
! dvdr_split(:,:,:,*) stores
! *=1 :  Coulomb
! *=2 :  LJ
! *=3 :  stretch
! *=4 :  bend

subroutine forces(r,v,dvdr,dvdr_split,nb,na,boxlxyz,z,virial,iopt)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Calculate the potential energy v and forces dvdr of the
  ! system.
  ! ------------------------------------------------------------------
  integer nb,na,k,j,iopt,rpmddft,nbdf3,rctdk
  real(8) r(3,na,nb),dvdr(3,na,nb),z(na),boxlxyz(3)
  real(8) dvdr_split(3,na,nb,4)
  real(8) v,dv,oo_eps,oo_sig,oo_gam,rcut,sig,cut,boxmin
  real(8) rgmax,vir(3,3),virial(3,3)
  real(8), allocatable :: rc(:,:)
  integer, allocatable :: point(:),list(:)
  common /oo_param/ oo_eps,oo_sig,oo_gam,rcut
  common /correct/ sig
  common /RPMDDFT/ rpmddft,nbdf3,rctdk
  ! Allocate

  allocate (point(na+3),list(maxnab*na))

  ! Rgmax = 0 gives a centroid cut-off in the potential
  ! for larger values it gives a bead cut-off

  rgmax = 5.d0

  v = 0.d0
  dvdr(:,:,:) = 0.d0
  dvdr_split(:,:,:,:) = 0.d0 ! if called with iopt!=0, some splits will remain 0
  vir(:,:) = 0.d0
  virial(:,:) = 0.d0

  ! Create centroid neighbour list for oxygens if required
  ! - used to speed LJ evaluations
  if ((iopt.eq.3).or.(iopt.le.1)) then
      list(1) = 0 ! nb=1 case
  endif

  ! Calculate the potential energy
  if (iopt.ne.5) then
      call potenl_opt(r(:,:,1),dv,dvdr(:,:,1),dvdr_split(:,:,1,:),vir(:,:), &
                        na,nb,boxlxyz,z,list,point,iopt,1) ! f_env_id = k
        v = v + dv
        virial(:,:) = virial(:,:) + vir(:,:)
  else
     ! Ewald Short Range Correction
      write(6,*) 'Error: forces called with iopt=5 in a FFM run. Aborting...'
      stop
  endif

  deallocate(point,list)
  
  return
end subroutine forces


subroutine potenl_opt(r,v,dvdr,dvdr_split,vir,na,nb,boxlxyz, &
                      z,list,point,iopt,bead) 
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Calculates the potential energy of a collection of nm = na/3
  ! flexible water molecules, with periodic boundary conditions.
  !
  ! Potential routine with iopt
  !  
  ! iopt = 0   - Full force evaluation
  ! iopt = 1   - Ewald and LJ force evaluation
  ! iopt = 2   - Ewald force evaluation
  ! iopt = 3   - LJ force evaluation and possibly EPSR correction
  ! iopt = 4   - Intramolecular force evaluation
  ! iopt = 5   - Ewald short range correction, doesn't call potenl_opt
  ! iopt = 8   - (Possibly EPSR correction)<-- not implemented
  ! iopt = 9   - RPMD-DFT force evaluation for AI-RPMD use
  ! ------------------------------------------------------------------
  integer na,nb,mol(na),nm,i,j,imol,ic,iopt,point(na+3),list(maxnab*na),bead,rpmde3b
  real(8) z(na),r(3,na),dvdr(3,na),vir(3,3),vir_ew(3,3),vir_oo(3,3),vir_epsr(3,3)
  real(8) dvdr_split(3,na,4)
  real(8) vir_int(3,3),boxlxyz(3)
  real(8) alpha,alpha2,wm,wh,ecut,voo,vepsr,oo_eps,oo_sig,oo_gam,rcut
  real(8) v,vint,vew,apot,alp,bpot,alpb,boxl
  real(8) f3B(na/3,4,3),vir_e3b(3,3),v_e3b,dvdr_e3b(3,na)
  real(8) box(3),pos(na/3,4,3) 
  real(8), allocatable :: ro(:,:),dvdroo(:,:),dvdrepsr(:,:)
  logical iamcub
  common /ew_param/ alpha,ecut,wm,wh
  common /oo_param/ oo_eps,oo_sig,oo_gam,rcut
  common /intra_param/ apot,bpot,alp,alpb
  common /symmetry/ iamcub
  logical epsr
  common /EPSR/ epsr
  common /E3B/ rpmde3b

  ! Zero the energy and forces

  v = 0.d0
  dvdr(:,:) = 0.d0
  dvdr_split(:,:,:) = 0.d0 ! if called with iopt!=0, some splits will remain 0
  vir(:,:) = 0.d0

  ! Check that na is divisible by 3!!!

  nm = na/3
  if (3*nm .ne. na) stop 'potenl 1'

  ! ** Decide whether to evaluate Ewald **

  if (iopt.le.2) then

     alpha2 = 0.5d0 * (1.d0 - alpha)

     ! Assign molecular identities to each atom.

     imol = 0
     do j = 1,na,3
        imol = imol+1
        mol(j) = imol
        mol(j+1) = imol
        mol(j+2) = imol
     enddo

     ! Store the oxygen coordinates and replace the oxygen
     ! coordinates with the M-coordinates.

     allocate (ro(3,na/3))

     ic = 0
     do i = 1, na, 3
        ic = ic + 1
        do j = 1, 3
           ro(j,ic) = r(j,i)
           r(j,i) = alpha * r(j,i) + alpha2*(r(j,i+1)+r(j,i+2))
        enddo
     enddo

     ! Calculate Coulomb energy and gradients

     if (iamcub) then
        boxl = boxlxyz(1)
        call ewald_cubic(na,mol,z,r,vew,vir_ew,dvdr_split(:,:,1),boxl)
     else          
        call ewald_nc(na,mol,z,r,vew,vir_ew,dvdr_split(:,:,1),boxlxyz)
     endif
     v = v + vew
     vir(:,:) = vir(:,:) + vir_ew(:,:)

     ! The Ewald routine has just calculated the forces on the
     ! m-site and the hydrogens. Use the chain rule to
     ! calculate the correct forces on the atoms. Also, replace
     ! the m-site with the original oxygen atoms.

     ic = 0
     do j = 1, na, 3
        ic = ic + 1
        do i = 1, 3
           dvdr_split(i,j+1,1) = dvdr_split(i,j+1,1) + alpha2 * dvdr_split(i,j,1)
           dvdr_split(i,j+2,1) = dvdr_split(i,j+2,1) + alpha2 * dvdr_split(i,j,1)
           dvdr_split(i,j,1) = alpha * dvdr_split(i,j,1)
           r(i,j) = ro(i,ic)
        enddo
     enddo
     deallocate (ro)

     dvdr(:,:) = dvdr(:,:) + dvdr_split(:,:,1) ! add ewald split

  endif

  !*** LJ calculation ***!

  if ((iopt.eq.0).or.(iopt.eq.1).or.(iopt.eq.3)) then
     ! OO interaction
     
     allocate (dvdroo(3,na))
     if (oo_gam.eq.0.d0) then
        call lj_driver(r,dvdroo,voo,vir_oo,list,point,na,boxlxyz,3)
     else
        call buck_driver(r,dvdroo,voo,vir_oo,list,point,na, &
                         boxlxyz,3)
     endif
     dvdr(:,:) = dvdr(:,:) + dvdroo(:,:)
     dvdr_split(:,:,2) = dvdroo(:,:)! the LJ split
     v = v + voo
     vir(:,:) = vir(:,:) + vir_oo(:,:)
     deallocate (dvdroo)

  endif
  
  ! **Decide whether to calculate high frequency parts**
  if ((iopt.eq.0).or.(iopt.eq.4)) then

     ! * Intramolecular Forces
     if (alp.eq.0.d0) then
         write(6,*)'Purely harmonic intramolecular potential not implemented for force matching. Aborting...'
         stop
     else
        if (alpb.ne.0.d0) then
           call intra_morse(r,dvdr,dvdr_split(:,:,3), dvdr_split(:,:,4),na,vint,vir_int)
        else
           call intra_morse_harm(r,dvdr,dvdr_split(:,:,3), dvdr_split(:,:,4),na,vint,vir_int)
        endif
     endif
     v = v + vint
     vir(:,:) = vir(:,:) + vir_int(:,:)

  endif

  !Corrections
  if (iopt.eq.8) then
      write(6,*) 'iopt=8 correction not implemented. Aborting...'
      stop
  endif

  return
end subroutine potenl_opt