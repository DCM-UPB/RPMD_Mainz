subroutine forces(r,v,dvdr,nb,na,boxlxyz,z,virial,iopt)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Calculate the potential energy v and forces dvdr of the
  ! system.
  ! ------------------------------------------------------------------
  integer nb,na,k,j,iopt,rpmddft,nbdf3,rctdk
  real(8) r(3,na,nb),dvdr(3,na,nb),z(na),boxlxyz(3)
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
  vir(:,:) = 0.d0
  virial(:,:) = 0.d0

  ! Create centroid neighbour list for oxygens if required
  ! - used to speed LJ evaluations
  if ((iopt.eq.3).or.(iopt.le.1)) then
      list(1) = 0 ! nb=1 case
  endif

  ! Calculate the potential energy
  if (iopt.ne.5) then
      call potenl_opt(r(:,:,1),dv,dvdr(:,:,1),vir(:,:), &
                        na,nb,boxlxyz,z,list,point,iopt,1) ! f_env_id = k
        v = v + dv
        virial(:,:) = virial(:,:) + vir(:,:)
  else
     ! Ewald Short Range Correction

     list(1) = 0 ! nb=1 case

     call ewc_driver_tem(r(1,1,1),z,dvdr(1,1,1),dv, &
                            vir,na,boxlxyz,list,point)
     v = v + dv
     virial(:,:) = virial(:,:) + vir(:,:)
  endif

  deallocate(point,list)
  
  return
end subroutine forces


subroutine potenl_opt(r,v,dvdr,vir,na,nb,boxlxyz, &
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
        call ewald_cubic(na,mol,z,r,vew,vir_ew,dvdr,boxl)
     else          
        call ewald_nc(na,mol,z,r,vew,vir_ew,dvdr,boxlxyz)
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
           dvdr(i,j+1) = dvdr(i,j+1) + alpha2 * dvdr(i,j)
           dvdr(i,j+2) = dvdr(i,j+2) + alpha2 * dvdr(i,j)
           dvdr(i,j) = alpha * dvdr(i,j)
           r(i,j) = ro(i,ic)
        enddo
     enddo
     deallocate (ro)
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
     v = v + voo
     vir(:,:) = vir(:,:) + vir_oo(:,:)
     deallocate (dvdroo)

    
    !e3b-Correction     
     if (rpmde3b.ne.0) then              
        call e3b_convert_pos (na,nm,r,boxlxyz,pos,box)
        call e3b_driver(nm,pos,f3B,vir_e3b,box,v_e3b)
        call e3b_convert_forces (nm,f3B,vir_e3b,v_e3b,dvdr_e3b)
        v = v + v_e3b
        vir(:,:)= vir(:,:) + vir_e3b(:,:)
        dvdr(:,:) = dvdr(:,:) - dvdr_e3b(:,:)     
     endif

  endif
  
  ! **Decide whether to calculate high frequency parts**
  if ((iopt.eq.0).or.(iopt.eq.4)) then

     ! * Intramolecular Forces
     if (alp.eq.0.d0) then
        call intra_harmonic(r,dvdr,na,vint,vir_int)
!!!        call intra_spcfw(r,dvdr,na,vint,vir_int)
     else
        if (alpb.ne.0.d0) then
           call intra_morse(r,dvdr,na,vint,vir_int)
        else
           call intra_morse_harm(r,dvdr,na,vint,vir_int)
        endif
     endif
     v = v + vint
     vir(:,:) = vir(:,:) + vir_int(:,:)


     if(epsr.eqv..true.) then
       allocate (dvdrepsr(3,na))
       ! Calculate epsr contribution and fix forces/virial/energy
       call epsr_driver(r,dvdrepsr,vepsr,vir_epsr,list,point,na,boxlxyz,3)
       dvdr(:,:) = dvdr(:,:) + dvdrepsr(:,:)
       v = v + vepsr
       vir(:,:) = vir(:,:) + vir_epsr(:,:)
       deallocate(dvdrepsr)
     endif

  endif

  !Corrections
  if (iopt.eq.8) then
    !e3b-Correction
    !if (rpmde3b.ne.0) then              
     !call e3b_convert_pos (na,nm,r,boxlxyz,pos,box)
     !call e3b_driver(nm,pos,f3B,vir_e3b,box,v_e3b)
     !call e3b_convert_forces (nm,f3B,vir_e3b,v_e3b,dvdr_e3b)       
     !dvdr(:,:) = -dvdr_e3b
     !v = v_e3b 
     !vir(:,:) = vir_e3b(:,:)
    !endif

  endif

  return
end subroutine potenl_opt


subroutine full_forces(r,na,nb,v,vew,voo,vint,vir,z,boxlxyz, &
                       dvdr,dvdr2)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Calculates the full system energy, force and virial
  ! ------------------------------------------------------------------
  integer na,nb,nbdf1,nbdf2,rpmddft,nbdf3,rctdk
  real(8) r(3,na,nb),dvdr(3,na,nb),dvdr2(3,na,nb),z(na),boxlxyz(3),rnm(3,na,nb),dvdrb(3,na,nbdf3),rb(3,na,nbdf3)
  real(8) vir(3,3),vir_oo(3,3),vir_ew(3,3)
  real(8) vir_itr(3,3),vir_ewc(3,3)
  real(8) v,vew,voo,vint,sig,ve,v_corr
  real(8) dvdr_corr(3,na,nb),vir_corr(3,3)
  real(8), allocatable :: dvdre(:,:,:),dvdrl(:,:,:)
  common /beaddiabatic/ nbdf1,nbdf2
  common /correct/ sig
  common /RPMDDFT/ rpmddft,nbdf3,rctdk

  ! Allocate

  allocate (dvdre(3,na,nb),dvdrl(3,na,nb))

  ! Zero constants and arrays

  vew = 0.d0
  voo = 0.d0
  vint = 0.d0

  dvdre(:,:,:) = 0.d0
  dvdrl(:,:,:) = 0.d0
  dvdr(:,:,:) = 0.d0                       

  vir(:,:) = 0.d0
  vir_itr(:,:) = 0.d0
  vir_oo(:,:) = 0.d0
  vir_ew(:,:) = 0.d0

  v_corr = 0.d0
  vir_corr(:,:) = 0.d0
  dvdr_corr(:,:,:) = 0.d0

    ! High Frequency forces

    call forces(r,vint,dvdr2,nb,na,boxlxyz,z,vir_itr,4)                     !dvdr2 is set to zero in forces

    ! Evaluate non beadd forces

    call forces(r,vew,dvdr,nb,na,boxlxyz,z,vir_ew,2)
    call forces(r,voo,dvdrl,nb,na,boxlxyz,z,vir_oo,3)

  !Corrections
   !call forces(r,v_corr,dvdr_corr,nb,na,boxlxyz,z,vir_corr,8)

    dvdr(:,:,:) = dvdr(:,:,:) + dvdr_corr(:,:,:)

    ! Sum ewald and O-O forces if needed

    dvdr(:,:,:) = dvdr(:,:,:) + dvdrl(:,:,:)

    ! Potential Energy
    voo = voo + v_corr
    v = vew + voo + vint

    ! Virial
    vir_oo(:,:) = vir_oo(:,:) + vir_corr(:,:)
    vir(:,:) = vir_oo(:,:) + vir_ew(:,:) + vir_itr(:,:)


    ! Intermolecular PE

  !!!!!  vint = vew + voo


  deallocate(dvdre,dvdrl)

  return
end subroutine full_forces
