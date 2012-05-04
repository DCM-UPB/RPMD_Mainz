subroutine forces(r,v,dvdr,nb,na,boxlxyz,z,virial,iopt)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Calculate the potential energy v and forces dvdr of the
  ! system.
  ! ------------------------------------------------------------------
  integer nb,na,k,j,iopt,rpmddft,nbdf3
  real(8) r(3,na,nb),dvdr(3,na,nb),z(na),boxlxyz(3)
  real(8) v,dv,oo_eps,oo_sig,oo_gam,rcut,sig,cut,boxmin
  real(8) rgmax,vir(3,3),virial(3,3)
  real(8), allocatable :: rc(:,:)
  integer, allocatable :: point(:),list(:)
  common /oo_param/ oo_eps,oo_sig,oo_gam,rcut
  common /correct/ sig
	common /RPMDDFT/ rpmddft,nbdf3
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
  ! - used for RPMD to speed LJ evaluations

  if ((iopt.eq.3).or.(iopt.le.1)) then
     if (nb.gt.3) then

        ! Oxygen centroids

        allocate (rc(3,na))
        rc(:,:) = 0.d0
        do j = 1,nb
           rc(:,:) = rc(:,:) + r(:,:,j)
        enddo
        rc(:,:) = rc(:,:)/dble(nb)
        cut = rcut
        boxmin = min(boxlxyz(1),boxlxyz(2),boxlxyz(3))
        cut = min(cut,0.5d0*boxmin)
        call nei_driver(rc,na,list,point,3,cut,boxlxyz)
        deallocate(rc)
     else
        list(1) = 0
     endif
  endif

  ! Calculate the potential energy for each bead system....

  if (iopt.ne.5) then
     do k = 1, nb
        call potenl_opt(r(:,:,k),dv,dvdr(:,:,k),vir(:,:), & 
                        na,nb,boxlxyz,z,list,point,iopt,k) ! f_env_id = k

				v = v + dv
        virial(:,:) = virial(:,:) + vir(:,:)
     enddo
  else

     ! Ewald Short Range Correction

     if (nb.gt.1) then

        ! Particle centroids

        allocate (rc(3,na))
        rc(:,:) = 0.d0
        do j = 1,nb
           rc(:,:) = rc(:,:) + r(:,:,j)
        enddo
        rc(:,:) = rc(:,:)/dble(nb)

        ! cut is set bigger than the cut-off

        cut = sig + rgmax
        boxmin = min(boxlxyz(1),boxlxyz(2),boxlxyz(3))
        if (cut.gt.0.5d0*boxmin) then
           write(6,*) 'ewc cutoff is too big for box'
           cut = min(cut,0.5d0*boxmin)
        endif
        call nei_driver(rc,na,list,point,1,cut,boxlxyz)
        deallocate (rc)
     else
        list(1) = 0
     endif

     do k = 1,nb
        call ewc_driver_tem(r(1,1,k),z,dvdr(1,1,k),dv, &
                            vir,na,boxlxyz,list,point)
        v = v + dv
        virial(:,:) = virial(:,:) + vir(:,:)
     enddo
  endif

  ! Average the PE over the bead systems....

  v = v / dble(nb)
  virial(:,:) = virial(:,:)/dble(nb)

  deallocate(point,list)
  
  return
end subroutine forces


subroutine potenl_opt(r,v,dvdr,vir,na,nb,boxlxyz, &
                      z,list,point,iopt,bead) ! da f_env_id (nb)
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
  ! iopt = 3   - LJ force evaluation
  ! iopt = 4   - Intramolecular force evaluation
	! iopt = 9   - RPMD-DFT force evaluation
  ! ------------------------------------------------------------------
  integer na,nb,mol(na),nm,i,j,imol,ic,iopt,point(na+3),list(maxnab*na),bead
  real(8) z(na),r(3,na),dvdr(3,na),vir(3,3),vir_ew(3,3),vir_oo(3,3)
  real(8) vir_int(3,3),boxlxyz(3)
  real(8) alpha,alpha2,wm,wh,ecut,voo,oo_eps,oo_sig,oo_gam,rcut
  real(8) v,vint,vew,apot,alp,bpot,alpb,boxl
  real(8), allocatable :: ro(:,:),dvdroo(:,:)
  logical iamcub
  common /ew_param/ alpha,ecut,wm,wh
  common /oo_param/ oo_eps,oo_sig,oo_gam,rcut
  common /intra_param/ apot,bpot,alp,alpb
  common /symmetry/ iamcub

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
  endif

#ifdef CP2K_BINDING
	!*** RPMD-DFT Force ***
	if (iopt.eq.9) then
				call RPMDDFT_force(r,dvdr,na,nb,v,vir,boxlxyz,bead) !!! übergebe auch boxlxyz, da NPT und bead, da PI

  endif
#endif

  return
end subroutine potenl_opt


subroutine full_forces(r,na,nb,v,vew,voo,vint,vir,z,boxlxyz, &
                       dvdr,dvdr2)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Calculates the full system energy, force and virial
  ! ------------------------------------------------------------------
  integer na,nb,nbdf1,nbdf2,rpmddft
  real(8) r(3,na,nb),dvdr(3,na,nb),dvdr2(3,na,nb),z(na),boxlxyz(3)
  real(8) vir(3,3),vir_oo(3,3),vir_ew(3,3)
  real(8) vir_itr(3,3),vir_ewc(3,3)
  real(8) v,vew,voo,vint,sig,ve
  real(8), allocatable :: dvdre(:,:,:),dvdrl(:,:,:)
  common /beaddiabatic/ nbdf1,nbdf2
  common /correct/ sig
	common /RPMDDFT/ rpmddft

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
	
	if (rpmddft.eq.0) then

  	! High Frequency forces

  	call forces(r,vint,dvdr2,nb,na,boxlxyz,z,vir_itr,4) 										!dvdr2 wird in forces = 0 gesetzt

  	! Evaluation of the low frequency forces using RP-contraction

  	if (nb.gt.1) then
    	 if ((nbdf1.gt.0).or.(nbdf2.gt.0)) then

    	    ! Entering normal mode representation

    	    call realft (r,3*na,nb,+1)

    	    ! Ewald

    	    if (nbdf1.gt.0) then
    	       call rp_contract_nm(r,vew,dvdr,nb,na,boxlxyz,z,vir_ew, &            !dvdr vorher = 0 gesetzt
                                nbdf1,2)

          	 if (sig.gt.0.d0) then

    	          ! Short range coulombic correction
	
     	         call rp_contract_nm(r,ve,dvdre,nb,na,boxlxyz,z, &
     	                              vir_ewc,nbdf1,5)
       	       vew = vew - ve
       	       vir_ew(:,:) = vir_ew(:,:) - vir_ewc(:,:)
       	       dvdr(:,:,:) = dvdr(:,:,:) - dvdre(:,:,:)
      	     endif
      	  endif

    	    ! Oxygen-Oxgen interaction

    	    if (nbdf2.gt.0) then
     	      call rp_contract_nm(r,voo,dvdrl,nb,na,boxlxyz,z,vir_oo, &
     	                           nbdf2,3)
     	   endif

     	   call realft (r,3*na,nb,-1)
    	 endif

   	  ! Transform dvdr back to bead representation (if needed)

  	   if ((nbdf1.gt.0).and.(nbdf2.gt.0)) then
  	      dvdr(:,:,:) = dvdr(:,:,:) + dvdrl(:,:,:)  
  	      call realft(dvdr,3*na,nb,-1)
		   else
		      if (nbdf1.gt.0) then
    	       call realft(dvdr,3*na,nb,-1)               
   	 	   else if (nbdf2.gt.0) then
     	      call realft(dvdrl,3*na,nb,-1)
     	   endif
    	 endif
	  endif

  	! Evaluate non beadd forces

  	if ((nb.eq.1).or.(nbdf1.eq.0)) then
 	    call forces(r,vew,dvdr,nb,na,boxlxyz,z,vir_ew,2) 
  	else if (sig.gt.0.d0) then
     
  	   ! Short range coulombic correction
     
  	   call forces(r,ve,dvdre,nb,na,boxlxyz,z,vir_ewc,5)
     
  	   vew = vew + ve
  	   vir_ew(:,:) = vir_ew(:,:) + vir_ewc(:,:)
   	  dvdr(:,:,:) = dvdr(:,:,:) + dvdre(:,:,:)							
 	 endif
  
 	 if ((nb.eq.1).or.(nbdf2.eq.0)) then
 	    call forces(r,voo,dvdrl,nb,na,boxlxyz,z,vir_oo,3)
 	 endif

  	! Sum ewald and O-O forces if needed

  if ((nb.eq.1).or.(nbdf1.le.0).or.(nbdf2.le.0)) then
     dvdr(:,:,:) = dvdr(:,:,:) + dvdrl(:,:,:)
  endif

 	 ! Potential Energy

  	v = vew + voo + vint

  	! Virial

  	vir(:,:) = vir_oo(:,:) + vir_ew(:,:) + vir_itr(:,:)

  	! Intermolecular PE

	!!!!!  vint = vew + voo
#ifdef CP2K_BINDING
  else

		! Calculate full force using CP2K
		call forces(r,v,dvdr,nb,na,boxlxyz,z,vir,9)
		dvdr2 = 0.d0
#endif
endif

  deallocate(dvdre,dvdrl)

  return
end subroutine full_forces


