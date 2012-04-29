subroutine evolve_pi(p,r,v,vew,vlj,vint,dvdr,dvdr2,dt,mass,na,nb, &
                     boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
  use thermostat
  use barostat
  use gle
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! RPMD/ACMD evolution using multiple time step method using normal
  ! mode primitive representation
  ! ------------------------------------------------------------------
  integer na,nb,irun,k,j,nbdf1,nbdf2,im,ic,mts,nbaro
  real(8) p(3,na,nb),r(3,na,nb),dvdr(3,na,nb),dvdr2(3,na,nb)
  real(8) mass(na),z(na),boxlxyz(3)
  real(8) vir(3,3),vir_lj(3,3),vir_ew(3,3),tvxyz(3),tvxyz_itr(3)
  real(8) vir_itr(3,3),vir_ewc(3,3),vir_lf(3,3),vir_hf(3,3)
  real(8) pprime,halfdtsmall,halfdt,om,gaussian
  real(8) dt,v,beta,dtsmall,vew,vlj,vint,sig,ve
  real(8) tv,tv_itr,tq1,tq2
  real(8), allocatable :: dvdre(:,:,:),rst(:,:,:),dvdrl(:,:,:)
  real(8), allocatable :: monod(:,:,:),delp(:,:)
  real(8) dheat, comx, comy, comz, mm !!GLE
  character*4 type
  external gaussian
  common /path_i/ om,type
  common /beaddiabatic/ nbdf1,nbdf2
  common /multiple_ts/ mts
  common /correct/ sig

  ! Allocate

  allocate (dvdre(3,na,nb),rst(3,na,nb),dvdrl(3,na,nb))
  allocate (monod(3,4,nb),delp(3,nb))

  ! Zero constants and arrays

  dvdre(:,:,:) = 0.d0
  monod(:,:,:) = 0.d0
  delp(:,:) = 0.d0

  vew = 0.d0
  vlj = 0.d0
  vint = 0.d0
  vir_itr(:,:) = 0.d0
  vir_lj(:,:) = 0.d0
  vir_ew(:,:) = 0.d0
  vir_hf(:,:) = 0.d0
  tv = 0.d0
  tvxyz(:) = 0.d0

  dtsmall = dt/dble(mts)
  halfdtsmall = 0.5d0*dtsmall
  halfdt = 0.5d0*dt

  if (type.eq.'RPMD') then
     call ring_rpmd (mass,monod,na,nb,dtsmall,beta)
  else if (type.eq.'ACMD') then
     call ring_acmd (mass,monod,na,nb,dtsmall,beta,om,delp)
  else
     write(6,*) 'evolve_mm : INVALID TYPE '
     stop
  endif

  ! Evolve the momenta under the low frequency force to half time step
  ! and under the high frequency to the first small time step

  if (type.eq.'RPMD') then
     if (therm.eq.'PRG') then
        call parinello_therm(p,mass,ttau,na,nb,halfdt,irun,beta)
     else if (therm.eq.'PRL') then
        call parinello_therm_loc(p,mass,ttau,na,nb, &
                                 dt,irun,beta)
     else if (therm.eq.'GLE') then     !! GLE
         call therm_gle(p,dheat,mass,na,nb,irun)
     endif
  endif
  p(:,:,:) = p(:,:,:) - halfdt*dvdr(:,:,:)
  p(:,:,:) = p(:,:,:) - halfdtsmall*dvdr2(:,:,:)

  ! Convert the positions and momenta to the normal mode representation
  
  call realft (p,3*na,nb,+1)
  call realft (r,3*na,nb,+1)

  ! Multiple time step under high frequency forces
  ! -----------------------------------------------

  do im = 1,mts

     ! Free ring polymer evolution in normal mode representation

     if (type.eq.'ACMD') then
        if (therm.eq.'PRA') then
           ! Thermostat the non-centroid modes in ACMD
           call parinello_therm_acmd(p,mass,ttau,delp,na,nb, &
                                     dtsmall,irun)
        endif
     elseif (type.eq.'RPMD' .and. therm.eq.'PRO') then  !!PRO
         call therm_pro(p,mass,beta,ttau,halfdtsmall,na,nb,irun)
     elseif (type.eq.'RPMD' .and. therm.eq.'PRM') then
         call therm_pro_global(p,mass,beta,ttau,halfdtsmall,na,nb,irun)
     endif
         
     do k = 1,nb
        do j = 1,na
           ic = mod(j,3)
           if (ic.eq.0) then
              ic = 3
           endif
           pprime = p(1,j,k)*monod(ic,1,k)+r(1,j,k)*monod(ic,2,k)
           r(1,j,k) = p(1,j,k)*monod(ic,3,k)+r(1,j,k)*monod(ic,4,k)
           p(1,j,k) = pprime
           pprime = p(2,j,k)*monod(ic,1,k)+r(2,j,k)*monod(ic,2,k)
           r(2,j,k) = p(2,j,k)*monod(ic,3,k)+r(2,j,k)*monod(ic,4,k)
           p(2,j,k) = pprime
           pprime = p(3,j,k)*monod(ic,1,k)+r(3,j,k)*monod(ic,2,k)
           r(3,j,k) = p(3,j,k)*monod(ic,3,k)+r(3,j,k)*monod(ic,4,k)
           p(3,j,k) = pprime
        enddo
     enddo
     
     if (type.eq.'ACMD') then
        if (therm.eq.'PRA') then
           ! Thermostat the non-centroid modes in ACMD
           call parinello_therm_acmd(p,mass,ttau,delp,na,nb, &
                                     dtsmall,irun)
        endif
     elseif (type.eq.'RPMD' .and. therm.eq.'PRO')  then !!PRO
         call therm_pro(p,mass,beta,ttau,halfdtsmall,na,nb,irun)
     elseif (type.eq.'RPMD' .and. therm.eq.'PRM')  then
         call therm_pro_global(p,mass,beta,ttau,halfdtsmall,na,nb,irun)
     endif
     
     ! Evaluate the high frequency force

     rst(:,:,:) = r(:,:,:)

     ! Form Bead Positions and evaluate forces
     
     call realft(rst,3*na,nb,-1)
     call forces(rst,vint,dvdr2,nb,na,boxlxyz,z,vir_itr,4)
     call virial_ke(rst,dvdr,dvdr2,tv_itr,tvxyz_itr, &
                    tq1,tq2,beta,na,nb,mass)
     vir_hf(:,:) = vir_hf(:,:) + vir_itr(:,:)
     tv = tv + tv_itr
     tvxyz(:) = tvxyz(:) + tvxyz_itr(:)

     ! FT forces to normal mode representation and
     ! Evolve normal mode momenta

     if (im.ne.mts) then
        call realft(dvdr2,3*na,nb,+1)
        p(:,:,:) = p(:,:,:) - dtsmall*dvdr2(:,:,:)
     endif
  enddo

  ! Average intramolecular-virial over multiple time steps
  ! (increases stability of method)
  
  tv = tv/dble(mts)
  tvxyz(:) = tvxyz(:)/dble(mts)
  vir_hf(:,:) = vir_hf(:,:)/dble(mts)
  vir(:,:) = vir_lf(:,:) + vir_hf(:,:)

  ! Barostat
  ! (note:// COMs scaled therefore do not need to recalculate
  !  intramolecular forces as they remain the same)

  if (nbaro.eq.1) then
     if (baro.eq.'BER') then
        call beren_driver(vir,tv,tvxyz,dt,rst,boxlxyz,na,nb)
        r(:,:,:) = rst(:,:,:)
        call realft(r,3*na,nb,+1)
     else if (baro.eq.'MCI') then
        call mc_baro(rst,dvdr,dvdr2,vir,v,z,beta,boxlxyz,na,nb,irun)
        r(:,:,:) = rst(:,:,:)
        call realft(r,3*na,nb,+1)
     else 
        write(6,*) ' ** Invalid barostat **'
        stop
     endif
  endif

  ! Evaluation of the low frequency forces using RP-contraction
  ! --------------------------------------------------------------

  ! Ewald

  if (nbdf1.gt.0) then
     call rp_contract_nm(r,vew,dvdr,nb,na,boxlxyz,z,vir_ew, &
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

  ! LJ

  if (nbdf2.gt.0) then
     call rp_contract_nm(r,vlj,dvdrl,nb,na,boxlxyz,z,vir_lj, &
                          nbdf2,3)
  endif

  ! Convert the positions and mometa back to the bead
  ! representation and perform the final momentum evolution
  ! with the high frequency forces

  call realft (p,3*na,nb,-1)

  ! Rather than transforming back the positions
  ! we copy the ones used to obtain the forces

  r(:,:,:) = rst(:,:,:)
  p(:,:,:) = p(:,:,:) - halfdtsmall*dvdr2(:,:,:)

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

  ! Evaluate non beadd forces

  if (nbdf1.le.0) then
     call forces(r,vew,dvdr,nb,na,boxlxyz,z,vir_ew,2)
  else if (sig.gt.0.d0) then

     ! Short range coulombic correction

     call forces(r,ve,dvdre,nb,na,boxlxyz,z,vir_ewc,5)

     vew = vew + ve
     vir_ew(:,:) = vir_ew(:,:) + vir_ewc(:,:)
     dvdr(:,:,:) = dvdr(:,:,:) + dvdre(:,:,:)
  endif

  if (nbdf2.le.0) then
     call forces(r,vlj,dvdrl,nb,na,boxlxyz,z,vir_lj,3)
  endif

  ! Sum ewald and LJ forces if needed

  if ((nbdf1.le.0).or.(nbdf2.le.0)) then
     dvdr(:,:,:) = dvdr(:,:,:) + dvdrl(:,:,:)
  endif

  ! Evolve the momenta under the low frequency forces

  p(:,:,:) = p(:,:,:) - halfdt*dvdr(:,:,:)
  
  if (type.eq.'RPMD') then
     if (therm.eq.'PRG') then
        call parinello_therm(p,mass,ttau,na,nb,halfdt,irun,beta)
     else if (therm.eq.'PRL') then
        call parinello_therm_loc(p,mass,ttau,na,nb, &
                                       dt,irun,beta)
     else if (therm.eq.'GLE') then     !! GLE
         call therm_gle(p,dheat,mass,na,nb,irun)
     endif
  endif
  ! Potential Energy

  v = vew + vlj + vint

  ! Virial

  vir_lf(:,:) = vir_lj(:,:) + vir_ew(:,:)
  vir(:,:) =  vir_lf(:,:) + vir_hf(:,:)

  deallocate(dvdre,rst,dvdrl)
  deallocate(monod,delp)

  !compute COM and outputs to make sure that we removed COM velocity correctly....
  comx=0.d0
  comy=0.d0
  comz=0.d0
  mm=0.0d0
  do j = 1,na
    do k = 1,nb
      comx=comx+mass(j)*r(1,j,k)
      comy=comy+mass(j)*r(2,j,k)
      comz=comz+mass(j)*r(3,j,k)
    enddo
    mm=mm+mass(j)
  enddo
  comx=comx/mm
  comy=comy/mm
  comz=comz/mm
!  write(*,*) "COM: ",comx,comy,comz
  return
end subroutine evolve_pi

subroutine rp_contract_nm(r,v,dvdr,nb,na,boxlxyz,z,vir, &
                           nbdf,iopt)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Calculates the RP contraction potential energy v and forces
  ! assumes you are using the NORMAL MODE primitive representation.
  ! i.e. r are normal mode coordinates and routine produces normal mode
  !      forces dvdr.
  ! ------------------------------------------------------------------
  integer nb,na,nbdf,iopt
  real(8) r(3,na,nb),dvdr(3,na,nb),vir(3,3),z(na),boxlxyz(3)
  real(8) v
  real(8), allocatable :: rb(:,:,:),dvdrb(:,:,:)

  ! nbdf = RP Contraction (beadibatic) factor
  ! (determines the number of normal modes to keep)

  ! Check that beadiabatic factor is odd

  if (mod(nbdf-1,2).ne.0) then
     write(6,*) '* Beadiabatic factor must be ODD'
     stop
  endif

  ! Allocate

  allocate (rb(3,na,nbdf),dvdrb(3,na,nbdf))

  ! Zero arrays

  dvdr(:,:,:) = 0.d0
  rb(:,:,:) = 0.d0
  dvdrb(:,:,:) = 0.d0
  vir(:,:) = 0.d0

  ! Form contracted RP coordinates

  call ring_contract(r,rb,na,nb,nbdf)

  ! Evaluate forces on the bead coordinates rb

  call forces(rb,v,dvdrb,nbdf,na,boxlxyz,z,vir,iopt)

  ! Create full bead forces from contracted coordinates

  call force_contract(dvdr,dvdrb,na,nb,nbdf)

  deallocate (rb,dvdrb)

  return
end subroutine rp_contract_nm

Subroutine force_contract(dvdr,dvdr_c,na,nb,nb_c)
  implicit none
  integer :: na, nb, nb_c
  integer :: j, k, i
  double precision :: dvdr(3,na,nb), dvdr_c(3,na,nb_c)
  double precision :: scale_f

  call realft(dvdr_c,3*na,nb_c,+1)

  scale_f  = dsqrt(dble(nb)/dble(nb_c))

  do j = 1,na
     do i = 1,3
        dvdr(i,j,1) = scale_f*dvdr_c(i,j,1)
     enddo
  enddo

  do k = 1,(nb_c-1)/2
     do j = 1,na
        do i = 1,3
           dvdr(i,j,k+1) = scale_f * dvdr_c(i,j,k+1)
           dvdr(i,j,nb-k+1) = scale_f * dvdr_c(i,j,nb_c-k+1)
        enddo
     enddo
  enddo

end Subroutine force_contract

Subroutine ring_contract(r,r_c,na,nb,nb_c)
  implicit none
  integer :: na, nb, nb_c
  integer :: j, k, i
  double precision :: r(3,na,nb), r_c(3,na,nb_c)
  double precision :: scale_f

  if (mod(nb_c-1,2).ne.0) then
     write(6,*) '* ring_contract : nc_c must be ODD'
     stop
  endif

  r_c(:,:,1) = r(:,:,1)
  do k = 1,(nb_c-1)/2
     do j = 1,na
        do i = 1,3
          r_c(i,j,k+1) = r(i,j,k+1)
          r_c(i,j,nb_c-k+1) = r(i,j,nb-k+1)
        enddo
     enddo
  enddo

  call realft (r_c,3*na,nb_c,-1)
  scale_f = dsqrt(dble(nb_c)/dble(nb))
  r_c(:,:,:) = scale_f*r_c(:,:,:)

end Subroutine ring_contract

