subroutine evolve_cl(p,r,v,v_lf,v_hf,dvdr,dvdr2,dt,mass,na,nb, &
                     boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
  use thermostat
  use barostat
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Classical evolution using multiple time step method
  ! ------------------------------------------------------------------
  integer na,nb,irun,i,mts,nbaro,reftraj
  real(8) p(3,na,nb), r(3,na,nb),dvdr(3,na,nb),dvdr2(3,na,nb)
  real(8) vir(3,3), vir_lf(3,3),vir_hf(3,3),vir_tmp(3,3)
  real(8) halfdtsmall,dt,boxlxyz(3),tvxyz(3),v,beta,dtsmall
  real(8) mass(na),z(na)
  real(8) halfdt,om,v_lf,v_hf,v_corr
  real(8) dvdr_corr(3,na,nb),vir_corr(3,3)
  real(8) tv,tq1,tq2
  !real(8), allocatable :: monod(:,:,:),delp(:,:)
  character(len=4) type
  common /multiple_ts/ mts
  common /path_i/ om,type
  common /reftraj/ reftraj

  real(8) boxlxyz_backup(3)
  real(8), allocatable :: r_backup(:,:,:)
  allocate (r_backup(3,na,nb))

  !allocate (monod(3,4,nb),delp(3,nb))
  !monod(:,:,:) = 0.d0
  !delp(:,:) = 0.d0
  vir_hf(:,:) = 0.d0

  v_corr = 0.d0
  vir_corr(:,:) = 0.d0
  dvdr_corr(:,:,:) = 0.d0

  halfdt = 0.5d0*dt
  dtsmall = dt/dble(mts)
  halfdtsmall = 0.5d0*dtsmall

  ! Evolve the momenta under the low frequency force to half time step

  if (therm.eq.'PRG') then
     call parinello_therm(p,mass,ttau,na,nb,halfdt,irun,beta)
  else if (therm.eq.'PRL') then
     call parinello_therm_loc(p,mass,ttau,na,nb,halfdt,irun,beta)
  endif
  p(:,:,:) = p(:,:,:) - halfdt*dvdr(:,:,:)

  ! Mutiple time step under high frequency forces

  do i = 1,mts

     if (therm.eq.'PRG') then
        call parinello_therm(p,mass,ttau,na,nb, &
                             halfdtsmall,irun,beta)
     else if (therm.eq.'PRL') then
        call parinello_therm_loc(p,mass,ttau,na,nb, &
                                 halfdtsmall,irun,beta)
     endif
     p(:,:,:) = p(:,:,:)- halfdtsmall*dvdr2(:,:,:)

     ! Free ring polymer evolution
     if (reftraj.ne.0) then
       r_backup(:,:,:) = r(:,:,:)
       boxlxyz_backup(:) = boxlxyz(:)
     endif

     if (type.eq.'ACMD') then
        call freerp_acmd (p,r,dtsmall,mass,na,nb,beta,irun,om)
     else if (type.eq.'RPMD') then
        call freerp_rpmd (p,r,dtsmall,mass,na,nb,beta)
     endif

     ! Evaluate the high frequency force
     if (reftraj.ne.0) then
       r(:,:,:) = r_backup(:,:,:)
       boxlxyz(:) = boxlxyz_backup(:)
     endif

     call forces(r,v_hf,dvdr2,nb,na,boxlxyz,z,vir_tmp,4)
     vir_hf(:,:) = vir_hf(:,:) + vir_tmp(:,:)

     p(:,:,:) = p(:,:,:)- halfdtsmall*dvdr2(:,:,:)
     if (reftraj.eq.0) then
       if (therm.eq.'PRG') then
          call parinello_therm(p,mass,ttau,na,nb, &
                               halfdtsmall,irun,beta)
       else if (therm.eq.'PRL') then
          call parinello_therm_loc(p,mass,ttau,na,nb, &
                                   halfdtsmall,irun,beta)
       endif
     endif
  enddo

  ! Average intramolecular-virial over multiple time steps
  ! (increases stability of method)

  vir_hf(:,:) = vir_hf(:,:) / dble(mts)
  vir(:,:) = vir_lf(:,:) + vir_hf(:,:)

  ! Barostat
  ! (note:// COMs scaled therefore do not need to recalculate
  ! intramolecular forces as they remain the same)

  if (reftraj.eq.0 .and. nbaro.eq.1) then
     if (baro.eq.'BER') then
        call virial_ke(r,dvdr,dvdr2,tv,tvxyz,tq1,tq2,beta,na,nb)
        call beren_driver(vir,tv,tvxyz,dt,r,boxlxyz,na,nb)
     else if (baro.eq.'MCI') then
        call mc_baro(r,dvdr,dvdr2,vir,v,z,beta,boxlxyz,na,nb,irun)
     else
        write(6,*) ' ** Invalid barostat **'
        stop
     endif
  endif

  ! Evaluatation of the low frequency forces

  call forces(r,v_lf,dvdr,nb,na,boxlxyz,z,vir_lf,1)

  !Corrections
   !call forces(r,v_corr,dvdr_corr,nb,na,boxlxyz,z,vir_corr,8)
  
  v_lf = v_lf + v_corr
  v = v_lf + v_hf

  ! Virial
  vir_lf(:,:)= vir_lf(:,:) + vir_corr(:,:)
  vir(:,:)= vir_lf(:,:) + vir_hf(:,:)

  ! Evolve the momenta
  dvdr(:,:,:) = dvdr(:,:,:) + dvdr_corr(:,:,:)
  p(:,:,:) = p(:,:,:)-halfdt*dvdr(:,:,:)
  if (therm.eq.'PRG') then
     call parinello_therm(p,mass,ttau,na,nb,halfdt,irun,beta)
  else if (therm.eq.'PRL') then
     call parinello_therm_loc(p,mass,ttau,na,nb,halfdt,irun,beta)
  endif

  !deallocate (monod,delp)
  
  return
end subroutine evolve_cl
