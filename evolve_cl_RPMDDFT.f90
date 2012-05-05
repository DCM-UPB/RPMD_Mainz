subroutine evolve_cl_RPMDDFT(p,r,v,v_lf,v_hf,dvdr,dvdr2,dt,mass,na,nb, &
                     boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
  use thermostat
  use barostat
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Classical evolution using RPMD-DFT-Method
  ! ------------------------------------------------------------------
  integer na,nb,irun,i,mts,nbaro,reftraj,ierr,rpmddft,iii1,jjj1
  real(8) p(3,na,nb), r(3,na,nb),dvdr(3,na,nb),dvdr2(3,na,nb)
  real(8) vir(3,3), vir_lf(3,3),vir_hf(3,3),vir_tmp(3,3)
  real(8) halfdtsmall,dt,boxlxyz(3),tvxyz(3),v,beta,dtsmall
  real(8) mass(na),z(na)
  real(8) halfdt,om,v_lf,v_hf
  real(8) tv,tq1,tq2
  character*4 type

  common /multiple_ts/ mts
  common /path_i/ om,type
  common /reftraj/ reftraj
  common /RPMDDFT/ rpmddft

  real(8) boxlxyz_backup(3)
  real(8), allocatable :: r_backup(:,:,:)
  allocate (r_backup(3,na,nb))


  vir_hf(:,:) = 0.d0  !vir_lf und vir gehen eh als 0.d0 rein, vir = volles virial


  halfdt = 0.5d0*dt

  ! Evolve the momenta to half time step
  if (therm.eq.'PRG') then                               ! hier wird themostat propagiert
     call parinello_therm(p,mass,ttau,na,nb,halfdt,irun,beta)
  else if (therm.eq.'PRL') then
     call parinello_therm_loc(p,mass,ttau,na,nb,halfdt,irun,beta)
  endif

  p(:,:,:) = p(:,:,:) - halfdt*(dvdr(:,:,:)+dvdr2(:,:,:)) !dvdr2 sollte nach erstem step = 0 sein, in erstem aber noch kraft aus äquilibrierung --> dvdr2 !=0
  
  !  get new coordinates
  call freerp_rpmd(p,r,dt,mass,na,nb,beta) ! Was ist mit ACMD???


!  ! Barostat
!  ! (note:// COMs scaled therefore do not need to recalculate
!  ! intramolecular forces as they remain the same)

  ! Barostat
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

  ! Evaluation of the t+timestep forces 
  call forces(r,v,dvdr,nb,na,boxlxyz,z,vir,9)  

  ! Set dvdr2 = 0 because we don't have a high frequency part
  dvdr2(:,:,:) = 0.d0

  ! Evolve the momenta to full timestep                        
  p(:,:,:) = p(:,:,:)-halfdt*dvdr(:,:,:)

  if (therm.eq.'PRG') then
     call parinello_therm(p,mass,ttau,na,nb,halfdt,irun,beta)
  else if (therm.eq.'PRL') then
     call parinello_therm_loc(p,mass,ttau,na,nb,halfdt,irun,beta)
  endif

  return
end subroutine evolve_cl_RPMDDFT
