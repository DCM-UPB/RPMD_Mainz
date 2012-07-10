subroutine evolve_pi_RPMDDFT(p,r,v,vew,vlj,vint,dvdr,dvdr2,dt,mass,na,nb, &
                     boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
  use thermostat
  use barostat
  use gle
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! RPMD/ACMD evolution using RPMDDFT method
  ! 
  ! ------------------------------------------------------------------
  integer na,nb,irun,k,j,i,nbdf1,nbdf2,im,ic,mts,nbaro,rpmddft
  real(8) p(3,na,nb),r(3,na,nb),dvdr(3,na,nb),dvdr2(3,na,nb)
  real(8) mass(na),z(na),boxlxyz(3)
  real(8) vir(3,3),vir_lj(3,3),vir_ew(3,3),tvxyz(3),tvxyz_itr(3)
  real(8) vir_itr(3,3),vir_ewc(3,3),vir_lf(3,3),vir_hf(3,3)
  real(8) pprime,halfdtsmall,halfdt,om,gaussian
  real(8) dt,v,beta,dtsmall,vew,vlj,vint,sig,ve
  real(8) tv,tv_itr,tq1,tq2,c1,c2,c3
  real(8) dheat, comx, comy, comz, mm !!GLE
  character(len=4) type
  external gaussian
  common /path_i/ om,type
  common /beaddiabatic/ nbdf1,nbdf2
  common /multiple_ts/ mts
  common /correct/ sig
  common /RPMDDFT/ rpmddft


  vew = 0.d0
  vlj = 0.d0
  vint = 0.d0
  vir_itr(:,:) = 0.d0
  vir_lj(:,:) = 0.d0
  vir_ew(:,:) = 0.d0
  vir_hf(:,:) = 0.d0
  tv = 0.d0 ! im ersten schritt falsch!!!!
  tvxyz(:) = 0.d0

  halfdt = 0.5d0*dt
                                   
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

  p(:,:,:) = p(:,:,:) - halfdt*(dvdr(:,:,:)+dvdr2(:,:,:)) !because dvdr2 is != 0 after equilibration

  !  get new coordinates
  call freerp_rpmd(p,r,dt,mass,na,nb,beta)

  ! Barostat
  ! (note:// COMs scaled therefore do not need to recalculate
  !  intramolecular forces as they remain the same)

  !get kinetic energy
  

  if (nbaro.eq.1) then
     if (baro.eq.'BER') then
				call virial_ke(r,dvdr,dvdr2,tv,tvxyz,tq1,tq2,beta,na,nb,mass)
        call beren_driver(vir,tv,tvxyz,dt,r,boxlxyz,na,nb) !tv != 0 sein, daher vorher berechnen
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
  return
end subroutine evolve_pi_RPMDDFT
