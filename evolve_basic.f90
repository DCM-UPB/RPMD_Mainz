subroutine evolve_basic(p,r,v,v_lf,v_hf,dvdr,dvdr2,dt,mass,na,nb, &
                        boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
  use thermostat
  use barostat
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Classical evolution - basic VV version
  ! ------------------------------------------------------------------
  integer na,nb,irun,i,j,k,nbaro
  real(8) p(3,na,nb), r(3,na,nb),dvdr(3,na,nb),dvdr2(3,na,nb)
  real(8) vir(3,3),vir_lf(3,3),vir_hf(3,3)
  real(8) dt,boxlxyz(3),tvxyz(3),v,beta
  real(8) mass(na),z(na),tau
  real(8) halfdt,v_lf,v_hf
  real(8) tv,tq1,tq2

  tau = 10.d0*tofs
  halfdt = 0.5d0*dt

  ! Parinello Thermostat step 1

  call parinello_therm(p,mass,tau,na,nb,halfdt,irun,beta)

  ! VV step 1

  p(:,:,:) = p(:,:,:) - halfdt*(dvdr(:,:,:)+dvdr2(:,:,:))

  do k = 1,nb
     do j = 1,na
        do i = 1,3
           r(i,j,k) = r(i,j,k) + p(i,j,k)*dt/mass(j)
        enddo
     enddo
  enddo

  ! Barostat

  if (nbaro.eq.1) then
     if (baro.eq.'BER') then
        call virial_ke(r,dvdr,dvdr2,tv,tvxyz,tq1,tq2,beta,na,nb)
        call beren_driver(vir,tv,tvxyz,dt,r,boxlxyz,na,nb)
     else if (baro.eq.'MCI') then
        call mc_baro(r,dvdr,dvdr2,vir,v,z,beta,boxlxyz, &
                     na,nb,irun)
     else
        write(6,*) ' ** Invalid barostat **'
        stop
     endif
  endif

  ! Evaluate the high frequency force

  call forces(r,v_hf,dvdr2,nb,na,boxlxyz,z,vir_hf,4)

  ! Evaluatation of the low frequency forces

  call forces(r,v_lf,dvdr,nb,na,boxlxyz,z,vir_lf,1)

  v = v_lf + v_hf

  ! Virial

  vir(:,:) = vir_lf(:,:) + vir_hf(:,:)

  ! VV step 2

  p(:,:,:) = p(:,:,:) - halfdt*(dvdr(:,:,:)+dvdr2(:,:,:))

  ! Parinello Thermostat step 2

  call parinello_therm(p,mass,tau,na,nb,halfdt,irun,beta)

  return
end subroutine evolve_basic
