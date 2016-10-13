subroutine freerp_acmd(p,r,dt,mass,na,nb,beta,irun,om)
  use thermostat
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Free harmonic ring-polymer evolution through a time interval dt.
  ! This version is for partially-adiabatic CMD, using Parinello's
  ! Langevin thermostat.
  ! ------------------------------------------------------------------      
  integer na,nb,i,j,k,irun,ic
  real(8) p(3,na,nb), r(3,na,nb),mass(na)
  real(8) dt,beta,om,monod(3,4,nb),delp(3,nb)
  real(8) pprime,gaussian
  external gaussian

  ! Zero the monodromy matrix and then calculate 
  ! the monodromy elements.

  monod(:,:,:) = 0.d0
  delp(:,:) = 0.d0
  
  call ring_acmd (mass,monod,na,nb,dt,beta,om,delp)
  
  if (nb .eq. 1) then

     ! Classical Evolution

     do j = 1,na
        ic = mod(j,3)
        if (ic.eq.0) then
           ic = 3
        endif
        do i = 1, 3
           r(i,j,1) = r(i,j,1) + p(i,j,1)*monod(ic,3,1)
        enddo
     enddo
  else

     ! ACMD Evolution

     call realft (p,3*na,nb,+1)
     call realft (r,3*na,nb,+1)

     ! Parinello's Langevin thermostat. ttau is originally
     ! set in the main.f routine. Note that setting ttau = infinity
     ! recovers deterministic dynamics.

     if (therm.eq.'PAR') then
        call parinello_therm_acmd(p,mass,ttau,delp,na,nb,dt,irun)
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
     
     ! Second part of Parrinello's thermostat

     if (therm.eq.'PAR') then
        call parinello_therm_acmd(p,mass,ttau,delp,na,nb,dt,irun)
     endif
     
     call realft (p,3*na,nb,-1)
     call realft (r,3*na,nb,-1)
  endif

  return
end subroutine freerp_acmd
     

subroutine ring_acmd (mass,monod,na,nb,dt,beta,om,delp)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Monodromy matrix elements for free ring-polymer evolution.
  ! This version is for partially-adiabatic CMD.
  ! ------------------------------------------------------------------
  integer na,nb,j,k
  real(8) beta, dt,pi,mass(na),monod(3,4,nb),delp(3,nb)
  real(8) pibyn,twown,wk,wm,cwt,swt,betan,om

  pi = dacos(-1.d0)
  betan = beta/nb
  twown = 2.d0/(betan*hbar)
  pibyn = pi/nb

  do j = 1, 3
     delp(j,1) = dsqrt(mass(j)/betan)
     monod(j,1,1) = 1.d0
     monod(j,2,1) = 0.d0
     monod(j,3,1) = dt/mass(j)
     monod(j,4,1) = 1.d0
     if (nb .gt. 1) then
        do k = 1,nb/2
           wk = twown*dsin(k*pibyn)
           wm = mass(j)*wk*wk/om
           cwt = dcos(om*dt)
           swt = dsin(om*dt)
           monod(j,1,k+1) = cwt
           monod(j,2,k+1) = -wm*swt
           monod(j,3,k+1) = swt/wm
           monod(j,4,k+1) = cwt
           delp(j,k+1)   = dsqrt(0.5d0*(wm/om)/betan)
        enddo
        if (2*(nb/2) .eq. nb) then
           delp(j,(nb/2)+1) = dsqrt(2.d0)*delp(j,(nb/2)+1)
        endif
        do k = 1,(nb-1)/2
           monod(j,1,nb-k+1) = monod(j,1,k+1)
           monod(j,2,nb-k+1) = monod(j,2,k+1)
           monod(j,3,nb-k+1) = monod(j,3,k+1)
           monod(j,4,nb-k+1) = monod(j,4,k+1)
           delp(j,nb-k+1)   = delp(j,k+1)
        enddo
     endif
  enddo
  
  return
end subroutine ring_acmd

subroutine freerp_rpmd (p,r,dt,mass,na,nb,beta)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Free-ring polymer evolution through a time-step dt.
  ! This version is for ring-polymer molecular dynamics.
  ! ------------------------------------------------------------------
  integer na,nb,j,i,k,ic
  real(8) p(3,na,nb), r(3,na,nb), mass(na)
  real(8) dt,beta,monod(3,4,nb)
  real(8) pprime
      
  ! Zero the monodromy matrix and then calculate 
  ! the monodromy elements:

  monod(:,:,:) = 0.d0
  
  call ring_rpmd (mass,monod,na,nb,dt,beta)
      
  if (nb .eq. 1) then
     do j = 1,na
        ic = mod(j,3)
        if (ic.eq.0) then
           ic = 3
        endif
        do i = 1, 3
           r(i,j,1) = r(i,j,1) + p(i,j,1)*monod(ic,3,1)
        enddo
     enddo
  else
     call realft (p,3*na,nb,+1)
     call realft (r,3*na,nb,+1)
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
     call realft (p,3*na,nb,-1)
     call realft (r,3*na,nb,-1)
  endif

  return
end subroutine freerp_rpmd

subroutine ring_rpmd (mass,monod,na,nb,dt,beta)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Monodromy matrix elements for free ring-polymer evolution.
  ! This version for ring-polymer molecular dynamics evolution.
  ! ------------------------------------------------------------------
  integer na,nb,j,k
  real(8) beta, dt,pi,mass(na),monod(3,4,nb)
  real(8) pibyn,twown,wk,wt,wm,cwt,swt,betan
      
  pi = dacos(-1.d0)
  do j = 1, 3
     monod(j,1,1) = 1.d0
     monod(j,2,1) = 0.d0
     monod(j,3,1) = dt/mass(j)
     monod(j,4,1) = 1.d0
     if (nb .gt. 1) then
        betan = beta/dble(nb)
        twown = 2.d0/(betan*hbar)
        pibyn = pi/dble(nb)
        do k = 1,nb/2
           wk = twown*dsin(k*pibyn)
           wt = wk*dt
           wm = wk*mass(j)
           cwt = dcos(wt)
           swt = dsin(wt)
           monod(j,1,k+1) = cwt
           monod(j,2,k+1) = -wm*swt
           monod(j,3,k+1) = swt/wm
           monod(j,4,k+1) = cwt
        enddo
        do k = 1,(nb-1)/2
           monod(j,1,nb-k+1) = monod(j,1,k+1)
           monod(j,2,nb-k+1) = monod(j,2,k+1)
           monod(j,3,nb-k+1) = monod(j,3,k+1)
           monod(j,4,nb-k+1) = monod(j,4,k+1)
        enddo
     endif
  enddo
  return
end subroutine ring_rpmd

