module thermostat
    implicit none
    real(8) ttau
    character*3 therm
contains
  subroutine parinello_therm(p,mass,tau,na,nb,dt,irun,beta)
    implicit none
    ! ------------------------------------------------------------------
    ! Global Parinello thermostat
    ! ------------------------------------------------------------------
    integer na,nb,irun,ndof,i,j,k,nctot,nbond
    real(8) p(3,na,nb),mass(na),dt,beta,ranno,ransq
    real(8) tau,rt,alpha,c,ktav,kt,snfm1,gaussian,fac,t1,t2,plusmin
    real(8) betan
    external gaussian
    common /constraint/ nctot,nbond
    
    betan = beta/dble(nb)
    ndof = 3*na*dble(nb)-nctot-3
    rt = gaussian(irun,1.d0)
    ktav = dble(ndof)/(2.d0*betan)
    snfm1 = 0.d0
    do i = 1,ndof-1
       ranno = gaussian(irun,1.d0)
       snfm1 = snfm1 + ranno*ranno
    enddo
    
    ! Current kinetic energy
    
    kt = 0.d0
    do k = 1,nb
       do j = 1,na
          do i = 1,3
             kt = kt + p(i,j,k)*p(i,j,k)/(2.d0*mass(j))
          enddo
       enddo
    enddo
    
    ranno = gaussian(irun,1.d0)
    ransq = ranno*ranno
    c = dexp(-dt/tau)
    fac = (1.d0-c)*ktav/(dble(ndof)*kt)
    t1 = fac*(snfm1+ransq)
    t2 = 2*ranno*dsqrt(c*fac)
    alpha = c + t1 + t2
    alpha = dsqrt(alpha)

    plusmin = ranno + dsqrt(c/fac)
    alpha = sign(alpha,plusmin)
    
    p(:,:,:) = alpha * p(:,:,:)

    return
  end subroutine parinello_therm
  
  subroutine parinello_therm_loc(p,mass,tau,na,nb,dt,irun,beta)
    implicit none
    ! ------------------------------------------------------------------
    ! Parinello Local Langevin Thermostat
    ! ------------------------------------------------------------------
    integer na,nb,irun,k,j,i
    real(8) p(3,na,nb),mass(na),tau,gamma,dt,beta
    real(8) c1,c2,c3,gaussian,betan
    external gaussian
    
    gamma = 1.d0/(2.d0*tau)
    betan = beta/dble(nb)
    
    c1 = dexp (-0.5d0 * gamma * dt)
    c2 = dsqrt(1.d0 - c1*c1)
    
    do k = 1,nb
       do j = 1, na
          c3 = c2 * dsqrt(mass(j)/betan)
          do i = 1, 3
             p(i,j,k) = c1 * p(i,j,k) + gaussian(irun,c3)
          enddo
       enddo
    enddo
    
    return
  end subroutine parinello_therm_loc

  subroutine parinello_therm_acmd(p,mass,tau,delp,na,nb,dt,irun)
    implicit none
    ! ------------------------------------------------------------------
    ! Parinello Local Langevin Thermostat for use in ACMD simulations
    ! ------------------------------------------------------------------
    integer na,nb,irun,k,j,i,ic
    real(8) p(3,na,nb),mass(na),delp(3,nb),tau,gamma,dt
    real(8) c1,c2,c3,gaussian
    external gaussian
    
    gamma = 1.d0/(2.d0*tau)
    
    c1 = dexp (-0.5d0 * gamma * dt)
    c2 = dsqrt(1.d0 - c1*c1)
    
    do k = 1,nb
       do j = 1, na
          ic = mod(j,3)
          if (ic.eq.0) then
             ic = 3
          endif
          c3 = c2 * delp(ic,k)
          do i = 1, 3
             p(i,j,k) = c1 * p(i,j,k) + gaussian(irun,c3)
          enddo
       enddo
    enddo
    
    return
  end subroutine parinello_therm_acmd

  subroutine sample(p,na,nb,mass,beta,irun,dt)
    implicit none
    ! ------------------------------------------------------------------
    ! Sampling of momenta
    ! ------------------------------------------------------------------
    integer na,nb,irun
    real(8) p(3,na,nb),mass(na),beta,dt,om
    character*4 type
    common /path_i/ om,type
    
    if ((type.eq.'RPMD').or.(nb.eq.1)) then
       call sample_rpmd (p,na,nb,mass,beta,irun)
    else if (type.eq.'ACMD') then
       call sample_acmd (p,na,nb,mass,beta,irun,dt,om)
    endif
    
    return
  end subroutine sample


  subroutine sample_acmd (p,na,nb,mass,beta,irun,dt,om)
    implicit none
    include 'globals.inc'
    ! ------------------------------------------------------------------
    ! Sample momenta for partially adiabatic CMD calculation.
    ! Note that there are some extra factors of two included which
    ! are neccessary to account for normalization in the FT routines.
    ! ------------------------------------------------------------------
    integer na,nb,j,k,l,irun
    real(8) p(3,na,nb),mass(na)
    real(8) dt,beta,wk,om,betan,gaussian,pibyn,twown,emk
    real(8), allocatable :: deltap(:,:)
    external gaussian

    ! Sample from fourier-transformed momenta
    
    betan = beta / dble(nb)
    twown = 2.d0 / (betan * hbar)
    pibyn = dacos(-1.d0) / dble(nb)
    
    allocate (deltap(na,nb))
    
    do j = 1, na
       deltap(j,1) = dsqrt(mass(j)/betan)
       if (nb.gt.1) then
          do k = 1, nb/2
             wk = twown * dsin(k * pibyn)
             emk = mass(j)*(wk/om)**2
             deltap(j,k+1) = dsqrt(0.5d0 * emk/betan)
          enddo
          if(2*(nb/2).eq.nb) then
             deltap(j,(nb/2)+1) = dsqrt(2.d0)*deltap(j,(nb/2)+1)
          endif
          do k = 1, (nb-1)/2
             deltap(j,nb-k+1) = deltap(j,k+1)
          enddo
       endif
    enddo
    
    do k = 1, nb
       do j = 1, na
          do l = 1,3
             p(l,j,k) = gaussian(irun,deltap(j,k))
          enddo
       enddo
    enddo
    
    ! Fourier transform back to real-space momenta
    
    if (nb.gt.1) then
       call realft(p,3*na,nb,-1)
    endif
    
    ! Shift the momenta such that the centre-of-mass momentum
    ! is zero.
    
    call pshift(p,mass,na,nb)
    
    deallocate (deltap)
    
    return
  end subroutine sample_acmd
  

  subroutine sample_rpmd (p,na,nb,mass,beta,irun)
    implicit none
    include 'globals.inc'
    ! ------------------------------------------------------------------
    ! Sample momenta from Maxwell distribution at inverse temperatures
    ! beta / nb, as required for an RPMD simulation.
    ! ------------------------------------------------------------------
    integer na,nb,irun,k,j,l
    real(8) mass(na),p(3,na,nb)
    real(8) beta,sigma,betan,gaussian
    external gaussian
    
    betan = beta/dble(nb)
    do k = 1,nb
       do j = 1,na
          do l = 1, 3
             sigma = dsqrt(mass(j)/betan)
             p(l,j,k) = gaussian(irun,sigma)
          enddo
       enddo
    enddo
    
    ! Shift the momenta so that the center-of-mass momentum is zero.
    
    call pshift (p,mass,na,nb)

    return
  end subroutine sample_rpmd

  
  subroutine pshift (p,mass,na,nb)
    implicit none
    ! ------------------------------------------------------------------
    ! Sets the centre-of-mass momentum of n particles to zero.
    ! ------------------------------------------------------------------
    integer j,k,na,nb
    real(8) p(3,na,nb),vx,vy,vz,mass(na),totmass
    
    vx = 0.d0
    vy = 0.d0
    vz = 0.d0
    totmass = 0.d0
    do j = 1,na
       totmass = totmass + mass(j)
    enddo
    totmass = dble(nb)*totmass
    do k = 1,nb
       do j = 1,na
          vx = vx + p(1,j,k)
          vy = vy + p(2,j,k)
          vz = vz + p(3,j,k)
       enddo
    enddo
    vx = vx/totmass
    vy = vy/totmass
    vz = vz/totmass
    do k = 1,nb
       do j = 1,na
          p(1,j,k) = p(1,j,k)-mass(j)*vx
          p(2,j,k) = p(2,j,k)-mass(j)*vy
          p(3,j,k) = p(3,j,k)-mass(j)*vz
       enddo
    enddo
    
    return
  end subroutine pshift
  
end module thermostat

