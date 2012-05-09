subroutine virial_ke(r,dvdr,dvdr2,tv,tvxyz,tq1,tq2,beta,na,nb)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Virial estimate of the total kinetic energy.
  ! ------------------------------------------------------------------
  integer na,nb,nc,nctot,nbond
  integer i,j,k
  real(8) tv, tvxyz(3),beta
  real(8) dr,tq1,tq2,tc,drdvr,drdvr2,tq1xyz(3),tq2xyz(3)
  real(8) r(3,na,nb),dvdr(3,na,nb),dvdr2(3,na,nb)
  real(8), allocatable :: rc(:,:)
  common /constraint/ nctot,nbond

  allocate (rc(3,na))
  
  rc(:,:) = 0.d0
  tq1 = 0.d0
  tq2 = 0.d0
  tq1xyz(:) = 0.d0
  tq2xyz(:) = 0.d0
  
  do k = 1,nb
     rc(:,:) = rc(:,:) + r(:,:,k)
  enddo
  rc(:,:) = rc(:,:)/dble(nb)
  
  do k = 1,nb
     do j = 1,na
        do i = 1,3
           dr = r(i,j,k)-rc(i,j)
           drdvr = dr*dvdr(i,j,k)
           drdvr2 = dr*dvdr2(i,j,k)
           tq1 = tq1 + drdvr
           tq2 = tq2 + drdvr2
           tq1xyz(i) = tq1xyz(i) + drdvr
           tq2xyz(i) = tq2xyz(i) + drdvr2
        enddo
     enddo
  enddo
      
  nc = 3                    ! NUMBER OF CONSTRAINTS !
  tq1 = 0.5d0*tq1/dble(nb)
  tq2 = 0.5d0*tq2/dble(nb)
  tq1xyz(:) = 0.5d0*tq1xyz(:)/dble(nb)
  tq2xyz(:) = 0.5d0*tq2xyz(:)/dble(nb)
  tc = 0.5d0*(3*na-nc-nctot)/beta
  tv = tc + tq1 + tq2
  tvxyz(:) = tc/3 + tq1xyz(:) + tq2xyz(:)
  
  deallocate (rc)

  return
end subroutine virial_ke

!subroutine virial_ke(r,dvdr,dvdr2,tv,tvxyz,tq1,tq2,beta,na,nb,mass)
!  implicit none
!  include 'globals.inc'
!  ! ------------------------------------------------------------------
!  ! Virial estimate of the total kinetic energy.
!  ! ------------------------------------------------------------------
!  integer na,nb,nc,nctot,nbond
!  integer i,j,k,nm,j1,j2,j3,ic
!  real(8) tv, tvxyz(3),beta,tq
!  real(8) dr,tq1,tq2,tc,drdvr,drdvr2,tq1xyz(3),tq2xyz(3)
!  real(8) r(3,na,nb),dvdr(3,na,nb),dvdr2(3,na,nb),mass(na),totmass,sum
!  real(8), allocatable :: rc(:,:),rcom(:,:,:),dvdrcom(:,:,:)
!  common /constraint/ nctot,nbond
!
!  nm = na/3
!
!  allocate (rc(3,nm))
!
!  rc(:,:) = 0.d0
!  tq1 = 0.d0
!  tq2 = 0.d0
!  tq1xyz(:) = 0.d0
!  tq2xyz(:) = 0.d0
!
!  totmass = mass(1) + mass(2) + mass(3)
!
!  ! For bead of each molecule, get the COM coordinates.
!  allocate( rcom(3,nm,nb) )
!  allocate( dvdrcom(3,nm,nb) )
!  dvdrcom(:,:,:) = 0.0 
!  ic = 0
!  do j = 1, na, 3
!     ic = ic + 1
!     do k = 1, nb
!        rcom(1,ic,k) = (mass(j)*r(1,j,k)+mass(j+1)*r(1,j+1,k)+mass(j+2)*r(1,j+2,k))/totmass 
!        rcom(2,ic,k) = (mass(j+1)*r(2,j,k)+mass(j+1)*r(2,j+1,k)+mass(j+2)*r(2,j+2,k))/totmass 
!        rcom(3,ic,k) = (mass(j+2)*r(3,j,k)+mass(j+1)*r(3,j+1,k)+mass(j+2)*r(3,j+2,k))/totmass
!
!        dvdrcom(1,ic,k) = dvdr(1,j,k) + dvdr(1,j+1,k) + dvdr(1,j+2,k) 
!        dvdrcom(2,ic,k) = dvdr(2,j,k)+dvdr(2,j+1,k)+dvdr(2,j+2,k)
!        dvdrcom(3,ic,k) = dvdr(3,j,k)+dvdr(3,j+1,k)+dvdr(3,j+2,k)
!
!        dvdrcom(1,ic,k) = dvdrcom(1,ic,k) + dvdr2(1,j,k)+dvdr2(1,j+1,k)+dvdr2(1,j+2,k)
!        dvdrcom(2,ic,k) = dvdrcom(2,ic,k) + dvdr2(2,j,k)+dvdr2(2,j+1,k)+dvdr2(2,j+2,k)
!        dvdrcom(3,ic,k) = dvdrcom(3,ic,k) + dvdr2(3,j,k)+dvdr2(3,j+1,k)+dvdr2(3,j+2,k)
!
!     enddo
!  enddo
!
!  ! Get the molecular centroid coordinates:
!  rc(:,:) = 0.0
!  do i = 1, nm
!     do k = 1, nb  
!        rc(1,i) = rc(1,i) + rcom(1,i,k)
!        rc(2,i) = rc(2,i) + rcom(2,i,k)
!        rc(3,i) = rc(3,i) + rcom(3,i,k)
!     enddo
!  enddo
!  rc(:,:) = rc(:,:) / dble(nb)
!
!  ! Sum the virial contributions.
!  do j = 1, nm
!     do k = 1,nb
!        do i = 1,3
!           tq = tq + (rcom(i,j,k)-rc(i,j))*dvdrcom(i,j,k)
!        enddo
!     enddo
!  enddo
!
!  nc = 3                    ! NUMBER OF CONSTRAINTS !
!  tq = 0.5d0*tq/dble(nb)
!  tc = 0.5d0*(3*nm)/beta
!  tv = tc + tq
!  tvxyz(:) = tc/3 + tq1xyz(:) + tq2xyz(:)
!
!  deallocate (rc)
!  deallocate( rcom )
!  deallocate( dvdrcom )
!
!  return
!end subroutine virial_ke


subroutine virial_ke_ee(r,dvdr,tv,tq,beta,na,nb)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Virial estimate of the total kinetic energy for exact estimator.
  ! ------------------------------------------------------------------
  integer na,nb,nc,nctot,nbond
  integer j,k,i
  real(8) r(3,na,nb),dvdr(3,na,nb)
  real(8) tv, beta,dr,tq,tc
  real(8), allocatable :: rc(:,:)
  common /constraint/ nctot,nbond

  allocate (rc(3,na))

  rc(:,:) = 0.d0
  tq = 0.d0

  do k = 1,nb
     rc(:,:) = rc(:,:)+r(:,:,k)
  enddo
  rc(:,:) = rc(:,:)/dble(nb)
  
  do k = 1,nb
     do j = 1,na
        do i = 1,3
           dr = r(i,j,k)-rc(i,j)
           tq = tq + dr*dvdr(i,j,k)
        enddo
     enddo
  enddo
  nc = 3                    ! NUMBER OF CONSTRAINTS !
  tq = 0.5d0*tq/dble(nb)
  tc = 0.5d0*(3*na-nc-nctot)/beta
  tv = tc + tq

  deallocate (rc)

  return
end subroutine virial_ke_ee

subroutine kinetic(p,na,tk,mass,nb,beta)
  implicit none
  ! ------------------------------------------------------------------
  ! Kinetic energy evaluation
  ! ------------------------------------------------------------------
  integer na,nb
  real(8) p(3,na,nb),mass(na),tk,beta,om
  character(len=4) type
  common /path_i/ om,type

  if ((type.eq.'RPMD').or.(nb.eq.1)) then
     call kinetic_rpmd(p,na,tk,mass,nb)
  else if (type.eq.'ACMD') then
     call kinetic_acmd(p,na,tk,mass,nb,beta,om)
  endif

  return
end subroutine kinetic

subroutine kinetic_acmd(p,na,tk,rmass,nb,beta,om)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Calculates the kinetic energy in a PACMD simulation.
  ! Note that this requires Fourier transformation of the momenta,
  ! and we also need to properly account for the normalization.
  ! ------------------------------------------------------------------
  integer na,nb,j,k,k1
  real(8) p(3,na,nb),rmass(na)
  real(8) om,wn,tk,pp,wk,rk,em,emk,pi,beta
  real(8) , allocatable :: ptemp(:,:,:)

  allocate (ptemp(3,na,nb))

  wn = dble(nb) / (beta * hbar)
  pi = dacos(-1.d0)
  ptemp(:,:,:) = p(:,:,:)
  call realft (ptemp,3*na,nb,+1)

  tk = 0.d0
  do j = 1, na
     do k = 1, nb
        k1 = k-1
        if (k1.eq.0) then
           em = rmass(j)
           pp = ptemp(1,j,k)**2+ptemp(2,j,k)**2+ptemp(3,j,k)**2
           tk = tk + pp / (2.d0*em)
        else if (k1.eq.(nb/2)) then
           em = rmass(j)
           wk = 2.d0 * wn * dsin(pi/2.d0)
           emk = em * (wk**2/OM**2)
           pp = ptemp(1,j,k)**2+ptemp(2,j,k)**2+ptemp(3,j,k)**2
           tk = tk + pp / (2.d0*emk)
        else
           em = rmass(j)
           rk = dble(k1)
           wk = 2.d0 * wn * dsin(rk*pi/dble(nb))
           emk = em * (wk**2/OM**2)
           pp = ptemp(1,j,k)**2+ptemp(2,j,k)**2+ptemp(3,j,k)**2
           tk = tk+ pp / emk
        endif
     enddo
  enddo
  
  return
end subroutine kinetic_acmd


subroutine kinetic_rpmd(p,na,tk,mass,nb)
  implicit none
  ! ------------------------------------------------------------------
  ! Calculates the kinetic energy in an RPMD simulation.
  ! ------------------------------------------------------------------
  integer na, i, nb,j
  real(8) p(3,na,nb), mass(na), tk

  tk = 0.d0
  do j = 1, nb
     do i = 1, na
        tk = tk + (1.d0/(2.d0*mass(i))) * &
            (p(1,i,j)**2 + p(2,i,j)**2 +p(3,i,j)**2)
     enddo
  enddo

  return
end subroutine kinetic_rpmd
