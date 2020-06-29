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
  real(8), allocatable :: dvdrb_split(:,:,:,:) ! dummy variable for forces

  ! nbdf = RP Contraction (beadibatic) factor
  ! (determines the number of normal modes to keep)

  ! Check that beadiabatic factor is odd

  if (mod(nbdf-1,2).ne.0) then
     write(6,*) '* Beadiabatic factor must be ODD'
     stop
  endif

  ! Allocate

  allocate (rb(3,na,nbdf),dvdrb(3,na,nbdf),dvdrb_split(3,na,nbdf,4))

  ! Zero arrays

  dvdr(:,:,:) = 0.d0
  rb(:,:,:) = 0.d0
  dvdrb(:,:,:) = 0.d0
  dvdrb_split(:,:,:,:) = 0.d0
  vir(:,:) = 0.d0

  ! Form contracted RP coordinates

  call ring_contract(r,rb,na,nb,nbdf)

  ! Evaluate forces on the bead coordinates rb

  call forces(rb,v,dvdrb,dvdrb_split,nbdf,na,boxlxyz,z,vir,iopt)

  ! Create full bead forces from contracted coordinates

  call force_contract(dvdr,dvdrb,na,nb,nbdf)

  deallocate (rb,dvdrb,dvdrb_split)

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
