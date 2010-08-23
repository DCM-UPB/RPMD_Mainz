subroutine print_structure(r,boxlxyz,nm,na,nb,iunit)
  implicit none
  ! ------------------------------------------------------------------
  ! Subroutine to print the current structure for use
  ! in restarts
  ! ------------------------------------------------------------------
  integer na,nm,nb,iunit
  real(8) r(3,na,nb),boxlxyz(3)

  ! Print the structure.

  write(iunit,*) nm,na,nb
  write(iunit,*) boxlxyz(1),boxlxyz(2),boxlxyz(3)
  write(iunit,*) r(:,:,:)

  return
end subroutine print_structure

subroutine print_rdf(ihoo,ihoh,ihhh,na,boxlxyz,ng,nb,nrdf)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Outputs RDF
  ! ------------------------------------------------------------------
  integer ihoo(imaxbin),ihoh(imaxbin)
  integer ihhh(imaxbin),na,ng,nb,nrdf
  integer no,nh,i
  real(8) delr,rlower,rupper,rh,vol,boxlxyz(3)
  real(8) const,rideal,gr,pi,boxmax

  pi = dacos(-1.d0)
  vol = boxlxyz(1)*boxlxyz(2)*boxlxyz(3)

  ! Open RDF output files.

  open (unit=50,file='g_oo.out')
  open (unit=51,file='g_oh.out')
  open (unit=52,file='g_hh.out')

  boxmax = max(boxlxyz(1),boxlxyz(2),boxlxyz(3))

  delr = dble( 0.5d0 * boxmax / imaxbin )
  no = na / 3
  nh = 2 * na / 3

  ! Print the RDFs

  do i = 1, imaxbin
     rlower = real( i - 1 ) * delr
     rupper = rlower + delr

     ! oxygen - oxygen distribution

     rh = dble(no) / vol
     const = 4.d0 * pi * rh / 3.d0
     rideal = const * ( rupper**3 - rlower**3)
     gr = real(ihoo(i)) / ( real( no * nb * nrdf) * rideal )
     write(50,*)(rlower+ delr / 2.d0)*ToA, gr

     ! oxygen - hydrogen distribution

     rh = dble(no) / vol
     const = 4.d0 * pi * rh / 3.d0
     rideal = const * ( rupper**3 - rlower**3)
     gr = real(ihoh(i)) / (real( nh*nb * nrdf) * rideal)
     write(51,*)(rlower + delr / 2.d0)*ToA, gr

     ! hydrogen - hydrogen distribution

     rh = dble(nh) / vol
     const = 4.d0 * pi * rh / 3.d0
     rideal = const * ( rupper**3 - rlower**3)
     gr = real(ihhh(i)) / (real( nh*nb * nrdf) * rideal)
     write(52,*) (rlower + delr / 2.d0)*ToA, gr
  enddo
  close (unit=50)
  close (unit=51)
  close (unit=52)
  return
end subroutine print_rdf

subroutine print_cvv(ct,nt,dtfs,dt)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Output the velocity autocorrelation function
  ! in (A/ps)**2 vs t in ps.
  ! Also calculates and prints the diffusion constant.
  !------------------------------------------------------------------
  integer nt,it
  real(8) ct(0:nt),cscale,dtfs,dtps,tops
  real(8) dc,tps,dt
  
  dtps = dtfs*1.d-3
  tops = dtps/dt
  open (unit=10,file='Cvv.out')
  cscale = (toA/tops)**2
  do it = 0,nt
     tps = it*dtps
     ct(it) = cscale*ct(it)
     write (10,11) tps,ct(it)
11   format(1x,2f16.8)
  enddo
  close (unit=10)
  
  dc = 0.5d0*ct(0)
  do it = 1,nt
     dc = dc+ct(it)
  enddo
  dc = dtps*(dc/3.d0)
  write (6,64) dc
64 format(/1x,'Diffusion constant  = ',f10.4,' A**2/ps')
  write(6,*)
  return
end subroutine print_cvv

subroutine print_dipole(dtfs,nt,rmtr,rmtv,dqq,dqx, &
                        dqy,dqz,boxlxyz,beta,dt)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Outputs the dipole and dipole derivative
  ! autocorrelation functions.
  ! ------------------------------------------------------------------
  integer nt,it
  real(8) rmtr(0:nt),rmtv(0:nt),dtfs,dtps,fac
  real(8) boxlxyz(3),beta,m2,vol,eps,dqq,dqx,dqy
  real(8) dqz,dt,tps,pi

  pi = dacos(-1.d0)
  open (unit=11,file='Dmmr.out')
  open (unit=13,file='Dmmv.out')
  dtps = dtfs*1.d-3
  do it = 0,nt,iskip
     tps = it*dtps
     write (11,*)tps, rmtr(it)
     write (13,*)tps, rmtv(it)
  enddo
  close (unit=11)
  close (unit=13)

  ! Print the calculated properties.

  write(6,*)'* Electric properties '
  write(6,*)
  write (6,*)'<Dipole**2> = ', dqq,' e**2 bohr**2 '
  write (6,*)'<Dipole(x)> = ', dqx,' e bohr '
  write (6,*)'<Dipole(y)> = ', dqy,' e bohr '
  write (6,*)'<Dipole(z)> = ', dqz,' e bohr '

  ! Dielectric constant :

  ! eps = 1 + ( (4*pi*beta)/(3*Vol) ) * (<M**2> - <M>**2)
  
  vol = boxlxyz(1)*boxlxyz(2)*boxlxyz(3)
  m2 = dqx**2 + dqy**2 + dqz**2
  fac = (4.d0*pi*beta) / (3.d0 * vol)
  eps = 1.d0 + fac * (dqq - m2)
  write(6,*)'Dielectric constant = ',eps
  write(6,*)
  
  return
end subroutine print_dipole

subroutine print_orientation(dtfs,nt,dmo1)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Outputs the orientational autocorrelation functions.
  ! ------------------------------------------------------------------
  integer nt,it
  real(8) dtfs,dmo1(6,0:nt),dtps,tps

  open (unit=14,file='Cori1.out')
  open (unit=15,file='Cori2.out')
  dtps = dtfs*1.d-3

  ! Note that we skip every iskip steps here.

  do it = 0,nt,iskip
     tps = it*dtps
     write (14,62)tps,dmo1(1,it),dmo1(2,it),dmo1(3,it)
     write (15,62)tps, 1.5d0*dmo1(4,it)-0.5d0, &
          1.5d0*dmo1(5,it)-0.5d0,1.5d0*dmo1(6,it)-0.5d0
  enddo
62 format(4f12.6)
  close (unit=14)
  close (unit=15)
end subroutine print_orientation

subroutine simp(n,h,fi,s)
  implicit none
  ! ------------------------------------------------------------------
  ! Subroutine for integration over f(x) with the Simpson rule.
  ! FI: integrand f(x); H: interval; S: integral.
  ! ------------------------------------------------------------------
  integer n,i
  real(8) fi(n),h,s,s0,s1,s2

  s  = 0.d0
  s0 = 0.d0
  s1 = 0.d0
  s2 = 0.d0
  do  i = 2, n-1, 2
     s1 = s1 + fi(i-1)
     s0 = s0 + fi(i)
     s2 = s2 + fi(i+1)
  enddo
  s = h*(s1 + 4.d0*s0 + s2)/3.d0
    
  ! If N is even, add the last slice separately
  
  if (mod(n,2).eq.0) s = s &
       + h*(5.d0*fi(n) + 8.d0*fi(n-1) - fi(n-2))/12.d0
  return
end subroutine simp

subroutine print_vmd_full(r,nb,na,nm,boxlxyz,nunit)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! VMD output
  ! ------------------------------------------------------------------
  integer na,nb,nm,i,j,k,ni,nunit
  real(8) r(3,na,nb),boxlxyz(3)
  real(8) xbox,ybox,zbox,onboxlx,onboxly,onboxlz
  real(8) dh1x,dh2x,dh1y,dh2y,dh1z,dh2z
  real(8), allocatable :: rc(:,:)

  logical use_traj
  common /reftraj/ use_traj

  allocate (rc(3,na))
  rc(:,:) = 0.d0

  ! Useful parameters

  xbox = boxlxyz(1)
  ybox = boxlxyz(2)
  zbox = boxlxyz(3)
  onboxlx = 1.d0 / xbox
  onboxly = 1.d0 / ybox
  onboxlz = 1.d0 / zbox

  ! Positions of centroids

  do i = 1,na
     do j = 1,3
        rc(j,i) = 0.d0
        do k = 1,nb
           rc(j,i) = rc(j,i) + r(j,i,k)
        enddo
        rc(j,i) = rc(j,i)/dble(nb)
     enddo
  enddo

  ! Apply Boundary Conditions to oxygen positions:
  ! when not in reftraj mode

  if (use_traj.eqv..false.) then
    do i = 0,nm-1
      ! Vector of Hydrogen relative to Oxygen:
      ni = 3*i+1
      dh1x = rc(1,ni+1) - rc(1,ni)
      dh2x = rc(1,ni+2) - rc(1,ni)
      dh1y = rc(2,ni+1) - rc(2,ni)
      dh2y = rc(2,ni+2) - rc(2,ni)
      dh1z = rc(3,ni+1) - rc(3,ni)
      dh2z = rc(3,ni+2) - rc(3,ni)
      ! Shift Oxygen atoms into box:
      rc(1,ni) = rc(1,ni) - xbox*(nint(rc(1,ni)*onboxlx)-0.5d0)
      rc(2,ni) = rc(2,ni) - ybox*(nint(rc(2,ni)*onboxly)-0.5d0)
      rc(3,ni) = rc(3,ni) - zbox*(nint(rc(3,ni)*onboxlz)-0.5d0)
      ! Move Hydrogens relative to Oxygen:
      rc(1,ni+1) = dh1x + rc(1,ni)
      rc(1,ni+2) = dh2x + rc(1,ni)
      rc(2,ni+1) = dh1y + rc(2,ni)
      rc(2,ni+2) = dh2y + rc(2,ni)
      rc(3,ni+1) = dh1z + rc(3,ni)
      rc(3,ni+2) = dh2z + rc(3,ni)
    enddo
  endif
  write(nunit,*) na
  write(nunit,*) "BOX ", toA*xbox, toA*ybox, toA*zbox          !!GLE: PRINT ALSO BOX PARS, THIS IS A GOOD PLACE!!!
  do k = 1,na,3
     write(nunit,*) 'O ', toA*rc(1,k), toA*rc(2,k),toA*rc(3,k)
     write(nunit,*) 'H ', toA*rc(1,k+1), toA*rc(2,k+1),toA*rc(3,k+1)
     write(nunit,*) 'H ', toA*rc(1,k+2), toA*rc(2,k+2),toA*rc(3,k+2)
  enddo
  deallocate (rc)
  return
end subroutine print_vmd_full

subroutine print_vmd_full_forces(dvdr,dvdr2,nb,na,boxlxyz,nunit)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! VMD output
  ! ------------------------------------------------------------------
  integer na,nb,i,j,k,ni,nunit
  real(8) dvdr(3,na,nb),dvdr2(3,na,nb),boxlxyz(3)
  real(8) xbox,ybox,zbox
  real(8), allocatable :: frc(:,:)

  allocate (frc(3,na))
  frc(:,:) = 0.d0

  ! Useful parameters

  xbox = boxlxyz(1)
  ybox = boxlxyz(2)
  zbox = boxlxyz(3)

  ! Forces on centroids

  do i = 1,na
     do j = 1,3
        frc(j,i) = 0.d0
        do k = 1,nb
           frc(j,i) = frc(j,i) + dvdr(j,i,k) + dvdr2(j,i,k)
        enddo
        frc(j,i) = frc(j,i)/dble(nb)     
     enddo
  enddo

  write(nunit,*) na
  write(nunit,*) "BOX ", toA*xbox, toA*ybox, toA*zbox          !!GLE: PRINT ALSO BOX PARS, THIS IS A GOOD PLACE!!!
  do k = 1,na,3
     write(nunit,*) 'O ', frc(1,k+0), frc(2,k+0),frc(3,k+0)
     write(nunit,*) 'H ', frc(1,k+1), frc(2,k+1),frc(3,k+1)
     write(nunit,*) 'H ', frc(1,k+2), frc(2,k+2),frc(3,k+2)
  enddo
  deallocate (frc)
  return
end subroutine print_vmd_full_forces

subroutine print_vmd_full_vels(p,mass,nb,na,boxlxyz,nunit)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! VMD output
  ! ------------------------------------------------------------------
  integer na,nb,i,j,k,ni,nunit
  real(8) p(3,na,nb),mass(na),boxlxyz(3)
  real(8) xbox,ybox,zbox
  real(8), allocatable :: vel(:,:)

  allocate (vel(3,na))
  vel(:,:) = 0.d0

  ! Useful parameters

  xbox = boxlxyz(1)
  ybox = boxlxyz(2)
  zbox = boxlxyz(3)

  ! Velocities on centroids

  do i = 1,na  
     do j = 1,3
        vel(j,i) = 0.d0
        do k = 1,nb
           vel(j,i) = vel(j,i) + p(j,i,k)/mass(i)
        enddo
        vel(j,i) = vel(j,i)/dble(nb)
     enddo
  enddo

  write(nunit,*) na
  write(nunit,*) "BOX ", toA*xbox, toA*ybox, toA*zbox          !!GLE: PRINT ALSO BOX PARS, THIS IS A GOOD PLACE!!!
  do k = 1,na,3
     write(nunit,*) 'O ', vel(1,k+0), vel(2,k+0),vel(3,k+0)
     write(nunit,*) 'H ', vel(1,k+1), vel(2,k+1),vel(3,k+1)
     write(nunit,*) 'H ', vel(1,k+2), vel(2,k+2),vel(3,k+2)
  enddo
  deallocate (vel)
  return
end subroutine print_vmd_full_vels

subroutine print_vmd_bead(r,nb,ib,na,nm,boxlxyz,nunit)
   implicit none
   include 'globals.inc'
   ! ------------------------------------------------------------------
   ! VMD output
   ! ------------------------------------------------------------------
   integer na,nb,ib,nm,i,j,k,ni,nunit
   real(8) r(3,na,nb),boxlxyz(3)
   real(8) xbox,ybox,zbox,onboxlx,onboxly,onboxlz
   real(8) dh1x,dh2x,dh1y,dh2y,dh1z,dh2z
   real(8), allocatable :: rb(:,:)

   logical use_traj
   common /reftraj/ use_traj

   allocate (rb(3,na))
   rb(:,:) = 0.d0

   ! Useful parameters

   xbox = boxlxyz(1)
   ybox = boxlxyz(2)
   zbox = boxlxyz(3)
   onboxlx = 1.d0 / xbox
   onboxly = 1.d0 / ybox
   onboxlz = 1.d0 / zbox

   do i = 1,na
      do j = 1,3
         rb(j,i) = 0.d0
         rb(j,i) = r(j,i,ib)
      enddo
   enddo

  ! Apply Boundary Conditions to oxygen positions:
  ! when not in reftraj mode

  if (use_traj.eqv..false.) then
    do i = 0,nm-1
      ! Vector of Hydrogen relative to Oxygen:
      ni = 3*i+1
      dh1x = rb(1,ni+1) - rb(1,ni)
      dh2x = rb(1,ni+2) - rb(1,ni)
      dh1y = rb(2,ni+1) - rb(2,ni)
      dh2y = rb(2,ni+2) - rb(2,ni)
      dh1z = rb(3,ni+1) - rb(3,ni)
      dh2z = rb(3,ni+2) - rb(3,ni)
      ! Shift Oxygen atoms into box:
      rb(1,ni) = rb(1,ni) - xbox*(nint(rb(1,ni)*onboxlx)-0.5d0)
      rb(2,ni) = rb(2,ni) - ybox*(nint(rb(2,ni)*onboxly)-0.5d0)
      rb(3,ni) = rb(3,ni) - zbox*(nint(rb(3,ni)*onboxlz)-0.5d0)
      ! Move Hydrogens relative to Oxygen:
      rb(1,ni+1) = dh1x + rb(1,ni)
      rb(1,ni+2) = dh2x + rb(1,ni)
      rb(2,ni+1) = dh1y + rb(2,ni)
      rb(2,ni+2) = dh2y + rb(2,ni)
      rb(3,ni+1) = dh1z + rb(3,ni)
      rb(3,ni+2) = dh2z + rb(3,ni)
    enddo
  endif
  write(nunit,*) na
  write(nunit,*) "BOX ", toA*xbox, toA*ybox, toA*zbox          !!GLE: PRINT ALSO BOX PARS, THIS IS A GOOD PLACE!!!
  do k = 1,na,3
     write(nunit,*) 'O ', toA*rb(1,k+0), toA*rb(2,k),toA*rb(3,k)
     write(nunit,*) 'H ', toA*rb(1,k+1), toA*rb(2,k+1),toA*rb(3,k+1)
     write(nunit,*) 'H ', toA*rb(1,k+2), toA*rb(2,k+2),toA*rb(3,k+2)
  enddo
  deallocate (rb)
  return
end subroutine print_vmd_bead

subroutine print_vmd_bead_forces(dvdr,dvdr2,nb,ib,na,boxlxyz,nunit)
   implicit none
   include 'globals.inc'
   ! ------------------------------------------------------------------
   ! VMD output
   ! ------------------------------------------------------------------
   integer na,nb,ib,i,j,k,ni,nunit
   real(8) dvdr(3,na,nb),dvdr2(3,na,nb),boxlxyz(3)
   real(8) xbox,ybox,zbox
   real(8), allocatable :: frc(:,:)

   allocate (frc(3,na))
   frc(:,:) = 0.d0

   ! Useful parameters

   xbox = boxlxyz(1)
   ybox = boxlxyz(2)
   zbox = boxlxyz(3)

   do i = 1,na
      do j = 1,3
         frc(j,i) = 0.d0
         frc(j,i) = dvdr(j,i,ib) + dvdr2(j,i,ib)
      enddo
   enddo

  write(nunit,*) na
  write(nunit,*) "BOX ", toA*xbox, toA*ybox, toA*zbox          !!GLE: PRINT ALSO BOX PARS, THIS IS A GOOD PLACE!!!
  do k = 1,na,3
     write(nunit,*) 'O ', frc(1,k+0), frc(2,k+0),frc(3,k+0)
     write(nunit,*) 'H ', frc(1,k+1), frc(2,k+1),frc(3,k+1)
     write(nunit,*) 'H ', frc(1,k+2), frc(2,k+2),frc(3,k+2)
  enddo
  deallocate (frc)
  return
end subroutine print_vmd_bead_forces

subroutine print_vmd_bead_vels(p,mass,nb,ib,na,boxlxyz,nunit)
   implicit none
   include 'globals.inc'
   ! ------------------------------------------------------------------
   ! VMD output
   ! ------------------------------------------------------------------
   integer na,nb,ib,i,j,k,ni,nunit
   real(8) p(3,na,nb),mass(na),boxlxyz(3)
   real(8) xbox,ybox,zbox
   real(8), allocatable :: vel(:,:)

   allocate (vel(3,na))
   vel(:,:) = 0.d0

   ! Useful parameters

   xbox = boxlxyz(1)
   ybox = boxlxyz(2)
   zbox = boxlxyz(3)

   do i = 1,na
      do j = 1,3
         vel(j,i) = 0.d0
         vel(j,i) = p(j,i,ib)/mass(i)
      enddo
   enddo

  write(nunit,*) na
  write(nunit,*) "BOX ", toA*xbox, toA*ybox, toA*zbox          !!GLE: PRINT ALSO BOX PARS, THIS IS A GOOD PLACE!!!
  do k = 1,na,3
     write(nunit,*) 'O ', vel(1,k+0), vel(2,k+0),vel(3,k+0)
     write(nunit,*) 'H ', vel(1,k+1), vel(2,k+1),vel(3,k+1)
     write(nunit,*) 'H ', vel(1,k+2), vel(2,k+2),vel(3,k+2)
  enddo
  deallocate (vel)
  return
end subroutine print_vmd_bead_vels
