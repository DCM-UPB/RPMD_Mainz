subroutine setup_positions(r,na,nb,nm,boxlxyz,qo,irun)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Performs setup of staring lattice - ice, fcc or bcc
  !
  ! If narg = 1, we start from a cubic lattice,
  ! if narg = 3, we read the coordinates from file.
  ! ------------------------------------------------------------------
  integer na,nb,ncellxyz(3),nm,i,k,nlat,irun,it,itmax
  real(8) r(3,na,nb),boxlxyz(3)
  real(8) dip,dipx,dipy,dipz,dip2,dipm,qh,qo
  real(8), allocatable :: z_tmp(:)
  common /lattice/ ncellxyz,nlat

  if (3*nm.ne.na) then
     write(6,*)'* Error in setup',nm,na
     stop
  endif

  ! Setup a new lattice of molecules
  ! ---------------------------------

  if (nlat.eq.1) then
     
     ! Cubic Lattice
 
     call cubic (r(1,1,1),nm,boxlxyz,ncellxyz)
  else if (nlat.eq.2) then

     ! Face Centred Cubic Lattice

     call fcc(r(1,1,1),nm,boxlxyz,ncellxyz)
  else if (nlat.eq.3) then
     
     ! Ice lattice (low dipole)

     allocate (z_tmp(na))
     z_tmp(:) = 0.d0
     qh = -0.5d0*qo

     do i = 1,na,3
        z_tmp(i) = qo
        z_tmp(i+1) = qh
        z_tmp(i+2) = qh
     enddo
     
     dip = 100.d0
     it = 0
     itmax = 50000
     do while ((dip.gt.2.3d0).and.(it.le.itmax))
        it = it + 1
        call ice_lattice(r(1,1,1),boxlxyz,ncellxyz,nm,irun,0)
        do k = 2,nb
           do i = 1,na
              r(1,i,k) = r(1,i,1)
              r(2,i,k) = r(2,i,1)
              r(3,i,k) = r(3,i,1)
           enddo
        enddo
        call dipole(r,dipx,dipy,dipz,dip2,dipm,z_tmp,na,nb)
        dip = dsqrt(dipx**2+dipy**2+dipz**2)/0.393456d0
     enddo
     if (it.ge.itmax) then
        write(6,*) 'setup_positions : TOO MANY INTERATIONS'
     endif
     write(6,*)
     write(6,*) 'Hexagonal Ice: Setup - Structure Found : '
     write(6,*) '-----------------------------------------'
     write(6,'(a,f8.3,a)') ' ice dipole   = ', dip ,' D'
     deallocate (z_tmp)
  endif

  ! Copy bead 1 coordinates to other beads
  
  do k = 2,nb
     do i = 1,na
        r(1,i,k) = r(1,i,1)
        r(2,i,k) = r(2,i,1)
        r(3,i,k) = r(3,i,1)
     enddo
  enddo

  return
end subroutine setup_positions
      

subroutine cubic(r,nm,boxlxyz,ncellxyz)
  implicit none
  ! ------------------------------------------------------------------
  ! Sets up a cubic lattice of water molecules with a specified
  ! geometry.
  ! ------------------------------------------------------------------
  integer i,j,k,ncellxyz(3),itel,nm
  real(8) r(9,nm),boxlxyz(3),del(3),theta,reoh
  real(8) gx,gy,gz,hx,hy,hz
  real(8) dx,dy,dz,ox,oy,oz
  common /geometry/ theta,reoh

  ! Water molecule geometry
  
  hx = reoh
  hy = 0.d0
  hz = 0.d0
  gx = reoh*dcos(theta)
  gy = reoh*dsin(theta)
  gz = 0.d0
  
  if (ncellxyz(1)*ncellxyz(2)*ncellxyz(3).ne.nm) then
     write(6,*)'* Problem in cubic setup! '
     stop
  endif
  
  del(:) = boxlxyz(:)/dble(ncellxyz(:))
  itel = 0
  dx = -del(1)
  do i = 1, ncellxyz(1)
     dx = dx + del(1)
     dy = -del(2)
     do j = 1, ncellxyz(2)
        dy = dy + del(2)
        dz = -del(3)
        do k = 1, ncellxyz(3)
           dz = dz + del(3)
           if (itel.lt.nm) then
              itel = itel + 1
              ox = dx
              oy = dy
              oz = dz
              r(1,itel) = ox
              r(2,itel) = oy
              r(3,itel) = oz
              r(4,itel) = ox + hx
              r(5,itel) = oy + hy
              r(6,itel) = oz + hz
              r(7,itel) = ox + gx
              r(8,itel) = oy + gy
              r(9,itel) = oz + gz
           endif
        enddo
     enddo
  enddo
  return
end subroutine cubic


subroutine fcc(r,nm,boxlxyz,ncellxyz)
  implicit none
  ! ------------------------------------------------------------------
  ! Sets up an fcc of water molecules with a specified geometry.
  ! ------------------------------------------------------------------
  integer i,j,k,l,ncellxyz(3),itel,nm
  real(8) r(9,nm),boxlxyz(3),theta,reoh
  real(8) gx,gy,gz,hx,hy,hz,lat(3),halfl(3),ox,oy,oz
  real(8) dx(4),dy(4),dz(4),x0,y0,z0
  common /geometry/ theta,reoh
  
  ! Water molecule geometry

  hx = reoh
  hy = 0.d0
  hz = 0.d0
  gx = reoh*dcos(theta)
  gy = reoh*dsin(theta)
  gz = 0.d0

  if (4*ncellxyz(1)*ncellxyz(2)*ncellxyz(3).ne.nm) then
     write(6,*)'* Problem in fcc setup! '
     write(6,*) ncellxyz(1),ncellxyz(2),ncellxyz(3),nm
     stop
  endif
  lat(:) = boxlxyz(:)/dble(ncellxyz(:))

  ! Unit cell postions

  halfl(1) = 0.5d0*lat(1)
  halfl(2) = 0.5d0*lat(2)
  halfl(3) = 0.5d0*lat(3)

  dx(1) = 0.d0
  dx(2) = halfl(1)
  dx(3) = halfl(2)
  dx(4) = 0.d0
  dy(1) = 0.d0
  dy(2) = halfl(2)
  dy(3) = 0.d0
  dy(4) = halfl(2)
  dz(1) = 0.d0
  dz(2) = 0.d0
  dz(3) = halfl(3)
  dz(4) = halfl(3)

  ! Loop over unit cells

  itel = 0
  do k = 1,ncellxyz(1)
     z0 = (k-1)*lat(1)
     do j = 1,ncellxyz(2)
        y0 = (j-1)*lat(2)
        do i = 1,ncellxyz(3)
           x0 = (i-1)*lat(3)

           ! Loop over the 4 water molecules per unit cell

           do l = 1,4
              
              ! Oxygens
              
              ox = x0 + dx(l)
              oy = y0 + dy(l)
              oz = z0 + dz(l)
              r(1,itel+l) = ox
              r(2,itel+l) = oy
              r(3,itel+l) = oz

              ! Hydrogens

              r(4,itel+l) = ox + hx
              r(5,itel+l) = oy + hy
              r(6,itel+l) = oz + hz
              r(7,itel+l) = ox + gx
              r(8,itel+l) = oy + gy
              r(9,itel+l) = oz + gz
           enddo
           itel = itel + 4
        enddo
     enddo
  enddo
  
  return
end subroutine fcc
