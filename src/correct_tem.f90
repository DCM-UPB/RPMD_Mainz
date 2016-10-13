subroutine ewc_driver_tem(r,z,dvdr,v,vir,na,boxlxyz,list,point)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Ewald PI Correction Driver - TEM Switched Method
  ! ------------------------------------------------------------------
  integer na,point(na+3),list(maxnab*na),ic,i,j
  real(8) boxlxyz(3),r(3,na),dvdr(3,na),vir(3,3),z(na)
  real(8) sig,v,alpha,alpha2,ecut,wm,wh,boxmax
  real(8), allocatable :: ro(:,:)
  common /correct/ sig
  common /ew_param/ alpha,ecut,wm,wh

  if ((wm.gt.0.d0).or.(wh.gt.0.d0)) then
     write(6,*) 'ewc_correction not yet &
                implemented for smeared charges'
     stop
  endif

  alpha2 = 0.5d0 * (1.d0 - alpha)
  if (abs(alpha-1.d0).gt.1d-4) then

     ! Store the oxygen coordinates and replace the oxygen
     ! coordinates with the M-coordinates.

     allocate (ro(3,na/3))

     ic = 0
     do i = 1,na,3
        ic = ic + 1
        do j = 1, 3
           ro(j,ic) = r(j,i)
           r(j,i) = alpha * r(j,i) + alpha2*(r(j,i+1)+r(j,i+2))
        enddo
     enddo
  endif

  if (list(1).ne.0) then
     call ewc_list_tem(r,z,dvdr,v,vir,na,boxlxyz,list,point)
  else
     boxmax = maxval(boxlxyz)
     if (3.d0*sig .gt. boxmax) then
        call ewc_basic_tem(r,z,dvdr,v,vir,na,boxlxyz)
     else
        
        ! Use linked cell list
        
        call ewc_cell_tem(r,z,dvdr,v,vir,na,boxlxyz)
     endif
  endif

  if (abs(alpha-1.d0).gt.1d-4) then

     ! Use the chain rule to calculate the correct forces on the atoms
     ! and replace the m-site with the original oxygen atoms.

     ic = 0
     do j = 1, na, 3
        ic = ic + 1
        do i = 1, 3
           dvdr(i,j+1) = dvdr(i,j+1) + alpha2 * dvdr(i,j)
           dvdr(i,j+2) = dvdr(i,j+2) + alpha2 * dvdr(i,j)
           dvdr(i,j) = alpha * dvdr(i,j)
           r(i,j) = ro(i,ic)
        enddo
     enddo
     deallocate (ro)
  endif

  return
end subroutine ewc_driver_tem

subroutine ewc_basic_tem(r,z,dvdr,v,vir,na,boxlxyz)
  implicit none
  ! ------------------------------------------------------------------
  ! Ewald correction using TEM switch - Basic Version
  ! ------------------------------------------------------------------
  integer na,j,i,imol,mol(na)
  real(8) r(3,na),z(na),dvdr(3,na),vir(3,3)
  real(8) v,dx,dy,dz,boxlxyz(3),sig
  real(8) drsq,dr,rcut,rcutsq,boxmin,boxlx,boxly,boxlz
  real(8) zij,dvr,dvx,dvy,dvz,dv,onr
  real(8) onboxx,onboxy,onboxz
  real(8) a1,a2,a3,dr3,dr4,vc,dvc,f,df
  common /correct/ sig

  ! TEM switch

  rcut = sig
  rcutsq = rcut*rcut
  a1 = 2.d0/sig
  a2 = 2.d0/(sig**3)
  a3 = 1.d0/(sig**4)

  boxlx = boxlxyz(1)
  boxly = boxlxyz(2)
  boxlz = boxlxyz(3)
  onboxx = 1.d0/boxlx
  onboxy = 1.d0/boxly
  onboxz = 1.d0/boxlz

  ! Check

  boxmin = max(boxlx,boxly,boxlz)
  if (rcut.gt.0.5d0*boxmin) then
     write(6,*) ' ewc_basic : rcut too large '
     write(6,*) rcut , boxmin
     stop
  endif

  ! Clear stores
  
  v = 0.d0
  vir(:,:) = 0.d0
  dvdr(:,:) = 0.d0

  ! Assign molecular identities to each atom.

  imol = 0
  do j = 1,na,3
     imol = imol+1
     mol(j) = imol
     mol(j+1) = imol
     mol(j+2) = imol
  enddo

  do j = 1,na
     do i = 1,j-1

        ! Check atoms are not in same molecule
        
        if (mol(i).ne.mol(j)) then
           dx = r(1,i)-r(1,j)
           dy = r(2,i)-r(2,j)
           dz = r(3,i)-r(3,j)
           dx = dx - boxlx*dble(nint(onboxx*dx))
           dy = dy - boxly*dble(nint(onboxy*dy))
           dz = dz - boxlz*dble(nint(onboxz*dz))
           drsq = dx*dx + dy*dy + dz*dz
           if (drsq.lt.rcutsq) then
              dr = dsqrt(drsq)
              onr = 1.d0/dr
              zij = z(i)*z(j)
              vc = zij*onr
              dvc = -vc*onr
              dr3 = drsq*dr
              dr4 = dr3*dr
              f = 1.d0 - a1*dr + a2*dr3 - a3*dr4
              df = -a1 + 3.d0*a2*drsq - 4.d0*a3*dr3
              dv = vc*f
              dvr = vc*df + dvc*f
              dvx = dvr * dx * onr
              dvy = dvr * dy * onr
              dvz = dvr * dz * onr
              v = v + dv
              dvdr(1,i) = dvdr(1,i) + dvx
              dvdr(2,i) = dvdr(2,i) + dvy
              dvdr(3,i) = dvdr(3,i) + dvz
              dvdr(1,j) = dvdr(1,j) - dvx
              dvdr(2,j) = dvdr(2,j) - dvy
              dvdr(3,j) = dvdr(3,j) - dvz
              vir(1,1) = vir(1,1) + dx * dvx
              vir(2,2) = vir(2,2) + dy * dvy
              vir(3,3) = vir(3,3) + dz * dvz
              vir(1,2) = vir(1,2) + dx * dvy
              vir(1,3) = vir(1,3) + dx * dvz
              vir(2,1) = vir(2,1) + dy * dvx
              vir(2,3) = vir(2,3) + dy * dvz
              vir(3,1) = vir(3,1) + dz * dvx
              vir(3,2) = vir(3,2) + dz * dvy
           endif
        endif
     enddo
  enddo
  
  return
end subroutine ewc_basic_tem

subroutine ewc_list_tem(r,z,dvdr,v,vir,na,boxlxyz,list,point)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! TEM switch : using Neighbour List
  ! ------------------------------------------------------------------
  integer na,j,i,imol,mol(na),point(na+3),list(maxnab*na)
  integer ibeg,iend,inab
  real(8) r(3,na),z(na),dvdr(3,na),vir(3,3)
  real(8) v,dx,dy,dz,boxlxyz(3),sig
  real(8) drsq,dr,rcut,rcutsq,boxmax,boxx,boxy,boxz
  real(8) zij,dvr,dvx,dvy,dvz,dv,onr
  real(8) onboxx,onboxy,onboxz,rxj,ryj,rzj
  real(8) dvdxj,dvdyj,dvdzj
  real(8) a1,a2,a3,dr3,dr4,vc,dvc,f,df
  common /correct/ sig

  ! TEM 2nd derivative correction method
  
  rcut = sig
  rcutsq = rcut*rcut
  a1 = 2.d0/sig
  a2 = 2.d0/(sig**3)
  a3 = 1.d0/(sig**4)

  boxx = boxlxyz(1)
  boxy = boxlxyz(2)
  boxz = boxlxyz(3)
  onboxx = 1.d0/boxx
  onboxy = 1.d0/boxy
  onboxz = 1.d0/boxz

  ! Check

  boxmax = max(boxx,boxy,boxz)
  if (rcut.gt.0.5d0*boxmax) then
     write(6,*) ' ewc_list : rcut too large '
     write(6,*) rcut , boxmax
     stop
  endif

  ! Clear stores

  v = 0.d0
  vir(:,:) = 0.d0
  dvdr(:,:) = 0.d0

  ! Assign molecular identities to each atom.

  imol = 0
  do j = 1,na,3
     imol = imol+1
     mol(j) = imol
     mol(j+1) = imol
     mol(j+2) = imol
  enddo

  do j = 1,na
     ibeg = point(j)
     iend = point(j+1) - 1
     if (ibeg.le.iend) then
        rxj = r(1,j)
        ryj = r(2,j)
        rzj = r(3,j)
        dvdxj = dvdr(1,j)
        dvdyj = dvdr(2,j)
        dvdzj = dvdr(3,j)
        do inab = ibeg,iend
           i = list(inab)
           if (mol(i).ne.mol(j)) then
              dx = r(1,i) - rxj
              dy = r(2,i) - ryj
              dz = r(3,i) - rzj
              dx = dx - boxx*dble(nint(onboxx*dx))
              dy = dy - boxy*dble(nint(onboxy*dy))
              dz = dz - boxz*dble(nint(onboxz*dz))
              drsq = dx*dx + dy*dy + dz*dz
              if (drsq .lt. rcutsq) then
                 dr = dsqrt(drsq)
                 onr = 1.d0/dr
                 zij = z(i)*z(j)
                 vc = zij*onr
                 dvc = -vc*onr
                 dr3 = drsq*dr
                 dr4 = dr3*dr
                 f = 1.d0 - a1*dr + a2*dr3 - a3*dr4
                 df = -a1 + 3.d0*a2*drsq - 4.d0*a3*dr3
                 dv = vc*f
                 dvr = vc*df + dvc*f
                 dvx = dvr*dx*onr
                 dvy = dvr*dy*onr
                 dvz = dvr*dz*onr
                 v = v + dv
                 dvdr(1,i) = dvdr(1,i) + dvx
                 dvdr(2,i) = dvdr(2,i) + dvy
                 dvdr(3,i) = dvdr(3,i) + dvz
                 dvdxj = dvdxj - dvx
                 dvdyj = dvdyj - dvy
                 dvdzj = dvdzj - dvz
                 vir(1,1) = vir(1,1) + dx * dvx
                 vir(2,2) = vir(2,2) + dy * dvy
                 vir(3,3) = vir(3,3) + dz * dvz
                 vir(1,2) = vir(1,2) + dx * dvy
                 vir(1,3) = vir(1,3) + dx * dvz
                 vir(2,1) = vir(2,1) + dy * dvx
                 vir(2,3) = vir(2,3) + dy * dvz
                 vir(3,1) = vir(3,1) + dz * dvx
                 vir(3,2) = vir(3,2) + dz * dvy
              endif
           endif
        enddo
        dvdr(1,j) = dvdxj
        dvdr(2,j) = dvdyj
        dvdr(3,j) = dvdzj
     endif
  enddo

  return
end subroutine ewc_list_tem

subroutine ewc_cell_tem(r,z,dvdr,v,vir,na,boxlxyz)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! TEM switch : using cell List
  ! -----------------------------------------------------------------
  integer na,j,i,k,imol,ix,iy,iz,jx,jy,jz,incx,incy
  integer incz,nn,list(na)
  integer mol(na),ncellxyz(3),mcellxyz(3)
  real(8) r(3,na),z(na),dvdr(3,na),boxlxyz(3),vir(3,3)
  real(8) v,dx,dy,dz,sig,drsq,dr,rcut,rcutsq
  real(8) zij,dvr,dvx,dvy,dvz,dv,onr
  real(8) a1,a2,a3,dr3,dr4,vc,dvc,f,df
  real(8) boxx,boxy,boxz,onboxx,onboxy,onboxz,boxmin
  integer, allocatable :: link(:,:,:),ink(:,:)
  integer, allocatable :: mapx(:,:),mapy(:,:),mapz(:,:)
  common /correct/ sig

  rcut = sig
  rcutsq = rcut*rcut

  ! generate linked cell list
  
  call cell_setup(boxlxyz,rcut,mcellxyz,ncellxyz,nn)

  allocate (link(ncellxyz(1),ncellxyz(2),ncellxyz(3)))
  allocate (ink(3,nn))
  allocate (mapx(ncellxyz(1),-mcellxyz(1):mcellxyz(1)))
  allocate (mapy(ncellxyz(2),-mcellxyz(2):mcellxyz(2)))
  allocate (mapz(ncellxyz(3),-mcellxyz(3):mcellxyz(3)))

  call cell_list(r,list,link,ink,mapx,mapy,mapz,nn, &
                 na,boxlxyz,ncellxyz,mcellxyz,1)

  ! TEM switch

  a1 = 2.d0/sig
  a2 = 2.d0/(sig**3)
  a3 = 1.d0/(sig**4)

  boxx = boxlxyz(1)
  boxy = boxlxyz(2)
  boxz = boxlxyz(3)

  onboxx = 1.d0/boxx
  onboxy = 1.d0/boxy
  onboxz = 1.d0/boxz

  boxmin = min(boxx,boxy,boxz)
  if (rcut.gt.0.5d0*boxmin) then
     write(6,*) ' ewc_basic : rcut too large '
     write(6,*) rcut , boxmin
     stop
  endif

  v = 0.d0
  vir(:,:) = 0.d0
  dvdr(:,:) = 0.d0

  ! Assign molecular identities to each atom.

  imol = 0
  do j = 1,na,3
     imol = imol+1
     mol(j) = imol
     mol(j+1) = imol
     mol(j+2) = imol
  enddo

  ! evaluate potential energy and forces using the linked cell list

  do iz = 1,ncellxyz(3)
     do iy = 1,ncellxyz(2)
        do ix = 1,ncellxyz(1)
           
           ! including contributions from within the same cell
           
           i = link(ix,iy,iz)
           do while (i .gt. 0)
              j = list(i)
              do while (j .gt. 0)
                 if (mol(i).ne.mol(j)) then
                    dx = r(1,i)-r(1,j)
                    dy = r(2,i)-r(2,j)
                    dz = r(3,i)-r(3,j)
                    dx = dx-boxx*dble(nint(onboxx*dx))
                    dy = dy-boxy*dble(nint(onboxy*dy))
                    dz = dz-boxz*dble(nint(onboxz*dz))
                    drsq = dx*dx + dy*dy + dz*dz
                    if (drsq .lt. rcutsq) then
                       dr = dsqrt(drsq)
                       onr = 1.d0/dr
                       zij = z(i)*z(j)
                       vc = zij*onr
                       dvc = -vc*onr
                       dr3 = drsq*dr
                       dr4 = dr3*dr
                       f = 1.d0 - a1*dr + a2*dr3 - a3*dr4
                       df = -a1 + 3.d0*a2*drsq - 4.d0*a3*dr3
                       dv = vc*f
                       dvr = vc*df + dvc*f
                       dvx = dvr*dx*onr
                       dvy = dvr*dy*onr
                       dvz = dvr*dz*onr
                       v = v + dv
                       dvdr(1,i) = dvdr(1,i) + dvx
                       dvdr(2,i) = dvdr(2,i) + dvy
                       dvdr(3,i) = dvdr(3,i) + dvz
                       dvdr(1,j) = dvdr(1,j) - dvx
                       dvdr(2,j) = dvdr(2,j) - dvy
                       dvdr(3,j) = dvdr(3,j) - dvz
                       vir(1,1) = vir(1,1) + dx * dvx
                       vir(2,2) = vir(2,2) + dy * dvy
                       vir(3,3) = vir(3,3) + dz * dvz
                       vir(1,2) = vir(1,2) + dx * dvy
                       vir(1,3) = vir(1,3) + dx * dvz
                       vir(2,1) = vir(2,1) + dy * dvx
                       vir(2,3) = vir(2,3) + dy * dvz
                       vir(3,1) = vir(3,1) + dz * dvx
                       vir(3,2) = vir(3,2) + dz * dvy
                    endif
                 endif
                 j = list(j)
              enddo
              i = list(i)
           enddo
        enddo
     enddo
  enddo
  
  ! and contributions from neighbouring cells

  do k = 1,nn
     incx = ink(1,k)
     incy = ink(2,k)
     incz = ink(3,k)
     do iz = 1,ncellxyz(3)
        jz = mapz(iz,incz)
        do iy = 1,ncellxyz(2)
           jy = mapy(iy,incy)
           do ix = 1,ncellxyz(1)
              jx = mapx(ix,incx)
              i = link(ix,iy,iz)
              do while (i .gt. 0)
                 j = link(jx,jy,jz)
                 do while (j .gt. 0)
                    if (mol(i).ne.mol(j)) then
                       dx = r(1,i)-r(1,j)
                       dy = r(2,i)-r(2,j)
                       dz = r(3,i)-r(3,j)
                       dx = dx-boxx*dble(nint(onboxx*dx))
                       dy = dy-boxy*dble(nint(onboxy*dy))
                       dz = dz-boxz*dble(nint(onboxz*dz))
                       drsq = dx*dx+dy*dy+dz*dz
                       if (drsq .lt. rcutsq) then
                          dr = dsqrt(drsq)
                          onr = 1.d0/dr
                          zij = z(i)*z(j)
                          vc = zij*onr
                          dvc = -vc*onr
                          dr3 = drsq*dr
                          dr4 = dr3*dr
                          f = 1.d0 - a1*dr + a2*dr3 - a3*dr4
                          df = -a1 + 3.d0*a2*drsq - 4.d0*a3*dr3
                          dv = vc*f
                          dvr = vc*df + dvc*f
                          dvx = dvr*dx*onr
                          dvy = dvr*dy*onr
                          dvz = dvr*dz*onr
                          v = v + dv
                          dvdr(1,i) = dvdr(1,i) + dvx
                          dvdr(2,i) = dvdr(2,i) + dvy
                          dvdr(3,i) = dvdr(3,i) + dvz
                          dvdr(1,j) = dvdr(1,j) - dvx
                          dvdr(2,j) = dvdr(2,j) - dvy
                          dvdr(3,j) = dvdr(3,j) - dvz
                          vir(1,1) = vir(1,1) + dx * dvx
                          vir(2,2) = vir(2,2) + dy * dvy
                          vir(3,3) = vir(3,3) + dz * dvz
                          vir(1,2) = vir(1,2) + dx * dvy
                          vir(1,3) = vir(1,3) + dx * dvz
                          vir(2,1) = vir(2,1) + dy * dvx
                          vir(2,3) = vir(2,3) + dy * dvz
                          vir(3,1) = vir(3,1) + dz * dvx
                          vir(3,2) = vir(3,2) + dz * dvy
                       endif
                    endif
                    j = list(j)
                 enddo
                 i = list(i)
              enddo
           enddo
        enddo
     enddo
  enddo

  deallocate (link,ink,mapx,mapy,mapz)
  
  return
end subroutine ewc_cell_tem
      
