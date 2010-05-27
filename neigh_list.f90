subroutine nei_driver(r,na,list,point,njump,rcut,boxlxyz)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Driver for neighbour list creation
  ! ------------------------------------------------------------------
  integer na,point(na+3),list(maxnab*na),njump
  real(8) r(3,na),boxlxyz(3),rcut,boxmax

  boxmax = maxval(boxlxyz)
  if (3.d0*rcut .gt. boxmax) then
     call nei_list(r,na,list,point,njump,rcut,boxlxyz)
  else

     ! Use linked cell list to make neighbour list

     call nei_cell(r,na,list,point,njump,rcut,boxlxyz)
  endif

  return
end subroutine nei_driver

subroutine nei_list(r,na,list,point,njump,rcut,boxlxyz)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Creates a neighbour list
  ! ------------------------------------------------------------------
  integer i,j,na,njump,point(na+3),list(maxnab*na),nlist,jend
  real(8) r(3,na),dx,dy,dz,boxlxyz(3)
  real(8) onboxx,onboxy,onboxz,boxx,boxy,boxz
  real(8) rxj,ryj,rzj,drsq,rcut,rcutsq

  rcutsq = rcut*rcut
  boxx = boxlxyz(1)
  boxy = boxlxyz(2)
  boxz = boxlxyz(3)
  onboxx = 1.d0/boxx
  onboxy = 1.d0/boxy
  onboxz = 1.d0/boxz

  point(:) = 0
  list(:) = 0

  nlist = 0
  jend = na-njump+1

  do j = 1,jend,njump
     point(j) = nlist + 1
     rxj = r(1,j)
     ryj = r(2,j)
     rzj = r(3,j)
     do i = j+njump,na,njump
        dx = r(1,i) - rxj
        dy = r(2,i) - ryj
        dz = r(3,i) - rzj
        dx = dx - boxx*nint(onboxx*dx)
        dy = dy - boxy*nint(onboxy*dy)
        dz = dz - boxz*nint(onboxz*dz)
        drsq = dx*dx + dy*dy + dz*dz
        if (drsq .lt. rcutsq) then
           nlist = nlist + 1
           if (nlist.eq.maxnab*na) stop 'neighbour list too small'
           list(nlist) = i
        endif
     enddo
  enddo

  do i = jend+njump,na+njump,njump
     nlist = nlist + 1
     point(i) = nlist
  enddo

  return
end subroutine nei_list

subroutine nei_cell(r,na,list,point,njump,rcut,boxlxyz)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Use a linked cell list to create a neighbour list
  ! ------------------------------------------------------------------
  integer i,j,k,na,njump,nlist,nn,jend
  integer incx,incy,incz,jx,jy,jz,ix,iy,iz
  integer mcellxyz(3),ncellxyz(3)
  integer clist(na),point(na+3),list(maxnab*na)
  real(8) r(3,na),boxlxyz(3),dx,dy,dz
  real(8) rxj,ryj,rzj,drsq,rcut,rcutsq
  real(8) boxx,boxy,boxz,onboxx,onboxy,onboxz,halfx,halfy,halfz
  real(8) cellx,celly,cellz
  integer, allocatable :: clink(:,:,:),ink(:,:)
  integer, allocatable :: mapx(:,:),mapy(:,:),mapz(:,:)

  ! generate linked cell list

  call cell_setup(boxlxyz,rcut,mcellxyz,ncellxyz,nn)

  allocate (clink(ncellxyz(1),ncellxyz(2),ncellxyz(3)))
  allocate (ink(3,nn))
  allocate (mapx(ncellxyz(1),-mcellxyz(1):mcellxyz(1)))
  allocate (mapy(ncellxyz(2),-mcellxyz(2):mcellxyz(2)))
  allocate (mapz(ncellxyz(3),-mcellxyz(3):mcellxyz(3)))

  call cell_list(r,clist,clink,ink,mapx,mapy,mapz,nn, &
                 na,boxlxyz,ncellxyz,mcellxyz,njump)

  ! use linked cell list to create neighbour list
  
  boxx = boxlxyz(1)
  boxy = boxlxyz(2)
  boxz = boxlxyz(3)

  onboxx = 1.d0/boxx
  onboxy = 1.d0/boxy
  onboxz = 1.d0/boxz

  halfx = 0.5d0*boxx
  halfy = 0.5d0*boxy
  halfz = 0.5d0*boxz

  rcutsq = rcut*rcut
  cellx = ncellxyz(1)/boxx
  celly = ncellxyz(2)/boxy
  cellz = ncellxyz(3)/boxz

  point(:) = 0
  list(:) = 0

  ! Loop over particles

  nlist = 0
  jend = na-njump+1

  do j = 1,jend,njump
     point(j) = nlist + 1
     rxj = r(1,j)
     ryj = r(2,j)
     rzj = r(3,j)

     ! Find cell particle is in:

     dx = r(1,j)-boxx*nint(onboxx*r(1,j))
     dy = r(2,j)-boxy*nint(onboxy*r(2,j))
     dz = r(3,j)-boxz*nint(onboxz*r(3,j))
     jx = min(ncellxyz(1),1+int((dx+halfx)*cellx))
     jy = min(ncellxyz(2),1+int((dy+halfy)*celly))
     jz = min(ncellxyz(3),1+int((dz+halfz)*cellz))

     ! Find it's place in the head of chain array

     i = clink(jx,jy,jz)
     do while (i.ne.j)
        i = clist(i)
     enddo

     ! Move to next particle in chain

     i = clist(i)

     ! Neigbours in same cell (below it in chain)

     do while (i.gt.0)
        dx = r(1,i) - rxj
        dy = r(2,i) - ryj
        dz = r(3,i) - rzj
        dx = dx - boxx*nint(onboxx*dx)
        dy = dy - boxy*nint(onboxy*dy)
        dz = dz - boxz*nint(onboxz*dz)
        drsq = dx*dx + dy*dy + dz*dz
        if (drsq .lt. rcutsq) then
           nlist = nlist + 1
           if (nlist.eq.maxnab*na) stop 'neighbour list too small'
           list(nlist) = i
        endif
        i = clist(i)
     enddo

     ! Neighbours in other cells

     do k = 1,nn

        ! Find a neighbouring cell

        incx = ink(1,k)
        incy = ink(2,k)
        incz = ink(3,k)
        ix = mapx(jx,incx)
        iy = mapy(jy,incy)
        iz = mapz(jz,incz)

        ! Loop through particles in neighbour cell

        i = clink(ix,iy,iz)
        do while (i .gt. 0)
           dx = r(1,i) - rxj
           dy = r(2,i) - ryj
           dz = r(3,i) - rzj
           dx = dx-boxx*nint(onboxx*dx)
           dy = dy-boxy*nint(onboxy*dy)
           dz = dz-boxz*nint(onboxz*dz)
           drsq = dx*dx+dy*dy+dz*dz
           if (drsq .lt. rcutsq) then
              nlist = nlist + 1
              if (nlist.eq.maxnab*na) stop 'neigh: list too small'
              list(nlist) = i
           endif
           i = clist(i)
        enddo
     enddo
  enddo

  do i = jend+njump,na+njump,njump
     nlist = nlist + 1
     point(i) = nlist
  enddo

  deallocate (clink,ink,mapx,mapy,mapz)
  
  return
end subroutine nei_cell
