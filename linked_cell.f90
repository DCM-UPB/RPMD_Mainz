subroutine cell_setup(boxlxyz,rcut,mcellxyz,ncellxyz,nn)
  implicit none
  ! ------------------------------------------------------------------
  ! Works out the cell and neighbour setup for linked cell list
  ! ------------------------------------------------------------------
  integer ncellxyz(3),mcellxyz(3),nn,i
  real(8) boxlxyz(3),rcut

  nn = 1
  do i = 1,3
     if (3.d0*rcut .gt. boxlxyz(i)) then

        ! No cell list in this direction
        
        ncellxyz(i) = 1
        mcellxyz(i) = 0
     else

        ! Cell list with cell size 0.5 rcut (or as close as possible)

        ncellxyz(i) = max(6,int(2.d0*boxlxyz(i)/rcut))
        mcellxyz(i) = 2
     endif
     nn = nn * (2*mcellxyz(i)+1)
  enddo

  ! nn is the number of neighbours around each cell to be checked

  nn = (nn-1)/2

  return
end subroutine cell_setup

subroutine cell_list(r,list,link,ink,mapx,mapy,mapz,nn, &
                     n,boxlxyz,ncellxyz,mcellxyz,ijump)
  implicit none
  ! ------------------------------------------------------------------
  ! Constructs a linked cell list - cell size is 0.5*rcut
  ! ------------------------------------------------------------------
  integer i,ix,iy,iz,ijump,n,nn
  integer ncellxyz(3),mcellxyz(3)
  real(8) r(3,n),boxlxyz(3)
  real(8) dx,dy,dz,halfx,halfy,halfz,cellx,celly,cellz
  real(8) boxx,boxy,boxz,onboxx,onboxy,onboxz
  integer link(ncellxyz(1),ncellxyz(2),ncellxyz(3)),list(n)
  integer ink(3,nn)
  integer mapx(ncellxyz(1),-mcellxyz(1):mcellxyz(1))
  integer mapy(ncellxyz(2),-mcellxyz(2):mcellxyz(2))
  integer mapz(ncellxyz(3),-mcellxyz(3):mcellxyz(3))

  boxx = boxlxyz(1)
  boxy = boxlxyz(2)
  boxz = boxlxyz(3)
  onboxx = 1.d0/boxx
  onboxy = 1.d0/boxy
  onboxz = 1.d0/boxz

  halfx = 0.5d0*boxx
  halfy = 0.5d0*boxy
  halfz = 0.5d0*boxz

  ! Assign particles to cells

  link(:,:,:) = 0
  list(:) = 0

  cellx = dble(ncellxyz(1))/boxx
  celly = dble(ncellxyz(2))/boxy
  cellz = dble(ncellxyz(3))/boxz

  do i = 1,n,ijump
     dx = r(1,i)-boxx*nint(onboxx*r(1,i))
     dy = r(2,i)-boxy*nint(onboxy*r(2,i))
     dz = r(3,i)-boxz*nint(onboxz*r(3,i))
     ix = min(ncellxyz(1),1+int((dx+halfx)*cellx))
     iy = min(ncellxyz(2),1+int((dy+halfy)*celly))
     iz = min(ncellxyz(3),1+int((dz+halfz)*cellz))
     list(i) = link(ix,iy,iz)
     link(ix,iy,iz) = i
  enddo

  ! Make a map of neigbouring cells

  if (nn.gt.0) then
     call cell_map_nc(mcellxyz,ncellxyz,mapx,mapy,mapz,ink,nn)
  endif

  return
end subroutine cell_list

subroutine cell_list_unit(r,list,link,ink,map,nn,n,ncell, &
                          mcell,ijump)
  implicit none
  ! ------------------------------------------------------------------
  ! Constructs a linked cell list for a box with all
  ! sides of unit length - cell size is 0.5*rcut
  ! ------------------------------------------------------------------
  integer i,ix,iy,iz,ncell,mcell,ijump,n,nn
  real(8) dx,dy,dz
  real(8) r(3,n),cell
  integer link(ncell,ncell,ncell),list(n)
  integer ink(3,62),map(ncell,-2:2)

  link(:,:,:) = 0
  
  cell = dble(ncell)
  do i = 1,n
     dx = r(1,i)-nint(r(1,i))
     dy = r(2,i)-nint(r(2,i))
     dz = r(3,i)-nint(r(3,i))
     ix = min(ncell,1+int((dx+0.5d0)*cell))
     iy = min(ncell,1+int((dy+0.5d0)*cell))
     iz = min(ncell,1+int((dz+0.5d0)*cell))
     list(i) = link(ix,iy,iz)
     link(ix,iy,iz) = i
  enddo

  ! Make a map of neigbouring cells

  nn = ((2*mcell+1)*(2*mcell+1)*(2*mcell+1)-1)/2
  if (mcell.gt.0) then
     call cell_map(mcell,ncell,map,ink,nn)
  else
     nn = 0
  endif

  return
end subroutine cell_list_unit

subroutine cell_map(mcell,ncell,map,ink,nn)
  implicit none
  ! ------------------------------------------------------------------
  ! Creates a neighbour map for the linked cell list
  ! (since cells are 0.5*rcut must check (5**3-1)/2 cells = 62)
  ! ------------------------------------------------------------------
  integer k,mcell,ncell,nn,ix,iy,iz,map(ncell,-2:2),ink(3,62)
  integer jx,kx

  k = 0
  do iz = 1,mcell
     do iy = -mcell,mcell
        do ix = -mcell,mcell
           k = k+1
           ink(1,k) = ix
           ink(2,k) = iy
           ink(3,k) = iz
        enddo
     enddo
  enddo
  iz = 0
  do iy = 1,mcell
     do ix = -mcell,mcell
        k = k+1
        ink(1,k) = ix
        ink(2,k) = iy
        ink(3,k) = iz
     enddo
  enddo
  iy = 0
  do ix = 1,mcell
     k = k+1
     ink(1,k) = ix
     ink(2,k) = iy
     ink(3,k) = iz
  enddo
  
  if (k.ne.nn) then
     write(6,*) 'cell_map : k is not nn',k, nn
     stop
  endif

  do ix = 1,ncell
     do jx = -mcell,mcell
        kx = ix+jx
        if (kx .lt. 1) then
           kx = kx + ncell
        else if (kx .gt. ncell) then
           kx = kx - ncell
        endif
        map(ix,jx) = kx
     enddo
  enddo
  
  return
end subroutine cell_map

subroutine cell_map_nc(mcellxyz,ncellxyz,mapx,mapy,mapz,ink,nn)
  implicit none
  ! ------------------------------------------------------------------
  ! Creates a neighbour map for the linked cell list
  ! (since cells are 0.5*rcut must check (5**3-1)/2 cells = 62)
  ! ------------------------------------------------------------------
  integer nn,ix,iy,iz,jx,jy,jz,kx,ky,kz
  integer k,mcellxyz(3),ncellxyz(3),ink(3,nn)
  integer mapx(ncellxyz(1),-mcellxyz(1):mcellxyz(1))
  integer mapy(ncellxyz(2),-mcellxyz(2):mcellxyz(2))
  integer mapz(ncellxyz(3),-mcellxyz(3):mcellxyz(3))

  ! Create the ink array which contains the relative positions
  ! of all neighbours
  
  k = 0
  do iz = 1,mcellxyz(3)
     do iy = -mcellxyz(2),mcellxyz(2)
        do ix = -mcellxyz(1),mcellxyz(1)
           k = k+1
           ink(1,k) = ix
           ink(2,k) = iy
           ink(3,k) = iz
        enddo
     enddo
  enddo
  iz = 0
  do iy = 1,mcellxyz(2)
     do ix = -mcellxyz(1),mcellxyz(1)
        k = k+1
        ink(1,k) = ix
        ink(2,k) = iy
        ink(3,k) = iz
     enddo
  enddo
  iy = 0
  do ix = 1,mcellxyz(1)
     k = k+1
     ink(1,k) = ix
     ink(2,k) = iy
     ink(3,k) = iz
  enddo

  if (k.ne.nn) then
     write(6,*) 'cell_map_nc : k is not nn',k, nn
     stop
  endif

  ! Create map - converts the relative positions provided by ink
  !               into cell location vectors.

  ! Map for x direction

  do ix = 1,ncellxyz(1)
     do jx = -mcellxyz(1),mcellxyz(1)
        kx = ix+jx
        if (kx .lt. 1) then
           kx = kx + ncellxyz(1)
        else if (kx .gt. ncellxyz(1)) then
           kx = kx - ncellxyz(1)
        endif
        mapx(ix,jx) = kx
     enddo
  enddo
  
  ! Map for y direction

  do iy = 1,ncellxyz(2)
     do jy = -mcellxyz(2),mcellxyz(2)
        ky = iy+jy
        if (ky .lt. 1) then
           ky = ky + ncellxyz(2)
        else if (ky .gt. ncellxyz(2)) then
           ky = ky - ncellxyz(2)
        endif
        mapy(iy,jy) = ky
     enddo
  enddo

  ! Map for z direction

  do iz = 1,ncellxyz(3)
     do jz = -mcellxyz(3),mcellxyz(3)
        kz = iz+jz
        if (kz .lt. 1) then
           kz = kz + ncellxyz(3)
        else if (kz .gt. ncellxyz(3)) then
           kz = kz - ncellxyz(3)
        endif
        mapz(iz,jz) = kz
     enddo
  enddo

  return
end subroutine cell_map_nc
