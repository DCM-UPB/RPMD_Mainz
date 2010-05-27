subroutine ice_lattice(r,boxlxyz,ncellxyz,nm,irun,iunit)
  implicit none
  ! -----------------------------------------------------------------
  ! Setup of a Hexagonal Ice Lattice
  ! -----------------------------------------------------------------
  integer nm,irun,ncellxyz(3),iunit
  real(8) r(9,nm),boxlxyz(3),alat,blat,clat
  
  ! Lattice parameters

  alat = boxlxyz(1)/dble(ncellxyz(1))
  blat = boxlxyz(2)/dble(ncellxyz(2))
  clat = boxlxyz(3)/dble(ncellxyz(3))

  ! Place Oxygen Atoms on Hexagonal Lattice Sites

  call ox_place(r,ncellxyz,boxlxyz,alat,blat,clat,nm)

  ! Now Place Hydrogens

  call h_place(r,nm,irun,boxlxyz,alat,blat,clat,iunit)
  
  return
end subroutine ice_lattice

subroutine h_place(r,nm,irun,boxlxyz,alat,blat,clat,iunit)
  implicit none
  ! -----------------------------------------------------------------
  ! Arrange the Hydrogen Atoms - Random Assignment and MC Improvement
  ! -----------------------------------------------------------------
  integer irun,nm,i,k,nn,nstore,iunit
  integer neigh(2*nm,2),ncb(nm)
  real(8) r(9,nm),boxlxyz(3),ran2
  real(8) boxlx,boxly,boxlz,alat,blat,clat,rb,tol,rb2
  real(8) dx,dy,dz,dr2,ran,onboxx,onboxy,onboxz
  external ran2

  ! * Method from : Buch et al. J.Phys.Chem. B 102,44,8641

  ! Arrays :
  ! ---------
  ! neigh - neigbour pairs : neigh(i,1) stores chemically bonded
  ! oxygen of bond i and neigh(i,2) stores h-bonded oxygen.

  ! ncb(i) -  gives the number of chemical bonds to oxygen i(target=2)

  ! Useful Parameters -

  boxlx = boxlxyz(1)
  boxly = boxlxyz(2)
  boxlz = boxlxyz(3)
  onboxx = 1.d0/boxlx
  onboxy = 1.d0/boxly
  onboxz = 1.d0/boxlz
  
  ! Shortest O-O distance (4-Oxygens located at this distance)

  rb = (alat/2.d0)**2+(blat/6.d0)**2+(clat/8.d0)**2
  rb = dsqrt(rb)

  ! Add a small tolerance to rb

  tol = 1.d-5
  rb = rb + tol
  rb2 = rb*rb

  ! Clear arrays

  ncb(:) = 0

  ! Create Oxygen nearest neighbour matrix

  ! - Search distances to find nearest oxygens (shorter than dr)

  nn = 0
  ! Loop over oxygen atoms
  do i = 1,nm
     do k = 1,i-1
        dx = r(1,i)-r(1,k)
        dy = r(2,i)-r(2,k)
        dz = r(3,i)-r(3,k)
        dx = dx - boxlx * nint(onboxx*dx)
        dy = dy - boxly * nint(onboxy*dy)
        dz = dz - boxlz * nint(onboxz*dz)
        dr2 = dx*dx + dy*dy + dz*dz
        if (dr2.lt.rb2) then
           if (i.ne.k) then
              ! These oxygen atoms are a nearest neighbour pair
              nn = nn + 1
              if (nn.gt.2*nm) then
                 write(6,*) 'h_place : Too many neighbours', nn,nm
              endif
              neigh(nn,1) = i
              neigh(nn,2) = k
           endif
        endif
     enddo
  enddo
  if (nn.ne.2*nm) then
     write(6,*) ' *** POSSIBLE ERROR : '
     write(6,*) 'Number of bonds not equal to 2*number of oxygens!'
  endif

  ! Now we have the neighbour pairs randomly assign h-chemical bonds
  
  do i = 1,nn
     ran = ran2(irun,0.d0,1.d0)
     if (ran.lt.0.5d0) then
        nstore = neigh(i,1)
        neigh(i,1) = neigh(i,2)
        neigh(i,2) = nstore
     endif
  enddo
  
  ! Count the number of chemical bonds to each atom and store in ncb

  do i = 1,nn
     ncb(neigh(i,1)) = ncb(neigh(i,1))+1
  enddo

  ! Use Monte Carlo to correct structure so that all oxygens have
  ! chemical bonds.

  call mc_h_correction(neigh,ncb,nm,irun,iunit)

  ! Convert neighbour bonding array into H-coordinates

  call h_coordinates(neigh,r,nm,boxlxyz)
  
  return
end subroutine h_place

subroutine mc_h_correction(neigh,ncb,nm,irun,iunit)
  implicit none
  ! -----------------------------------------------------------------
  ! Performs Monte Carlo Moves to Correct the Ice Bonding
  ! -----------------------------------------------------------------
  integer nm,nhappy,irun,nit,i,jj,kk,iunit,neigh(2*nm,2),ncb(nm)
  integer nstore,nch,nch1,nch2,nmv
  real(8) ran2,test
  logical accept
  external ran2

  ! nhappy - stores the number of water molecules that have correct
  !          chemical bonding. When {nhappy = nm} structure is done.

  nhappy = 0
  do i = 1,nm
     if (ncb(i).eq.2) then
        nhappy = nhappy+1
     endif
  enddo

  ! Begin correction moves

  nit = 0
  do while (nhappy.lt.nm)
     nit = nit+1
     ! Bond to move
     nmv = 1 + int(2.d0*dble(nm)*ran2(irun,0.d0,1.d0))
     
     ! Decide whether to switch bond

     jj = neigh(nmv,1)
     kk = neigh(nmv,2)
     ! Difference in no. of chemical bonds before
     nch1 = abs(ncb(jj)-ncb(kk))
     ! Difference in no. of chemical bonds after
     nch2 = abs((ncb(jj)-1)-(ncb(kk)+1))
     ! Change in bonding difference
     nch = nch2 - nch1

     ! Acceptance conditions

     if (nch.lt.0) then
        ! If bonding difference decreased then accept
        accept = .true.
     else if (nch.gt.1) then
        ! If bonding difference increased then reject
        accept = .false.
     else if (nch.eq.0) then
        ! If bonding difference unchanged then 50:50 acceptance
        test = ran2(irun,0.d0,1.d0)
        if (test.lt.0.5) then
           accept = .true.
        else
           accept = .false.
        endif
     endif

     if (accept) then
        ! Swtich bond

        nstore = neigh(nmv,1)
        neigh(nmv,1) = neigh(nmv,2)
        neigh(nmv,2) = nstore

        ! Recompute happyness

        if (ncb(jj).eq.2) then
           nhappy = nhappy-1
        endif
        if (ncb(kk).eq.2) then
           nhappy = nhappy-1
        endif
        ncb(jj) = ncb(jj)-1
        ncb(kk) = ncb(kk)+1
        if (ncb(jj).eq.2) then
           nhappy = nhappy+1
        endif
        if (ncb(kk).eq.2) then
           nhappy = nhappy+1
        endif
     endif
  enddo
  if (iunit.ne.0) then
     write(iunit,*) '-------------------------------------------'
     write(iunit,*) ' Hexagonal Ice: Setup - Structure Found  '
     write(iunit,*) '    ', nit , '  MC moves required'
     write(iunit,*) '-------------------------------------------'
     write(iunit,*)
  endif

  return
end subroutine mc_h_correction

subroutine h_coordinates(neigh,r,nm,boxlxyz)
  implicit none
  !-----------------------------------------------------------------
  ! H-coordinates from the neighbour list
  !-----------------------------------------------------------------
  common /geometry/ theta,reoh
  integer nm,neigh(2*nm,2),nh,k,i,jj,kk
  real(8) r(9,nm),vec(3),theta,reoh,boxlxyz(3)
  real(8) dx,dy,dz,dr,onboxx,onboxy,onboxz
  real(8) boxlx,boxly,boxlz

  ! Useful Parameters

  boxlx = boxlxyz(1)
  boxly = boxlxyz(2)
  boxlz = boxlxyz(3)
  onboxx = 1.d0/boxlx
  onboxy = 1.d0/boxly
  onboxz = 1.d0/boxlz

  ! Place Oxygen at correct distance on vector to coordinates oxygen

  do i = 1,2*nm
     jj = neigh(i,1)
     kk = neigh(i,2)

     ! Normalized vector between oxygens

     dx = r(1,kk)-r(1,jj)
     dy = r(2,kk)-r(2,jj)
     dz = r(3,kk)-r(3,jj)
     dx = dx - boxlx*nint(onboxx*dx)
     dy = dy - boxly*nint(onboxy*dy)
     dz = dz - boxlz*nint(onboxz*dz)
     dr = dx*dx + dy*dy + dz*dz
     dr = dsqrt(dr)
     vec(1) = dx/dr
     vec(2) = dy/dr
     vec(3) = dz/dr

     ! Place Hydrogens at correct bond length along that vector

     if (r(4,jj).lt.1d-10) then
        nh = 3
     else
        nh = 6
     endif

     do k = 1,3
        r(nh+k,jj) = r(k,jj) + reoh*vec(k)
     enddo
  enddo

  ! Correct the bond angle to fit the water model chosen

  call correct_angle(r,nm)
  
  return
end subroutine h_coordinates

subroutine correct_angle(r,nm)
  implicit none
  !-----------------------------------------------------------------
  !Performs Axis Angle Rotation to Give correct bond angle
  !-----------------------------------------------------------------
  integer nm,i,k
  real(8) r(9,nm),v1(3),v2(3),vn(3)
  real(8) theta,reoh,dot,sz,sz1,sz2,angle,ach
  common /geometry/ theta,reoh
  
  ! Loop over oxygen atoms

  do i = 1,nm

     ! Vector to oxygen atoms

     do k = 1,3
        v1(k) = r(3+k,i) - r(k,i)
        v2(k) = r(6+k,i) - r(k,i)
     enddo
     
     ! Calculate bond angle between hydrogens -

     call dotp(v1,v2,dot)
     call dotp(v1,v1,sz1)
     call dotp(v2,v2,sz2)
     angle = dacos(dot/(dsqrt(sz1)*dsqrt(sz2)))

     ! Rotation of Hydrogens required -

     ach = 0.5d0*(theta-angle)

     ! Rotation axis (perpendicular to plane of H20)

     call xproduct(v1,v2,vn)

     ! Normalize vn
     
     call dotp(vn,vn,sz)
     sz = dsqrt(sz)
     do k = 1,3
        vn(k) = vn(k)/sz
     enddo

     ! Axis Angle rotation by angle ach about axis vn

     call axis_angle(v1,vn,ach)
     call axis_angle(v2,vn,-ach)

     ! Check new bond angle

     call dotp(v1,v2,dot)
     call dotp(v1,v1,sz1)
     call dotp(v2,v2,sz2)
     angle = dacos(dot/(dsqrt(sz1)*dsqrt(sz2)))

     ! Write new Hydrogen atoms vectors

     do k = 1,3
        r(3+k,i) = v1(k) + r(k,i)
        r(6+k,i) = v2(k) + r(k,i)
     enddo
  enddo
  return
end subroutine correct_angle

subroutine ox_place(r,ncellxyz,boxlxyz,alat,blat,clat,nm)
  implicit none
  ! -----------------------------------------------------------------
  ! Arrange the Oxygen Atoms
  ! -----------------------------------------------------------------
  integer i,j,k,ns,ia,nm,ncellxyz(3)
  real(8) r(9,nm),site(3,8),boxlxyz(3)
  real(8) x0,y0,z0,alat,blat,clat

  ! Clear positions

  r(:,:) = 0.d0

  ! Site Positions in Unit cell

  site(1,1) = 1.d0/4.d0*alat
  site(2,1) = 1.d0/6.d0*blat
  site(3,1) = 3.d0/16.d0*clat

  site(1,2) = 3.d0/4.d0*alat
  site(2,2) = 2.d0/6.d0*blat
  site(3,2) = 5.d0/16.d0*clat

  site(1,3) = 1.d0/4.d0*alat
  site(2,3) = 1.d0/6.d0*blat
  site(3,3) = 13.d0/16.d0*clat

  site(1,4) = 3.d0/4.d0*alat
  site(2,4) = 2.d0/6.d0*blat
  site(3,4) = 11.d0/16.d0*clat
  
  site(1,5) = 3.d0/4.d0*alat
  site(2,5) = 4.d0/6.d0*blat
  site(3,5) = 3.d0/16.d0*clat
  
  site(1,6) = 1.d0/4.d0*alat
  site(2,6) = 5.d0/6.d0*blat
  site(3,6) = 5.d0/16.d0*clat

  site(1,7) = 3.d0/4.d0*alat
  site(2,7) = 4.d0/6.d0*blat
  site(3,7) = 13.d0/16.d0*clat

  site(1,8) = 1.d0/4.d0*alat
  site(2,8) = 5.d0/6.d0*blat
  site(3,8) = 11.d0/16.d0*clat

  ! Periodic Repetitions

  ia = 0
  do k = 1,ncellxyz(3)
     z0 = (k-1)*clat
     do j = 1,ncellxyz(2)
        y0 = (j-1)*blat
        do i = 1,ncellxyz(1)
           x0 = (i-1)*alat
           do ns = 1,8
              r(1,ia+ns) = x0 + site(1,ns)
              r(2,ia+ns) = y0 + site(2,ns)
              r(3,ia+ns) = z0 + site(3,ns)
           enddo
           ia = ia+8
        enddo
     enddo
  enddo

  return
end subroutine ox_place

subroutine axis_angle(r,axis,ang)
  implicit none
  ! -----------------------------------------------------------------
  ! Axis angle rotation of vector r about axis by ang
  ! -----------------------------------------------------------------
  integer k
  real(8) r(3),axis(3),xr(3),ang,sang,cang,fac,dot

  cang = dcos(ang)
  sang = dsin(ang)
  call xproduct(r,axis,xr)
  call dotp(r,axis,dot)

  fac = dot*(1.d0-cang)
  
  do k = 1,3
     r(k) = r(k)*cang + axis(k)*fac + xr(k)*sang
  enddo

  return
end subroutine axis_angle

subroutine dotp(r1,r2,dot)
  implicit none
  ! -----------------------------------------------------------------
  ! Dot product of two vectors r1 and r2
  ! -----------------------------------------------------------------
  real(8) r1(3),r2(3),dot
  
  dot = r1(1)*r2(1) + r1(2)*r2(2) + r1(3)*r2(3)
  
  return
end subroutine dotp

subroutine xproduct(r1,r2,xr)
  implicit none
  ! -----------------------------------------------------------------
  ! Cross product of two vectors: r1 x r2 = xr
  ! -----------------------------------------------------------------
  real(8) r1(3),r2(3),xr(3)
  
  xr(1) = r1(2)*r2(3) - r2(2)*r1(3)
  xr(2) = r2(1)*r1(3) - r1(1)*r2(3)
  xr(3) = r1(1)*r2(2) - r2(1)*r1(2)
  
  return
end subroutine xproduct




