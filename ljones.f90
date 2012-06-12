subroutine lj_driver(r,dvdr,v,vir,list,point,na,boxlxyz,njump)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Lennard-Jones Driver
  ! ------------------------------------------------------------------
  integer na,point(na+3),list(maxnab*na),njump
  real(8) r(3,na),dvdr(3,na),vir(3,3),boxlxyz(3)
  real(8) v,oo_eps,oo_sig,oo_gam,rcut,boxmax
  common /oo_param/ oo_eps,oo_sig,oo_gam,rcut

  ! Cubic

  if (list(1).ne.0) then
     call lj_list(r,dvdr,v,vir,na,boxlxyz,list,point,njump)
  else
     boxmax = maxval(boxlxyz)
     if (3.d0*rcut .gt. boxmax) then
        call lj_basic(r,dvdr,v,vir,na,boxlxyz,njump)
     else
        
        ! Use linked cell list
        
        call lj_cell(r,v,vir,dvdr,na,boxlxyz,njump)
     endif
  endif
  
  return
end subroutine lj_driver


subroutine lj_basic(r,dvdr,v,vir,na,boxlxyz,njump)
  implicit none
  ! ----------------------------------------------------------------
  ! Lennard Jones - General Version
  ! ----------------------------------------------------------------
  integer na,i,j,njump
  real(8) r(3,na),dvdr(3,na),vir(3,3),boxlxyz(3),v
  real(8) onboxx,onboxy,onboxz,oo_eps,oo_sig,oo_gam,rcut
  real(8) sigsq,rcutsq,boxx,boxy,boxz,vij
  real(8) drsq,onr2,fij,dfx,dfy,dfz,sr2,sr6,wij
  real(8) dx,dy,dz,vscale,dscale
  common /oo_param/ oo_eps,oo_sig,oo_gam,rcut

  sigsq = oo_sig*oo_sig
  rcutsq = rcut*rcut
  boxx = boxlxyz(1)
  boxy = boxlxyz(2)
  boxz = boxlxyz(3)
  onboxx = 1.d0/boxx
  onboxy = 1.d0/boxy
  onboxz = 1.d0/boxz

  v = 0.d0
  vir(:,:) = 0.d0
  dvdr(:,:) = 0.d0

  do j = 1+njump,na,njump
     do i = 1,j-njump,njump
        dx = r(1,i)-r(1,j)
        dy = r(2,i)-r(2,j)
        dz = r(3,i)-r(3,j)
        dx = dx - boxx*dble(nint(onboxx*dx))
        dy = dy - boxy*dble(nint(onboxy*dy))
        dz = dz - boxz*dble(nint(onboxz*dz))
        drsq = dx*dx + dy*dy + dz*dz
        if (drsq .lt. rcutsq) then
           onr2 = 1.d0/drsq
           sr2 = sigsq * onr2
           sr6 = sr2 * sr2 * sr2
           vij = sr6 * (sr6-1.d0)
           v = v + vij
           wij = sr6 * (sr6-0.5d0)
           fij = wij * onr2
           dfx = fij * dx
           dfy = fij * dy
           dfz = fij * dz
           dvdr(1,i) = dvdr(1,i) - dfx
           dvdr(2,i) = dvdr(2,i) - dfy
           dvdr(3,i) = dvdr(3,i) - dfz
           dvdr(1,j) = dvdr(1,j) + dfx
           dvdr(2,j) = dvdr(2,j) + dfy
           dvdr(3,j) = dvdr(3,j) + dfz
           vir(1,1) = vir(1,1) - dx * dfx
           vir(2,2) = vir(2,2) - dy * dfy
           vir(3,3) = vir(3,3) - dz * dfz
           vir(1,2) = vir(1,2) - dx * dfy
           vir(1,3) = vir(1,3) - dx * dfz
           vir(2,1) = vir(2,1) - dy * dfx
           vir(2,3) = vir(2,3) - dy * dfz
           vir(3,1) = vir(3,1) - dz * dfx
           vir(3,2) = vir(3,2) - dz * dfy
        endif
     enddo
  enddo

  vscale = 4.d0*oo_eps
  dscale = 48.d0*oo_eps

  v = vscale * v
  vir(:,:) = dscale * vir(:,:)
  do i = 1,na,njump
     dvdr(1,i) = dscale*dvdr(1,i)
     dvdr(2,i) = dscale*dvdr(2,i)
     dvdr(3,i) = dscale*dvdr(3,i)
  enddo

  return
end subroutine lj_basic

subroutine lj_list(r,dvdr,v,vir,na,boxlxyz,list,point,njump)
  implicit none
  include 'globals.inc'
  ! ----------------------------------------------------------------
  ! Lennard Jones - Neighbour list version
  ! ----------------------------------------------------------------
  integer na,i,j,point(na+3),list(maxnab*na),ibeg,iend,inab,njump
  real(8) r(3,na),dvdr(3,na),vir(3,3),boxlxyz(3),v
  real(8) oo_eps,oo_sig,oo_gam,rcut
  real(8) sigsq,rcutsq,vij
  real(8) drsq,onr2,fij,dfx,dfy,dfz,sr2,sr6,wij
  real(8) dx,dy,dz,vscale,dscale
  real(8) rxj,ryj,rzj,dvdxj,dvdyj,dvdzj
  real(8) onboxx,onboxy,onboxz,boxx,boxy,boxz
  common /oo_param/ oo_eps,oo_sig,oo_gam,rcut

  ! clear arrays

  v = 0.d0
  vir(:,:) = 0.d0
  dvdr(:,:) = 0.d0

  rcutsq = rcut*rcut
  boxx = boxlxyz(1)
  boxy = boxlxyz(2)
  boxz = boxlxyz(3)
  onboxx = 1.d0/boxx
  onboxy = 1.d0/boxy
  onboxz = 1.d0/boxz

  ! use the neighbour list to calculate energy and forces

  sigsq = oo_sig*oo_sig
  do j = 1,na,njump
     ibeg = point(j)
     iend = point(j+njump) - 1
     if (ibeg.le.iend) then
        rxj = r(1,j)
        ryj = r(2,j)
        rzj = r(3,j)
        dvdxj = dvdr(1,j)
        dvdyj = dvdr(2,j)
        dvdzj = dvdr(3,j)
        do inab = ibeg,iend
           i = list(inab)
           dx = r(1,i) - rxj
           dy = r(2,i) - ryj
           dz = r(3,i) - rzj
           dx = dx - boxx*dble(nint(onboxx*dx))
           dy = dy - boxy*dble(nint(onboxy*dy))
           dz = dz - boxz*dble(nint(onboxz*dz))
           drsq = dx*dx + dy*dy + dz*dz
           if (drsq .lt. rcutsq) then
              onr2 = 1.d0/drsq
              sr2 = sigsq * onr2
              sr6 = sr2 * sr2 * sr2
              vij = sr6 * (sr6-1.d0)
              v = v + vij
              wij = sr6 * (sr6-0.5d0)
              fij = wij * onr2
              dfx = fij * dx
              dfy = fij * dy
              dfz = fij * dz
              dvdr(1,i) = dvdr(1,i) - dfx
              dvdr(2,i) = dvdr(2,i) - dfy
              dvdr(3,i) = dvdr(3,i) - dfz
              dvdxj = dvdxj + dfx
              dvdyj = dvdyj + dfy
              dvdzj = dvdzj + dfz
              vir(1,1) = vir(1,1) - dx * dfx
              vir(2,2) = vir(2,2) - dy * dfy
              vir(3,3) = vir(3,3) - dz * dfz
              vir(1,2) = vir(1,2) - dx * dfy
              vir(1,3) = vir(1,3) - dx * dfz
              vir(2,1) = vir(2,1) - dy * dfx
              vir(2,3) = vir(2,3) - dy * dfz
              vir(3,1) = vir(3,1) - dz * dfx
              vir(3,2) = vir(3,2) - dz * dfy
           endif
        enddo
        dvdr(1,j) = dvdxj
        dvdr(2,j) = dvdyj
        dvdr(3,j) = dvdzj
     endif
  enddo

  vscale = 4.d0*oo_eps
  dscale = 48.d0*oo_eps

  v = vscale * v
  vir(:,:) = dscale * vir(:,:)
  do i = 1,na,njump
     dvdr(1,i) = dscale*dvdr(1,i)
     dvdr(2,i) = dscale*dvdr(2,i)
     dvdr(3,i) = dscale*dvdr(3,i)
  enddo

  return
end subroutine lj_list

subroutine lj_cell(r,v,vir,dvdr,na,boxlxyz,njump)
  implicit none
  ! ------------------------------------------------------------------
  ! Lennard Jones - Linked Cell Version
  ! ------------------------------------------------------------------
  integer na,ix,iy,iz,k,i,j,incx,incy,incz,njump
  integer jx,jy,jz,nn,ncellxyz(3),mcellxyz(3)
  integer list(na)
  real(8) r(3,na),dvdr(3,na),vir(3,3),boxlxyz(3)
  real(8) oo_eps,oo_sig,oo_gam,rcut,rcutsq,vij,sigsq,v
  real(8) drsq,onr2,fij,dfx,dfy,dfz,sr2,sr6,wij
  real(8) dx,dy,dz,vscale,dscale
  real(8) boxx,boxy,boxz,onboxx,onboxy,onboxz
  integer, allocatable :: link(:,:,:),ink(:,:)
  integer, allocatable :: mapx(:,:),mapy(:,:),mapz(:,:)
  common /oo_param/ oo_eps,oo_sig,oo_gam,rcut

  ! generate linked cell list

  call cell_setup(boxlxyz,rcut,mcellxyz,ncellxyz,nn)

  allocate (link(ncellxyz(1),ncellxyz(2),ncellxyz(3)))
  allocate (ink(3,nn))
  allocate (mapx(ncellxyz(1),-mcellxyz(1):mcellxyz(1)))
  allocate (mapy(ncellxyz(2),-mcellxyz(2):mcellxyz(2)))
  allocate (mapz(ncellxyz(3),-mcellxyz(3):mcellxyz(3)))

  call cell_list(r,list,link,ink,mapx,mapy,mapz,nn, &
                 na,boxlxyz,ncellxyz,mcellxyz,njump)

  ! and evaluate the potential energy and forces

  boxx = boxlxyz(1)
  boxy = boxlxyz(2)
  boxz = boxlxyz(3)

  onboxx = 1.d0/boxx
  onboxy = 1.d0/boxy
  onboxz = 1.d0/boxz

  v = 0.d0
  vir(:,:) = 0.d0
  dvdr(:,:) = 0.d0

  rcutsq = rcut*rcut
  sigsq = oo_sig*oo_sig

  do iz = 1,ncellxyz(3)
     do iy = 1,ncellxyz(2)
        do ix = 1,ncellxyz(1)
           
           ! including contributions from within the same cell
           
           i = link(ix,iy,iz)
           do while (i .gt. 0)
              j = list(i)
              do while (j .gt. 0)
                 dx = r(1,i)-r(1,j)
                 dy = r(2,i)-r(2,j)
                 dz = r(3,i)-r(3,j)
                 dx = dx-boxx*dble(nint(onboxx*dx))
                 dy = dy-boxy*dble(nint(onboxy*dy))
                 dz = dz-boxz*dble(nint(onboxz*dz))
                 drsq = dx*dx + dy*dy + dz*dz
                 if (drsq .lt. rcutsq) then
                    onr2 = 1.d0/drsq
                    sr2 = sigsq*onr2
                    sr6 = sr2*sr2*sr2
                    vij = sr6*(sr6-1.d0)
                    v = v + vij
                    wij = sr6*(sr6-0.5d0)
                    fij = wij*onr2
                    dfx = fij*dx
                    dfy = fij*dy
                    dfz = fij*dz
                    dvdr(1,i) = dvdr(1,i) - dfx
                    dvdr(2,i) = dvdr(2,i) - dfy
                    dvdr(3,i) = dvdr(3,i) - dfz
                    dvdr(1,j) = dvdr(1,j) + dfx
                    dvdr(2,j) = dvdr(2,j) + dfy
                    dvdr(3,j) = dvdr(3,j) + dfz
                    vir(1,1) = vir(1,1) - dx * dfx
                    vir(2,2) = vir(2,2) - dy * dfy
                    vir(3,3) = vir(3,3) - dz * dfz
                    vir(1,2) = vir(1,2) - dx * dfy
                    vir(1,3) = vir(1,3) - dx * dfz
                    vir(2,1) = vir(2,1) - dy * dfx
                    vir(2,3) = vir(2,3) - dy * dfz
                    vir(3,1) = vir(3,1) - dz * dfx
                    vir(3,2) = vir(3,2) - dz * dfy
                 endif
                 j = list(j)
              enddo
              i = list(i)
           enddo
           
           ! neigbouring cells

           do k = 1,nn
              incx = ink(1,k)
              incy = ink(2,k)
              incz = ink(3,k)
              jx = mapx(ix,incx)
              jy = mapy(iy,incy)
              jz = mapz(iz,incz)

              i = link(ix,iy,iz)
              do while (i .gt. 0)
                 j = link(jx,jy,jz)
                 do while (j .gt. 0)
                    dx = r(1,i)-r(1,j)
                    dy = r(2,i)-r(2,j)
                    dz = r(3,i)-r(3,j)
                    dx = dx-boxx*dble(nint(onboxx*dx))
                    dy = dy-boxy*dble(nint(onboxy*dy))
                    dz = dz-boxz*dble(nint(onboxz*dz))
                    drsq = dx*dx+dy*dy+dz*dz
                    if (drsq .lt. rcutsq) then
                       onr2 = 1.d0/drsq
                       sr2 = sigsq*onr2
                       sr6 = sr2*sr2*sr2
                       vij = sr6*(sr6-1.d0)
                       v = v + vij
                       wij = sr6*(sr6-0.5d0)
                       fij = wij*onr2
                       dfx = fij*dx
                       dfy = fij*dy
                       dfz = fij*dz
                       dvdr(1,i) = dvdr(1,i) - dfx
                       dvdr(2,i) = dvdr(2,i) - dfy
                       dvdr(3,i) = dvdr(3,i) - dfz
                       dvdr(1,j) = dvdr(1,j) + dfx
                       dvdr(2,j) = dvdr(2,j) + dfy
                       dvdr(3,j) = dvdr(3,j) + dfz
                       vir(1,1) = vir(1,1) - dx * dfx
                       vir(2,2) = vir(2,2) - dy * dfy
                       vir(3,3) = vir(3,3) - dz * dfz
                       vir(1,2) = vir(1,2) - dx * dfy
                       vir(1,3) = vir(1,3) - dx * dfz
                       vir(2,1) = vir(2,1) - dy * dfx
                       vir(2,3) = vir(2,3) - dy * dfz
                       vir(3,1) = vir(3,1) - dz * dfx
                       vir(3,2) = vir(3,2) - dz * dfy
                    endif
                    j = list(j)
                 enddo
                 i = list(i)
              enddo
           enddo
        enddo
     enddo
  enddo
  
  ! common factors

  vscale = 4.d0*oo_eps
  dscale = 48.d0*oo_eps
  v = vscale*v
  vir(:,:) = dscale*vir(:,:)
  do i = 1,na,njump
     dvdr(1,i) = dscale*dvdr(1,i)
     dvdr(2,i) = dscale*dvdr(2,i)
     dvdr(3,i) = dscale*dvdr(3,i)
  enddo

  deallocate (link,ink,mapx,mapy,mapz)

  return
end subroutine lj_cell
