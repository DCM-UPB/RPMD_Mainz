subroutine buck_driver(r,dvdr,v,vir,list,point,na,boxlxyz,njump)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Buckingham exp 6 driver
  ! ------------------------------------------------------------------
  integer na,point(na+3),list(maxnab*na),njump
  real(8) r(3,na),dvdr(3,na),vir(3,3),boxlxyz(3)
  real(8) v,oo_eps,oo_sig,oo_gam,rcut,boxmax
  common /oo_param/ oo_eps,oo_sig,oo_gam,rcut

  ! Cubic

  if (list(1).ne.0) then
     call buck_list(r,dvdr,v,vir,na,boxlxyz,list,point,njump)
  else
     boxmax = maxval(boxlxyz)
     if (3.d0*rcut .gt. boxmax) then
        call buck_basic(r,dvdr,v,vir,na,boxlxyz,njump)
     else

        ! Use linked cell list

        call buck_cell(r,v,vir,dvdr,na,boxlxyz,njump)
     endif
  endif

  return
end subroutine buck_driver

subroutine buck_basic(r,dvdr,v,vir,na,boxlxyz,njump)
  implicit none
  ! ----------------------------------------------------------------
  ! Buckingham 6 potential - General Version
  ! ----------------------------------------------------------------
  integer na,i,j,njump
  real(8) r(3,na),dvdr(3,na),vir(3,3),boxlxyz(3),v
  real(8) onboxx,onboxy,onboxz,oo_eps,oo_sig,oo_gam,rcut
  real(8) rcutsq,boxx,boxy,boxz,vij,drsq,fij,dfx,dfy,dfz
  real(8) a,b,c,br6,exp6,dr,onr
  real(8) dx,dy,dz
  common /oo_param/ oo_eps,oo_sig,oo_gam,rcut

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

  a = dexp(oo_gam)*6.0d0*oo_eps/(oo_gam-6.d0)
  b = oo_eps*(oo_sig**6)/(1.d0-(6.d0/oo_gam))
  c = -oo_gam/oo_sig

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
           br6 = b/(drsq*drsq*drsq)
           dr = dsqrt(drsq)
           exp6 = a*dexp(c*dr)
           vij = exp6 - br6
           v = v + vij
           onr = 1.d0/dr
           fij = c*exp6 + 6.d0*br6*onr
           dfx = fij * dx * onr
           dfy = fij * dy * onr
           dfz = fij * dz * onr
           dvdr(1,i) = dvdr(1,i) + dfx
           dvdr(2,i) = dvdr(2,i) + dfy
           dvdr(3,i) = dvdr(3,i) + dfz
           dvdr(1,j) = dvdr(1,j) - dfx
           dvdr(2,j) = dvdr(2,j) - dfy
           dvdr(3,j) = dvdr(3,j) - dfz
           vir(1,1) = vir(1,1) + dx * dfx
           vir(2,2) = vir(2,2) + dy * dfy
           vir(3,3) = vir(3,3) + dz * dfz
           vir(1,2) = vir(1,2) + dx * dfy
           vir(1,3) = vir(1,3) + dx * dfz
           vir(2,1) = vir(2,1) + dy * dfx
           vir(2,3) = vir(2,3) + dy * dfz
           vir(3,1) = vir(3,1) + dz * dfx
           vir(3,2) = vir(3,2) + dz * dfy
        endif
     enddo
  enddo

  return
end subroutine buck_basic


subroutine buck_list(r,dvdr,v,vir,na,boxlxyz,list,point,njump)
  implicit none
  include 'globals.inc'
  ! ----------------------------------------------------------------
  ! Buckingham exp6 potential - Neighbour list version
  ! ----------------------------------------------------------------
  integer na,i,j,point(na+3),list(maxnab*na),ibeg,iend,inab,njump
  real(8) r(3,na),dvdr(3,na),vir(3,3),boxlxyz(3),v
  real(8) oo_eps,oo_sig,oo_gam,rcut
  real(8) rcutsq,vij,drsq,fij,dfx,dfy,dfz
  real(8) dx,dy,dz
  real(8) a,b,c,br6,exp6,dr,onr
  real(8) rxj,ryj,rzj,dvdxj,dvdyj,dvdzj
  real(8) onboxx,onboxy,onboxz,boxx,boxy,boxz
  common /oo_param/ oo_eps,oo_sig,oo_gam,rcut

  ! clear arrays

  v = 0.d0
  vir(:,:) = 0.d0
  dvdr(:,:) = 0.d0

  a = dexp(oo_gam)*6.0d0*oo_eps/(oo_gam-6.d0)
  b = oo_eps*(oo_sig**6)/(1.d0-(6.d0/oo_gam))
  c = -oo_gam/oo_sig
  rcutsq = rcut*rcut

  boxx = boxlxyz(1)
  boxy = boxlxyz(2)
  boxz = boxlxyz(3)
  onboxx = 1.d0/boxx
  onboxy = 1.d0/boxy
  onboxz = 1.d0/boxz

  ! use the neighbour list to calculate energy and forces

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
              br6 = b/(drsq*drsq*drsq)
              dr = dsqrt(drsq)
              exp6 = a*dexp(c*dr)
              vij = exp6 - br6
              v = v + vij
              onr = 1.d0/dr
              fij = c*exp6 + 6.d0*br6*onr
              dfx = fij * dx * onr
              dfy = fij * dy * onr
              dfz = fij * dz * onr
              dvdr(1,i) = dvdr(1,i) + dfx
              dvdr(2,i) = dvdr(2,i) + dfy
              dvdr(3,i) = dvdr(3,i) + dfz
              dvdxj = dvdxj - dfx
              dvdyj = dvdyj - dfy
              dvdzj = dvdzj - dfz
              vir(1,1) = vir(1,1) + dx * dfx
              vir(2,2) = vir(2,2) + dy * dfy
              vir(3,3) = vir(3,3) + dz * dfz
              vir(1,2) = vir(1,2) + dx * dfy
              vir(1,3) = vir(1,3) + dx * dfz
              vir(2,1) = vir(2,1) + dy * dfx
              vir(2,3) = vir(2,3) + dy * dfz
              vir(3,1) = vir(3,1) + dz * dfx
              vir(3,2) = vir(3,2) + dz * dfy
           endif
        enddo
        dvdr(1,j) = dvdxj
        dvdr(2,j) = dvdyj
        dvdr(3,j) = dvdzj
     endif
  enddo
  
  return
end subroutine buck_list


subroutine buck_cell(r,v,vir,dvdr,na,boxlxyz,njump)
  implicit none
  ! ------------------------------------------------------------------
  ! Buckingham exp6 potential - Linked Cell Version
  ! ------------------------------------------------------------------
  integer na,ix,iy,iz,k,i,j,incx,incy,incz,njump
  integer jx,jy,jz,nn,ncellxyz(3),mcellxyz(3)
  integer list(na)
  real(8) r(3,na),dvdr(3,na),vir(3,3),boxlxyz(3)
  real(8) oo_eps,oo_sig,oo_gam,rcut,rcutsq,vij,v
  real(8) drsq,dx,dy,dz,fij,dfx,dfy,dfz
  real(8) a,b,c,br6,exp6,dr,onr
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

  a = dexp(oo_gam)*6.0d0*oo_eps/(oo_gam-6.d0)
  b = oo_eps*(oo_sig**6)/(1.d0-(6.d0/oo_gam))
  c = -oo_gam/oo_sig
  rcutsq = rcut*rcut

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
                    br6 = b/(drsq*drsq*drsq)
                    dr = dsqrt(drsq)
                    exp6 = a*dexp(c*dr)
                    vij = exp6 - br6
                    v = v + vij
                    onr = 1.d0/dr
                    fij = c*exp6 + 6.d0*br6*onr
                    dfx = fij * dx * onr
                    dfy = fij * dy * onr
                    dfz = fij * dz * onr
                    dvdr(1,i) = dvdr(1,i) + dfx
                    dvdr(2,i) = dvdr(2,i) + dfy
                    dvdr(3,i) = dvdr(3,i) + dfz
                    dvdr(1,j) = dvdr(1,j) - dfx
                    dvdr(2,j) = dvdr(2,j) - dfy
                    dvdr(3,j) = dvdr(3,j) - dfz
                    vir(1,1) = vir(1,1) + dx * dfx
                    vir(2,2) = vir(2,2) + dy * dfy
                    vir(3,3) = vir(3,3) + dz * dfz
                    vir(1,2) = vir(1,2) + dx * dfy
                    vir(1,3) = vir(1,3) + dx * dfz
                    vir(2,1) = vir(2,1) + dy * dfx
                    vir(2,3) = vir(2,3) + dy * dfz
                    vir(3,1) = vir(3,1) + dz * dfx
                    vir(3,2) = vir(3,2) + dz * dfy
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
                       br6 = b/(drsq*drsq*drsq)
                       dr = dsqrt(drsq)
                       exp6 = a*dexp(c*dr)
                       vij = exp6 - br6
                       v = v + vij
                       onr = 1.d0/dr
                       fij = c*exp6 + 6.d0*br6*onr
                       dfx = fij * dx * onr
                       dfy = fij * dy * onr
                       dfz = fij * dz * onr
                       dvdr(1,i) = dvdr(1,i) + dfx
                       dvdr(2,i) = dvdr(2,i) + dfy
                       dvdr(3,i) = dvdr(3,i) + dfz
                       dvdr(1,j) = dvdr(1,j) - dfx
                       dvdr(2,j) = dvdr(2,j) - dfy
                       dvdr(3,j) = dvdr(3,j) - dfz
                       vir(1,1) = vir(1,1) + dx * dfx
                       vir(2,2) = vir(2,2) + dy * dfy
                       vir(3,3) = vir(3,3) + dz * dfz
                       vir(1,2) = vir(1,2) + dx * dfy
                       vir(1,3) = vir(1,3) + dx * dfz
                       vir(2,1) = vir(2,1) + dy * dfx
                       vir(2,3) = vir(2,3) + dy * dfz
                       vir(3,1) = vir(3,1) + dz * dfx
                       vir(3,2) = vir(3,2) + dz * dfy
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
end subroutine buck_cell
