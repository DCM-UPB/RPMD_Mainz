subroutine rwald_smeared(n,mol,q,r,v,vir,dvdr,rcut,alpha, &
                         mcell,ncell,wid)
  implicit none
  ! ------------------------------------------------------------------
  ! Cubic real space ewald sum for SMEARED CHARGES
  ! ------------------------------------------------------------------
  integer n,mcell,ncell,nn,i,j,k,ix,iy,iz,jx,jy,jz,incx,incy,incz
  integer mol(n),link(ncell,ncell,ncell),list(n),matm(n)
  integer ink(3,62),map(ncell,-2:2)
  real(8) q(n),r(3,n),dvdr(3,n),vir(3,3),wid(3)
  real(8) pi,rfac,sfac,rcut,rcutsq,alpha,dx,dy,dz,drsq,dr
  real(8) x,e,t,erfc,du,dur,qij,dv,dvr,dvx,dvy,dvz,v
  real(8) sfacs,es,wt,xs,ts,erfcs
  real(8) p,a1,a2,a3,a4,a5

  ! parameters in the approximation to erfc(x) [A&S, 7.1.26]

  parameter (p = 3.0525860d0)
  parameter (a1 = 0.254829592d0)
  parameter (a2 = -0.284496736d0)
  parameter (a3 = 1.421413741d0)
  parameter (a4 = -1.453152027d0)
  parameter (a5 = 1.061405429d0)

  ! Create atom list

  do i = 1,n,3
     matm(i) = 1
     matm(i+1) = 2
     matm(i+2) = 2
  enddo

  ! construct the linked cell list

  call cell_list_unit(r,list,link,ink,map,nn,n,ncell,mcell,1)

  ! and evaluate the real space sum

  pi = dacos(-1.d0)
  rfac = -2.d0/dsqrt(pi)
  sfac = alpha*rfac
  rcutsq = rcut*rcut

  ! including contributions from within the same cell

  do iz = 1,ncell
     do iy = 1,ncell
        do ix = 1,ncell
           i = link(ix,iy,iz)
           do while (i .gt. 0)
              j = list(i)
              do while (j .gt. 0)
                 dx = r(1,i)-r(1,j)
                 dy = r(2,i)-r(2,j)
                 dz = r(3,i)-r(3,j)
                 dx = dx-nint(dx)
                 dy = dy-nint(dy)
                 dz = dz-nint(dz)
                 drsq = dx*dx+dy*dy+dz*dz
                 if (drsq .lt. rcutsq) then
                    dr = dsqrt(drsq)
                    x = alpha*dr
                    e = dexp(-x*x)
                    t = p/(p+x)
                    erfc = e*t*(a1+t*(a2+t*(a3+t*(a4+t*a5))))
                    if (mol(i) .eq. mol(j)) then
                       
                       ! No intramolecular interactions
                       
                       erfc = erfc-1.d0
                       sfacs = 0.d0
                       es = 0.d0
                    else

                       ! Smeared charge

                       if (matm(i).eq.matm(j)) then
                          ! M-M or H-H interaction
                          wt = wid(matm(i))
                          xs = dr * wt
                       else
                          ! Cross Interaction (H-M)
                          wt = wid(3)
                          xs = dr * wt
                       endif
                       sfacs = wt*rfac
                       es = dexp(-xs*xs)
                       ts = p/(p+xs)
                       erfcs = es*ts*(a1+ts*(a2+ts*(a3+ts* &
                                         (a4+ts*a5))))
                       erfc = erfc - erfcs
                    endif
                    du = erfc/dr
                    dur = (sfac*e-sfacs*es-du)/drsq
                    qij = q(i)*q(j)
                    dv = qij*du
                    dvr = qij*dur
                    dvx = dvr*dx
                    dvy = dvr*dy
                    dvz = dvr*dz
                    v = v + dv
                    dvdr(1,i) = dvdr(1,i)+dvx
                    dvdr(2,i) = dvdr(2,i)+dvy
                    dvdr(3,i) = dvdr(3,i)+dvz
                    dvdr(1,j) = dvdr(1,j)-dvx
                    dvdr(2,j) = dvdr(2,j)-dvy
                    dvdr(3,j) = dvdr(3,j)-dvz
                    vir(1,1) = vir(1,1) + dx * dvx
                    vir(1,2) = vir(1,2) + dx * dvy
                    vir(1,3) = vir(1,3) + dx * dvz
                    vir(2,1) = vir(2,1) + dy * dvx
                    vir(2,2) = vir(2,2) + dy * dvy
                    vir(2,3) = vir(2,3) + dy * dvz
                    vir(3,1) = vir(3,1) + dz * dvx
                    vir(3,2) = vir(3,2) + dz * dvy
                    vir(3,3) = vir(3,3) + dz * dvz
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
     do iz = 1,ncell
        jz = map(iz,incz)
        do iy = 1,ncell
           jy = map(iy,incy)
           do ix = 1,ncell
              jx = map(ix,incx)
              i = link(ix,iy,iz)
              do while (i .gt. 0)
                 j = link(jx,jy,jz)
                 do while (j .gt. 0)
                    dx = r(1,i)-r(1,j)
                    dy = r(2,i)-r(2,j)
                    dz = r(3,i)-r(3,j)
                    dx = dx-nint(dx)
                    dy = dy-nint(dy)
                    dz = dz-nint(dz)
                    drsq = dx*dx+dy*dy+dz*dz
                    if (drsq .lt. rcutsq) then
                       dr = dsqrt(drsq)
                       x = alpha*dr
                       e = dexp(-x*x)
                       t = p/(p+x)
                       erfc = e*t*(a1+t*(a2+t*(a3+t*(a4+t*a5))))
                       if (mol(i) .eq. mol(j)) then

                          ! No intramolecular interactions
                          
                          erfc = erfc-1.d0
                          sfacs = 0.d0
                          es = 0.d0
                       else

                          ! Smeared charge

                          if (matm(i).eq.matm(j)) then
                             ! M-M or H-H interaction
                             wt = wid(matm(i))
                             xs = dr * wt
                          else
                             ! Cross Interaction (H-M)
                             wt = wid(3)
                             xs = dr * wt
                          endif
                          sfacs = wt*rfac
                          es = dexp(-xs*xs)
                          ts = p/(p+xs)
                          erfcs = es*ts*(a1+ts*(a2+ts*(a3+ts* &
                                 (a4+ts*a5))))
                          erfc = erfc - erfcs
                       endif
                       du = erfc/dr
                       dur = (sfac*e-sfacs*es-du)/drsq
                       qij = q(i)*q(j)
                       dv = qij*du
                       dvr = qij*dur
                       dvx = dvr*dx
                       dvy = dvr*dy
                       dvz = dvr*dz
                       v = v + dv
                       dvdr(1,i) = dvdr(1,i)+dvx
                       dvdr(2,i) = dvdr(2,i)+dvy
                       dvdr(3,i) = dvdr(3,i)+dvz
                       dvdr(1,j) = dvdr(1,j)-dvx
                       dvdr(2,j) = dvdr(2,j)-dvy
                       dvdr(3,j) = dvdr(3,j)-dvz
                       vir(1,1) = vir(1,1) + dx * dvx
                       vir(1,2) = vir(1,2) + dx * dvy
                       vir(1,3) = vir(1,3) + dx * dvz
                       vir(2,1) = vir(2,1) + dy * dvx
                       vir(2,2) = vir(2,2) + dy * dvy
                       vir(2,3) = vir(2,3) + dy * dvz
                       vir(3,1) = vir(3,1) + dz * dvx
                       vir(3,2) = vir(3,2) + dz * dvy
                       vir(3,3) = vir(3,3) + dz * dvz
                    endif
                    j = list(j)
                 enddo
                 i = list(i)
              enddo
           enddo
        enddo
     enddo
  enddo

  return
end subroutine rwald_smeared

subroutine rwald_smeared_nc(n,mol,q,r,v,vir,dvdr,rcut,alpha, &
                            boxlxyz,wid)
  implicit none
  ! ------------------------------------------------------------------
  ! Non-cubic real space ewald sum for SMEARED CHARGES
  ! ------------------------------------------------------------------
  integer n,nn,i,j,k,ix,iy,iz,jx,jy,jz,incx,incy,incz
  integer mcellxyz(3),ncellxyz(3),mol(n),list(n),matm(n)
  real(8) q(n),r(3,n),dvdr(3,n),vir(3,3),boxlxyz(3),wid(3)
  real(8) pi,rfac,sfac,rcut,rcutsq,alpha,dx,dy,dz,drsq,dr
  real(8) x,e,t,erfc,du,dur,qij,dv,dvr,dvx,dvy,dvz,v
  real(8) p,a1,a2,a3,a4,a5
  real(8) sfacs,es,wt,xs,ts,erfcs
  real(8) xbox,ybox,zbox,onboxx,onboxy,onboxz
  integer, allocatable :: link(:,:,:),ink(:,:)
  integer, allocatable :: mapx(:,:),mapy(:,:),mapz(:,:)

  ! parameters in the approximation to erfc(x) [A&S, 7.1.26]

  parameter (p = 3.0525860d0)
  parameter (a1 = 0.254829592d0)
  parameter (a2 = -0.284496736d0)
  parameter (a3 = 1.421413741d0)
  parameter (a4 = -1.453152027d0)
  parameter (a5 = 1.061405429d0)

  ! Create atom list

  do i = 1,n,3
     matm(i) = 1
     matm(i+1) = 2
     matm(i+2) = 2
  enddo

  ! construct the linked cell list

  call cell_setup(boxlxyz,rcut,mcellxyz,ncellxyz,nn)

  allocate (link(ncellxyz(1),ncellxyz(2),ncellxyz(3)))
  allocate (ink(3,nn))
  allocate (mapx(ncellxyz(1),-mcellxyz(1):mcellxyz(1)))
  allocate (mapy(ncellxyz(2),-mcellxyz(2):mcellxyz(2)))
  allocate (mapz(ncellxyz(3),-mcellxyz(3):mcellxyz(3)))

  call cell_list(r,list,link,ink,mapx,mapy,mapz,nn, &
                 n,boxlxyz,ncellxyz,mcellxyz,1)

  ! and evaluate the real space sum

  pi = dacos(-1.d0)
  rfac = -2.d0/dsqrt(pi)
  sfac = alpha*rfac
  rcutsq = rcut*rcut

  xbox = boxlxyz(1)
  ybox = boxlxyz(2)
  zbox = boxlxyz(3)
  onboxx = 1.d0/xbox
  onboxy = 1.d0/ybox
  onboxz = 1.d0/zbox

  ! including contributions from within the same cell

  do iz = 1,ncellxyz(3)
     do iy = 1,ncellxyz(2)
        do ix = 1,ncellxyz(1)
           i = link(ix,iy,iz)
           do while (i .gt. 0)
              j = list(i)
              do while (j .gt. 0)
                 dx = r(1,i)-r(1,j)
                 dy = r(2,i)-r(2,j)
                 dz = r(3,i)-r(3,j)
                 dx = dx-xbox*nint(dx*onboxx)
                 dy = dy-ybox*nint(dy*onboxy)
                 dz = dz-zbox*nint(dz*onboxz)
                 drsq = dx*dx+dy*dy+dz*dz
                 if (drsq .lt. rcutsq) then
                    dr = dsqrt(drsq)
                    x = alpha*dr
                    e = dexp(-x*x)
                    t = p/(p+x)
                    erfc = e*t*(a1+t*(a2+t*(a3+t*(a4+t*a5))))
                    if (mol(i) .eq. mol(j)) then

                       ! No intramolecular interactions

                       erfc = erfc-1.d0
                       sfacs = 0.d0
                       es = 0.d0
                    else

                       ! Smeared charge

                       if (matm(i).eq.matm(j)) then
                          ! M-M or H-H interaction
                          wt = wid(matm(i))
                          xs = dr * wt
                       else
                          ! Cross Interaction (H-M)
                          wt = wid(3)
                          xs = dr * wt
                       endif
                       sfacs = wt*rfac
                       es = dexp(-xs*xs)
                       ts = p/(p+xs)
                       erfcs = es*ts*(a1+ts*(a2+ts*(a3+ts* &
                                         (a4+ts*a5))))
                       erfc = erfc - erfcs
                    endif
                    du = erfc/dr
                    dur = (sfac*e-sfacs*es-du)/drsq
                    qij = q(i)*q(j)
                    dv = qij*du
                    dvr = qij*dur
                    dvx = dvr*dx
                    dvy = dvr*dy
                    dvz = dvr*dz
                    v = v+dv
                    dvdr(1,i) = dvdr(1,i)+dvx
                    dvdr(2,i) = dvdr(2,i)+dvy
                    dvdr(3,i) = dvdr(3,i)+dvz
                    dvdr(1,j) = dvdr(1,j)-dvx
                    dvdr(2,j) = dvdr(2,j)-dvy
                    dvdr(3,j) = dvdr(3,j)-dvz
                    vir(1,1) = vir(1,1) + dx * dvx
                    vir(1,2) = vir(1,2) + dx * dvy
                    vir(1,3) = vir(1,3) + dx * dvz
                    vir(2,1) = vir(2,1) + dy * dvx
                    vir(2,2) = vir(2,2) + dy * dvy
                    vir(2,3) = vir(2,3) + dy * dvz
                    vir(3,1) = vir(3,1) + dz * dvx
                    vir(3,2) = vir(3,2) + dz * dvy
                    vir(3,3) = vir(3,3) + dz * dvz
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
                    dx = r(1,i)-r(1,j)
                    dy = r(2,i)-r(2,j)
                    dz = r(3,i)-r(3,j)
                    dx = dx-xbox*nint(dx*onboxx)
                    dy = dy-ybox*nint(dy*onboxy)
                    dz = dz-zbox*nint(dz*onboxz)
                    drsq = dx*dx+dy*dy+dz*dz
                    if (drsq .lt. rcutsq) then
                       dr = dsqrt(drsq)
                       x = alpha*dr
                       e = dexp(-x*x)
                       t = p/(p+x)
                       erfc = e*t*(a1+t*(a2+t*(a3+t*(a4+t*a5))))
                       if (mol(i) .eq. mol(j)) then

                          ! No intramolecular interactions

                          erfc = erfc-1.d0
                          sfacs = 0.d0
                          es = 0.d0
                       else

                          ! Smeared charge
                          
                          if (matm(i).eq.matm(j)) then
                             ! M-M or H-H interaction
                             wt = wid(matm(i))
                             xs = dr * wt
                          else
                             ! Cross Interaction (H-M)
                             wt = wid(3)
                             xs = dr * wt
                          endif
                          sfacs = wt*rfac
                          es = dexp(-xs*xs)
                          ts = p/(p+xs)
                          erfcs = es*ts*(a1+ts*(a2+ts*(a3+ts* &
                                           (a4+ts*a5))))
                          erfc = erfc - erfcs
                       endif
                       du = erfc/dr
                       dur = (sfac*e-sfacs*es-du)/drsq
                       qij = q(i)*q(j)
                       dv = qij*du
                       dvr = qij*dur
                       dvx = dvr*dx
                       dvy = dvr*dy
                       dvz = dvr*dz
                       v = v+dv
                       dvdr(1,i) = dvdr(1,i)+dvx
                       dvdr(2,i) = dvdr(2,i)+dvy
                       dvdr(3,i) = dvdr(3,i)+dvz
                       dvdr(1,j) = dvdr(1,j)-dvx
                       dvdr(2,j) = dvdr(2,j)-dvy
                       dvdr(3,j) = dvdr(3,j)-dvz
                       vir(1,1) = vir(1,1) + dx * dvx
                       vir(1,2) = vir(1,2) + dx * dvy
                       vir(1,3) = vir(1,3) + dx * dvz
                       vir(2,1) = vir(2,1) + dy * dvx
                       vir(2,2) = vir(2,2) + dy * dvy
                       vir(2,3) = vir(2,3) + dy * dvz
                       vir(3,1) = vir(3,1) + dz * dvx
                       vir(3,2) = vir(3,2) + dz * dvy
                       vir(3,3) = vir(3,3) + dz * dvz
                    endif
                    j = list(j)
                 enddo
                 i = list(i)
              enddo
           enddo
        enddo
     enddo
  enddo

  return
end subroutine rwald_smeared_nc
