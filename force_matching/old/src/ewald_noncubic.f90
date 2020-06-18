subroutine ewald_nc(n,mol,q,r,v,vir,dvdr,boxlxyz)
  implicit none
  ! ------------------------------------------------------------------
  ! Non-Cubic Ewald sum
  ! ------------------------------------------------------------------
  common /ew_param/ alpha,ecut,wm,wh
  common /ew_setup/ wrcut,walpha,rkmax,kmax
  integer n,kmax,mol(n)
  real(8) boxlxyz(3),dvdr(3,n),r(3,n),q(n),vir(3,3),v
  real(8) wm,wh,ecut,alpha,wrcut,walpha,rkmax,rkmax2,boxmax
  real(8) wid(3),sm,sh

  v = 0.d0
  vir(:,:) = 0.d0

  ! Real space Ewald Sum

  if ((wm.gt.0.d0).or.(wh.gt.0.d0)) then

     ! Smeared Charge Ewald sum

     ! Values of sigma and width list

     sm = dsqrt(1.d0/(2.d0*wm))
     sh = dsqrt(1.d0/(2.d0*wh))

     wid(1) = 1.d0/(2.d0*sm)
     wid(2) = 1.d0/(2.d0*sh)
     wid(3) = 1.d0/(dsqrt(2.d0*(sm*sm + sh*sh)))

     call rwald_smeared_nc(n,mol,q,r,v,vir,dvdr,wrcut,walpha, &
                           boxlxyz,wid)

  else

     ! Point Charge Ewald sum

     boxmax = maxval(boxlxyz)
     if (3.d0*wrcut .gt. boxmax) then
        call rwald_basic(n,mol,q,r,v,vir,dvdr,wrcut,walpha,boxlxyz)
     else
        call rwald_cell(n,mol,q,r,v,vir,dvdr,wrcut,walpha,boxlxyz)
     endif
  endif

  ! Reciprocal Space Ewald Sum

  rkmax2 = rkmax*rkmax
  call kwald_nc(r,q,n,v,vir,dvdr,walpha,rkmax2,kmax,boxlxyz)

  return
end subroutine ewald_nc

subroutine rwald_cell(n,mol,q,r,v,vir,dvdr,rcut,alpha,boxlxyz)
  implicit none
  ! ------------------------------------------------------------------
  ! Real space part of the Ewald sum for non-cubic systems
  ! Low precision version with an approximation to erfc(x)
  ! ------------------------------------------------------------------
  integer n,nn,i,j,k,ix,iy,iz,jx,jy,jz,incx,incy,incz
  integer mcellxyz(3),ncellxyz(3),mol(n),list(n)
  real(8) q(n),r(3,n),dvdr(3,n),vir(3,3),boxlxyz(3)
  real(8) pi,rfac,sfac,rcut,rcutsq,alpha,dx,dy,dz,drsq,dr
  real(8) x,e,t,erfc,du,dur,qij,dv,dvr,dvx,dvy,dvz,v
  real(8) p,a1,a2,a3,a4,a5
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
  rfac = alpha/dsqrt(pi)
  sfac = -2.d0*rfac
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
                 dx = dx-xbox*dble(nint(dx*onboxx))
                 dy = dy-ybox*dble(nint(dy*onboxy))
                 dz = dz-zbox*dble(nint(dz*onboxz))
                 drsq = dx*dx+dy*dy+dz*dz
                 if (drsq .lt. rcutsq) then
                    dr = dsqrt(drsq)
                    x = alpha*dr
                    e = dexp(-x*x)
                    t = p/(p+x)
                    erfc = e*t*(a1+t*(a2+t*(a3+t*(a4+t*a5))))
                    if (mol(i) .eq. mol(j)) erfc = erfc-1.d0
                    du = erfc/dr
                    dur = (sfac*e-du)/drsq
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
                    dx = dx-xbox*dble(nint(dx*onboxx))
                    dy = dy-ybox*dble(nint(dy*onboxy))
                    dz = dz-zbox*dble(nint(dz*onboxz))
                    drsq = dx*dx+dy*dy+dz*dz
                    if (drsq .lt. rcutsq) then
                       dr = dsqrt(drsq)
                       x = alpha*dr
                       e = dexp(-x*x)
                       t = p/(p+x)
                       erfc = e*t*(a1+t*(a2+t*(a3+t*(a4+t*a5))))
                       if (mol(i) .eq. mol(j)) erfc = erfc-1.d0
                       du = erfc/dr
                       dur = (sfac*e-du)/drsq
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
end subroutine rwald_cell

subroutine rwald_basic(n,mol,q,r,v,vir,dvdr,rcut,alpha,boxlxyz)
  implicit none
  ! ------------------------------------------------------------------
  ! Real space part of the Ewald sum for non-cubic systems
  ! Low precision version with an approximation to erfc(x)
  ! ------------------------------------------------------------------
  integer n,i,j
  integer mol(n)
  real(8) q(n),r(3,n),dvdr(3,n),vir(3,3),boxlxyz(3)
  real(8) xbox,ybox,zbox,onboxx,onboxy,onboxz
  real(8) pi,rfac,sfac,rcut,rcutsq,alpha,dx,dy,dz,drsq,dr
  real(8) x,e,t,erfc,du,dur,qij,dv,dvr,dvx,dvy,dvz,v
  real(8) p,a1,a2,a3,a4,a5

  ! parameters in the approximation to erfc(x) [A&S, 7.1.26]

  parameter (p = 3.0525860d0)
  parameter (a1 = 0.254829592d0)
  parameter (a2 = -0.284496736d0)
  parameter (a3 = 1.421413741d0)
  parameter (a4 = -1.453152027d0)
  parameter (a5 = 1.061405429d0)

  ! and evaluate the real space sum

  pi = dacos(-1.d0)
  rfac = alpha/dsqrt(pi)
  sfac = -2.d0*rfac
  rcutsq = rcut*rcut

  xbox = boxlxyz(1)
  ybox = boxlxyz(2)
  zbox = boxlxyz(3)
  onboxx = 1.d0/xbox
  onboxy = 1.d0/ybox
  onboxz = 1.d0/zbox

  ! including contributions from within the same cell

  do i = 1,n
     do j = 1,i-1
        dx = r(1,i)-r(1,j)
        dy = r(2,i)-r(2,j)
        dz = r(3,i)-r(3,j)
        dx = dx-xbox*dble(nint(dx*onboxx))
        dy = dy-ybox*dble(nint(dy*onboxy))
        dz = dz-zbox*dble(nint(dz*onboxz))
        drsq = dx*dx+dy*dy+dz*dz
        if (drsq .lt. rcutsq) then
           dr = dsqrt(drsq)
           x = alpha*dr
           e = dexp(-x*x)
           t = p/(p+x)
           erfc = e*t*(a1+t*(a2+t*(a3+t*(a4+t*a5))))
           if (mol(i) .eq. mol(j)) erfc = erfc-1.d0
           du = erfc/dr
           dur = (sfac*e-du)/drsq
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
     enddo
  enddo  
  return
end subroutine rwald_basic

subroutine kwald_nc(r,z,n,v,vir,dvdr,alpha,rkmax2,kmax,boxlxyz)
  implicit none
  ! ------------------------------------------------------------------
  ! Reciprocal space part of the Ewald sum for non-cubic systems
  ! ------------------------------------------------------------------
  integer n,kmax,i,kx,ky,kz
  real(8) r(3,n),z(n),dvdr(3,n),vir(3,3),boxlxyz(3)
  real(8) ckx(n,0:kmax),skx(n,0:kmax)
  real(8) cky(n,-kmax:kmax),sky(n,-kmax:kmax)
  real(8) ckz(n,-kmax:kmax),skz(n,-kmax:kmax)
  real(8) cxy(n),sxy(n),dsdr(6,n)
  real(8) v,alpha,rkmax2,pi,rfac,twopi,xbox,ybox,zbox
  real(8) xlat,ylat,zlat,xi,yi,zi,b,f,rkx,rky,rkz,vs
  real(8) c,s,rk2,sr,si,tr,ti,w,et,xx,yy,zz,xy,xz,yz,term

  pi = dacos(-1.d0)
  rfac = alpha/dsqrt(pi)

  twopi = 2.d0*pi
  xbox = boxlxyz(1)
  ybox = boxlxyz(2)
  zbox = boxlxyz(3)

  xlat = twopi/xbox
  ylat = twopi/ybox
  zlat = twopi/zbox

  ! setup the trigonometric arrays

  do i = 1,n
     ckx(i,0) = z(i)
     skx(i,0) = 0.d0
     xi = xlat*r(1,i)
     c = dcos(xi)
     s = dsin(xi)
     do kx = 1,kmax
        ckx(i,kx) = c*ckx(i,kx-1) - s*skx(i,kx-1)
        skx(i,kx) = s*ckx(i,kx-1) + c*skx(i,kx-1)
     enddo
  enddo
  do i = 1,n
     cky(i,0) = 1.d0
     sky(i,0) = 0.d0
     yi = ylat*r(2,i)
     c = dcos(yi)
     s = dsin(yi)
     do ky = 1,kmax
        cky(i,ky) = c*cky(i,ky-1) - s*sky(i,ky-1)
        sky(i,ky) = s*cky(i,ky-1) + c*sky(i,ky-1)
        cky(i,-ky) = cky(i,ky)
        sky(i,-ky) = -sky(i,ky)
     enddo
  enddo
  do i = 1,n
     ckz(i,0) = 1.d0
     skz(i,0) = 0.d0
     zi = zlat*r(3,i)
     c = dcos(zi)
     s = dsin(zi)
     do kz = 1,kmax
        ckz(i,kz) = c*ckz(i,kz-1) - s*skz(i,kz-1)
        skz(i,kz) = s*ckz(i,kz-1) + c*skz(i,kz-1)
        ckz(i,-kz) = ckz(i,kz)
        skz(i,-kz) = -skz(i,kz)
     enddo
  enddo

  ! and evaluate the reciprocal space sum

  b = 0.25d0/(alpha*alpha)
  f = twopi/(xbox*ybox*zbox)      ! 2pi / V !
  do kx = 0,kmax
     if (kx .eq. 1) f = 2.d0*f
     rkx = xlat*kx
     do ky = -kmax,kmax
        rky = ylat*ky
        do i = 1,n
           cxy(i) = ckx(i,kx)*cky(i,ky)-skx(i,kx)*sky(i,ky)
           sxy(i) = ckx(i,kx)*sky(i,ky)+skx(i,kx)*cky(i,ky)
        enddo
        do kz = -kmax,kmax
           rkz = zlat*kz
           rk2 = rkx*rkx + rky*rky + rkz*rkz
           if (rk2.lt.rkmax2 .and. rk2.ne.0.d0) then
              sr = 0.d0
              si = 0.d0
              do i = 1,n
                 tr = cxy(i)*ckz(i,kz) - sxy(i)*skz(i,kz)
                 ti = cxy(i)*skz(i,kz) + sxy(i)*ckz(i,kz)
                 sr = sr + tr
                 si = si + ti
                 dsdr(1,i) = -rkx*ti
                 dsdr(2,i) =  rkx*tr
                 dsdr(3,i) = -rky*ti
                 dsdr(4,i) =  rky*tr
                 dsdr(5,i) = -rkz*ti
                 dsdr(6,i) =  rkz*tr
              enddo
              w = (f/rk2)*dexp(-b*rk2)
              et = w*(sr*sr+si*si)
              v = v + et
              w = 2.d0*w
              sr = w*sr
              si = w*si
              do i = 1,n
                 dvdr(1,i) = dvdr(1,i) + sr*dsdr(1,i) + si*dsdr(2,i)
                 dvdr(2,i) = dvdr(2,i) + sr*dsdr(3,i) + si*dsdr(4,i)
                 dvdr(3,i) = dvdr(3,i) + sr*dsdr(5,i) + si*dsdr(6,i)
              enddo
              
              ! Virial

              term = 2.d0*(1.d0/rk2 + b)
              xx = et * (term*rkx*rkx-1.d0)
              xy = et * (term*rkx*rky)
              xz = et * (term*rkx*rkz)
              yy = et * (term*rky*rky-1.d0)
              yz = et * (term*rky*rkz)
              zz = et * (term*rkz*rkz-1.d0)
              vir(1,1) = vir(1,1) + xx
              vir(1,2) = vir(1,2) + xy
              vir(1,3) = vir(1,3) + xz
              vir(2,1) = vir(2,1) + xy
              vir(2,2) = vir(2,2) + yy
              vir(2,3) = vir(2,3) + yz
              vir(3,1) = vir(3,1) + xz
              vir(3,2) = vir(3,2) + yz
              vir(3,3) = vir(3,3) + zz
           endif
        enddo
     enddo
  enddo

  ! ...minus the self term

  vs = 0.d0
  do i = 1,n
     vs = vs+z(i)*z(i)
  enddo
  vs = rfac*vs
  v = v-vs

  return
end subroutine kwald_nc



