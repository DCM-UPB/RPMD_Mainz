subroutine ewald_cubic(n,mol,q,r,v,vir,dvdr,boxl)
  implicit none
  ! ------------------------------------------------------------------
  ! Cubic Ewald sum with derivatives with respect to r.
  !
  ! The interactions between charges within the same molecule are
  ! excluded by subtracting the Coulomb potential between q(i)
  ! and q(j) when mol(i) = mol(j).
  ! ------------------------------------------------------------------
  common /ew_param/ alpha,ecut,wm,wh
  common /ew_setup/ wrcut,walpha,rkmax,kmax
  integer n,kmax,mcell,ncell,mol(n)
  real(8) q(n),r(3,n),dvdr(3,n),vir(3,3),wid(3),sm,sh
  real(8) wrcut,walpha,rkmax,rkmax2,v,boxl,onboxl,onbox2
  real(8) alpha,ecut,wm,wh
  real(8), allocatable :: rs(:,:)

  ! Scale box to unit length

  allocate (rs(3,n))
  onboxl = 1.d0/boxl
  rs(:,:) = onboxl*r(:,:)

  v = 0.d0
  vir(:,:) = 0.d0
  rkmax2 = rkmax*rkmax

  ! Linked cell list

  if (3.d0*wrcut .gt. 1.d0) then
     ncell = 1
     mcell = 0
  else
     ncell = max(6,int(2.d0/wrcut))
     mcell = 2
  endif
  
  ! Real space part

  if ((wm.gt.0.d0).or.(wh.gt.0.d0)) then

     ! Smeared Charge Ewald sum

     ! Values of sigma and width list

     sm = dsqrt(1.d0/(2.d0*wm))*onboxl
     sh = dsqrt(1.d0/(2.d0*wh))*onboxl

     wid(1) = 1.d0/(2.d0*sm)
     wid(2) = 1.d0/(2.d0*sh)
     wid(3) = 1.d0/(dsqrt(2.d0*(sm*sm + sh*sh)))

     call rwald_smeared(n,mol,q,rs,v,vir,dvdr,wrcut,walpha, &
                        mcell,ncell,wid)
  else

     ! Point Charge Ewald sum

     call rwald (n,mol,q,rs,v,vir,dvdr,wrcut,walpha,mcell,ncell)
  endif

  ! Reciprocal space part
  
  call kwald (n,q,rs,v,vir,dvdr,walpha,rkmax2,kmax)

  ! Scaling

  v = onboxl*v
  vir(:,:) = onboxl*vir(:,:)
  onbox2 = onboxl*onboxl
  dvdr(:,:) = onbox2*dvdr(:,:)

  deallocate (rs)
  
  return
end subroutine ewald_cubic

subroutine rwald (n,mol,q,r,v,vir,dvdr,rcut,alpha,mcell,ncell)
  implicit none
  ! ------------------------------------------------------------------
  ! Real space part of the Ewald sum, using a linked cell list.
  ! Low precision version with an approximation to erfc(x).
  ! ------------------------------------------------------------------
  integer n,mcell,ncell,nn,i,j,k,ix,iy,iz,jx,jy,jz,incx,incy,incz
  integer mol(n),link(ncell,ncell,ncell),list(n)
  integer ink(3,62),map(ncell,-2:2)
  real(8) q(n),r(3,n),dvdr(3,n),vir(3,3)
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

  ! construct the linked cell list

  call cell_list_unit(r,list,link,ink,map,nn,n,ncell,mcell,1)

  ! and evaluate the real space sum

  pi = dacos(-1.d0)
  rfac = alpha/dsqrt(pi)
  sfac = -2.d0*rfac
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
end subroutine rwald

subroutine kwald (n,q,r,v,vir,dvdr,alpha,rkmax2,kmax)
  implicit none
  ! ------------------------------------------------------------------
  ! Reciprocal space part of the Ewald sum.
  ! ------------------------------------------------------------------
  integer i,k,n,kmax,kx,ky,kz
  real(8) q(n),r(3,n),dvdr(3,n),vir(3,3)
  real(8) cx(2,n,0:kmax),cy(2,n,0:kmax),cz(2,n,0:kmax)
  real(8) cj(4,n),ck(8,n)
  real(8) alpha,rkmax2,rfac,vs,pi,rlat,arg,c,s,f,b,rkx,rk2
  real(8) v,s1,s2,s3,s4,s5,s6,s7,s8,cc,cca,csa,cs,ss,sc,sca
  real(8) w,wkx,wky,wkz,t2,t4,t6,t8,rky,rkz,et,term
  real(8) xx,xy,xz,yy,yz,zz
  real(8) ssa,ccb,csb,scb,ssb

  ! setup the trigonometric arrays
  
  pi = dacos(-1.d0)
  rlat = 2.d0*pi
  do i = 1,n
     cx(1,i,0) = q(i)
     cx(2,i,0) = 0.d0
     arg = rlat*r(1,i)
     c = dcos(arg)
     s = dsin(arg)
     do k = 1,kmax
        cx(1,i,k) = c*cx(1,i,k-1)-s*cx(2,i,k-1)
        cx(2,i,k) = s*cx(1,i,k-1)+c*cx(2,i,k-1)
     enddo
  enddo
  do i = 1,n
     cy(1,i,0) = 1.d0
     cy(2,i,0) = 0.d0
     arg = rlat*r(2,i)
     c = dcos(arg)
     s = dsin(arg)
     do k = 1,kmax
        cy(1,i,k) = c*cy(1,i,k-1)-s*cy(2,i,k-1)
        cy(2,i,k) = s*cy(1,i,k-1)+c*cy(2,i,k-1)
     enddo
  enddo
  do i = 1,n
     cz(1,i,0) = 1.d0
     cz(2,i,0) = 0.d0
     arg = rlat*r(3,i)
     c = dcos(arg)
     s = dsin(arg)
     do k = 1,kmax
        cz(1,i,k) = c*cz(1,i,k-1)-s*cz(2,i,k-1)
        cz(2,i,k) = s*cz(1,i,k-1)+c*cz(2,i,k-1)
     enddo
  enddo

  ! and evaluate the reciprocal space sum

  b = 0.25d0/(alpha*alpha)
  f = rlat
  do kx = 0,kmax
     if (kx .eq. 1) f = 2.d0*f
     rkx = rlat*kx
     rky = 0.d0
     rkz = 0.d0
     rk2 = rkx*rkx

     ! including contributions from ky = 0 and kz = 0

     if (rk2 .gt. 0.d0) then
        s1 = 0.d0
        s2 = 0.d0
        do i = 1,n
           s1 = s1+cx(1,i,kx)
           s2 = s2+cx(2,i,kx)
        enddo
        w = (f/rk2)*dexp(-b*rk2)
        et = w*(s1*s1+s2*s2)
        v = v+et

        ! Virial - missing terms are trivially zero

        term = 2.d0*(1.d0/rk2 + b)
        xx = et * (term*rkx*rkx-1.d0)
        yy = et * (-1.d0)
        zz = et * (-1.d0)
        vir(1,1) = vir(1,1) + xx
        vir(2,2) = vir(2,2) + yy
        vir(3,3) = vir(3,3) + zz

        w = w+w
        wkx = w*rkx
        do i = 1,n
           t2 = s2*cx(1,i,kx)-s1*cx(2,i,kx)
           dvdr(1,i) = dvdr(1,i) + wkx*t2
        enddo
     endif

     ! ky = 0 and |kz| > 0

     do kz = 1,kmax
        rky = 0.d0
        rkz = rlat*kz
        rk2 = rkx*rkx+rkz*rkz
        if (rk2 .gt. rkmax2) exit
        s1 = 0.d0
        s2 = 0.d0
        s3 = 0.d0
        s4 = 0.d0
        do i = 1,n
           cc = cx(1,i,kx)*cz(1,i,kz)
           cs = cx(1,i,kx)*cz(2,i,kz)
           sc = cx(2,i,kx)*cz(1,i,kz)
           ss = cx(2,i,kx)*cz(2,i,kz)
           cj(1,i) = cc-ss
           cj(2,i) = sc+cs
           cj(3,i) = cc+ss
           cj(4,i) = sc-cs
           s1 = s1+cj(1,i)
           s2 = s2+cj(2,i)
           s3 = s3+cj(3,i)
           s4 = s4+cj(4,i)
        enddo
        w = (f/rk2)*dexp(-b*rk2)
        et = w*(s1*s1+s2*s2+s3*s3+s4*s4)
        v = v + et

        ! Virial - missing terms are trivially zero

        term = 2.d0*(1.d0/rk2 + b)
        xx = et * (term*rkx*rkx-1.d0)
        xz = et * (term*rkx*rkz)
        yy = et * (-1.d0)
        zz = et * (term*rkz*rkz-1.d0)
        vir(1,1) = vir(1,1) + xx
        vir(1,3) = vir(1,3) + xz
        vir(2,2) = vir(2,2) + yy
        vir(3,1) = vir(3,1) + xz
        vir(3,3) = vir(3,3) + zz

        w = w+w
        wkx = w*rkx
        wkz = w*rkz
        do i = 1,n
           t2 = s2*cj(1,i)-s1*cj(2,i)
           t4 = s4*cj(3,i)-s3*cj(4,i)
           dvdr(1,i) = dvdr(1,i) + wkx*(t2+t4)
           dvdr(3,i) = dvdr(3,i) + wkz*(t2-t4)
        enddo
     enddo

     ! |ky| > 0 and kz = 0

     do ky = 1,kmax
        rkz = 0.d0
        rky = rlat*ky
        rk2 = rkx*rkx+rky*rky
        if (rk2 .gt. rkmax2) exit
        s1 = 0.d0
        s2 = 0.d0
        s3 = 0.d0
        s4 = 0.d0
        do i = 1,n
           cc = cx(1,i,kx)*cy(1,i,ky)
           cs = cx(1,i,kx)*cy(2,i,ky)
           sc = cx(2,i,kx)*cy(1,i,ky)
           ss = cx(2,i,kx)*cy(2,i,ky)
           cj(1,i) = cc-ss
           cj(2,i) = sc+cs
           cj(3,i) = cc+ss
           cj(4,i) = sc-cs
           s1 = s1+cj(1,i)
           s2 = s2+cj(2,i)
           s3 = s3+cj(3,i)
           s4 = s4+cj(4,i)
        enddo
        w = (f/rk2)*dexp(-b*rk2)
        et = w*(s1*s1+s2*s2+s3*s3+s4*s4)
        v = v + et

        ! Virial - missing terms are trivially zero

        term = 2.d0*(1.d0/rk2 + b)
        xx = et * (term*rkx*rkx-1.d0)
        xy = et * (term*rkx*rky)
        yy = et * (term*rky*rky-1.d0)
        zz = et * (-1.d0)
        vir(1,1) = vir(1,1) + xx
        vir(1,2) = vir(1,2) + xy
        vir(1,3) = vir(1,3) + xz
        vir(2,1) = vir(2,1) + xy
        vir(2,2) = vir(2,2) + yy
        vir(2,3) = vir(2,3) + yz
        vir(3,1) = vir(3,1) + xz
        vir(3,2) = vir(3,2) + yz
        vir(3,3) = vir(3,3) + zz

        w = w+w
        wkx = w*rkx
        wky = w*rky
        do i = 1,n
           t2 = s2*cj(1,i)-s1*cj(2,i)
           t4 = s4*cj(3,i)-s3*cj(4,i)
           dvdr(1,i) = dvdr(1,i) + wkx*(t2+t4)
           dvdr(2,i) = dvdr(2,i) + wky*(t2-t4)
        enddo

        ! and |ky| > 0 and |kz| > 0

        do kz = 1,kmax
           rkz = rlat*kz
           rk2 = rkx*rkx+rky*rky+rkz*rkz
           if (rk2 .gt. rkmax2) exit
           s1 = 0.d0
           s2 = 0.d0
           s3 = 0.d0
           s4 = 0.d0
           s5 = 0.d0
           s6 = 0.d0
           s7 = 0.d0
           s8 = 0.d0
           do i = 1,n
              cca = cj(1,i)*cz(1,i,kz)
              csa = cj(1,i)*cz(2,i,kz)
              sca = cj(2,i)*cz(1,i,kz)
              ssa = cj(2,i)*cz(2,i,kz)
              ccb = cj(3,i)*cz(1,i,kz)
              csb = cj(3,i)*cz(2,i,kz)
              scb = cj(4,i)*cz(1,i,kz)
              ssb = cj(4,i)*cz(2,i,kz)
              ck(1,i) = cca-ssa
              ck(2,i) = sca+csa
              ck(3,i) = cca+ssa
              ck(4,i) = sca-csa
              ck(5,i) = ccb-ssb
              ck(6,i) = scb+csb
              ck(7,i) = ccb+ssb
              ck(8,i) = scb-csb
              s1 = s1+ck(1,i)
              s2 = s2+ck(2,i)
              s3 = s3+ck(3,i)
              s4 = s4+ck(4,i)
              s5 = s5+ck(5,i)
              s6 = s6+ck(6,i)
              s7 = s7+ck(7,i)
              s8 = s8+ck(8,i)
           enddo
           w = (f/rk2)*dexp(-b*rk2)
           et = w*(s1*s1+s2*s2+s3*s3+s4*s4+s5*s5+s6*s6+s7*s7+s8*s8)
           v = v + et

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

           w = w+w
           wkx = w*rkx
           wky = w*rky
           wkz = w*rkz
           do i = 1,n
              t2 = s2*ck(1,i)-s1*ck(2,i)
              t4 = s4*ck(3,i)-s3*ck(4,i)
              t6 = s6*ck(5,i)-s5*ck(6,i)
              t8 = s8*ck(7,i)-s7*ck(8,i)
              dvdr(1,i) = dvdr(1,i) + wkx*(t2+t4+t6+t8)
              dvdr(2,i) = dvdr(2,i) + wky*(t2+t4-t6-t8)
              dvdr(3,i) = dvdr(3,i) + wkz*(t2-t4+t6-t8)
           enddo
        enddo
     enddo
  enddo

  ! subtract the self term

  vs = 0.d0
  do i = 1,n
     vs = vs + q(i)*q(i)
  enddo
  rfac = alpha/dsqrt(pi)
  vs = rfac*vs
  v = v-vs

  return
end subroutine kwald
