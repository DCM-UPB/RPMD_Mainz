subroutine intra_morse(r,dvdr,dvdr_stretch,dvdr_bend,na,v,vir)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Intramolecular potential energy function based on Morse
  ! potentials for the two stretches and the bending motion of the
  ! monomer.
  !
  ! A quartic expansion of the morse potential is used rather than
  ! the actual Morse potential, in order to prevent dissociation.
  !
  ! Edited version for force field matching:
  ! dvdr_stretch : stretch part of intram. force
  ! dvdr_bend    : bend part of intram. force
  ! -------------------------------------------------------------------
  integer na,j
  real(8) r(3,na),dvdr(3,na),vir(3,3),dvdr_stretch(3,na),dvdr_bend(3,na)
  real(8) dr1,dr2,dr3,theta,reoh
  real(8) dx1,dy1,dz1,v,dx2,dy2,dz2,dr1sq,dr2sq,dr3dx3
  real(8) dx3,dy3,dz3,dr3sq,dr1dx1,dr1dy1,dr1dz1
  real(8) dr2dx2,dr2dy2,dr2dz2,dr3dy3,dr3dz3,u,vv,v2,arg,ang
  real(8) dang,uprime,vprime,grad
  real(8) dvdx1,dvdx2,dvdx3,dvdy1,dvdy2,dvdy3
  real(8) dvdz1,dvdz2,dvdz3,dvdr1,dvdr2,dvdr3
  real(8) dthetadr1,dthetadr2,dthetadr3
  real(8) darg,dr1a,dr2a,dr3a,apot,bpot
  real(8) xx,xy,xz,yy,yz,zz
  real(8) de,alp,alp2,alp3,alp4,drasq,drbsq
  real(8) alpb,alpb2,alpb3,alpb4,f1,f2,deb,a1,a2,a3
  common /geometry/ theta, reoh
  common /intra_param/ apot,bpot,alp,alpb

  de = apot
  deb = bpot
  alpb2 = alpb*alpb
  alpb3 = alpb2*alpb
  alpb4 = alpb2 * alpb2
  alp2  = alp*alp
  alp3  = alp*alp2
  alp4  = alp3*alp
  f1 = 7.d0 / 12.d0
  f2 = 7.d0 / 3.d0

  v = 0.d0
  vir(:,:) = 0.d0

  ! Loop over molecules

  do j = 1,na,3
     dx1 = r(1,j+1)-r(1,j)
     dy1 = r(2,j+1)-r(2,j)
     dz1 = r(3,j+1)-r(3,j)
     dr1sq = dx1*dx1 + dy1*dy1 + dz1*dz1
     dr1 = dsqrt(dr1sq)
     dr1a = dr1
     dr1dx1 = dx1/dr1
     dr1dy1 = dy1/dr1
     dr1dz1 = dz1/dr1
     dr1 = dr1-reoh
     drasq=dr1*dr1
     
     dx2 = r(1,j+2)-r(1,j)
     dy2 = r(2,j+2)-r(2,j)
     dz2 = r(3,j+2)-r(3,j)
     dr2sq = dx2*dx2 + dy2*dy2 + dz2*dz2
     dr2 = dsqrt(dr2sq)
     dr2a = dr2
     dr2dx2 = dx2/dr2
     dr2dy2 = dy2/dr2
     dr2dz2 = dz2/dr2
     dr2 = dr2-reoh
     drbsq=dr2*dr2
     
     dx3 = r(1,j+2)-r(1,j+1)
     dy3 = r(2,j+2)-r(2,j+1)
     dz3 = r(3,j+2)-r(3,j+1)
     dr3sq = dx3*dx3 + dy3*dy3 + dz3*dz3
     dr3 = dsqrt(dr3sq)
     dr3a = dr3
     dr3dx3 = dx3/dr3
     dr3dy3 = dy3/dr3
     dr3dz3 = dz3/dr3

     u = (dr1sq + dr2sq - dr3sq)
     vv = (2.d0 * dr1a * dr2a)         
     v2 = 1.d0/(vv * vv)
     arg = u / vv
     darg = -1.d0/dsqrt(1.d0 - arg*arg)
     ang = dacos( arg )
     dang = (ang - theta)
     uprime = 2.d0 * dr1a
     vprime = 2.d0 * dr2a
     grad = (uprime*vv - vprime*u) * v2
     dthetadr1= darg * grad
     uprime = 2.d0 * dr2a
     vprime = 2.d0 * dr1a
     grad = (uprime*vv - vprime*u) * v2
     dthetadr2= darg * grad
     uprime = -2.d0*dr3a
     grad = (uprime*vv) * v2
     dthetadr3= darg * grad

     v = v+de*(alp2*drasq-alp3*dr1*drasq+f1*alp4*drasq*drasq)
     v = v+de*(alp2*drbsq-alp3*dr2*drbsq+f1*alp4*drbsq*drbsq)
     v = v+deb*(alpb2*dang**2-alpb3*dang**3+f1*alpb4*dang**4)

     a1 = de*(2.d0*alp2*dr1-3.d0*alp3*drasq+f2*alp4*dr1*drasq)
     a2 = de*(2.d0*alp2*dr2-3.d0*alp3*drbsq+f2*alp4*dr2*drbsq)
     a3 = deb*(2.d0*alpb2*dang-3.d0*alpb3*dang**2+f2*alpb4*dang**3)
     dvdr1 = a1+a3*dthetadr1
     dvdr2 = a2+a3*dthetadr2
     dvdr3 = a3*dthetadr3
     dvdx1 = dvdr1*dr1dx1
     dvdy1 = dvdr1*dr1dy1
     dvdz1 = dvdr1*dr1dz1 
     dvdx2 = dvdr2*dr2dx2
     dvdy2 = dvdr2*dr2dy2
     dvdz2 = dvdr2*dr2dz2
     dvdx3 = dvdr3*dr3dx3
     dvdy3 = dvdr3*dr3dy3
     dvdz3 = dvdr3*dr3dz3
     dvdr(1,j) = dvdr(1,j) - dvdx1 - dvdx2
     dvdr(2,j) = dvdr(2,j) - dvdy1 - dvdy2
     dvdr(3,j) = dvdr(3,j) - dvdz1 - dvdz2
     dvdr(1,j+1) = dvdr(1,j+1) + dvdx1 - dvdx3
     dvdr(2,j+1) = dvdr(2,j+1) + dvdy1 - dvdy3
     dvdr(3,j+1) = dvdr(3,j+1) + dvdz1 - dvdz3
     dvdr(1,j+2) = dvdr(1,j+2) + dvdx2 + dvdx3
     dvdr(2,j+2) = dvdr(2,j+2) + dvdy2 + dvdy3
     dvdr(3,j+2) = dvdr(3,j+2) + dvdz2 + dvdz3

     ! hacked-in force split

     dvdr_stretch(1,j) = - a1*dr1dx1 - a2*dr2dx2
     dvdr_stretch(2,j) = - a1*dr1dy1 - a2*dr2dy2
     dvdr_stretch(3,j) = - a1*dr1dz1 - a2*dr2dz2
     dvdr_stretch(1,j+1) = a1*dr1dx1
     dvdr_stretch(2,j+1) = a1*dr1dy1
     dvdr_stretch(3,j+1) = a1*dr1dz1
     dvdr_stretch(1,j+2) = a2*dr2dx2
     dvdr_stretch(2,j+2) = a2*dr2dy2
     dvdr_stretch(3,j+2) = a2*dr2dz2

     dvdr_bend(1,j) = - a3*dthetadr1*dr1dx1 - a3*dthetadr2*dr2dx2
     dvdr_bend(2,j) = - a3*dthetadr1*dr1dy1 - a3*dthetadr2*dr2dy2
     dvdr_bend(3,j) = - a3*dthetadr1*dr1dz1 - a3*dthetadr2*dr2dz2
     dvdr_bend(1,j+1) = a3*dthetadr1*dr1dx1 - a3*dthetadr3*dr3dx3
     dvdr_bend(2,j+1) = a3*dthetadr1*dr1dy1 - a3*dthetadr3*dr3dy3
     dvdr_bend(3,j+1) = a3*dthetadr1*dr1dz1 - a3*dthetadr3*dr3dz3
     dvdr_bend(1,j+2) = a3*dthetadr2*dr2dx2 + a3*dthetadr3*dr3dx3
     dvdr_bend(2,j+2) = a3*dthetadr2*dr2dy2 + a3*dthetadr3*dr3dy3
     dvdr_bend(3,j+2) = a3*dthetadr2*dr2dz2 + a3*dthetadr3*dr3dz3

     ! Virial

     xx = dx1*dvdx1 + dx2*dvdx2 + dx3*dvdx3
     xy = dx1*dvdy1 + dx2*dvdy2 + dx3*dvdy3
     xz = dx1*dvdz1 + dx2*dvdz2 + dx3*dvdz3
     yy = dy1*dvdy1 + dy2*dvdy2 + dy3*dvdy3
     yz = dy1*dvdz1 + dy2*dvdz2 + dy3*dvdz3
     zz = dz1*dvdz1 + dz2*dvdz2 + dz3*dvdz3
     vir(1,1) = vir(1,1) + xx
     vir(1,2) = vir(1,2) + xy
     vir(1,3) = vir(1,3) + xz
     vir(2,1) = vir(2,1) + xy
     vir(2,2) = vir(2,2) + yy
     vir(2,3) = vir(2,3) + yz
     vir(3,1) = vir(3,1) + xz
     vir(3,2) = vir(3,2) + yz
     vir(3,3) = vir(3,3) + zz
     
  enddo
  
  return 
end subroutine intra_morse
