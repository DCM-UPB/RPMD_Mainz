subroutine intra_harmonic(r,dvdr,na,v,vir)
  implicit none
  ! ------------------------------------------------------------------
  ! This subroutine calculates the intramolecular potential energy
  ! using Toukan and Rahman's SPC/F harmonic potential.
  !  
  ! Atomic units used throughout.
  ! ------------------------------------------------------------------
  integer na,j
  real(8) r(3,na),dvdr(3,na),vir(3,3)
  real(8) dr1,dr2,dr3,theta,reoh
  real(8) dx1,dy1,dz1,v
  real(8) dx2,dy2,dz2,dr3dx3
  real(8) dx3,dy3,dz3,dr3dy3,dr3dz3
  real(8) dr1dx1,dr1dy1,dr1dz1,dr2dx2,dr2dy2,dr2dz2
  real(8) dvdx1,dvdx2,dvdx3,dvdy1,dvdy2,dvdy3
  real(8) dvdz1,dvdz2,dvdz3,rehh,dvdr1,dvdr2,dvdr3
  real(8) apot,bpot,cpot,dpot
  common /geometry/ theta,reoh

  ! SPC/F parameters from Toukan and Rahman.

  apot = 0.5993d0
  bpot = 0.1466d0
  cpot = -0.09435d0
  dpot = 0.04984d0

  rehh = 2.d0*dsin(0.5d0*theta)*reoh

  v = 0.d0
  vir(:,:) = 0.d0

  do j = 1,na,3
     dx1 = r(1,j+1)-r(1,j)
     dy1 = r(2,j+1)-r(2,j)
     dz1 = r(3,j+1)-r(3,j)
     dr1 = dsqrt(dx1*dx1+dy1*dy1+dz1*dz1)
     dr1dx1 = dx1/dr1
     dr1dy1 = dy1/dr1
     dr1dz1 = dz1/dr1
     dr1 = dr1-reoh
     dx2 = r(1,j+2)-r(1,j)
     dy2 = r(2,j+2)-r(2,j)
     dz2 = r(3,j+2)-r(3,j)
     dr2 = dsqrt(dx2*dx2+dy2*dy2+dz2*dz2)
     dr2dx2 = dx2/dr2
     dr2dy2 = dy2/dr2
     dr2dz2 = dz2/dr2
     dr2 = dr2-reoh
     dx3 = r(1,j+2)-r(1,j+1)
     dy3 = r(2,j+2)-r(2,j+1)
     dz3 = r(3,j+2)-r(3,j+1)
     dr3 = dsqrt(dx3*dx3+dy3*dy3+dz3*dz3)
     dr3dx3 = dx3/dr3
     dr3dy3 = dy3/dr3
     dr3dz3 = dz3/dr3
     dr3 = dr3-rehh
     v = v+0.5d0*apot*(dr1*dr1+dr2*dr2)
     v = v+0.5d0*bpot*dr3*dr3
     v = v+cpot*(dr1+dr2)*dr3
     v = v+dpot*dr1*dr2
     dvdr1 = apot*dr1+cpot*dr3+dpot*dr2
     dvdr2 = apot*dr2+cpot*dr3+dpot*dr1
     dvdr3 = bpot*dr3+cpot*(dr1+dr2)
     dvdx1 = dvdr1*dr1dx1
     dvdy1 = dvdr1*dr1dy1
     dvdz1 = dvdr1*dr1dz1 
     dvdx2 = dvdr2*dr2dx2
     dvdy2 = dvdr2*dr2dy2
     dvdz2 = dvdr2*dr2dz2
     dvdx3 = dvdr3*dr3dx3
     dvdy3 = dvdr3*dr3dy3
     dvdz3 = dvdr3*dr3dz3
     dvdr(1,j) = dvdr(1,j)-dvdx1-dvdx2
     dvdr(2,j) = dvdr(2,j)-dvdy1-dvdy2
     dvdr(3,j) = dvdr(3,j)-dvdz1-dvdz2
     dvdr(1,j+1) = dvdr(1,j+1)+dvdx1-dvdx3
     dvdr(2,j+1) = dvdr(2,j+1)+dvdy1-dvdy3
     dvdr(3,j+1) = dvdr(3,j+1)+dvdz1-dvdz3
     dvdr(1,j+2) = dvdr(1,j+2)+dvdx2+dvdx3
     dvdr(2,j+2) = dvdr(2,j+2)+dvdy2+dvdy3
     dvdr(3,j+2) = dvdr(3,j+2)+dvdz2+dvdz3

     ! Virial.
     
     vir(1,1) = vir(1,1) + dx1*dvdx1 + dx2*dvdx2 + dx3*dvdx3
     vir(2,2) = vir(2,2) + dy1*dvdy1 + dy2*dvdy2 + dy3*dvdy3
     vir(3,3) = vir(3,3) + dz1*dvdz1 + dz2*dvdz2 + dz3*dvdz3
     vir(1,2) = vir(1,2) + dx1*dvdy1 + dx2*dvdy2 + dx3*dvdy3
     vir(1,3) = vir(1,3) + dx1*dvdz1 + dx2*dvdz2 + dx3*dvdz3
     vir(2,3) = vir(2,3) + dy1*dvdz1 + dy2*dvdz2 + dy3*dvdz3
     vir(3,2) = vir(3,2) + dz1*dvdy1 + dz2*dvdy2 + dz3*dvdy3
     vir(2,1) = vir(2,1) + dy1*dvdx1 + dy2*dvdx2 + dy3*dvdx3
     vir(3,1) = vir(3,1) + dz1*dvdx1 + dz2*dvdx2 + dz3*dvdx3
  enddo

  return
end subroutine intra_harmonic
