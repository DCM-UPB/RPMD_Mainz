subroutine com_velocity(px,py,pz,mass,p,na,nb)
  implicit none
  ! ------------------------------------------------------------------
  ! Subroutine to calculate system COM velocity
  ! ------------------------------------------------------------------
  integer na,nb,j,k
  real(8) px,py,pz,p(3,na,nb),mass(na),totmass
  
  px = 0.d0
  py = 0.d0
  pz = 0.d0
  
  do k = 1,nb
     do j = 1,na
        px = px + p(1,j,k)
        py = py + p(2,j,k)
        pz = pz + p(3,j,k)
     enddo
  enddo
  
  do j = 1,na
     totmass = totmass + mass(j)
  enddo

  totmass = totmass*dble(nb)

  px = px/totmass
  py = py/totmass
  pz = pz/totmass

  return
end subroutine com_velocity


subroutine molprop (r,avang,avoh,nm,nb,z,na)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Subroutine to calculate the angle and bond lengths in water
  ! ------------------------------------------------------------------
  integer nm,nb,k,j,ii,jj,kk,na
  real(8) r(3,na,nb),avang,avoh,z(na)
  real(8) dx1,dy1,dz1,dx2,dy2,dz2,d1,d2,pi

  pi = dacos(-1.d0)
  avang = 0.d0
  avoh = 0.d0

  do k = 1, nb
     do j = 1, nm

        ii = 3*(j-1) + 1
        jj = 3*(j-1) + 2
        kk = 3*(j-1) + 3
        dx1 = r(1,ii,k) - r(1,jj,k)
        dy1 = r(2,ii,k) - r(2,jj,k)
        dz1 = r(3,ii,k) - r(3,jj,k)
        d1 = dsqrt(dx1**2+dy1**2+dz1**2)
        avoh = avoh + d1
        
        dx2 = r(1,ii,k) - r(1,kk,k)
        dy2 = r(2,ii,k) - r(2,kk,k)
        dz2 = r(3,ii,k) - r(3,kk,k)
        d2 = dsqrt(dx2**2+dy2**2+dz2**2)
        avoh = avoh + d2
        avang=avang+dacos((1.d0/(d1*d2))*(dx1*dx2+dy1*dy2+dz1*dz2))
     enddo
  enddo
  avoh = avoh / 2.d0
  avang = avang * (180.d0/pi)

  return
end subroutine molprop

subroutine dipole(r,dipx,dipy,dipz,dip2,dipm,z,na,nb)
  implicit none
  ! ------------------------------------------------------------------
  ! Calculate the value of the total system dipole moment.
  ! dipx, dipy and dipz are the three components of
  ! the total system dipole.
  !
  ! dip2 collects the total system dipole moment squared
  ! dipm collects the total dipole moment per molecule
  ! ------------------------------------------------------------------
  integer na,nb,k,j,nm
  real(8) r(3,na,nb),z(na),dipx,dipy,dipz,dip2
  real(8) dipm,dipmx,dipmy,dipmz
  real(8) alpha,ecut,alpha2,fac1,fac2,wm,wh
  common /ew_param/ alpha,ecut,wm,wh

  alpha2 = 0.5d0 * (1.d0 - alpha)
  dipx = 0.d0
  dipy = 0.d0
  dipz = 0.d0
  dipm = 0.d0
  dip2 = 0.d0

  do k = 1, nb
     do j = 1,na,3
        fac1 = (z(j+1)+alpha2*z(j))
        fac2 = (z(j+2)+alpha2*z(j))
        dipmx = alpha*z(j)*r(1,j,k) + fac1 * r(1,j+1,k) &
                   + fac2 * r(1,j+2,k)
        dipmy = alpha*z(j)*r(2,j,k) + fac1 * r(2,j+1,k) &
                   + fac2 * r(2,j+2,k)
        dipmz = alpha*z(j)*r(3,j,k) + fac1 * r(3,j+1,k) &
                   + fac2 * r(3,j+2,k)
        dipx = dipx + dipmx
        dipy = dipy + dipmy
        dipz = dipz + dipmz
        dipm = dipm + dsqrt(dipmx*dipmx + dipmy*dipmy + dipmz*dipmz)
     enddo
  enddo
  nm = na/3
  dipm = dipm / dble(nb*nm)
  dipx = dipx / dble(nb)
  dipy = dipy / dble(nb)
  dipz = dipz / dble(nb)
  dip2 = dipx*dipx + dipy*dipy + dipz*dipz

  return
end subroutine dipole

subroutine rad_gyr(rgo,rgh,rgcm,r,na,nb,mass)
  implicit none
  ! ------------------------------------------------------------------
  ! Calculates the Radius of Gyration
  ! ------------------------------------------------------------------
  integer na,nb,nm,i,j,k,m
  real(8) r(3,na,nb),rmean(3),rcm(3,na/3,nb),mass(na),massw
  real(8) rgo,rgh,rgcm,dx,dy,dz,dr,sum

  ! Radius of Gyration of Oxygen Ring Polymer

  nm = na/3
  rgo = 0.d0
  do i = 0,nm-1
     do j = 1,3
        rmean(j) = 0.d0
        do k = 1,nb
           rmean(j) = rmean(j) + r(j,3*i+1,k)
        enddo
        rmean(j) = rmean(j)/nb
     enddo
     sum = 0.d0
     do k = 1,nb
        dx = r(1,3*i+1,k)-rmean(1)
        dy = r(2,3*i+1,k)-rmean(2)
        dz = r(3,3*i+1,k)-rmean(3)
        dr =  dx*dx+dy*dy+dz*dz
        sum = sum + dr
     enddo
     rgo = rgo+sum
  enddo
  rgo = rgo/(nm*nb)

  ! Radius of Gyration Hydrogen

  rgh = 0.d0
  do m = 2,3
     do i = 0,nm-1
        do j = 1,3
           rmean(j) = 0.d0
           do k = 1,nb
              rmean(j) = rmean(j) + r(j,3*i+m,k)
           enddo
           rmean(j) = rmean(j)/dble(nb)
        enddo
        sum = 0.d0
        do k = 1,nb
           dx = r(1,3*i+m,k)-rmean(1)
           dy = r(2,3*i+m,k)-rmean(2)
           dz = r(3,3*i+m,k)-rmean(3)
           dr = dx*dx+dy*dy+dz*dz
           sum = sum + dr
        enddo
        rgh = rgh+sum
     enddo
  enddo
  rgh = rgh/(2.d0*nm*nb)

  ! Radius of Gyration of the Centre of Mass

  ! Calculate COM for each molecule
  
  massw = mass(1)+mass(2)+mass(3)
  do i = 0,nm-1
     do j = 1,3
        do k =1,nb
           rcm(j,i+1,k) = mass(1)*r(j,3*i+1,k)+ &
                       mass(2)*r(j,3*i+2,k)+ mass(3)*r(j,3*i+3,k)
           rcm(j,i+1,k) = rcm(j,i+1,k)/massw
        enddo
     enddo
  enddo
  
  rgcm = 0.d0
  do i = 1,nm
     do j = 1,3
        rmean(j) = 0.d0
        do k = 1,nb
           rmean(j) = rmean(j) + rcm(j,i,k)
        enddo
        rmean(j) = rmean(j)/nb
     enddo
     sum = 0.d0
     do k = 1,nb
        dx = rcm(1,i,k)-rmean(1)
        dy = rcm(2,i,k)-rmean(2)
        dz = rcm(3,i,k)-rmean(3)
        dr = dx*dx+dy*dy+dz*dz
        sum = sum + dr
     enddo
     rgcm = rgcm+sum
  enddo
  rgcm = rgcm/(nm*nb)

  return
end subroutine rad_gyr

subroutine orientation (r,ax1,ay1,az1,ax2,ay2,az2,it,na,nb,nm,nt,z)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Subroutine to calculate the orientation vectors of each water
  ! molecule.
  !  
  ! In this version, the three vectors simply correspond to:
  ! 1 ) The dipole axis,
  ! 2 ) The inter-hydrogen vector,
  ! 3 ) The the vector perpendicular to the two above.
  ! ------------------------------------------------------------------
  integer na,nm,nb,k,ic,i,it,nt
  real(8) r(3*na,nb),dxx(3),dyy(3),dzz(3)
  real(8) ax1(3,nm,0:(2*nt/iskip))
  real(8) ay1(3,nm,0:(2*nt/iskip))
  real(8) az1(3,nm,0:(2*nt/iskip))
  real(8) ax2(6,nm,0:(2*nt/iskip))
  real(8) ay2(6,nm,0:(2*nt/iskip))
  real(8) az2(6,nm,0:(2*nt/iskip))
  real(8) dx,dy,dz,size,z(na),alpha,alpha2,rmx,rmy,rmz,ecut,wt
  real(8) wm,wh
  common /ew_param/ alpha,ecut,wm,wh

  ! Calculate the vectors.

  alpha2 = 0.5d0 * (1.d0 - alpha)
  do k = 1, nb
     ic = 0
     do i = 1, nm
        rmx = alpha * r(ic+1,k) + alpha2*(r(ic+4,k)+r(ic+7,k))
        rmy = alpha * r(ic+2,k) + alpha2*(r(ic+5,k)+r(ic+8,k))
        rmz = alpha * r(ic+3,k) + alpha2*(r(ic+6,k)+r(ic+9,k))
        dx=z(3*i-2)*rmx+z(3*i-1)*r(ic+4,k) &
                   +z(3*i)*r(ic+7,k)
        dy=z(3*i-2)*rmy+z(3*i-1)*r(ic+5,k) &
                   +z(3*i)*r(ic+8,k)
        dz=z(3*i-2)*rmz+z(3*i-1)*r(ic+6,k) &
                   +z(3*i)*r(ic+9,k)
        size =dsqrt(dx**2+dy**2+dz**2) 
        dxx(1)=-dx/size
        dyy(1)=-dy/size
        dzz(1)=-dz/size            
        dx= r(ic+7,k) - r(ic+4,k)
        dy= r(ic+8,k) - r(ic+5,k)
        dz= r(ic+9,k) - r(ic+6,k)
        size =dsqrt(dx**2+dy**2+dz**2)
        dxx(2)=dx/size
        dyy(2)=dy/size
        dzz(2)=dz/size
        dx = dyy(1) * dzz(2) - dzz(1) * dyy(2) 
        dy = dzz(1) * dxx(2) - dxx(1) * dzz(2) 
        dz = dxx(1) * dyy(2) - dyy(1) * dxx(2)
        size =dsqrt(dx**2+dy**2+dz**2) 
        dxx(3)=dx/size
        dyy(3)=dy/size
        dzz(3)=dz/size
        ax1(1,i,it) = ax1(1,i,it) + dxx(1)
        ax1(2,i,it) = ax1(2,i,it) + dyy(1)
        ax1(3,i,it) = ax1(3,i,it) + dzz(1)
        ay1(1,i,it) = ay1(1,i,it) + dxx(2)
        ay1(2,i,it) = ay1(2,i,it) + dyy(2)
        ay1(3,i,it) = ay1(3,i,it) + dzz(2)
        az1(1,i,it) = az1(1,i,it) + dxx(3)
        az1(2,i,it) = az1(2,i,it) + dyy(3)
        az1(3,i,it) = az1(3,i,it) + dzz(3)

        ! Calculate matrices for second-order Legendre polynomials.
        
        ax2(1,i,it) = ax2(1,i,it) + dxx(1) * dxx(1)
        ax2(2,i,it) = ax2(2,i,it) + dyy(1) * dyy(1)
        ax2(3,i,it) = ax2(3,i,it) + dzz(1) * dzz(1)
        ax2(4,i,it) = ax2(4,i,it) + dxx(1) * dyy(1)
        ax2(5,i,it) = ax2(5,i,it) + dxx(1) * dzz(1)
        ax2(6,i,it) = ax2(6,i,it) + dyy(1) * dzz(1)
        
        ay2(1,i,it) = ay2(1,i,it) + dxx(2) * dxx(2)
        ay2(2,i,it) = ay2(2,i,it) + dyy(2) * dyy(2)
        ay2(3,i,it) = ay2(3,i,it) + dzz(2) * dzz(2)
        ay2(4,i,it) = ay2(4,i,it) + dxx(2) * dyy(2)
        ay2(5,i,it) = ay2(5,i,it) + dxx(2) * dzz(2)
        ay2(6,i,it) = ay2(6,i,it) + dyy(2) * dzz(2)
        
        az2(1,i,it) = az2(1,i,it) + dxx(3) * dxx(3)
        az2(2,i,it) = az2(2,i,it) + dyy(3) * dyy(3)
        az2(3,i,it) = az2(3,i,it) + dzz(3) * dzz(3)
        az2(4,i,it) = az2(4,i,it) + dxx(3) * dyy(3)
        az2(5,i,it) = az2(5,i,it) + dxx(3) * dzz(3)
        az2(6,i,it) = az2(6,i,it) + dyy(3) * dzz(3)
        
        ic = ic + 9
     enddo
  enddo

  ! Finally, average the 3 x 3 matrix values for the inertia
  ! tensors over the beads and create a 9 x nm matrix which
  ! is returned to "solve".

  wt = 1/dble(nb)
  do k = 1, nm
     do i = 1, 3
        ax1(i,k,it) = ax1(i,k,it) * wt
        ay1(i,k,it) = ay1(i,k,it) * wt
        az1(i,k,it) = az1(i,k,it) * wt
     enddo
     do i = 1, 6
        ax2(i,k,it) = ax2(i,k,it) * wt
        ay2(i,k,it) = ay2(i,k,it) * wt
        az2(i,k,it) = az2(i,k,it) * wt
     enddo
     do i = 4, 6
        ax2(i,k,it) = dsqrt(2.d0)*ax2(i,k,it)
        ay2(i,k,it) = dsqrt(2.d0)*ay2(i,k,it)
        az2(i,k,it) = dsqrt(2.d0)*az2(i,k,it)
     enddo
  enddo
  
  return
end subroutine orientation

subroutine center_atoms (p,pc,na,nb)
  implicit none
  ! ------------------------------------------------------------------
  ! Calculates ring-polymer coordinate centroids qc(j).
  ! ------------------------------------------------------------------
  integer na,nb,k
  real(8) p(3,na,nb), pc(3,na),wt

  pc(:,:) = 0.d0

  wt = 1.d0/dble(nb)

  do k = 1,nb
     pc(:,:) = pc(:,:)+p(:,:,k)
  enddo
  pc(:,:) = wt*pc(:,:)

  return
end subroutine center_atoms

subroutine center_water (p,pc,nm,nb)
  implicit none
  ! ------------------------------------------------------------------
  ! Molecular momentum centroids for nm = na/3 water molecules.
  ! ------------------------------------------------------------------
  integer nm, nb,j,k
  real(8) p(9,nm,nb),pc(3,nm),wt

  pc(:,:) = 0.d0
  wt = 1/dble(nb)

  do j = 1,nm
     do k = 1,nb
        pc(1,j) = pc(1,j) + p(1,j,k) + p(4,j,k) + p(7,j,k)
        pc(2,j) = pc(2,j) + p(2,j,k) + p(5,j,k) + p(8,j,k)
        pc(3,j) = pc(3,j) + p(3,j,k) + p(6,j,k) + p(9,j,k)
     enddo
     pc(1,j) = pc(1,j) * wt
     pc(2,j) = pc(2,j) * wt
     pc(3,j) = pc(3,j) * wt
  enddo
  
  return
end subroutine center_water

