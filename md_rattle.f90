subroutine rattle_setup(nm)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Setup of rattle routine
  ! ------------------------------------------------------------------
  integer nm, nctot, bond(nbondmax,2),i,nbond
  real(8) rcb(nbondmax)
  real(8) reoh,rehh,theta
  common /rattle/ rcb,bond
  common /constraint/ nctot,nbond
  common /geometry/ theta,reoh

  rehh = 2.d0*dsin(0.5d0*theta)*reoh
  nctot = 0
  nbond = 3
  do i = 1,nm
     rcb(1+nctot) = reoh * reoh
     bond(1+nctot,1) = 3*i - 2
     bond(1+nctot,2) = 3*i - 1
     rcb(2+nctot) = reoh * reoh
     bond(2+nctot,1) = 3*i - 2
     bond(2+nctot,2) = 3*i
     rcb(3+nctot) = rehh * rehh
     bond(3+nctot,1) = 3*i - 1
     bond(3+nctot,2) = 3*i
     nctot = nctot + nbond
  enddo

  return
end subroutine rattle_setup

subroutine rattle_s1(r,rold,p,nm,mass,dt,dvdr,vir,na)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Shake/Rattle Routine : Stage 1
  !
  ! nam     number of atoms per molecule
  ! nbond   number of bonds per molecule
  ! itmax   maximum number of shake iterations
  ! tol     fractional tolerance of constraints.
  ! rcb     squared constrained bond lengths.
  ! bond    bonded atom matrix
  ! nctot   total number of bonds
  ! nbtot   total number of bonds - local counter
  ! ------------------------------------------------------------------
  integer na,nm,nctot,nbtot,bond(nbondmax,2),itmax
  integer ii,jj,it,i,nbond,nam,noff,ik
  real(8) r(3,na),p(3,na),rold(3,na),dvdr(3,na)
  real(8) mass(na),rcb(nbondmax),vir(3,3)
  real(8) dx1,dy1,dz1,d2,diff,rm1,rm2
  real(8) dx2,dy2,dz2,dot,gamma,dmax,gdx,gdy,gdz
  real(8) tol,dt
  logical done
  common /rattle/ rcb,bond
  common /constraint/ nctot,nbond

  vir(:,:) = 0.d0
  tol = 1d-6
  it = 0
  itmax = 1000
  nam = 3
  noff = 0
  nbtot = 0

  do ik = 1, nm
     done = .false.
     it = 0
     do while ((.not.done).and.(it.lt.itmax))
        done = .true.
        dmax = -1d10
        do i = 1, nbond
           ii = bond(i+nbtot,1)
           jj = bond(i+nbtot,2)
           dx1 = r(1,jj) - r(1,ii)
           dy1 = r(2,jj) - r(2,ii)
           dz1 = r(3,jj) - r(3,ii)
           d2 = dx1*dx1 + dy1*dy1 + dz1*dz1
           diff = rcb(i+nbtot) - d2
           if (dabs(diff).gt.dmax)then
              dmax = dabs(diff)
           endif
           if (dabs(diff) .gt. rcb(i+nbtot)*tol) then
              rm1 = 1.d0 / mass(jj)
              rm2 = 1.d0 / mass(ii)
              dx2 = rold(1,jj) - rold(1,ii)
              dy2 = rold(2,jj) - rold(2,ii)
              dz2 = rold(3,jj) - rold(3,ii)
              dot = dx1*dx2 + dy1*dy2 + dz1*dz2
              gamma = diff / (2.d0*(rm1+rm2)*dot)

              ! Positions

              gdx = gamma * dx2
              gdy = gamma * dy2
              gdz = gamma * dz2
              r(1,jj) = r(1,jj) + rm1 * gdx
              r(2,jj) = r(2,jj) + rm1 * gdy
              r(3,jj) = r(3,jj) + rm1 * gdz
              r(1,ii) = r(1,ii) - rm2 * gdx
              r(2,ii) = r(2,ii) - rm2 * gdy
              r(3,ii) = r(3,ii) - rm2 * gdz

              ! Momenta

              gdx = gdx / dt
              gdy = gdy / dt
              gdz = gdz / dt
              p(1,jj) = p(1,jj) + gdx
              p(2,jj) = p(2,jj) + gdy
              p(3,jj) = p(3,jj) + gdz
              p(1,ii) = p(1,ii) - gdx
              p(2,ii) = p(2,ii) - gdy
              p(3,ii) = p(3,ii) - gdz

              ! Virial

              gdx = 2.d0 * gdx
              gdy = 2.d0 * gdy
              gdz = 2.d0 * gdz
              vir(1,1) = vir(1,1) - gdx * dx2
              vir(1,2) = vir(1,2) - gdx * dy2
              vir(1,3) = vir(1,2) - gdx * dz2
              vir(2,3) = vir(2,3) - gdy * dz2
              vir(2,2) = vir(2,2) - gdy * dy2
              vir(3,3) = vir(3,3) - gdz * dz2

              ! Forces

              dvdr(1,jj) = dvdr(1,jj) + gdx
              dvdr(2,jj) = dvdr(2,jj) + gdy
              dvdr(3,jj) = dvdr(3,jj) + gdz
              dvdr(1,ii) = dvdr(1,ii) - gdx
              dvdr(2,ii) = dvdr(2,ii) - gdy
              dvdr(3,ii) = dvdr(3,ii) - gdz

              done = .false.
           endif
        enddo
        it = it + 1
     enddo
     noff = noff + nam
     nbtot = nbtot + nbond
  enddo

  ! Stop if the constraints have not converged.
  
  if (it.ge.itmax) then
     stop 'max iterations exceeded in rattle stage 1'
  endif

  if (nctot.ne.nbtot) then
     stop 'number of bonds inconsistent in shake/rattle routine'
  endif

  vir(2,1) = vir(1,2)
  vir(3,1) = vir(1,3)
  vir(3,2) = vir(2,3)
  
  return
end subroutine rattle_s1

subroutine rattle_s2(r,p,nm,mass,dt,dvdr,vir,na)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Shake/Rattle Routine : Stage 2
  ! ------------------------------------------------------------------
  integer na,nm,nctot,nbtot,bond(nbondmax,2),itmax
  integer ii,jj,it,i,nbond,nam,noff,ik
  real(8) r(3,na),p(3,na),dvdr(3,na)
  real(8) mass(na),rcb(nbondmax),vir(3,3)
  real(8) tol,dt,halfdt
  real(8) dx1,dy1,dz1,dvx1,dvy1,dvz1,rm1,rm2
  real(8) dot,gamma,dmax,gdx,gdy,gdz
  logical done
  common /rattle/ rcb,bond
  common /constraint/ nctot,nbond

  vir(:,:) = 0.d0
  tol = 1d-6
  it = 0
  itmax = 1000
  nam = 3
  noff = 0
  nbtot = 0
  halfdt = 0.5d0 * dt

  do ik = 1,nm
     done = .false.
     it = 0
     do while ((.not.done).and.(it.le.itmax))
        done = .true.
        dmax = -1d10
        do i = 1, nbond
           ii = bond(i+nbtot,1)
           jj = bond(i+nbtot,2)
           dx1 = r(1,jj) - r(1,ii)
           dy1 = r(2,jj) - r(2,ii)
           dz1 = r(3,jj) - r(3,ii)
           rm1 = 1.d0 / mass(jj)
           rm2 = 1.d0 / mass(ii)
           dvx1 = p(1,jj)*rm1 - p(1,ii)*rm2
           dvy1 = p(2,jj)*rm1 - p(2,ii)*rm2
           dvz1 = p(3,jj)*rm1 - p(3,ii)*rm2
           dot = dx1*dvx1 + dy1*dvy1 + dz1*dvz1
           gamma = -dot / ((rm1+rm2)*rcb(i+nbtot))
           if (dabs(dot).gt.dmax) then
              dmax = dabs(dot)
           endif
           if (dabs(gamma) .gt. tol) then
              
              ! Momenta

              gdx = gamma * dx1
              gdy = gamma * dy1
              gdz = gamma * dz1
              p(1,jj) = p(1,jj) + gdx
              p(2,jj) = p(2,jj) + gdy
              p(3,jj) = p(3,jj) + gdz
              p(1,ii) = p(1,ii) - gdx
              p(2,ii) = p(2,ii) - gdy
              p(3,ii) = p(3,ii) - gdz

              ! Virial

              gdx = gdx/halfdt
              gdy = gdy/halfdt
              gdz = gdz/halfdt
              vir(1,1) = vir(1,1) - gdx * dx1
              vir(1,2) = vir(1,2) - gdx * dy1
              vir(1,3) = vir(1,2) - gdx * dz1
              vir(2,3) = vir(2,3) - gdy * dz1
              vir(2,2) = vir(2,2) - gdy * dy1
              vir(3,3) = vir(3,3) - gdz * dz1

              ! Forces

              dvdr(1,jj) = dvdr(1,jj) + gdx
              dvdr(2,jj) = dvdr(2,jj) + gdy
              dvdr(3,jj) = dvdr(3,jj) + gdz
              dvdr(1,ii) = dvdr(1,ii) - gdx
              dvdr(2,ii) = dvdr(2,ii) - gdy
              dvdr(3,ii) = dvdr(3,ii) - gdz

              done = .false.
           endif
        enddo
        it = it + 1
     enddo
     nbtot = nbtot + nbond
     noff = noff + nam
  enddo

  ! Stop if the constraints have not converged.

  if (it.ge.itmax) then
     stop 'max iterations exceeded in rattle stage 2'
  endif

  if (nctot.ne.nbtot) then
     stop 'number of bonds inconsistent in shake/rattle routine'
  endif

  vir(2,1) = vir(1,2)
  vir(3,1) = vir(1,3)
  vir(3,2) = vir(2,3)

  return
end subroutine rattle_s2
