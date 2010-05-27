subroutine md_melt(ne,p,r,dvdr,dvdr2,na,nb,boxlxyz,z,beta, &
                   dt,mass,irun,nm_ice)
  use thermostat
  use barostat
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Melting/Freezing Routine
  ! ------------------------------------------------------------------
  integer na,nb,ne,irun,jout,je,k,k1,k2,i,nm,nprog,nbaro
  integer jq,nq,nqprog,imax,nctot,nbond,nm_ice,j
  real(8) boxlxyz(3),beta,v,v1,v2,v3,dt
  real(8) p(3,na,nb),dvdr(3,na,nb),dvdr2(3,na,nb),z(na),mass(na)
  real(8) r(3,na,nb),vir(3,3),vir_lf(3,3),tvxyz(3)
  real(8) omegan,dx,dy,dz,dist
  real(8) vrp,vol,tps,tfs,tk,thresh,den,wmass,tq1,tq2
  real(8) dtfs,dtqt,tv,pres,ran2,gaussian,temp
  real(8) pcx,pcy,pcz,fcx,fcy,fcz
  real(8), allocatable :: rst(:,:)
  integer, allocatable :: id(:)
  external gaussian, ran2
  common /constraint/ nctot,nbond

  nbaro = 0

  vir(:,:) = 0.d0
  vir_lf(:,:) = 0.d0

  imax = 100
  allocate (rst(3,na/3),id(imax))

  ! Water starting centroid positions

  nm = na / 3
  call center_water(r,rst,nm,nb)

  ! Define some useful local constants

  dtfs = dt/tofs
  nprog = max(ne/10,1)
  jout = max(ne/100000,1)
  vol = boxlxyz(1)*boxlxyz(2)*boxlxyz(3)

  ! NPT equilibration (max 10000 time steps)
  ! very heavy thermostatting used to avoid melting

  thresh = 0.25d0
  open (unit=45,file='E_pre_eq.out')

  nq = min(50000,ne/10)
  nqprog = max(nq/10,1)
  dtqt = 1.d-3*dble(nq)*dtfs
  write(6,'(a,f6.2,a)') ' * NPT pre-equilibration for: ', &
                            dtqt, ' ps'

  ! Initial forces and momenta

  call full_forces(r,na,nb,v,v1,v2,v3,vir,z,boxlxyz,dvdr,dvdr2)
  call sample(p,na,nb,mass,beta,irun,dt)

  do jq = 1,nq

     ! After 1/4 of pre equlibration turn on barostat

     if(jq.gt.nq/4) then
       nbaro = 1
     endif

     ! Evolve the system in time

     call evolve(p,r,v,v1,v2,v3,dvdr,dvdr2,dt,mass,na,nb, &
                 boxlxyz,z,beta,vir,vir_lf,irun,nbaro)

     tfs = jq*dtfs
     tps = 1.d-3*tfs

     ! Check energy conservation

     if (mod(jq,10).eq.0) then
        call kinetic(p,na,tk,mass,nb,beta)
        vrp = 0.d0
        omegan = dble(nb) / (beta * hbar)
        if (nb.gt.0) then
           do i = 1, na
              do k = 1, nb
                 if (k.eq.1) then
                    k1 = k
                    k2 = nb
                 else
                    k1 = k
                    k2 = k - 1
                 endif
                 dx = r(1,i,k1)-r(1,i,k2)
                 dy = r(2,i,k1)-r(2,i,k2)
                 dz = r(3,i,k1)-r(3,i,k2)
                 dist = (dx*dx+dy*dy+dz*dz)
                 vrp = vrp + 0.5d0*mass(i)*omegan**2*dist
              enddo
           enddo
        endif
        write (45,'(5(1x,f13.6))') tps,tokjmol*(nb*v+tk+vrp)/nm, &
                   tokjmol*(vrp/nm),tokjmol*(nb*v)/nm,tokjmol*tk/nm
     endif

     ! time for a collision with the heat bath?

     if ((therm.eq.'AND').or.(therm.eq.'PRA')) then
        if (ran2(irun,0.d0,1.d0) .lt. thresh) then
           call sample(p,na,nb,mass,beta,irun,dt)
        endif
     endif

     if (mod(jq,nqprog).eq.0) then
        write(6,*) 10*(jq/nqprog), ' %'
     endif
  enddo
  close (unit=45)

  open (unit=12,file='vmd_postpre_eq.xyz')
  call print_vmd_full(r,nb,na,nm,boxlxyz,12)
  close (unit=12)

  call sample(p,na,nb,mass,beta,irun,dt)

  ! Melting
  ! --------

  ! Open output files

  open (unit=20,file='temperature.out')
  open (unit=35,file='density.out')
  open (unit=37,file='box_len.out')
  open (unit=40,file='E.out')
  open (unit=50,file='pressure.out')

  ! Initial Density as a function of melting axis:

  id(:) = 0
  open (unit=56,file='density_profile_start.out')
  call den_axis(r,boxlxyz,na,nb,id,imax,1,56)
  close (unit=56)

  write(6,*) '* Beginning Melting/Freezing'
  do je = 1,ne
     if (je.gt.50000) then
        thresh = 0.0025d0
     endif

     ! Evolve the system by dt

     call evolve(p,r,v,v1,v2,v3,dvdr,dvdr2,dt,mass,na,nb, &
                 boxlxyz,z,beta,vir,vir_lf,irun,nbaro)

     tps = 1.d-3*je*dtfs

     ! Pressure and Virial Kinetic Energy
     write(6,*)'BEFORE VIRIAL KE'
     call virial_ke(r,dvdr,dvdr2,tv,tvxyz,tq1,tq2,beta,na,nb,mass)
     call pressure(pres,vir,tv,na,boxlxyz)

     ! Density

     vol = boxlxyz(1)*boxlxyz(2)*boxlxyz(3)
     wmass = mass(1) + mass(2) + mass(3)
     den = tokgcm3*(wmass*dble(nm)/vol)

     ! Density as a function of melting axis:

     if (mod(je,500).eq.0) then
        open (unit=56,file='density_profile.out')
        call den_axis(r,boxlxyz,na,nb,id,imax,1,56)
        close (unit=56)
     else
        call den_axis(r,boxlxyz,na,nb,id,imax,1,0)
     endif

     ! Current Density

     if (mod(je,100).eq.0) then
        write(35,*) tps , den
        write(37,'(f10.3,3f10.5)') tps, boxlxyz(1), &
                                   boxlxyz(2),boxlxyz(3)
     endif

     ! Pressure

     if (mod(je,100).eq.0) then
        write(50,*)tps,pres*tobar
     endif

     if (mod(je,jout) .eq. 0) then

        ! Kinetic Energy

        call kinetic(p,na,tk,mass,nb,beta)

        ! Calculate the contribution of the ring-polymer harmonic
        ! springs to the total Hamiltonian

        vrp = 0.d0
        omegan = dble(nb) / (beta * hbar)
        if (nb.gt.0) then
           do i = 1, na
              do k = 1, nb
                 if (k.eq.1) then
                    k1 = k
                    k2 = nb
                 else
                    k1 = k
                    k2 = k - 1
                 endif
                 dx = r(1,i,k1)-r(1,i,k2)
                 dy = r(2,i,k1)-r(2,i,k2)
                 dz = r(3,i,k1)-r(3,i,k2)
                 dist = (dx*dx+dy*dy+dz*dz)
                 vrp = vrp + 0.5d0*mass(i)*omegan**2*dist
              enddo
           enddo
        endif

        temp = toK*(2.d0*tk/(dble(3*na-nctot)*dble(nb)))/dble(nb)
        write (20,*)tps,temp
        write (40,'(5(1x,f13.6))') tps,tokjmol*(nb*v+tk+vrp)/nm, &
                tokjmol*(vrp/nm),tokjmol*(nb*v)/nm,tokjmol*tk/nm

     endif

     ! Every 10000 time steps output the current configuration

     if (mod(je,10000).eq.0) then
        open (unit=46,file='vmd_current.xyz')
        call print_vmd_full(r,nb,na,nm,boxlxyz,46)
        close (unit=46)
     endif

     ! time for a collision with the heat bath?

     if ((therm.eq.'AND').or.(therm.eq.'PRA')) then
        if (ran2(irun,0.d0,1.d0) .lt. thresh) then
           call sample(p,na,nb,mass,beta,irun,dt)
        endif
     endif

     ! Progress indicator

     if (mod(je,nprog).eq.0) then
        write(6,*) 10*(je/nprog), ' %'
     endif

  enddo
  close (unit=20)
  close (unit=35)
  close (unit=37)
  close (unit=40)
  close (unit=50)
  close (unit=56)
  
  write (6,*)'* Melting Simulation complete. '
  call flush (6)

  deallocate (rst)
  
  return
end subroutine md_melt

subroutine den_axis(r,boxlxyz,na,nb,id,imax,iaxis,iout)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Density as a function of axis
  ! ------------------------------------------------------------------
  integer na,nb,k,j,j1,i,iaxis,ibin,iout,imax,id(imax),norm
  real(8) r(3,na,nb),boxlxyz(3),onbox,box,delr
  real(8) rupper,rlower,wt
  real(8), allocatable :: ro(:)

  allocate (ro(na/3))
  ro(:) = 0.d0

  ! Oxygen centroids along required axis

  do k =1,nb
     do j = 1,na,3
        j1 = (j-1)/3 + 1
        ro(j1) = ro(j1) + r(iaxis,j,k)
     enddo
  enddo
  ro(:) = ro(:)/dble(nb)

  ! Periodic boundary conditions and histogram

  box = boxlxyz(iaxis)
  onbox = 1.d0/box
  delr = (box/dble(imax))

  do j = 1,na/3
     ro(j) = ro(j) - box*(nint(ro(j)*onbox)-0.5d0)
     ibin = int(ro(j)/delr)+1
     if (ibin.le.imax) then
        id(ibin) = id(ibin)+1
     endif
  enddo

  deallocate (ro)

  ! If iout is not equal to 0 then output current density profile

  if (iout.ne.0) then
     norm = 0
     do i = 1,imax
        norm = norm + id(i)
     enddo
     wt = 1.d0/dble(norm)
     
     do i = 1,imax
        rlower = dble(i-1)*delr
        rupper = rlower + delr
        write(iout,*) (rupper)*toA , dble(id(i))*wt
     enddo
     
     ! Clear histogram array

     id(:) = 0
  endif

  return
end subroutine den_axis


