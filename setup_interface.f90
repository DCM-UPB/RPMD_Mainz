subroutine setup_int_size(rho_ice,rho_wat,box_ice,box_wat, &
                          nc_ice,nc_wat,nm_ice,nm_wat,wmass)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Calculates number of water molecules in interface boxes
  ! ------------------------------------------------------------------
  integer nc_ice(3),nc_wat(3),nm_ice,nm_wat
  real(8) rho_ice,rho_wat,box_ice(3),box_wat(3),wmass
  real(8) vbox,vcell,rlat,fac,denmol,denlen

  ! Size of ice lattice
  ! --------------------

  nm_ice = 8*nc_ice(1)*nc_ice(2)*nc_ice(3)

  vbox = (nm_ice*emass*wmass)/(rho_ice*1000.d0*(ToA*1d-10)**3)
  vcell = vbox / (nc_ice(1)*nc_ice(2)*nc_ice(3))
  rlat = ((3.d0*dsqrt(3.d0)/64.d0) * vcell)**(1.d0/3.d0)
  box_ice(1) = dsqrt(8.d0/3.d0) * rlat * nc_ice(1)
  box_ice(2) = dsqrt(8.d0) * rlat * nc_ice(2)
  box_ice(3) = (8.d0/3.d0) * rlat * nc_ice(3)

  ! Setup a water lattice to match up (y-z) plane
  ! ----------------------------------------------

  box_wat(2) = box_ice(2)
  box_wat(3) = box_ice(3)

  ! Box_water length 2 is chosen to give the required number of
  ! water molecules

  ! fac converts g/cm^3 to molecules/ao^3

  fac = 1000.d0*(ToA*1d-10)**3/(wmass*emass)
  denmol = rho_wat * fac                ! Molecules per ao^3
  denlen = denmol**(1.d0/3.d0)          ! Molecules per ao

  ! Calculate the number of water molecules over sides
  ! of interface region

  nc_wat(2) =  nint(box_wat(2)*denlen)
  nc_wat(3) =  nint(box_wat(3)*denlen)

  ! Calculate number of water molecules in the third direction
  ! such that they roughly equals number of ice ones.

  nc_wat(1) = nint(dble(nm_ice)/dble(nc_wat(2)*nc_wat(3)))

  ! Now calculate the third box length such that the starting
  ! density of the water is equal to rho_water

  nm_wat = nc_wat(1)*nc_wat(2)*nc_wat(3)
  vbox = (nm_wat*emass*wmass)/(rho_wat*1000.d0*(ToA*1d-10)**3)
  box_wat(1) = vbox/(box_wat(2)*box_wat(3))

  write(6,*)
  write(6,*) '-------------------------------------------'
  write(6,*) 'Interface Setup : '
  write(6,*)  '-------------------------------------------'
  write(6,'(a,i8)') ' nm_ice    = ', nm_ice
  write(6,'(a,i8)') ' nm_wat    = ', nm_wat

  return
end subroutine setup_int_size

subroutine setup_interface(r,na,ne,irun,nm_ice,nm_wat,box_ice, &
                           box_wat,boxlxyz,nc_ice,nc_wat, &
                           wmass,omass,hmass,qo,beta,dt,nb)
  implicit none
  include 'globals.inc'
  !------------------------------------------------------------------
  !Setup of an ice-water interface for melting simulations
  !
  !Interface is the (1120) plane - J.Cryst.Growth. 283 242-256(2005)
  !------------------------------------------------------------------
  integer i,j,k,irun,nb,nc_ice(3),nc_wat(3)
  integer na,nm_ice,na_ice,nm_wat,na_wat,ne,neq
  real(8) r(3,na,nb),box_ice(3),box_wat(3),boxlxyz(3)
  real(8) wmass,omass,hmass,qo,beta,dt,theta,reoh
  real(8) rcut,dx,dy,dz,boxmin,vol,den,dip,dipx,dipy,dipz,dip2,dipm
  real(8) oo_eps,oo_sig,oo_gam
  logical iamrigid
  real(8), allocatable :: r_ice(:,:,:),r_wat(:,:,:),z(:)
  common /geometry/ theta,reoh
  common /oo_param/ oo_eps,oo_sig,oo_gam,rcut
  common /structure/ iamrigid

  na_ice = 3*nm_ice
  na_wat = 3*nm_wat

  if (na_ice+na_wat.ne.na) then
     write(6,*) '** Interface: na_ice + na_wat not = na'
     stop
  endif

  allocate (r_ice(3,na_ice,nb))
  r_ice(:,:,:) = 0.d0

  ! Generate an ice lattice

  allocate (z(na_ice))
  do j = 1,na_ice,3
     z(j) = qo
     z(j+1) = -0.5d0*qo
     z(j+2) = -0.5d0*qo
  enddo

  ! Find low dipole ice configuration

  dip = 100.d0
  do while (dip.gt.2.3d0)
     call ice_lattice(r_ice(1,1,1),box_ice,nc_ice,nm_ice,irun,0)
     do k = 2,nb
        do j = 1,na_ice
           r_ice(1,j,k) = r_ice(1,j,1)
           r_ice(2,j,k) = r_ice(2,j,1)
           r_ice(3,j,k) = r_ice(3,j,1)
        enddo
     enddo
     call dipole(r_ice,dipx,dipy,dipz,dip2,dipm,z,na_ice,nb)
     dip = dsqrt(dipx**2+dipy**2+dipz**2)/0.393456d0
  enddo
  write(6,'(a,f8.3,a)') ' ice dip   = ', dip ,' D'
  write(6,*)
  deallocate (z)

  ! Equilibrate the ice lattice

  open (unit=26,file='E_if_ice.out')
  open (unit=27,file='den_if_ice.out')

  call setup_ewald(na_ice,box_ice)
  if (iamrigid) then
     call rattle_setup(nm_ice)
  endif
  boxmin = min(box_ice(1),box_ice(2),box_ice(3))
  rcut = min(rcut,0.5d0*boxmin)

  neq = min(ne/20,20000)
  neq = max(neq,500)
  write(6,'(a,f8.3,a)') ' Pre-eq Ice   = ', 1d-3*dt*neq/tofs, ' ps'
  call quick_eq(neq,r_ice,na_ice,nb,box_ice,beta,dt, &
                irun,omass,hmass,qo,26,1)
  close (unit=26)
  close (unit=27)

  open (unit=12,file='vmd_ice_start.xyz')
  call print_vmd_full(r_ice,nb,na_ice,nm_ice,box_ice,12)
  close (unit=12)

  dipx = 0.d0
  dipy = 0.d0
  dipz = 0.d0
  do j = 1,na_ice,3
     dipx = dipx + qo*r_ice(1,j,1)
     dipy = dipy + qo*r_ice(2,j,1)
     dipz = dipz + qo*r_ice(3,j,1)
     dipx = dipx - 0.5d0*qo*(r_ice(1,j+1,1)+r_ice(1,j+2,1))
     dipy = dipy - 0.5d0*qo*(r_ice(2,j+1,1)+r_ice(2,j+2,1))
     dipz = dipz - 0.5d0*qo*(r_ice(3,j+1,1)+r_ice(3,j+2,1))
  enddo
  dip = dsqrt(dipx**2+dipy**2+dipz**2)/0.393456d0
  write(6,*)
  write(6,'(a,f8.3,a)') ' ice dipeq = ', dip , ' D'

  vol = box_ice(1)*box_ice(2)*box_ice(3)
  den = tokgcm3*(wmass*dble(na_ice/3)/vol)

  write(6,'(a,f8.3)')   ' box_ice_x = ', box_ice(1)
  write(6,'(a,f8.3)')   ' box_ice_y = ', box_ice(2)
  write(6,'(a,f8.3)')   ' box_ice_z = ', box_ice(3)
  write(6,'(a,f8.3,a)') ' den_ice   = ', den , ' gcm**3'

  box_wat(2) = box_ice(2)
  box_wat(3) = box_ice(3)
  vol = box_wat(1)*box_wat(2)*box_wat(3)
  den = tokgcm3*(wmass*dble(na_wat/3)/vol)

  write(6,'(a,f8.3)')   ' box_wat_x = ', box_wat(1)
  write(6,'(a,f8.3)')   ' box_wat_y = ', box_wat(2)
  write(6,'(a,f8.3)')   ' box_wat_z = ', box_wat(3)
  write(6,'(a,f8.3,a)') ' den_wat   = ', den , ' gcm**3'
  write(6,*)

  ! Remove two bonds length from each side
  ! (to avoid overlap when putting boxes together)

  box_wat(1) = box_wat(1) - 2*reoh

  ! Place water molecules in water on a cubic starting lattice

  allocate (r_wat(3,na_wat,nb))
  r_wat(:,:,:) = 0.d0
  call cubic(r_wat(1,1,1),nm_wat,box_wat,nc_wat)
  do k = 2,nb
     do j = 1,na_wat
        r_wat(1,j,k) = r_wat(1,j,1)
        r_wat(2,j,k) = r_wat(2,j,1)
        r_wat(3,j,k) = r_wat(3,j,1)
     enddo
  enddo

  ! Equilibrate the water lattice

  open (unit=26,file='E_if_wat.out')

  call setup_ewald(na_wat,box_wat)
  if (iamrigid) then
     call rattle_setup(nm_wat)
  endif
  boxmin = min(box_wat(1),box_wat(2),box_wat(3))
  rcut = min(rcut,0.5d0*boxmin)

  neq = min(ne/10,20000)
  neq = max(neq,500)
  write(6,'(a,f8.3,a)') ' Pre-eq Water = ', 1d-3*dt*neq/tofs, ' ps'
  call quick_eq(neq,r_wat,na_wat,nb,box_wat,beta,dt, &
                irun,omass,hmass,qo,26,0)
  close (unit=26)

  open (unit=12,file='vmd_wat_start.xyz')
  call print_vmd_full(r_wat,nb,na_wat,nm_wat,box_ice,12)
  close (unit=12)

  ! Combine the boxes such that we have an interface
  ! -------------------------------------------------

  r(:,:,:) = 0.d0
  boxlxyz(1) = box_ice(1) + box_wat(1) + 2*reoh
  boxlxyz(2) = box_ice(2)
  boxlxyz(3) = box_ice(3)

  ! Move all ice particles back into box in y direction -
  ! **** MAKE SURE O and H on same side! ****

  do j =1,na_ice,3
     dx = box_ice(1)*nint(r_ice(1,j,1)/box_ice(1)-0.5d0)
     dy = box_ice(2)*nint(r_ice(2,j,1)/box_ice(2)-0.5d0)
     dz = box_ice(3)*nint(r_ice(3,j,1)/box_ice(3)-0.5d0)
     do k = 1,nb
        r_ice(1,j,k) = r_ice(1,j,k) - dx
        r_ice(1,j+1,k) = r_ice(1,j+1,k) - dx
        r_ice(1,j+2,k) = r_ice(1,j+2,k) - dx
        r_ice(2,j,k) = r_ice(2,j,k) - dy
        r_ice(2,j+1,k) = r_ice(2,j+1,k) - dy
        r_ice(2,j+2,k) = r_ice(2,j+2,k) - dy
        r_ice(3,j,k) = r_ice(3,j,k) - dz
        r_ice(3,j+1,k) = r_ice(3,j+1,k) - dz
        r_ice(3,j+2,k) = r_ice(3,j+2,k) - dz
     enddo
  enddo

  do j =1,na_wat,3
     dx = box_wat(1)*(nint(r_wat(1,j,1)/box_wat(1))-0.5d0)
     dy = box_wat(2)*(nint(r_wat(2,j,1)/box_wat(2))-0.5d0)
     dz = box_wat(3)*(nint(r_wat(3,j,1)/box_wat(3))-0.5d0)
     do k = 1,nb
        r_wat(1,j,k) = r_wat(1,j,k) - dx
        r_wat(1,j+1,k) = r_wat(1,j+1,k) - dx
        r_wat(1,j+2,k) = r_wat(1,j+2,k) - dx
        r_wat(2,j,k) = r_wat(2,j,k) - dy
        r_wat(2,j+1,k) = r_wat(2,j+1,k) - dy
        r_wat(2,j+2,k) = r_wat(2,j+2,k) - dy
        r_wat(3,j,k) = r_wat(3,j,k) - dz
        r_wat(3,j+1,k) = r_wat(3,j+1,k) - dz
        r_wat(3,j+2,k) = r_wat(3,j+2,k) - dz
     enddo
  enddo

  ! Add a boxlength to coordinates of water y axis
  ! also add one bond length to avoid overlap

  do k = 1,nb
     do j = 1,na_wat
        r_wat(1,j,k) = r_wat(1,j,k) + box_ice(1) + reoh
     enddo
  enddo

  ! Combine coordinates into one array

  do k = 1,nb
     do j = 1,na_ice
        r(1,j,k) = r_ice(1,j,k)
        r(2,j,k) = r_ice(2,j,k)
        r(3,j,k) = r_ice(3,j,k)
     enddo
     do j = 1,na_wat
        i = na_ice + j
        r(1,i,k) = r_wat(1,j,k)
        r(2,i,k) = r_wat(2,j,k)
        r(3,i,k) = r_wat(3,j,k)
     enddo
  enddo

  ! Deallocate temporary arrays

  deallocate (r_ice,r_wat)

  return
end subroutine setup_interface

subroutine quick_eq(nqe,r,na,nb,box,beta,dt,irun,omass,hmass,qo, &
                    nunit,nbaro)
  use thermostat
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Quick equilibration of water and ice structrues
  ! ------------------------------------------------------------------
  integer j,je,i,k,k1,k2,irun,na,nb,nqe,nunit,nbaro,ibaro
  real(8) r(3,na,nb),box(3),vir(3,3),vir_lf(3,3)
  real(8) beta,dt,omass,hmass,qo,qh
  real(8) v,v1,v2,v3,ran2,vrp,tk,dx,dy,dz,dist,omegan,dtfs,tps
  real(8) den,vol,wmass
  external ran2
  real(8), allocatable :: p(:,:,:),mass(:),z(:)
  real(8), allocatable :: dvdr(:,:,:),dvdr2(:,:,:)

  allocate (p(3,na,nb),mass(na),z(na))
  allocate (dvdr(3,na,nb),dvdr2(3,na,nb))

  p(:,:,:) = 0.d0
  dvdr(:,:,:) = 0.d0
  dvdr2(:,:,:) = 0.d0
  vir(:,:) = 0.d0
  vir_lf(:,:) = 0.d0
  mass(:) = 0.d0
  z(:) = 0.d0

  qh = -0.5d0*qo
  do j = 1,na,3
     mass(j) = omass      ! Oxygen
     mass(j+1) = hmass    ! Hydrogen
     mass(j+2) = hmass    ! Hydrogen
     z(j) = qo
     z(j+1) = qh
     z(j+2) = qh
  enddo

  call sample(p,na,nb,mass,beta,irun,dt)

  do je = 1,nqe

     if (je.gt.(nqe/10)) then
        ibaro = nbaro
     else
        ibaro = 0
     endif

     call evolve(p,r,v,v1,v2,v3,dvdr,dvdr2,dt,mass,na,nb, &
                 box,z,beta,vir,vir_lf,irun,ibaro)

     dtfs = dt/tofs
     tps = 1.d-3*je*dtfs

     ! Density

     if (ibaro.eq.1) then
        vol = box(1)*box(2)*box(3)
        wmass = mass(1) + mass(2) + mass(3)
        den = tokgcm3*(wmass*dble(na/3)/vol)
        write(nunit+1,*) tps,den
        write(nunit+2,'(f10.3,3f10.5)') tps,box(1),box(2),box(3)
     endif

     ! Check energy conservation

     if (mod(je,10).eq.0) then
        call kinetic(p,na,tk,mass,nb,beta)
        vrp = 0.d0
        omegan = dble(nb) / (beta)
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
                 dist = dx*dx + dy*dy + dz*dz
                 vrp = vrp + 0.5d0*mass(i)*omegan**2*dist
              enddo
           enddo
        endif
        write (nunit,'(5(1x,f13.8))') tps,(nb*v+tk+vrp)/na, &
              (vrp/na),(nb*v)/na,tk/na

     endif

     ! Momentum resampling

     if ((therm.eq.'AND').or.(therm.eq.'PRA')) then
        if (ran2(irun,0.d0,1.d0) .lt. 0.01d0) then
           call sample(p,na,nb,mass,beta,irun,dt)
        endif
     endif
  enddo

  deallocate(p,mass,z,dvdr,dvdr2)

  return
end subroutine quick_eq
