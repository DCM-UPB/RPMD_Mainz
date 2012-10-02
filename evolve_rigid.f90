subroutine evolve_rig_cl(p,r,v,v_lf,v_hf,dvdr,dvdr2,dt,mass,na,nb, &
                         boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
  use thermostat
  use barostat
  implicit none
  ! ------------------------------------------------------------------
  ! Rigid body classical evolution by time step dt
  ! ------------------------------------------------------------------
  integer na,nb,nm,irun,nbaro,i,j,k
  real(8) p(3,na,nb),r(3,na,nb),dvdr(3,na,nb),dvdr2(3,na,nb)
  real(8) mass(na),z(na),vir(3,3),vir_lf(3,3),vir_hf(3,3)
  real(8) beta,tv,tq1,tq2,tvxyz(3)
  real(8) dt,halfdt,boxlxyz(3),v,v_lf,v_hf
  real(8), allocatable :: rold(:,:)

  allocate (rold(3,na))

  dvdr2(:,:,:) = 0.d0
  halfdt = 0.5d0*dt
  nm = na/3

  ! Evolve : Velocity Verlet Stage 1
  
  if (therm.eq.'PRG') then
     call parinello_therm(p,mass,ttau,na,nb,halfdt,irun,beta)
  else if (therm.eq.'PRL') then
     call parinello_therm_loc(p,mass,ttau,na,nb,halfdt,irun,beta)
  endif

  do k = 1,nb
     rold(:,:) = r(:,:,k)
     do j = 1,na
        do i = 1,3
           p(i,j,k) = p(i,j,k) - halfdt*dvdr(i,j,k)
           r(i,j,k) = r(i,j,k) + p(i,j,k)*dt/mass(j)
        enddo
     enddo

     ! Rattle Stage 1

     call rattle_s1(r(1,1,k),rold,p(1,1,k),nm,mass,dt, &
                    dvdr(1,1,k),vir_hf,na)
  enddo
  
  ! Intermolecular Forces

  call forces(r,v,dvdr,nb,na,boxlxyz,z,vir_lf,1)

  ! Evolve : Velocity Verlet Stage 2

  p(:,:,:) = p(:,:,:) - halfdt*dvdr(:,:,:)

  ! Rattle Stage 2

  do k = 1,nb
     call rattle_s2(r(1,1,k),p(1,1,k),nm,mass,dt, &
                    dvdr(1,1,k),vir_hf,na)
  enddo
  
  if (therm.eq.'PRG') then
     call parinello_therm(p,mass,ttau,na,nb,halfdt,irun,beta)
  else if (therm.eq.'PRL') then
     call parinello_therm_loc(p,mass,ttau,na,nb,halfdt,irun,beta)
  endif

  ! Barostat
  ! (note:// COMs scaled therefore does not affect RATTLE)

  vir(:,:) = vir_lf(:,:) + vir_hf(:,:)
  if (nbaro.eq.1) then
     if (baro.eq.'BER') then
        call virial_ke(r,dvdr,dvdr2,tv,tvxyz,tq1,tq2,beta,na,nb)
        call beren_driver(vir,tv,tvxyz,dt,r,boxlxyz,na,nb)
     else if (baro.eq.'MCI') then
        call mc_baro(r,dvdr,dvdr2,vir,v,z,beta,boxlxyz,na,nb,irun)
     else
        write(6,*) ' ** Invalid barostat **'
        stop
     endif
  endif

  vir(:,:) = vir_lf(:,:) + vir_hf(:,:)
  v_lf = v
  v_hf = 0.d0

  deallocate (rold)

  return
end subroutine evolve_rig_cl

subroutine evolve_rig_pi(p,r,v,v_lf,v_hf,dvdr,dvdr2,dt,mass,na,nb, &
                         boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
  use thermostat
  use barostat
  implicit none
  ! ------------------------------------------------------------------
  ! Rigid Body RPMD evolution by time step dt
  ! ------------------------------------------------------------------
  integer na,nb,nm,irun,nbaro,nstep,i,j,k,is
  real(8) p(3,na,nb),r(3,na,nb),dvdr(3,na,nb),dvdr2(3,na,nb)
  real(8) mass(na),z(na),boxlxyz(3),vir(3,3),vir_lf(3,3),vir_hf(3,3)
  real(8) vir_tmp(3,3),dt,v,v_lf,v_hf,beta,tv,tvxyz(3),tq1,tq2
  real(8) halfdt,smalldt,halfsmalldt
  real(8), allocatable :: rold(:,:)

  allocate (rold(3,na))

  halfdt = 0.5d0 * dt
  nstep = nb
  smalldt = dt/dble(nstep)
  halfsmalldt = 0.5d0*smalldt
  nm = na / 3

  if (therm.eq.'PRG') then
     call parinello_therm(p,mass,ttau,na,nb,halfdt,irun,beta)
  else if (therm.eq.'PRL') then
     call parinello_therm_loc(p,mass,ttau,na,nb,halfdt,irun,beta)
  endif

  ! Update momenta using intermolecular forces

  p(:,:,:) = p(:,:,:) - halfdt * dvdr(:,:,:)

  ! Loop over the small time steps....

  call bead_forces(r,na,nb,dvdr2,beta,mass)
  do is = 1,nstep

  ! Evolve : Velocity Verlet Stage 1

     do k = 1,nb
        do j = 1,na
           do i = 1,3
              rold(i,j) = r(i,j,k)
              p(i,j,k) = p(i,j,k) - halfsmalldt*dvdr2(i,j,k)
              r(i,j,k) = r(i,j,k) + p(i,j,k)*smalldt/mass(j)
           enddo
        enddo
        
        ! Rattle Stage 1

        call rattle_s1(r(1,1,k),rold,p(1,1,k),nm,mass, &
                       smalldt,dvdr2(1,1,k),vir_hf,na)
     enddo

     ! Bead forces

     call bead_forces(r,na,nb,dvdr2,beta,mass)

     ! Evolve : Velocity Verlet Stage 2

     p(:,:,:) = p(:,:,:) - halfsmalldt*dvdr2(:,:,:)
     
     ! Rattle Stage 2

     do k = 1, nb
        call rattle_s2(r(1,1,k),p(1,1,k),nm,mass, &
                       smalldt,dvdr2(1,1,k),vir_hf,na)
     enddo
  enddo

  ! Intermolecular Forces

  call forces(r,v,dvdr,nb,na,boxlxyz,z,vir_lf,1)

  ! Velocity Verlet Stage 2 under intermolecular forces

  p(:,:,:) = p(:,:,:) - halfdt * dvdr(:,:,:)

  ! Rattle stage 2

  vir_hf(:,:) = 0.d0
  do k = 1,nb
     call rattle_s2(r(1,1,k),p(1,1,k),nm,mass, &
                    dt,dvdr(1,1,k),vir_tmp,na)
     vir_hf(:,:) = vir_hf(:,:) + vir_tmp(:,:)
  enddo

  if (therm.eq.'PRG') then
     call parinello_therm(p,mass,ttau,na,nb,halfdt,irun,beta)
  else if (therm.eq.'PRL') then
     call parinello_therm_loc(p,mass,ttau,na,nb,halfdt,irun,beta)
  endif

  vir_hf(:,:) = vir_hf(:,:)/dble(nb)
  vir(:,:) = vir_lf(:,:) + vir_hf(:,:)

  ! Barostat

  dvdr2(:,:,:) = 0.d0
  if (nbaro.eq.1) then
     if (baro.eq.'BER') then
        call virial_ke(r,dvdr,dvdr2,tv,tvxyz,tq1,tq2,beta,na,nb)
        call beren_driver(vir,tv,tvxyz,dt,r,boxlxyz,na,nb)
     else if (baro.eq.'MCI') then
        call mc_baro(r,dvdr,dvdr2,vir,v,z,beta,boxlxyz,na,nb,irun)
     else
        write(6,*) ' ** Invalid barostat **'
        stop
     endif
  endif

  v_lf = v
  v_hf = 0.d0

  deallocate (rold)

  return
end subroutine evolve_rig_pi

subroutine bead_forces(r,na,nb,dvdrb,beta,mass)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Calculates forces between the beads
  ! ------------------------------------------------------------------
  integer na, nb,i,k,k1,k2,k3
  real(8) r(3,na,nb),dvdrb(3,na,nb),mass(na)
  real(8) omegan,beta,dx1,dx2,dx3

  omegan = nb / (beta * hbar)
  omegan = omegan * omegan
  do k = 1, nb
     do i = 1, na
        do k1 = 1, 3
           dvdrb(k1,i,k) = 0.d0
        enddo
     enddo
  enddo

  if (nb.gt.1) then
     do i = 1, na
        do k = 1, nb
           if (k.eq.1) then
              k1 = k
              k2 = nb
              k3 = 2
           else if (k.eq.nb) then
              k1 = k
              k2 = k - 1
              k3 = 1
           else
              k1 = k
              k2 = k - 1
              k3 = k + 1
           endif
           dx1 = r(1,i,k1)
           dx2 = r(1,i,k2)
           dx3 = r(1,i,k3)
           dvdrb(1,i,k) = mass(i)*omegan*(2.d0*dx1-dx2-dx3)
           dx1 = r(2,i,k1)
           dx2 = r(2,i,k2)
           dx3 = r(2,i,k3)
           dvdrb(2,i,k) = mass(i)*omegan*(2.d0*dx1-dx2-dx3)
           dx1 = r(3,i,k1)
           dx2 = r(3,i,k2)
           dx3 = r(3,i,k3)
           dvdrb(3,i,k) = mass(i)*omegan*(2.d0*dx1-dx2-dx3)
        enddo
     enddo
  endif
  return
end subroutine bead_forces

subroutine evolve_rigid_pi_RPMDDFT(p,r,v,v_lf,v_hf,dvdr,dvdr2,dt,mass,na,nb, &
                         boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
  use thermostat
  use barostat
#ifdef PARALLEL_BINDING
  include 'mpif.h' !parallel
#endif
  implicit none
  ! ------------------------------------------------------------------
  ! Rigid Body RPMD evolution by time step dt
  ! ------------------------------------------------------------------
  integer na,nb,nm,irun,nbaro,nstep,i,j,k,is
  real(8) p(3,na,nb),r(3,na,nb),dvdr(3,na,nb),dvdr2(3,na,nb)
  real(8) mass(na),z(na),boxlxyz(3),vir(3,3),vir_lf(3,3),vir_hf(3,3)
  real(8) vir_tmp(3,3),dt,v,v_lf,v_hf,beta,tv,tvxyz(3),tq1,tq2
  real(8) halfdt,smalldt,halfsmalldt
  real(8), allocatable :: rold(:,:)

  halfdt = 0.5d0 * dt
  nstep = nb
  smalldt = dt/dble(nstep)
  halfsmalldt = 0.5d0*smalldt
  nm = na / 3

#ifdef PARALLEL_BINDING
	call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr)
!write(*,*) "myid iin evolve:",myid
#endif
#ifdef PARALLEL_BINDING
if(myid.eq.0) then     
#endif  
  allocate (rold(3,na))
  if (therm.eq.'PRG') then
     call parinello_therm(p,mass,ttau,na,nb,halfdt,irun,beta)
  else if (therm.eq.'PRL') then
     call parinello_therm_loc(p,mass,ttau,na,nb,halfdt,irun,beta)
  endif

  ! Update momenta using intermolecular forces

  p(:,:,:) = p(:,:,:) - halfdt * dvdr(:,:,:)

  ! Loop over the small time steps....

  call bead_forces(r,na,nb,dvdr2,beta,mass)
  do is = 1,nstep

  ! Evolve : Velocity Verlet Stage 1

     do k = 1,nb
        do j = 1,na
           do i = 1,3
              rold(i,j) = r(i,j,k)
              p(i,j,k) = p(i,j,k) - halfsmalldt*dvdr2(i,j,k)
              r(i,j,k) = r(i,j,k) + p(i,j,k)*smalldt/mass(j)
           enddo
        enddo
        
        ! Rattle Stage 1

        call rattle_s1(r(1,1,k),rold,p(1,1,k),nm,mass, &
                       smalldt,dvdr2(1,1,k),vir_hf,na)
     enddo

     ! Bead forces

     call bead_forces(r,na,nb,dvdr2,beta,mass)

     ! Evolve : Velocity Verlet Stage 2

     p(:,:,:) = p(:,:,:) - halfsmalldt*dvdr2(:,:,:)
     
     ! Rattle Stage 2

     do k = 1, nb
        call rattle_s2(r(1,1,k),p(1,1,k),nm,mass, &
                       smalldt,dvdr2(1,1,k),vir_hf,na)
     enddo
  enddo

#ifdef PARALLEL_BINDING
endif
	call MPI_bcast(r,SIZE(r),MPI_real8,0,MPI_COMM_WORLD,ierr)
	call MPI_bcast(p,SIZE(p),MPI_real8,0,MPI_COMM_WORLD,ierr)
	call MPI_bcast(dvdr,SIZE(dvdr),MPI_real8,0,MPI_COMM_WORLD,ierr)
	call MPI_bcast(dvdr2,SIZE(dvdr2),MPI_real8,0,MPI_COMM_WORLD,ierr)
	call MPI_bcast(vir_lf,SIZE(vir_lf),MPI_real8,0,MPI_COMM_WORLD,ierr)
#endif

  ! Intermolecular Forces

  call forces(r,v,dvdr,nb,na,boxlxyz,z,vir_lf,9)

#ifdef PARALLEL_BINDING
if(myid.eq.0) then
#endif
  ! Velocity Verlet Stage 2 under intermolecular forces

  p(:,:,:) = p(:,:,:) - halfdt * dvdr(:,:,:)

  ! Rattle stage 2

  vir_hf(:,:) = 0.d0
  do k = 1,nb
     call rattle_s2(r(1,1,k),p(1,1,k),nm,mass, &
                    dt,dvdr(1,1,k),vir_tmp,na)
     vir_hf(:,:) = vir_hf(:,:) + vir_tmp(:,:)
  enddo

  if (therm.eq.'PRG') then
     call parinello_therm(p,mass,ttau,na,nb,halfdt,irun,beta)
  else if (therm.eq.'PRL') then
     call parinello_therm_loc(p,mass,ttau,na,nb,halfdt,irun,beta)
  endif

  vir_hf(:,:) = vir_hf(:,:)/dble(nb)
  vir(:,:) = vir_lf(:,:) + vir_hf(:,:)

  ! Barostat

  dvdr2(:,:,:) = 0.d0
  if (nbaro.eq.1) then
     if (baro.eq.'BER') then
        call virial_ke(r,dvdr,dvdr2,tv,tvxyz,tq1,tq2,beta,na,nb)
        call beren_driver(vir,tv,tvxyz,dt,r,boxlxyz,na,nb)
     else if (baro.eq.'MCI') then
        call mc_baro(r,dvdr,dvdr2,vir,v,z,beta,boxlxyz,na,nb,irun)
     else
        write(6,*) ' ** Invalid barostat **'
        stop
     endif
  endif

  v_lf = v
  v_hf = 0.d0

  deallocate (rold)
#ifdef PARALLEL_BINDING
endif
#endif
  return
end subroutine evolve_rigid_pi_RPMDDFT
