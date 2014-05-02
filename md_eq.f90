subroutine md_eq(ne,p,r,dvdr,dvdr2,na,nb,boxlxyz,z,beta, &
                 dt,mass,irun)
  use thermostat
  use barostat
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Equilibration Routine
  ! ------------------------------------------------------------------
  integer na,nb,ne,npre_eq,irun,jout,je,k,k1,k2,i,nm,nprog,nbaro
  integer jq,nq,nqprog,nden,nctot,nbond
  real(8) boxlxyz(3),beta,v,v1,v2,v3,dt,px,py,pz
  real(8) p(3,na,nb),dvdr(3,na,nb),dvdr2(3,na,nb),z(na),mass(na)
  real(8) r(3,na,nb),vir(3,3),vir_lf(3,3),tvxyz(3)
  real(8) omegan,dx,dy,dz,dist,denav
  real(8) vrp,vol,tps,tk,thresh,den,wmass,tq1,tq2
  real(8) dtfs,dtqt,tv,pres,ran2,gaussian,dr,rmsq,temp,ttaufs
  real(8), allocatable :: rst(:,:),rcm(:,:)
  logical iamcub
  external gaussian, ran2
  character(len=3) ens
  common /ensemble/ ens
  common /symmetry/ iamcub
  common /constraint/ nctot,nbond
  common /inp/ npre_eq
  common /thinp/ ttaufs

  nbaro = 0
  if (ens.eq.'NPT') then
     nbaro = 1
  endif

  vir(:,:) = 0.d0
  vir_lf(:,:) = 0.d0
  denav = 0.d0
  nden = 0

  allocate (rst(3,na/3),rcm(3,na/3))

  ! Water starting centroid positions

  nm = na / 3
  call center_water(r,rst,nm,nb)

  ! Define some useful local constants

  dtfs = dt/tofs
  nprog = max(ne/10,1)
  jout = max(ne/100000,1)
  vol = boxlxyz(1)*boxlxyz(2)*boxlxyz(3)

  ! Need to thermostat NPT calculations more strongly

  if (ens.eq.'NVT') then
    if (ttaufs.gt.0.d0 .and. (therm.eq.'AND' .or. therm.eq.'PRA')) then
        thresh = dtfs/ttaufs
    else
        thresh = 1.d0/dsqrt(dble(ne))
        thresh = max(0.01d0,thresh)
    end if
  else
     thresh = 0.01d0
  endif

  if (ens.ne.'NVT') then

     ! Quick NVT equilibration

     open (unit=45,file='E_pre_eq.out')

     nq = min(20000,ne/10)
     if (npre_eq.ne.0) then
        write(6,*) 'Using npre_eq pre-equilibration steps from input:', npre_eq
        nq = npre_eq
     endif
     nqprog = max(nq/10,1)
     dtqt = 1.d-3*dble(nq)*dtfs
     write(6,'(a,f6.2,a)') ' * NVT pre-equilibration for: ', &
                               dtqt, ' ps'
     do jq = 1,nq
        istep=jq

        ! Evolve system

        call evolve(p,r,v,v1,v2,v3,dvdr,dvdr2,dt,mass,na,nb, &
                    boxlxyz,z,beta,vir,vir_lf,irun,0)

        tps = 1.d-3*jq*dtfs

        ! Check energy conservation

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
        write (45,'(5(1x,f13.8))') tps,(nb*v+tk+vrp)/na, &
              (vrp/na),(nb*v)/na,tk/na

        ! time for a collision with the heat bath?

        if ((therm.eq.'AND').or.(therm.eq.'PRA')) then
           if (ran2(irun,0.d0,1.d0).lt.thresh) then
              call sample(p,na,nb,mass,beta,irun,dt)
           endif
        endif

        if (mod(jq,nqprog).eq.0) then
           write(6,*) 10*(jq/nqprog), ' %'
        endif
     enddo

     close (unit=45)

  endif

  ! Full equilibration

  ! Open output files

  open (unit=10,file='V_eq.out')
  open (unit=20,file='temp_eq.out')
  open (unit=30,file='RMS_eq.out')
  if (ens.eq.'NPT') then
     if (iamcub) then
        open (unit=35,file='density_eq.out')
        open (unit=36,file='density_av_eq.out')
     else
        open (unit=35,file='density_eq.out')
        open (unit=36,file='density_av_eq.out')
        open (unit=37,file='box_len_eq.out')
     endif
  endif
  open (unit=40,file='E_eq.out')
  open (unit=50,file='pressure_eq.out')
  open (unit=55,file='virial_ke_eq.out')

  write(6,*) '* Beginning Equilibration'
  do je = 1,ne
     istep = je

     ! Evolve the system by dt

     call evolve(p,r,v,v1,v2,v3,dvdr,dvdr2,dt,mass,na,nb, &
                 boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
     
     tps = 1.d-3*je*dtfs

     ! Pressure and Virial Kinetic Energy
     call virial_ke(r,dvdr,dvdr2,tv,tvxyz,tq1,tq2,beta,na,nb,mass) !!!!!
     call pressure(pres,vir,tv,na,boxlxyz)

     ! Density

     vol = boxlxyz(1)*boxlxyz(2)*boxlxyz(3)
     wmass = mass(1) + mass(2) + mass(3)
     den = tokgcm3*(wmass*dble(nm)/vol)
     
     if (ens.eq.'NPT') then
        
        ! Current Density

        if (mod(je,10).eq.0) then
           if (iamcub) then
              write(35,*) tps , den
           else
              write(35,*) tps , den 
              write(37,'(f10.3,3f10.5)') tps, boxlxyz(1), &
                                         boxlxyz(2),boxlxyz(3)
           endif
        endif
        
        ! Running average density

        if ((je.gt.(ne/10)).or.(je.gt.100000)) then
           denav = denav + den
           nden = nden + 1
           if (mod(je,1000).eq.0) then
              write(36,*) tps , denav/dble(nden)
           endif
        endif
     endif
     
     ! Pressure

     if (mod(je,10).eq.0) then
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
        
        temp = toK*(2.d0*tk/(dble(3*na-nctot-3)*dble(nb)))/dble(nb)
        write (20,*) tps,temp
        write (10,'(3f14.6)') tps,v,toKjmol*(v/na)
        write (40,'(5(1x,f13.8))') tps,(nb*v+tk+vrp)/na, &
                (vrp/na),(nb*v)/na,tk/na
        write(55,*)tps,tv*tokJmol/dble(nm)

     endif

     if (mod(je,20).eq.0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!falsch, da center_water nicht mit Massen gewichtet, da für Impulse gemacht
        ! Output RMS distance moved by centroid COM

        rmsq = 0.d0
        call center_water(r,rcm,nm,nb)
        
        do i = 1,nm
           do k = 1,3
              dr = rcm(k,i)-rst(k,i)
              rmsq = rmsq + dr*dr
           enddo
        enddo
        rmsq = rmsq / dble(nm)
        write (30,'(2f12.6)') tps, dsqrt(rmsq)
     endif

     ! Every 100 time steps output the current configuration  !!GLE 10000->100, append rather than overwrite

     if (mod(je,10000).eq.0) then
        open (unit=46,file='vmd_current.xyz', POSITION='APPEND')
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
  close (unit=10)
  close (unit=20)
  close (unit=30)
  if (ens.eq.'NPT') then
     if (iamcub) then
        close (unit=35)
        close (unit=36)
     else
        close (unit=35)
        close (unit=36)
        close (unit=37)
     endif
  endif
  close (unit=40)
  close (unit=50)
  close (unit=55)

  write (6,*)'* Equilibration complete. '
  if (ens.eq.'NPT') then
     write (6,'(a,f8.4,a)') 'Average density = ', &
                            denav/dble(nden),' g cm^-3'
  endif
  call flush (6)

  deallocate (rst,rcm)

  return
end subroutine md_eq

