program qmd
    use thermostat
    use barostat
    use gle
    use f_env
    implicit none
    include 'globals.inc'
#ifdef PARALLEL_BINDING
    include 'mpif.h' !parallel
#endif
    !--------------------------------------------------------------------
    !Classical/RPMD/PACMD Simulation Program for Flexible water
    !--------------------------------------------------------------------
    integer nt,ne,nb,m,ng,npre_eq,pt,pb,irun,nm,na,narg,iargc,nbdf1,nbdf2,nbdf3,rctdk,i,j,k,ptd
    integer nbr,mts,nlat,itcf(3),itst(3),ncellxyz(3),print(4),ntherm,iskip!,sizeMPI(9)

    ! used for reftrj and RPMD-DFT
    integer reftraj,rpmddft,ierr,rpmddfthelp,myid,numid,epsr_update,rpmde3b
    character(len=35) CP2K_path
    ! r_traj(xyz,molecules,nb,reftraj)
    real(8), allocatable :: r_traj(:,:,:,:),boxlxyz_traj(:,:)
    character(len=240) line
    common /reftraj/ reftraj,line
  
    integer nc_ice(3),nc_wat(3),nm_ice,nm_wat,nctot,nbond,vacfac,na_expected
    integer natoms ! Used in reftraj mode for obtaining read in atoms for consitency check
    real(8) temp,rho,dtfs,ecut,test,beta,dt,dtps,boxmin,pres
    real(8) teqm,tsim,trdf,gaussian,wmass,rcut,om,ttaufs
    real(8) taufs,v,vew,vlj,vint,pi,vir(3,3),cell(3,3)
    real(8) qo,qh,alpha,oo_sig,oo_eps,oo_gam,theta,reoh,thetad
    real(8) apot,bpot,alp,alpb,wm,wh,omass,hmass,dmass,sig,boxlxyz(3),vdum
    real(8) box_ice(3),box_wat(3),rcut_old
    real(8) vave
    real(8), allocatable :: mass(:),z(:),r(:,:,:)
    real(8), allocatable :: p(:,:,:),dvdr(:,:,:),dvdr2(:,:,:)
    character(len=25) filename
    character(len=4) type
    character(len=3) lattice,ens,therm_backup
    character(len=3) isotope
    logical iamcub,iamrigid
    logical epsr
    logical H2O,HDO,D2O
    external gaussian

    namelist/input/ens,isotope,temp,pres,rho,lattice,vacfac,iamcub,dtfs, &
    ecut,nt,ne,npre_eq,ntherm,nb,m,ng,print,reftraj,iskip,pt,pb, &
    ncellxyz,irun,itcf,itst,rcut, &
    type,therm,ttaufs,baro,taufs,mts,om,nbdf1,nbdf2,sig,rpmddft, &
    CP2K_path,nbdf3,rctdk,epsr,epsr_update,rpmde3b
    namelist/param/ wmass,omass,hmass,dmass,qo,alpha,oo_sig,oo_eps,oo_gam, &
    thetad,reoh,apot,bpot,alp,alpb,wm,wh

    common /path_i/ om,type
    common /beaddiabatic/ nbdf1,nbdf2
    common /geometry/ theta,reoh
    common /multiple_ts/ mts
    common /oo_param/ oo_eps,oo_sig,oo_gam,rcut
    common /ew_param/ alpha,ecut,wm,wh
    common /intra_param/ apot,bpot,alp,alpb
    common /correct/ sig
    common /lattice/ ncellxyz,nlat
    common /symmetry/ iamcub
    common /structure/ iamrigid
    common /ensemble/ ens
    common /constraint/ nctot,nbond
    common /inp/ npre_eq
    common /thinp/ ttaufs
    common /RPMDDFT/ rpmddft,nbdf3,rctdk
    common /EPSR/ epsr, epsr_update
    common /E3B/ rpmde3b
    common /isotope/ H2O,D2O,HDO
    common /global_input/ print,pb,pt,ng

  
    vacfac = 1
    ntherm = 0
    ttaufs = 0.d0
    iskip = 1
    rpmddft = 0
    rpmde3b = 0
    reftraj = 0
    npre_eq = 0
    cell(:,:) = 0.d0
    CP2K_path =""
    nbdf3 = 0
    rctdk = 0
    itst(3) = 0 ! itst 3 argument is for RMS calculation
    epsr = .false.
    epsr_update = 1
    isotope = "H2O"
    dmass = 0.d0
    print(4) = 0 !print 4 argument is for printing mol-dipole
    H2O = .false.
    D2O = .false.
    HDO = .false. 

		! Using CP2K parallel?
#ifdef PARALLEL_BINDING
#ifdef CP2K_BINDING
	call cp_init_cp2k(1,ierr)
	!call MPI_init(ierr) 
	!call cp_init_cp2k(0,ierr)
#endif
#ifndef CP2K_BINDING
	write(*,*) "Parallel Binding using CP2K seriell"
	call MPI_init(ierr) 		
#endif
		call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr)
 		call MPI_COMM_SIZE( MPI_COMM_WORLD, numid, ierr)
#endif

#ifndef PARALLEL_BINDING
#ifdef CP2K_BINDING
		call cp_init_cp2k(0,ierr)
#endif
#endif

#ifdef PARALLEL_BINDING
if(myid.eq.0) then
#endif
    write(6,*)  '-------------------------------------------'
    write(6,*)  '            Flexible Water Code            '
    write(6,*)  '-------------------------------------------'
    write(6,*)

    ! Read the input file, given by the first argument to the program

    narg = iargc()
    if ( narg .gt. 3 .or. narg .lt. 2 ) then
        write(6,*)'* Program usage:'
        write(6,*)' qmd.x <input file> <param file> <eq file> '
        stop
    endif

    ! Input file

    call getarg (1,filename)
    open(60,file=filename)
    read(60,input)
    close (unit=60)

    ! Parameters File

    call getarg (2,filename)
    open(61,file=filename)
    read(61,param)
    close (unit=61)

    ! Rigid Simulation?
    ! ------------------

    if ((apot.eq.0.d0).and.(bpot.eq.0.d0)) then
        iamrigid = .true.
        if ((type.eq.'ACMD').and.(nb.ne.1)) then
            write(6,*) '* main : incompatible options: ACMD and rigid'
            stop
        endif
    else
        iamrigid = .false.
    endif

    ! Temperature, bond angles and times
    ! ------------------------------------

    pi = dacos(-1.d0)
    patm = pres/tobar
    theta = thetad*(pi/180.d0)
    beta = 1.d0/(3.166829d-6*temp)   ! (1/k_b (Hartree) *temp)
    dt = tofs*dtfs
    dtps = 1.d-3*dtfs
    teqm = ne*dtps
    tsim = nt*dtps
    trdf = ng*dtps
    rcut = rcut / ToA

    write(6,60) isotope,ens,temp,pres,dtfs,teqm,tsim,trdf,m,mts,irun, &
    rho,ncellxyz(1),ncellxyz(2),ncellxyz(3)
60  format (' Simulation Parameters:' /1x &
    '-----------------------' /1x &
    'isotope =',a9 /1x &
    'ens    = ',a9 /1x &
    'temp   = ',f9.2,' K'/1x, &
    'pres   = ',f9.2,' bar'/1x &
    'dtfs   = ',f9.3,' fs'/1x, &
    'teqm   = ',f9.3,' ps'/1x, &
    'tsim   = ',f9.3,' ps'/1x, &
    'trdf   = ',f9.3,' ps'/1x, &
    'm      = ',i9,' trajectories'/1x, &
    'mts    = ',i9,' multiple time steps' /1x, &
    'irun   = ',i9 /1x, &
    'rho    = ',f9.3,' g/cm**3'/1x,  &
    'ncellx = ',i9,' unit cells'/1x, &
    'ncelly = ',i9,' unit cells'/1x, &
    'ncellz = ',i9,' unit cells') 
    

    ! ACMD frequencies
    ! -----------------
    ! Set omega, the ring-polymer harmonic frequencies
    ! for ACMD simulation.

    if (nb.gt.1) then
        if ((type.eq.'ACMD').and.(nb.gt.1)) then
            if (om.eq.0.d0) then
                om = dble(nb)**(dble(nb)/(dble(nb-1)))/beta
            endif
        else
            om = 0.d0
        endif
    else
        om = 0.d0
    endif
 

    ! Thermostat setup
    ! -----------------

    if (therm.eq.'PRA') then
        ! Parinello Local Thermostat for ACMD
        write(6,*) 'therm  = PARINELLO-A '
        ttau = 1.d0/(2.d0*om)
        ttau = max(41.341373d0*ttaufs,ttau)
        ttaufs = ttau/41.341373d0
        if ((type.ne.'ACMD').or.(nb.eq.1)) then
            write(6,*) '*** PARINELLO-A ONLY SUITABLE FOR ACMD ***'
        endif
        write(6,56) ttaufs
    else if (therm.eq.'PRG') then
        write(6,*) 'therm  = PARINELLO-G '
        ttau = 41.341373d0*ttaufs
        write(6,56) ttaufs
    else if (therm.eq.'PRL') then
        write(6,*) 'therm  = PARINELLO-L '
        ttau = 41.341373d0*ttaufs
        write(6,56) ttaufs
    else if (therm.eq.'AND') then
        write(6,*) 'therm  =  ANDERSEN '
    else if (therm.eq.'PRO') then
        write(6,*) 'therm  = PARINELLO-OPTIMIZED '       !!PRO
        ttau = 41.341373d0*ttaufs
        write(6,56) ttaufs
    else if (therm.eq.'PRM') then
        write(6,*) 'therm  = PARINELLO-OPTIMIZED GLOBAL'
        ttau = 41.341373d0*ttaufs
        write(6,56) ttaufs
    else if (therm.eq.'GLE') then                       !!GLE
        write(6,*) 'therm  =  GENERALIZED LANGEVIN '    !!GLE
        ttau = 41.341373d0*ttaufs
        write(6,56) ttaufs
    else
        write(6,*) 'therm  =      NONE '
    endif
56  format (' ttau   = ',f9.3' fs')

    ! Barostat setup
    ! ---------------

    if (ens.eq.'NPT') then
        if (baro.eq.'BER') then

            ! Berendsen Barostat

            taub = 41.341373d0*taufs
            write(6,*) 'baro   = BERENDSEN '
            write(6,59) 1d-3*taufs
        else if (baro.eq.'MCI') then

            ! Stochastic MC Barostat

            taub = taufs
            dv = 200.d0
            nacc = 0
            nmove = 0
            write(6,*) 'baro   =        MC '
            write(6,58) dv
        else
            write(6,*) 'Invalid choice of barostat', baro
            stop
        endif
    endif
59  format (' taufs  = ',f9.3,' ps')
58  format (' dv     = ',f9.3)

    if (epsr.eqv..true.) then
        write (6,*) "Using EPSR"
#ifdef EPSR_STARTUP_READ
        write (6,*) "Reading EPSR on STARTUP"
#else
        write (6,*) "NOT reading EPSR on STARTUP"
#endif
    else
        write (6,*) "NOT using EPSR"
    endif

    ! Initialize the random number generator:
    ! ---------------------------------------

    irun = -irun
    test = gaussian(irun,1.d0)

    ! PI Options
    ! ------------
    if (nb.gt.1) then
        call setup_pi(nb,beta,nbdf1,nbdf2,sig)
    endif

    ! Number of molecules and starting positions :
    ! --------------------------------------------

    if (narg.eq.3) then

        ! Read from file

        call getarg(3,filename)
        open (61, file = filename)

        if (reftraj.eq.0) then
            ! Third argument is equilibrium file
            read(61,*) nm,na,nbr

            allocate(r(3,na,nb))

            r(:,:,:) = 0.d0
            if (nb.eq.nbr) then
                read(61,*) boxlxyz(1),boxlxyz(2),boxlxyz(3)
                read(61,*) r(:,:,:)
                write(6,*)'* Initialized from file: ',filename
            else

                ! This allows reading of classical equilibrated
                ! configurations to start PI trajectories.

                write(6,*) ' Note:// all beads will start at centroid'

                ! Read first bead coordinate

                read(61,*) boxlxyz(1),boxlxyz(2),boxlxyz(3)
                !do j = 1,na
                !    read(61,*) r(1,j,1),r(2,j,1),r(3,j,1)
                !enddo
                read(61,*) r(:,:,1)
                ! Copy to all other beads

                do k = 2,nb
                    do j = 1,na
                        r(1,j,k) = r(1,j,1)
                        r(2,j,k) = r(2,j,1)
                        r(3,j,k) = r(3,j,1)
                    enddo
                enddo
            endif
            close (unit=61)
            na_expected = 3.0*ncellxyz(1)*ncellxyz(2)*ncellxyz(3)
            if (na.ne.na_expected) then
                write(6,*) "The restart file contains a different number of atoms than expected:", na, "vs.", na_expected
                stop
            endif
        else
            ! Thrid argument is trajectories file
            write(6,*)
            write(6,*) 'The current reftraj ID is: ', reftraj
            write(6,*) 'And is loaded from file: ', filename
            write(6,*)

            !call setup_box_size(lattice,rho,nm,boxlxyz,wmass)
            nm = ncellxyz(1)*ncellxyz(2)*ncellxyz(3)
            na = 3*nm

            allocate(r(3,na,nb))
            r(:,:,:) = 0.d0

            allocate(r_traj(3,na,nb,reftraj))
            r_traj(:,:,:,:) = 0.d0
            allocate(boxlxyz_traj(3,reftraj))

            do k = 1, nb
                do i = 1, reftraj
                    read(61,*) natoms
                    write (6,*) "ATOMS", natoms
                    if (natoms.gt.na) then
                        write(6,*) "natoms from reftraj file too big! You must set a bigger cell in the input file !!! Aborting..."
                        stop
                    endif
                    if (natoms.lt.na) then
                        write(6,*) "WARNING: overwriting na=", na, "in reftraj mode with", natoms, "!"
                        na = natoms
                        nm = na/3
                    endif

                    ! setup_ewald is not needed, ecut, wrcut, walpha, rkmax and kmax
                    !   are all independent of boxlxyz, so just set that
                    read(61,*) line,boxlxyz_traj(:,i)
                    boxlxyz_traj(:,i) = boxlxyz_traj(:,i) * (1.0d0/toA)
                    do j = 1, na
                        ! line will get the kind (O or H)
                        read(61,*) line, r_traj(:,j,k,i)
                        r_traj(:,j,k,i) = r_traj(:,j,k,i) * (1.0d0/toA)
                    enddo
                enddo
            enddo
            ! use physical default values
            boxlxyz(:) = boxlxyz_traj(:,1)
            r(:,:,:)   = r_traj(:,:,:,1)
        endif

    else
        if (lattice.eq.'INT') then

            ! Setup an ice-water interface

            nc_ice(:) = ncellxyz(:)
            iamcub = .false.
            call setup_int_size(0.920d0,0.990d0,box_ice,box_wat, &
            nc_ice,nc_wat,nm_ice,nm_wat,wmass)
            nm = nm_wat + nm_ice
            na = 3*nm
            rcut_old = rcut
            allocate(r(3,na,nb))
            r(:,:,:) = 0.d0

            therm_backup=therm
            therm='AND'

            call setup_interface(r,na,ne,irun,nm_ice,nm_wat,box_ice, &
            box_wat,boxlxyz,nc_ice,nc_wat, &
            wmass,omass,hmass,dmass,qo,beta,dt,nb)

            therm=therm_backup
            rcut = rcut_old
            write(6,*)  '-------------------------------------------'
        !------Vacuum preparation----------
        else if (lattice.eq.'VAC') then
            write(6,*)
            write(6,*) 'Beginning Vacuum preparation'
            iamcub=.true.

            call setup_box_size('CUB',rho,nm,boxlxyz,wmass)
            na = 3*nm
            allocate(r(3,na,nb))
            r(:,:,:) = 0.d0
            call setup_positions(r,na,nb,nm,boxlxyz,qo,irun)

            if (vacfac.gt.1) then
                boxlxyz(3)=vacfac*boxlxyz(3)
                r(3,:,:)=r(3,:,:)+0.5d0*boxlxyz(3)
                iamcub=.false.
                write(6,*)'* Vacuum preparation done'
            else if (vacfac.eq.1) then
                write(6,*) 'Vacuum not enabled (vacfac=1)'
                lattice='CUB'
            else
                write(6,*) 'Invalid vacfac choice (lt 1)'
                stop
            endif
            write(6,*)
        else

            ! Use lattice

            call setup_box_size(lattice,rho,nm,boxlxyz,wmass)
            na = 3*nm
            allocate(r(3,na,nb))
            r(:,:,:) = 0.d0

            call setup_positions(r,na,nb,nm,boxlxyz,qo,irun)
        endif
    endif
#ifdef PARALLEL_BINDING
endif


	call MPI_bcast(na,1,MPI_integer,0,MPI_COMM_WORLD,ierr) 
	call MPI_bcast(nb,1,MPI_integer,0,MPI_COMM_WORLD,ierr) 
	call MPI_bcast(nbdf1,1,MPI_integer,0,MPI_COMM_WORLD,ierr) 
	call MPI_bcast(nbdf2,1,MPI_integer,0,MPI_COMM_WORLD,ierr) 
	call MPI_bcast(nbdf3,1,MPI_integer,0,MPI_COMM_WORLD,ierr) 
	call MPI_bcast(rpmddft,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
	call MPI_bcast(rctdk,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
	call MPI_bcast(CP2K_path,34,MPI_character,0,MPI_COMM_WORLD,ierr) 
	call MPI_bcast(ng,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
	call MPI_bcast(pt,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
	call MPI_bcast(pb,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
	call MPI_bcast(nt,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
	call MPI_bcast(m,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
	call MPI_bcast(print,3,MPI_integer,0,MPI_COMM_WORLD,ierr)
	call MPI_bcast(irun,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
	call MPI_bcast(itst,3,MPI_integer,0,MPI_COMM_WORLD,ierr)
	call MPI_bcast(itcf,3,MPI_integer,0,MPI_COMM_WORLD,ierr)
	call MPI_bcast(ntherm,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
	call MPI_bcast(vacfac,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
	call MPI_bcast(iskip,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
	call MPI_bcast(beta,1,MPI_real8,0,MPI_COMM_WORLD,ierr)	
	call MPI_bcast(dt,1,MPI_real8,0,MPI_COMM_WORLD,ierr)
	call MPI_bcast(mts,1,MPI_integer,0,MPI_COMM_WORLD,ierr) 
	call MPI_bcast(ens,3,MPI_character,0,MPI_COMM_WORLD,ierr) 
	call MPI_bcast(type,4,MPI_character,0,MPI_COMM_WORLD,ierr) 
#endif

#ifdef PARALLEL_BINDING
if(myid.gt.0) then
	allocate(r(3,na,nb))
endif
	call MPI_bcast(r,SIZE(r),MPI_real8,0,MPI_COMM_WORLD,ierr)

if(myid.eq.0) then
#endif
		if (lattice.ne.'VAC' .and. vacfac.ne.1) then
        vacfac = 1
    end if
    ! Rigid body setup
    ! -----------------

    if (iamrigid) then
        call rattle_setup(nm)
    else
        nctot = 0
        nbond = 0
    endif

    ina = na
    inb = nb
    inm = nm
    ! Assign masses and charges :
    ! -----------------------------

    allocate (mass(na),z(na))
    z(:) = 0.d0
    mass(:) = 0.d0
   
   if(isotope.ne.'H2O' .and. dmass.eq. 0) then 
        write(6,*)
        write(6,*) 'Please specify mass of deuterium (dmass) in parameter-input!'
        write(6,*)
        stop
   
 
     else if(isotope.eq.'H2O')then                !Isotope selection H20, D2O or HDO resp. DHO.
        H2O = .true.
        do i = 1,na,3                             !D/H resp. H/D Mixture mode H2O/D2O not yet implemented.
          mass(i) = omass      ! Oxygen
          mass(i+1) = hmass    ! Hydrogen
          mass(i+2) = hmass    ! Hydrogen
        enddo
    
     else if(isotope.eq.'D2O') then
       D2O = .true.
       do i = 1,na,3
          mass(i) = omass      ! Oxygen
          mass(i+1) = dmass    ! Deuterium
          mass(i+2) = dmass    ! Deuterium
       enddo
    
     else if(isotope.eq.'HDO'.or.isotope.eq.'DHO') then
       HDO = .true.
       do i = 1,na,3
          mass(i) = omass      ! Oxygen
          mass(i+1) = dmass    ! Deuterium
          mass(i+2) = hmass    ! Hydrogen
       enddo
    
     !else if(isotope.eq.'H/D'.or.isotope.eq.'D/H') then        Doesn't work for uneven
       !do i = 1,na,6                                           number of molecules.
          !mass(i) = omass      ! Oxygen
          !mass(i+1) = dmass    ! Hydrogen
          !mass(i+2) = dmass    ! Hydrogen
       !enddo

      !do i = 2,na,6
          !mass(i) = omass      ! Oxygen
          !mass(i+1) = dmass    ! Deuterium
          !mass(i+2) = dmass    ! Deuterium
      !enddo

      else
        write(6,*)
        write(6,*) 'Please choose isotope = H2O, D2O or HDO!'
        write(6,*)
        stop
   endif   

    qh = -0.5d0*qo
    do i = 1,na,3
        z(i) = qo
        z(i+1) = qh
        z(i+2) = qh
    enddo
#ifdef PARALLEL_BINDING
endif

if(myid.gt.0) then
    allocate (mass(na),z(na))
endif
	call MPI_bcast(mass,na,MPI_real8,0,MPI_COMM_WORLD,ierr)	
	call MPI_bcast(z,na,MPI_real8,0,MPI_COMM_WORLD,ierr)
#endif

#ifdef PARALLEL_BINDING
if(myid.eq.0) then
#endif
    ! Write a vmd output of starting structure
    if (reftraj.eq.0) then
        open (unit=12,file='vmd_start.xyz')
        call print_vmd_full(r,nb,na,nm,boxlxyz,12)
        close (unit=12)
    end if

    ! structure is ready, we can initialize GLE thermostat
    if (therm.eq.'GLE') then                       !!GLE
        call therm_gle_init(ttau,na,nb,0.5*dt,irun,beta)      !!GLE
    end if

    ! Check Cut-off for LJ interactions
    ! ----------------------------

    boxmin = min(boxlxyz(1),boxlxyz(2),boxlxyz(3))
    rcut = min(rcut,0.5d0*boxmin)

    write (6,61) na,nm,boxlxyz(1),boxlxyz(2),boxlxyz(3),rcut
61  format( /1x, 'System setup : ' /1x, &
    '---------------' /1x, &
    'na     = ',i9,' atoms'/1x, &
    'nm     = ',i9,' molecules'/1x, &
    'boxlx  = ',f9.3,' bohr'/1x, &
    'boxly  = ',f9.3,' bohr'/1x, &
    'boxlz  = ',f9.3,' bohr'/1x, &
    'rcut   = ',f9.3,' bohr'/1x)

    ! Ewald parameters
    ! ------------------

    call setup_ewald(na,boxlxyz)

    if (reftraj.eq.0) then
        write(6,*)'Operations to be performed : '
        write(6,*)'-----------------------------'

        ! Static properties to be calculated
        ! --------------------------------------

        if (itst(1).eq.1) then
            write(6,*)'* RDF WILL be calculated... '
        else
            write(6,*)'* RDF WILL NOT be calculated... '
        endif
        if ((nbdf1.gt.0).or.(nbdf2.gt.0)) then
            if (itst(2).eq.1) then
                write(6,*)'* Exact Estimators WILL be calculated... '
            else
                write(6,*)'* Exact Estimators WILL NOT be calculated... '
            endif
        endif
        if (itst(3).eq.1) then
            write(6,*)'* MSD WILL be calculated... '
        else
            write(6,*)'* MSD WILL NOT be calculated... '
        endif

        ! Dynamics properties to be calculated
        ! --------------------------------------

        if (itcf(1).eq.1) then
            write(6,*)'* Cvv WILL be calculated... '
        else
            write(6,*)'* Cvv WILL NOT be calculated... '
        endif
        if (itcf(2).eq.1) then
            write(6,*)'* Dipole spectra WILL be calculated... '
        else
            write(6,*)'* Dipole spectra WILL NOT be calculated... '
        endif
        if (itcf(3).eq.1) then
            write(6,*)'* Orientational CFs WILL be calculated... '
        else
            write(6,*)'* Orientational CFs WILL NOT be calculated... '
        endif
          
        !E3B Correction to be calculated
        !.......................................

        if (rpmde3b.ne.0) then
           write(6,*) '* Three body correction WILL be calculated...' 
        else 
           write(6,*) '* Three body correction WILL NOT be calculated...' 
        endif
    

    !pt Check for dipole printing 
    if (ng .gt. 0) then
      if((print(4).eq.1).and.(mod(pt,10).ne.0)) then !dipolemoment is calculated every 10th step in md_static  
        if(pt.lt.10) then                            !therefore this Warning
          ptd = 10
        else 
          ptd = pt-mod(pt,10)
        endif
      write(6,*)
      write(6,*) 'WARNING:The molecular dipole moment is calculated every 10th step.'  
      write(6,*) 'Therefore the dipole moment will be calculated every', ptd, 'step.'
      write(6,*) 'Instead of every', pt, 'as the other parameters (position etc.).'
      else
       ptd = pt
      endif
    endif

    write(6,*)

        if (epsr.eqv..true.) then
            write (6,*) "EPSR WILL be used..."
#ifdef EPSR_STARTUP_READ
            write (6,*) "EPSR WILL be read on startup"
#else
            write (6,*) "EPSR WILL NOT be read on startup"
#endif
        else
            write (6,*) "EPSR WILL NOT be used..."
        endif

    !isotope Warninng
    if (D2O .or. HDO) then
      write(6,*)
      write(6,*) '     !!RPMD will use deuterium for the simulation!!'
      write(6,*) 'If you want to use some calculated properties like the radial distribution function.'
      write(6,*) 'Please check if the routines in this program are handling deuterium properly.'
      write(6,*)   
    endif


    endif
    ! Water Parameters used
    ! ----------------------

    open (unit=10,file='parameters.out')
    write(10,*) 'Parameters : '
    write(10,*) '-------------'
    write(10,*) ' Water mass              = ', wmass
    write(10,*) ' Oxygen mass             = ', omass
    write(10,*) ' Hydrogen mass           = ', hmass
    write(10,*) ' Deuterium mass          = ', dmass
    write(10,*) ' Charge on Oxygen        = ', qo
    write(10,*) ' Alpha (M-site)          = ', alpha
    write(10,*) ' LJ Sigma                = ', oo_sig
    write(10,*) ' LJ Epsilon              = ', oo_eps
    write(10,*) ' Buckingham gamma        = ', oo_gam
    write(10,*) ' Equilibrium Bond Angle  = ', thetad
    write(10,*) ' Equilibrium Bond Length = ', reoh
    write(10,*) ' Bond Stretch De         = ', apot
    write(10,*) ' Bond Bend De            = ', bpot
    write(10,*) ' Bond Stretch Alpha      = ', alp
    write(10,*) ' Bond Bend Alpha         = ', alpb
    write(10,*) ' Gaussian Width M-site   = ', wm
    write(10,*) ' Gaussian Width Hydrogen = ', wh
    close (unit=10)
#ifdef PARALLEL_BINDING
endif
#endif

    ! Allocate the evolution arrays
    ! ------------------------------

    allocate (p(3,na,nb),dvdr(3,na,nb),dvdr2(3,na,nb))
    p(:,:,:) = 0.d0
    dvdr(:,:,:) = 0.d0
    dvdr2(:,:,:) = 0.d0


    ! Allocate force enviroment id:
    ! -----------------------
    if (nbdf3.eq.0) then
        allocate(f_env_id(nb))
    else
        allocate(f_env_id(nbdf3))
    endif

#ifdef PARALLEL_BINDING
if(myid.eq.0) then
#endif
    ! set rpmddft = 0 because classical equilibration, set it back to old value later
    rpmddfthelp = rpmddft
    rpmddft = 0


    if (reftraj.eq.0) then
        ! Initial forces and momenta
        ! ---------------------------

        call full_forces(r,na,nb,v,vew,vlj,vint,vir,z,boxlxyz,dvdr,dvdr2)
        write(6,'(a,f10.5,a)') ' * Initial energy =', v/dble(na), ' E_h per atom'
        write(6,*)
        call sample(p,na,nb,mass,beta,irun,dt)

		if (iamrigid) then
				write(*,*) "Rigid Simulation"
				write(*,*) ""
		endif

        ! Equilibration or Interface Melting
        ! -----------------------------------

        if (lattice.ne.'INT') then
            call md_eq(ne,p,r,dvdr,dvdr2,na,nb,boxlxyz,z,beta, &
            dt,mass,irun)
        else
            call md_melt(ne,p,r,dvdr,dvdr2,na,nb,boxlxyz,z,beta, &
            dt,mass,irun,nm_ice)
        endif

        ! Write Equilibration file to use for restarts
        ! ----------------------------------------------

        open(10, file = 'eq.xyz')
        call print_structure(r,boxlxyz,nm,na,nb,10)
        close (unit=10)
        open (unit=12,file='vmd_eq.xyz')
        call print_vmd_full(r,nb,na,nm,boxlxyz,12)
        close (unit=12)
    endif


    ! set rpmddft back to old value
    rpmddft = rpmddfthelp
#ifdef PARALLEL_BINDING
endif
		
		call MPI_bcast(boxlxyz,3,MPI_real8,0,MPI_COMM_WORLD,ierr)
#endif
#ifdef CP2K_BINDING
      !  RPMD-DFT
      ! -----------------------
			if (rpmddft.eq.1) then
			  if (nbdf3.eq.0) then
            do i = 1, nb 
                call cp_create_fenv(f_env_id(i),CP2K_path,"out.out",ierr)
                ! set cell in CP2K
                cell(1,1) = boxlxyz(1)
                cell(2,2) = boxlxyz(2)
                cell(3,3) = boxlxyz(3)
                call cp_set_cell(f_env_id(i),cell,ierr)
                if (ierr.ne.0) STOP "set_cell"
								!write(*,*) "f_env_id:",f_env_id(i)
            enddo
        else
            do i = 1, nbdf3 ! if nbdf3 >= 1 set up nbdf3 versions of cp2k f_env_id = i
                call cp_create_fenv(f_env_id(i),CP2K_path,"out.out",ierr)
                !write(*,*) "f_env_id:",f_env_id(i)
                !set cell in CP2K
                cell(1,1) = boxlxyz(1)
                cell(2,2) = boxlxyz(2)
                cell(3,3) = boxlxyz(3)
                call cp_set_cell(f_env_id(i),cell,ierr)
                if (ierr.ne.0) STOP "set_cell"
            enddo
        endif
			endif
#endif




      ! Static Properties
      ! ------------------

    if (reftraj.eq.0) then
#ifdef PARALLEL_BINDING
						if(myid.eq.0) then
#endif
        if (ng .gt. 0) then
            write(6,*) '* Sampling Static Properties '


	            call sample(p,na,nb,mass,beta,irun,dt)

#ifdef PARALLEL_BINDING
						endif
						call MPI_bcast(r,SIZE(r),MPI_real,0,MPI_COMM_WORLD,ierr)
						call MPI_bcast(p,SIZE(p),MPI_real,0,MPI_COMM_WORLD,ierr)
						call MPI_bcast(dvdr,SIZE(p),MPI_real,0,MPI_COMM_WORLD,ierr)
						call MPI_bcast(dvdr2,SIZE(p),MPI_real,0,MPI_COMM_WORLD,ierr)
#endif

						call md_static(ng,p,r,dvdr,dvdr2,na,nb,boxlxyz,z,beta,&
            dt,mass,irun,itst,pt,pb,print,vave,ptd)
        endif
    else
        if (reftraj.lt.0) then
#ifdef ENABLE_IPI_DRIVER
            ! act as i-pi driver
            call run_driver(p,r,dvdr,dvdr2,na,nb,boxlxyz,z,beta,&
            dt,mass,irun,itst,pt,pb,print)
#else
            write(6,*) "Error: ipi-driver not compiled. Compiling with flag -DENABLE_IPI_DRIVER is required."
            stop
#endif
        else
				do i = 1, reftraj
            boxlxyz(:) = boxlxyz_traj(:,i)
            r(:,:,:) = r_traj(:,:,:,i)
            call md_static(1,p,r,dvdr,dvdr2,na,nb,boxlxyz,z,beta,&
            dt,mass,irun,itst,pt,pb,print,vave,ptd)
            write (6,*) "potential of snapshot", i, "is", vave*toKjmol, "kJ/mol"
        enddo
        endif
    endif


    if (reftraj.eq.0) then

        ! Dynamical Properties
        ! ----------------------

        if (nt.gt.0.and.m.gt.0) then
#ifdef PARALLEL_BINDING
						if(myid.eq.0) then
#endif
            write(6,*) '* Sampling Dynamical Properties'


	            call sample(p,na,nb,mass,beta,irun,dt)
#ifdef PARALLEL_BINDING
						endif
						call MPI_bcast(r,SIZE(r),MPI_real,0,MPI_COMM_WORLD,ierr)
						call MPI_bcast(p,SIZE(p),MPI_real,0,MPI_COMM_WORLD,ierr)
						call MPI_bcast(dvdr,SIZE(p),MPI_real,0,MPI_COMM_WORLD,ierr)
						call MPI_bcast(dvdr2,SIZE(p),MPI_real,0,MPI_COMM_WORLD,ierr)




#endif
						
            call dynamics(nt,m,p,r,dvdr,dvdr2,na,nm,nb,boxlxyz,z,beta, &
            dt,mass,irun,itcf,itst,pt,pb,print,iskip,ntherm,vacfac)
        endif
    endif
    ! Write a vmd output of ending structure
    if (reftraj.eq.0) then
        open (unit=12,file='vmd_end.xyz')
        call print_vmd_full(r,nb,na,nm,boxlxyz,12)
        close (unit=12)
        open(10, file = 'end.xyz')
        call print_structure(r,boxlxyz,nm,na,nb,10)
        close (unit=10)
    end if

#ifdef PARALLEL_BINDING
if(myid.eq.0) then
#endif
    close (unit=61)
    deallocate(r,mass,z,p,dvdr,dvdr2)
#ifdef PARALLEL_BINDING
endif
if(myid.gt.0) then
    deallocate(r,p,dvdr,dvdr2)
endif
#endif
    if (reftraj.gt.0) then
        deallocate(r_traj)
    endif

#ifdef CP2K_BINDING
#ifdef PARALLEL_BINDING
		call cp_finalize_cp2k(1,ierr)
		!call cp_finalize_cp2k(0,ierr)		
		!call MPI_finalize(ierr) 
#endif
#ifndef PARALLEL_BINDING
		call cp_finalize_cp2k(0,ierr)
#endif
#endif

#ifdef PARALLEL_BINDING
#ifndef CP2K_BINDING
		call MPI_finalize(ierr) 		
#endif
#endif



    return
end program