program wff_forces

    use f_env
    implicit none
    include 'globals.inc'

    !--------------------------------------------------------------------
    !Classical/RPMD/PACMD Simulation Program for Flexible water
    !--------------------------------------------------------------------
    integer nt,ne,nb,m,ng,npre_eq,pt,pb,irun,nm,na,narg,iargc,nbdf1,nbdf2,nbdf3,rctdk,i,j,k
    integer mts,nlat,itcf(3),itst(3),ncellxyz(3),print(3),ntherm,iskip!,sizeMPI(9)

    ! used for reftrj and RPMD-DFT
    integer reftraj,rpmddft
    character(len=35) CP2K_path
    ! r_traj(xyz,molecules,nb,reftraj)
    real(8), allocatable :: r_traj(:,:,:,:),boxlxyz_traj(:,:)
    character(len=240) line
    common /reftraj/ reftraj,line

    integer nctot,nbond,vacfac
    real(8) temp,rho,dtfs,ecut,test,beta,dt,dtps,boxmin,pres
    real(8) teqm,tsim,trdf,gaussian,wmass,rcut,om,ttaufs
    real(8) taufs,v,pi,vir(3,3),cell(3,3)
    real(8) qo,qh,alpha,oo_sig,oo_eps,oo_gam,theta,reoh,thetad
    real(8) apot,bpot,alp,alpb,wm,wh,omass,hmass,sig,boxlxyz(3)
    real(8) virial
    real(8), allocatable :: mass(:),z(:),r(:,:,:)
    real(8), allocatable :: dvdr(:,:,:),dvdr2(:,:,:),dvdr_split(:,:,:,:)
    character(len=25) filename
    character(len=4) type
    character(len=3) lattice,ens
    logical iamcub,iamrigid
    external gaussian

    namelist/input/iamcub,dtfs,ecut,nb,reftraj,ncellxyz,irun,rcut,nbdf1,nbdf2,sig
    namelist/param/ wmass,omass,hmass,qo,alpha,oo_sig,oo_eps,oo_gam, &
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

    vacfac = 1
    ntherm = 0
    ttaufs = 0.d0
    iskip = 1
    rpmddft = 0
    npre_eq = 0
    cell(:,:) = 0.d0
    CP2K_path =""
    nbdf3 = 0
    rctdk = 0
    itst(3) = 0 ! itst 3 argument is for RMS calculation

    !FFM defaults

    m=1
    mts=1
    ne=0
    ng=0
    irun=1
    temp=300.0d0
    pres=1.0d0
    dtfs=1.0d0
    ens='NVE'
    lattice='CUB'
    type='RPMD'
    print = (/ 0 , 1 , 0 /)
    pt=1
    pb=1
    itcf(:)=0
    itst(:)=0
    om=0.d0
    reftraj=1

    write(6,*)  '-------------------------------------------'
    write(6,*)  '            Flexible Water Code            '
    write(6,*)  '-------------------------------------------'
    write(6,*)

    ! Read the input file, given by the first argument to the program

    narg = iargc()
    if ( narg .ne. 3 ) then
        write(6,*)'* Program usage:'
        write(6,*)' wff.x <input file> <param file> <eq file> '
        stop
    endif

    ! Input file

    call getarg (1,filename)
    open(60,file=filename)
    read(60,input)
    close (unit=60)

    if (reftraj.le.0) then
        write(6,*) 'reftraj must not be zero! Stop.'
        stop
    end if
    nt=reftraj

    ! Parameters File

    call getarg (2,filename)
    open(61,file=filename)
    read(61,param)
    close (unit=61)

    ! Rigid Simulation?
    ! ------------------

    if ((apot.eq.0.d0).and.(bpot.eq.0.d0)) then
        iamrigid = .true.

        write(6,*) 'Warning: Detected rigid water model -> dvdr_split(:,:,:,3:4)=0'

    else
        iamrigid = .false.
    endif

    ! Temperature, bond angles and times
    ! ------------------------------------

    pi = dacos(-1.d0)
    theta = thetad*(pi/180.d0)
    beta = 1.d0/(3.166829d-6*temp)   ! (1/k_b (Hartree) *temp)
    dt = tofs*dtfs
    dtps = 1.d-3*dtfs
    teqm = ne*dtps
    tsim = nt*dtps
    trdf = ng*dtps
    rcut = rcut / ToA

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

    ! Read from file

    call getarg(3,filename)
    open (61, file = filename)

    ! Thrid argument is trajectories file
    write(6,*)
    write(6,*) 'The current reftraj ID is: ', reftraj
    write(6,*) 'And is loaded from file: ', filename
    write(6,*)

    !call setup_box_size(lattice,rho,nm,boxlxyz,wmass)
    nm = ncellxyz(1)*ncellxyz(2)*ncellxyz(3)
    na = 3*nm
    write(6,*) ' Note: Based on the ncellxyz setting we asume',na,'atoms.'// &
               ' Please check for correctness!'

    allocate(r(3,na,nb))
    r(:,:,:) = 0.d0

    allocate(r_traj(3,na,nb,reftraj))
    r_traj(:,:,:,:) = 0.d0
    allocate(boxlxyz_traj(3,reftraj))

    do k = 1, nb
        do i = 1, reftraj
            read(61,*) line

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
    rho=18.d0*nm*1.660538921d0/dsqrt(dot_product(boxlxyz,boxlxyz))

    ! Rigid body setup
    ! -----------------

    if (iamrigid) then
        write(6,*) 'Rigid not implemented... Stop.'
        stop
    else
        nctot = 0
        nbond = 0
    endif
    ! Assign masses and charges :
    ! -----------------------------

    allocate (mass(na),z(na))
    z(:) = 0.d0
    mass(:) = 0.d0

    do i = 1,na,3
        mass(i) = omass      ! Oxygen
        mass(i+1) = hmass    ! Hydrogen
        mass(i+2) = hmass    ! Hydrogen
    enddo

    qh = -0.5d0*qo
    do i = 1,na,3
        z(i) = qo
        z(i+1) = qh
        z(i+2) = qh
    enddo

    write(6,60) 'FFM',dtfs,tsim, &
        rho,ncellxyz(1),ncellxyz(2),ncellxyz(3)
60  format (' Simulation Parameters:' /1x &
        '-----------------------' /1x &
        'ens    = ',a9 /1x &
        'dtfs   = ',f9.3,' fs'/1x, &
        'tsim   = ',f9.3,' ps'/1x, &
        'rho    = ',f9.3,' g/cm**3'/1x,  &
        'ncellx = ',i9,' unit cells'/1x, &
        'ncelly = ',i9,' unit cells'/1x, &
        'ncellz = ',i9,' unit cells')

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

    ! Water Parameters used
    ! ----------------------

    open (unit=10,file='parameters.out')
    write(10,*) 'Parameters : '
    write(10,*) '-------------'
    write(10,*) ' Water mass              = ', wmass
    write(10,*) ' Oxygen mass             = ', omass
    write(10,*) ' Hydrogen mass           = ', hmass
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


    ! Allocate the evolution arrays
    ! ------------------------------

    allocate (dvdr(3,na,nb),dvdr2(3,na,nb),dvdr_split(3,na,nb,4))

    dvdr(:,:,:) = 0.d0
    dvdr2(:,:,:) = 0.d0
    dvdr_split(:,:,:,:) = 0.d0


    ! Allocate force enviroment id:
    ! -----------------------
    if (nbdf3.eq.0) then
        allocate(f_env_id(nb))
    else
        allocate(f_env_id(nbdf3))
    endif

    ! Forces
    ! ------------------
    if (nb.ne.1) then
        write(6,*) 'forces are called directly, switch to full_forces for nb>1! Stop.'
        stop
    end if

    dvdr2(:,:,:)=0.d0
    open(1250,file='vmd.frc')
    open(1251,file='vmd.frc1')
    open(1252,file='vmd.frc2')
    open(1253,file='vmd.frc3')
    open(1254,file='vmd.frc4')
    open(1255,file='vmd.frcT')
    do i = 1, reftraj
        boxlxyz(:) = boxlxyz_traj(:,i)
        r(:,:,:) = r_traj(:,:,:,i)
        call forces(r,v,dvdr,dvdr_split,nb,na,boxlxyz,z,virial,0)
        call print_vmd_full_forces(dvdr,dvdr2,nb,na,boxlxyz,1250)
        do j = 1, 4
            k=1250+j
            call print_vmd_full_forces(dvdr_split(:,:,:,j),dvdr2,nb,na,boxlxyz,k)
        enddo
        call print_vmd_full_forces(sum(dvdr_split(:,:,:,:),4),dvdr2,nb,na,boxlxyz,1255)
    enddo

    close (unit=1250)
    close (unit=1251)
    close (unit=1252)
    close (unit=1253)
    close (unit=1254)
    close (unit=1255)
    close (unit=61)

    deallocate(r,mass,z,dvdr,dvdr2)

    deallocate(r_traj)

    return
end program wff_forces


