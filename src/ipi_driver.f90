#ifdef ENABLE_IPI_DRIVER

!
! DOES NOT COMPILE with GCC10 (even with std=legacy) as the code relies
! on non-standard behavior. It is not clear how to easily fix the problem.
!
! *****************************************************************************
!> \brief Driver mode - To communicate with i-PI Python wrapper
!> \par History
!>      none
!> \author Michele Ceriotti 03.2012
! *****************************************************************************
SUBROUTINE run_driver (p,r,dvdr,dvdr2,na,nb,boxlxyz,z,beta,&
        dt,mass,irun,itst,pt,pb,print)
    implicit none
    integer nt,ne,nb,m,ng,npre_eq,pt,pb,irun,nm,na,narg,iargc,nbdf1,nbdf2,nbdf3,rctdk,i,j,k
    integer nbr,mts,nlat,itcf(3),itst(3),ncellxyz(3),print(4),ntherm,iskip!,sizeMPI(9)

    ! used for reftrj and RPMD-DFT
    integer reftraj,rpmddft,ierr,rpmddfthelp,myid,numid
    character(len=35) CP2K_path
    ! r_traj(xyz,molecules,nb,reftraj)
    real(8), allocatable :: r_traj(:,:,:,:),boxlxyz_traj(:,:)
    character(len=240) line
    common /reftraj/ reftraj,line

    integer nc_ice(3),nc_wat(3),nm_ice,nm_wat,nctot,nbond,vacfac
    real(8) temp,rho,dtfs,ecut,test,beta,dt,dtps,boxmin,pres
    real(8) teqm,tsim,trdf,gaussian,wmass,rcut,om,ttaufs
    real(8) taufs,v,vew,voo,vlj,vint,pi,vir(3,3),send_vir(3,3),cell(3,3)
    real(8) qo,qh,alpha,oo_sig,oo_eps,oo_gam,theta,reoh,thetad
    real(8) apot,bpot,alp,alpb,wm,wh,omass,hmass,sig,boxlxyz(3),iboxlxyz(3,3),vdum
    real(8) box_ice(3),box_wat(3),rcut_old
    !real(8), allocatable :: mass(:),z(:),r(:,:,:)
    !real(8), allocatable :: p(:,:,:),dvdr(:,:,:),dvdr2(:,:,:)
    real(8) :: mass(na),z(na),r(3,na,nb)
    real(8) :: p(3,na,nb),dvdr(3,na,nb),dvdr2(3,na,nb)
    character(len=25) filename
    character(len=4) type
    character(len=3) lattice,ens,therm_backup
    logical iamcub,iamrigid
    external gaussian

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
    !TYPE(force_env_type), POINTER            :: force_env
    !TYPE(global_environment_type), POINTER   :: globenv
    !TYPE(cp_error_type), INTENT(inout)       :: error

    !TYPE(section_vals_type), POINTER         :: drv_section, motion_section
    !TYPE(virial_type), POINTER               :: virial

    INTEGER, PARAMETER :: MSGLEN=12
    !CHARACTER(len=MSGLEN)       :: drv_hostname, c_hostname
    !INTEGER                                  :: drv_port
    !LOGICAL                                  :: drv_unix

    !MIK DRIVER
    ! server address parsing
    CHARACTER(1024) :: serveraddr
    CHARACTER(MSGLEN) :: host
    INTEGER socket, port, inet, nread, readbuffer, slock, swait, uwait
    ! buffers and temporaries for communication
    LOGICAL :: isinit=.true., hasdata=.false., ionode=.false.
    CHARACTER(12) :: header
    CHARACTER(1024) :: parbuffer
    INTEGER nat
    !REAL(8) :: cellh(3,3), cellih(3,3), pot
    REAL(8), ALLOCATABLE :: combuf(:)
    REAL(8) :: sigma(3,3)
    ! access cp2k structures
    !TYPE(cp_subsys_type), POINTER            :: subsys
    integer ip, ii, idir
    !TYPE(cell_type), POINTER            :: cpcell
    logical should_stop
    should_stop=.false.
    ionode=.true.!(default_para_env%source==default_para_env%mepos)

    ! reads driver parameters from input
    !motion_section => section_vals_get_subs_vals(force_env%root_section,"MOTION",error=error)
    !drv_section     => section_vals_get_subs_vals(motion_section,"DRIVER",error=error)

    !CALL section_vals_val_get(drv_section,"HOST",c_val=drv_hostname,error=error)
    !CALL section_vals_val_get(drv_section,"PORT",i_val=drv_port,error=error)
    !CALL section_vals_val_get(drv_section,"UNIX",l_val=drv_unix,error=error)


    if (nb.ne.1) then
        write(*,*) "i-PI mode must have nb = 1!"
        stop
    endif

    ! opens the socket
    socket=0
    inet=1
    if (ionode) then
        write(*,*) "@ i-PI DRIVER BEING LOADED"
        open(12,file="serverdata")
        read(12,*) serveraddr
        close(12)
        write(*,*) "SERVERDATA ", serveraddr

        host=serveraddr(1:INDEX(serveraddr,':',back=.true.)-1)//achar(0)
        read(serveraddr(INDEX(serveraddr,':',back=.true.)+1:),*) port
        IF (trim(serveraddr(1:INDEX(serveraddr,':')-1)).eq.trim("unix")) inet=0
        host=serveraddr(INDEX(serveraddr,':')+1:INDEX(serveraddr,':',back=.true.)-1)//achar(0)
        write(*,*) "INET: ",inet, " HOST: ", trim(host), " PORT: ", port
        CALL open_socket(socket, inet, port, trim(host))

        !write(*,*) "@ INPUT DATA: ", TRIM(drv_hostname), drv_port, drv_unix
        !c_hostname=TRIM(drv_hostname)//achar(0)
        !CALL open_socket(socket, .not. drv_unix, drv_port, c_hostname)
    endif

    !now we have a socket, so we can initialize the CP2K environments.
    !NULLIFY(cpcell)
    !call cell_create(cpcell,error=error)
    ! CELL ALREADY CREATED BASED ON OWN INPUT FILE
    uwait=10000  ! number of MICROseconds to be waited in filesystem lock
    driver_loop: DO
        ! do communication on master node only...
        header = ""

        ! this syncs the processes, possibly (see sockets.c) without calling MPI_Barrier,
        ! which is nice as MPI_barrier eats up a lot of CPU for nothing
        !inet=slock(default_para_env%source, default_para_env%mepos)
        !CALL mp_sync(default_para_env%group)

        if (ionode) nread=readbuffer(socket, header, MSGLEN)
        if (ionode)  write(0,*) "returned from readbuffer"

        !inet=swait(uwait, default_para_env%source, default_para_env%mepos)
        !CALL mp_sync(default_para_env%group)

        !call mp_bcast(nread,default_para_env%source, default_para_env%group)
        if (nread .eq. 0) then
            if (ionode) write(*,*) " @ DRIVER MODE: Could not read from socket, exiting now.", nread
            exit
        endif

        !call mp_bcast(header,default_para_env%source, default_para_env%group)

        if (ionode) write(*,*) " @ DRIVER MODE: Message from server: ", trim(header)
        if (trim(header) == "STATUS") then

            !inet=slock(default_para_env%source, default_para_env%mepos)
            !CALL mp_sync(default_para_env%group)
            if (ionode) then  ! does not  need init (well, maybe it should, just to check atom numbers and the like... )
                if (hasdata) then
                    call writebuffer(socket,"HAVEDATA    ",MSGLEN)
                else
                    call writebuffer(socket,"READY       ",MSGLEN)
                endif
            endif
            !inet=swait(uwait,default_para_env%source, default_para_env%mepos)
            !CALL mp_sync(default_para_env%group)
        else if (trim(header) == "POSDATA") then
            if (ionode) then
                nread=readbuffer(socket, iboxlxyz, 9*8)
                iboxlxyz=transpose(iboxlxyz)
                do ii=1,3
                    boxlxyz(ii) = iboxlxyz(ii,ii)
                enddo
                nread=readbuffer(socket, iboxlxyz, 9*8)
                nread=readbuffer(socket, nat, 4)
            endif
            !call mp_bcast(cellh,default_para_env%source, default_para_env%group)
            !call mp_bcast(cellih,default_para_env%source, default_para_env%group)
            !call mp_bcast(nat,default_para_env%source, default_para_env%group)
            if (.not.allocated(combuf)) allocate(combuf(3*nat))
            if (ionode) nread=readbuffer(socket, combuf, nat*3*8)
            !call mp_bcast(combuf,default_para_env%source, default_para_env%group)

            !CALL force_env_get(force_env,subsys=subsys,error=error)
            if (nat/=na) WRITE(*,*) " @DRIVER MODE: Uh-oh! Particle number mismatch between i-pi and rpmd input!"
            ii=0
            DO ip=1,na
                DO idir=1,3
                    ii=ii+1
                    r(idir,ip,1)=combuf(ii)
                END DO
                !write(6,*) "GOT ", r(:,ip,1)
            END DO
            ! TODO rpmd equivalent
            !CALL init_cell(cpcell, hmat=cellh)
            !CALL force_env_set_cell(force_env,cell=cpcell,error=error)

            !CALL force_env_calc_energy_force(force_env,calc_force=.TRUE. ,error=error)
            !call md_static(1,p,r,dvdr,dvdr2,na,nb,boxlxyz,z,beta,&
            !                dt,mass,irun,itst,pt,pb,print)
            call full_forces(r,na,nb,v,vew,voo,vint,vir,z,boxlxyz, &
                    dvdr,dvdr2)

            if (ionode) write(*,*) " @ DRIVER MODE: Received positions "

            combuf=0
            ii=0
            DO ip=1,na
                DO idir=1,3
                    ii=ii+1
                    combuf(ii)=-dvdr(idir,ip,1)-dvdr2(idir,ip,1)
                END DO
                !write(6,*) "SEND ", -dvdr(:,ip,1)-dvdr2(:,ip,1)
                !write(6,*) "SEND2 ", combuf(ii-2), combuf(ii-1), combuf(ii)
            END DO
            !CALL force_env_get(force_env, potential_energy=pot, error=error)
            !CALL force_env_get(force_env,cell=cpcell, virial=virial, error=error)
            !vir = transpose(virial%pv_virial)
            send_vir = transpose(vir)
            send_vir = -1.0*send_vir

            !CALL external_control(should_stop,"DPI",globenv=globenv,error=error)
            IF (should_stop) EXIT

            hasdata=.true.
        else if (trim(header)=="GETFORCE") then
            if (ionode) write(*,*) " @ DRIVER MODE: Returning v,forces,stress "
            if (ionode) then
                call writebuffer(socket,"FORCEREADY  ",MSGLEN)
                call writebuffer(socket,v,8)
                call writebuffer(socket,nat,4)
                call writebuffer(socket,combuf,3*nat*8)
                call writebuffer(socket,send_vir,9*8)

                ! i-pi can also receive an arbitrary string, that will be printed out to the "extra"
                ! trajectory file. this is useful if you want to return additional information, e.g.
                ! atomic charges, wannier centres, etc. one must return the number of characters, then
                ! the string. here we just send back zero characters.
                nat=0
                call writebuffer(socket,nat,4)  ! writes out zero for the length of the "extra" field (not implemented yet!)
            endif
            hasdata=.false.
        else
            if (ionode) write(*,*) " @DRIVER MODE:  Socket disconnected, time to exit. "
            exit
        endif
    ENDDO driver_loop

END SUBROUTINE run_driver

#endif