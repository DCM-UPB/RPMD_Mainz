module static_files
use FileMod, only : FileHandle
implicit none
public
character(LEN=4),Dimension(3),parameter :: suffix=(/'.xyz','.frc','.vel'/)
type(FileHandle),Dimension(3) :: vmd_traj_set
type(FileHandle),Dimension(:,:),allocatable :: vmd_bead_set
type(FileHandle) :: fh_mol_Dipole,fh_vt,fh_dielec_conv,fh_NPT_prop
type(FileHandle) :: fh_density,fh_pressure,fh_MSD,fh_dipole_sys
logical print_traj,print_beads
logical,Dimension(4):: print_set_element
end module
!----------------------------------------------------------------------
!<summary>
!   open all files neccessary and requested for a static calculation
!</summary>
!<param name="nb">the number of beads</param>
!<param name="print_traj">activates the output of trajectory files</param>
!<param name="print_beads">activates the output of bead files</param>
!<param name="print">
!<param name="reftraj">information about restarts</param>
!   an array that controls the printing of the varios output species
!</param>
!<param name=printNPT>should files with NPT output be opened</param>
!<param name=printMSD>should the MSD be printed </param>
!----------------------------------------------------------------------
subroutine md_static_prepare_traj(nb,pt,pb,print,reftraj,printNPT,printMSD)
    use thermostat
    use barostat
    use useful, only : toChar
    use FileMod, only : FileHandle
    use static_files, only : vmd_traj_set,vmd_bead_set,print_traj,print_beads
    use static_files, only : print_set_element,suffix,fh_mol_dipole
    use static_files, only : fh_vt,fh_dielec_conv,fh_NPT_prop
    use static_files, only : fh_density,fh_pressure,fh_MSD,fh_dipole_sys

    implicit none
    include 'globals.inc'
    character(*),parameter :: routineN="md_static_prepare_traj"
    character(*),parameter :: routineP="__NONE__"//":"//routineN
    integer,intent(in)  :: nb,pt,pb,reftraj
    integer,intent(in)  :: print(4)
    logical worked
    integer ib
    character(len=128) file_name
    character(LEN=25) file_position,file_status
    integer iter,j
    logical printNPT,printMSD

    print_set_element(:) = print(:).eq.1
    print_traj=pt.gt.0
    print_beads=pb.gt.0
    ! init all file handlers
    !suffix=(/'.xyz','.frc','.vel'/)
    allocate(vmd_bead_set(nb,3))
    do iter=lbound(vmd_bead_set,2),ubound(vmd_bead_set,2),1
        vmd_traj_set(iter)=FileHandle('vmd_traj'//suffix(iter)) 
        do j=lbound(vmd_bead_set,1),ubound(vmd_bead_set,1),1
            vmd_bead_set(j,iter)=FileHandle('vmd_bead-'//toChar(j)//&
            suffix(iter))
        end do
    end do
    fh_mol_dipole=FileHandle('mol-dipole_traj.txt')

    fh_vt=FileHandle('vt_st.out') !probably for E_pot and virial E_kin
    fh_dipole_sys=FileHandle('dipole_sys_st.out')
    fh_dielec_conv=FileHandle('dielec_con_st.out')
    fh_NPT_prop=FileHandle('NPT_prop_st.out')
    fh_density=FileHandle('density_st.out')
    fh_pressure=FileHandle('pressure_st.out')
    fh_MSD=FileHandle('MSD_st.out')

    !---------------------------------------------
    !set appropriate options for restarted run or 
    !new run
    if (reftraj.ne.0) then
        file_position="append"
        file_status="old"
    else 
        file_position="asis"
        file_status="replace"
    end if
    worked = fh_vt%open(STATUS=file_status,ACTION='write',&
    POSITION=file_position)
    worked = fh_dipole_sys%open(STATUS=file_status,ACTION='write',&
    POSITION=file_position)
    worked = fh_dielec_conv%open(STATUS=file_status,ACTION='write',&
    POSITION=file_position)
    worked = fh_pressure%open(STATUS=file_status,ACTION='write',&
    POSITION=file_position)
    if(printNPT)then
        worked = fh_NPT_prop%open(STATUS=file_status,ACTION='write',&
        POSITION=file_position)
        worked = fh_density%open(STATUS=file_status,ACTION='write',&
        POSITION=file_position)
    end if
    if(printMSD)then
        worked = fh_MSD%open(STATUS=file_status,ACTION='write',&
        POSITION=file_position)
    end if
    if (print_traj) then
        do iter=lbound(vmd_traj_set,1),ubound(vmd_traj_set,1),1
            if(print_set_element(iter))then
                worked=vmd_traj_set(iter)%open(STATUS=file_status,ACTION='write',&
                POSITION=file_position)
            end if
        end do
        if (print_set_element(4))then
            worked = fh_mol_dipole%open(STATUS=file_status,ACTION='write',&
            POSITION='file_position')
        endif
    endif
    !-----------------------------------------------
    !now make sure the vmd_bead files are replaced
    if (print_beads) then
        !bead set ioUnits=127+ , 100127+, 200127+
        do iter=lbound(vmd_bead_set,2),ubound(vmd_bead_set,2),1
            if(print_set_element(iter))then
                do ib=1,nb
                    worked=vmd_bead_set(ib,iter)%open(STATUS=file_status,&
                    ACTION='write',POSITION=file_position)
                end do
            end if
        end do
    endif
end subroutine md_static_prepare_traj
!------------------------------------------------------------------
!<summary>
!   Routine to calculate the static properties
!</summary>
!<param name="p">the array of all momenta of all particles</param>
!<param name="r">the array of all coordinates of all particles</param>
!<param name="pt">controls the printing of the trajectory output</param>
!<param name="pb">controls the printing of the bead output</param>
!------------------------------------------------------------------
subroutine md_static(ng,p,r,dvdr,dvdr2,na,nb,boxlxyz,z,beta, &
                     dt,mass,irun,itst,pt,pb,print,vave,ptd)
    use thermostat
    use barostat
    use useful, only : get_Unit
    use static_files, only : vmd_traj_set,vmd_bead_set,print_traj,print_beads
    use static_files, only : print_set_element,suffix,fh_mol_dipole
    use static_files, only : fh_vt,fh_dielec_conv,fh_NPT_prop
    use static_files, only : fh_density,fh_pressure,fh_MSD,fh_dipole_sys
    use RDFMod, only : calc_rdf
    implicit none
    include 'globals.inc'
#ifdef PARALLEL_BINDING
    include 'mpif.h' !parallel
    integer myid
#endif
#if RDF_IN_THIS_ROUTINE == 1
    integer jj,ibin
    real(8),rij,dx,dy,dz
#endif
    character(*),parameter :: routineN='md_static'
    character(*),parameter :: routineP="__NONE__"//":"//routineN
    integer na,nb,ng,irun,nrdf,nbdf1,nbdf2,nbaro,pt,pb,ierr,ptd
    integer no,i,je,k,j,ii,nm,itst(3),ib,print(4)
    real(8), allocatable :: ihhh(:),ihoo(:),ihoh(:)
    real(8) boxlxyz(3),beta,v,dt
    real(8) p(3,na,nb), dvdr(3,na,nb),dvdr2(3,na,nb),ran2,z(na)
    real(8) r(3,na,nb),mass(na),vir(3,3),vir_lf(3,3),rcm(3,na/3),rst(3,na/3),r0(3,na,nb)
    real(8) delr,thresh,dconv,dconv2,pi,fac,vol
    real(8) avang,avoh,tavang,tavoh,msd,msdO
    real(8) ttaufs
    real(8) vave,tave,tq1ave,tq2ave,tv,tvxyz(3),tavee,tq1aee,tq2aee
    real(8) tavdip,tavdipx,tavdipy,tavdipz,tavdip2,tavdipsq,eps
    real(8) dipx,dipy,dipz,dip2,dipm,dtps,dr
    real(8) boxmax
    real(8) rgh,rgo,rgcm,rghav,rgoav,rgcmav,tq1,tq2,tqe
    real(8) v1,v2,v3,v1ave,v2ave,v1eeav,v2eeav,v3ave,v3eeav,veeav,denc
    real(8) pres,volav,pav,pvav,xav,yav,zav,denav,wmass,v_corr
    real(8) dvdr_corr(3,na,nb),vir_corr(3,3)
    real(8), allocatable :: dist(:,:), dvdre(:,:,:)
    real(8), allocatable :: m_dipole(:,:)
    character(len=3) ens
    character(LEN=2),Dimension(:),allocatable :: element
    integer                         :: iter
    logical                         :: worked
    external ran2
    common /ensemble/ ens
    common /beaddiabatic/ nbdf1,nbdf2
    integer reftraj,rpmddft
    common /reftraj/ reftraj
    common /RPMDDFT/ rpmddft
    common /thinp/ ttaufs

    nbaro = 0
    if (ens.eq.'NPT') then
        nbaro = 1
    endif
#ifdef PARALLEL_BINDING
    call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr)
#endif
    allocate (dist(na,na),dvdre(3,na,nb))
    allocate (ihhh(imaxbin),ihoo(imaxbin),ihoh(imaxbin))
    allocate (m_dipole(3,na/3))

    ! Define some useful local constants and zero-out arrays
    v = 0.d0
    v1 = 0.d0
    v2 = 0.d0
    v3 = 0.d0
    pi = dacos(-1.d0)
    dtps = 1d-3*dt/tofs
    dconv = (ToA * 1d-10 * echarge) / ToDebye
    dconv2 = dconv*dconv

    v_corr = 0.d0
    vir_corr(:,:) = 0.d0
    dvdr_corr(:,:,:) = 0.d0

    if (ttaufs.gt.0.d0 .and. (therm.eq.'AND' .or. therm.eq.'PRA')) then
        thresh = 1d-3*dtps/ttaufs !dtfs is not used in md_static
    else
        thresh = 1.d0/dsqrt(dble(ng))
        thresh = max(0.005d0,thresh)
    end if

    nm = na/3
    no = na/3
    wmass = mass(1)+mass(2)+mass(3)
  

    ihoo(:) = 0.0d0
    ihoh(:) = 0.0d0
    ihhh(:) = 0.0d0
    dist(:,:) = 0.d0
    vir(:,:) = 0.d0
    vir_lf(:,:) = 0.d0

    tavang = 0.d0
    tavoh  = 0.d0

    tavdip = 0.d0
    tavdipx = 0.d0
    tavdipy = 0.d0
    tavdipz = 0.d0
    tavdip2 = 0.d0
    
    vave   = 0.d0
    v1ave  = 0.d0
    v2ave  = 0.d0
    v3ave  = 0.d0
    v1eeav = 0.d0
    v2eeav = 0.d0
    v3eeav = 0.d0
    veeav  = 0.d0
    
    tave   = 0.d0
    tq1ave = 0.d0
    tq2ave = 0.d0
    tavee  = 0.d0
    tq1aee = 0.d0
    tq2aee = 0.d0
    
    rghav = 0.d0
    rgoav = 0.d0
    rgcmav = 0.d0
    volav = 0.d0
    xav = 0.d0
    yav = 0.d0
    zav = 0.d0
    pav = 0.d0
    pvav = 0.d0
    nrdf = 0

#ifdef PARALLEL_BINDING
    if(myid.eq.0) then
#endif
    !zero time COM and Coordinates
    ! Calculate COM for each molecule
    do i = 0,nm-1
        do j = 1,3
            rst(j,i+1) = mass(1)*r(j,3*i+1,1)+ &
            mass(2)*r(j,3*i+2,1)+ mass(3)*r(j,3*i+3,1)
            rst(j,i+1) = rst(j,i+1)/wmass
        end do
    end do
    r0 = r

    call md_static_prepare_traj(nb,pt,pb,print,reftraj,ens.eq.'NPT',&
    itst(3).eq.1.AND.nb.eq.1)

#ifdef PARALLEL_BINDING
    endif
#endif
    ! Evolve for ng steps, calculating properties
    do je = 1,ng
        istep = je
#ifdef PARALLEL_BINDING
        call MPI_bcast(r,SIZE(r),MPI_real8,0,MPI_COMM_WORLD,ierr)
        call MPI_bcast(p,SIZE(p),MPI_real8,0,MPI_COMM_WORLD,ierr)
        call MPI_bcast(dvdr,SIZE(dvdr),MPI_real8,0,MPI_COMM_WORLD,ierr)
        call MPI_bcast(dvdr2,SIZE(dvdr2),MPI_real8,0,MPI_COMM_WORLD,ierr)  
        call MPI_bcast(v,1,MPI_real8,0,MPI_COMM_WORLD,ierr)   
        call MPI_bcast(v1,1,MPI_real8,0,MPI_COMM_WORLD,ierr)   
        call MPI_bcast(v2,1,MPI_real8,0,MPI_COMM_WORLD,ierr)   
        call MPI_bcast(v3,1,MPI_real8,0,MPI_COMM_WORLD,ierr)   
        call MPI_bcast(vir,SIZE(vir),MPI_real8,0,MPI_COMM_WORLD,ierr)   
        call MPI_bcast(vir_lf,SIZE(vir_lf),MPI_real8,0,MPI_COMM_WORLD,ierr)   
#endif

        call evolve(p,r,v,v1,v2,v3,dvdr,dvdr2,dt,mass,na,nb, &
        boxlxyz,z,beta,vir,vir_lf,irun,nbaro)

#ifdef PARALLEL_BINDING
        if (myid.eq.0) then
#endif
        !now calculate all properties
        vave = vave + v
        v1ave = v1ave + v1
        v2ave = v2ave + v2
        if (nb.ne.1) v3ave = v3ave + v3
        if ((print_traj).and.(mod(je,pt).eq.0)) then
            if (print_set_element(1)) then
                call print_vmd_full(r,nb,na,nm,boxlxyz,&
                vmd_traj_set(1)%Unit())
            endif
            if (print_set_element(2)) then
                call print_vmd_full_forces(dvdr,dvdr2,nb,na,boxlxyz,&
                vmd_traj_set(2)%Unit())
            endif
            if (print_set_element(3)) then
                call print_vmd_full_vels(p,mass,nb,na,boxlxyz,&
                vmd_traj_set(3)%Unit())
            end if
        end if
        if (print_beads .AND. (mod(je,pb).eq.0)) then
            if(print_set_element(1))then
                do ib = lbound(vmd_bead_set,1),ubound(vmd_bead_set,1),1
                    call print_vmd_bead(r,nb,ib,na,nm,boxlxyz,&
                    vmd_bead_set(ib,1)%Unit())
                end do
            end if
            if(print_set_element(2))then
                do ib = lbound(vmd_bead_set,1),ubound(vmd_bead_set,1),1
                    call print_vmd_bead_forces(dvdr,dvdr2,nb,ib,na,boxlxyz,&
                    vmd_bead_set(ib,2)%Unit())
                end do
            end if
            if(print_set_element(3))then
                do ib = lbound(vmd_bead_set,1),ubound(vmd_bead_set,1),1
                    call print_vmd_bead_vels(p,mass,nb,ib,na,boxlxyz,&
                    vmd_bead_set(ib,3)%Unit())
                end do
            end if
        end if

        if (reftraj.eq.0) then
            if (mod(je,10).eq.0) then
                nrdf = nrdf + 1
                ! Potential and Virial KE
                call virial_ke(r,dvdr,dvdr2,tv,tvxyz,tq1,tq2,beta,na,nb,mass)  !!!! CHANGED
                tave = tave + tv
                if(rpmddft.eq.0) then
                    tq1ave = tq1ave + tq1
                    tq2ave = tq2ave + tq2
                end if
                if (mod(je,100).eq.0) then
                    write(fh_vt%unit(),*) je*dtps,v,tv/nm
                end if
                ! Calculate pressure
                vol = boxlxyz(1)*boxlxyz(2)*boxlxyz(3)
                call pressure(pres,vir,tv,na,boxlxyz)
                pav = pav + pres
                if (mod(je,100).eq.0) then
                    write(fh_pressure%unit(),*) je*dtps,tobar*pres
                end if
                ! In NPT simulations calculate volume and PV contribution
                if (ens.eq.'NPT') then
                    volav = volav + vol
                    xav = xav + boxlxyz(1)
                    yav = yav + boxlxyz(2)
                    zav = zav + boxlxyz(3)
                    pvav = pvav + pres*vol
                    if (mod(je,100).eq.0) then
                        write(fh_NPT_prop%unit(),*) je*dtps,tokjmol*pres*vol,vol
                        denc = (dble(nm)*wmass*tokgcm3)/vol
                        write(fh_density%unit(),*) je*dtps,denc
                    end if
                end if
                if ((itst(2).eq.1).and.(rpmddft.eq.0)) then   ! not usable for AI-RPMD use
                    if ((nbdf1.gt.0).or.(nbdf2.gt.0)) then
                        ! Exact Estimators for Beadiabatic
                        ! Ewald
                        call forces(r,v1,dvdre,nb,na,boxlxyz,z,vir,2)
                        call virial_ke_ee(r,dvdre,tv,tqe,beta,na,nb)
                        ! LJ
                        call forces(r,v2,dvdre,nb,na,boxlxyz,z,vir,3)
                        !Corrections
                        !call forces(r,v_corr,dvdr_corr,nb,na,boxlxyz,z,vir_corr,8)
                        v2 = v2 + v_corr
                        vir(:,:) = vir(:,:) + vir_corr(:,:)
                        dvdre(:,:,:) = dvdre(:,:,:) + dvdr_corr(:,:,:)              
                        call virial_ke_ee(r,dvdre,tv,tq1,beta,na,nb)
                        tq1 = tq1 + tqe
                 
                        ! Intramolecular
                        call forces(r,v3,dvdre,nb,na,boxlxyz,z,vir,4)
                        call virial_ke_ee(r,dvdre,tv,tq2,beta,na,nb)
                        tv = tv + tq1
                        tavee = tavee + tv
                        tq1aee = tq1aee + tq1
                        tq2aee = tq2aee + tq2
                        v1eeav = v1eeav + v1
                        v2eeav = v2eeav + v2
                        v3eeav = v3eeav + v3
                    end if
                end if
                ! Calculate molecular properties.
                call molprop (r,avang,avoh,nm,nb,z,na)
                tavang = tavang + avang
                tavoh = tavoh + avoh
                ! radius of gyration
                if (nb.gt.1) then
                    call rad_gyr(rgo,rgh,rgcm,r,na,nb,mass)
                    rgoav = rgoav + rgo
                    rghav = rghav + rgh
                    rgcmav = rgcmav + rgcm
                end if
                ! System dipole and dipole squared
                call dipole(r,dipx,dipy,dipz,dip2,dipm,z,na,nb,m_dipole)
                tavdip = tavdip + dipm
                tavdipx = tavdipx + dipx
                tavdipy = tavdipy + dipy
                tavdipz = tavdipz + dipz
                tavdip2 = tavdip2 + dip2
                if ((ptd.gt.0).and.(mod(je,ptd).eq.0)) then
                    if (print(4).eq.1) then
                        call print_mol_dipole(na,m_dipole, 41)
                    end if
                end if
                ! Check convergence of average system dipole moment
                if (mod(je,500).eq.0) then
                    write(fh_dipole_sys%unit(),'(4f12.5)')&
                    je*dtps,dconv*tavdipx/dble(nrdf),dconv*tavdipy/dble(nrdf),&
                    dconv*tavdipz/dble(nrdf)
                end if
                ! Check convergence of dielectric constant
                if (mod(je,500).eq.0) then
                    vol = boxlxyz(1)*boxlxyz(2)*boxlxyz(3)
                    fac = (4.d0*pi*beta)/(3.d0*vol*nrdf)
                    tavdipsq=(tavdipx**2+tavdipy**2+tavdipz**2)/dble(nrdf)
                    eps = 1.d0 + fac*(tavdip2-tavdipsq)
                    write(fh_dielec_conv%unit(),*) je*dtps,eps
                end if
                ! RDF
                !old method of calculating RDF can be found at the bottom of
                !the file
                if(itst(1).eq.1)then
                    ! for the RDF, create an array that contains the element 
                    ! of the water molecules in  OHH order
                    if(.not.allocated(element))then
                        allocate(element(size(r,2)))
                        do i=lbound(r,2),ubound(r,2),3
                            element(i)='O'
                            element(i+1:i+2)='H'
                        end do
                    end if
                    boxmax = max(boxlxyz(1),boxlxyz(2),boxlxyz(3))
                    delr = 0.5d0*boxmax
                    call calc_rdf(r,element,"O ","O ",boxlxyz,0.0d0,delr,ihoo)
                    call calc_rdf(r,element,"O ","H ",boxlxyz,0.0d0,delr,ihoh)
                    call calc_rdf(r,element,"H ","H ",boxlxyz,0.0d0,delr,ihhh)
                end if
            end if
            ! time for a collision with the heat bath?
            if ((therm.eq.'AND').or.(therm.eq.'PRA')) then
                if (ran2(irun,0.d0,1.d0) .lt. thresh) then
                    call sample(p,na,nb,mass,beta,irun,dt)
                endif
            endif
    
            ! Progress indicator
         
            if (ng.gt.1 .and.  mod(je,(ng/10)).eq.0) then
                write(6,*) 10*(je/(ng/10)), ' %'
            endif
        endif
    
        ! Output MSD moved by centroid COM and O molecules in Angstrom^2/ps 
    
        !***** Careful, there is still some mistake here ***** !
        if(itst(3).eq.1 .and. nb.eq.1) then
            if (mod(je,10).eq.0) then
                msd = 0.d0
                !Calculate COM for each molecule
                do i = 0,nm-1
                    do j = 1,3
                        rcm(j,i+1) = mass(1)*r(j,3*i+1,1)+ &
                        mass(2)*r(j,3*i+2,1)+ mass(3)*r(j,3*i+3,1)
                        rcm(j,i+1) = rcm(j,i+1)/wmass
                    end do
                end do
                do i = 1,nm
                    do k = 1,3
                        dr = rcm(k,i)-rst(k,i)
                        msd = msd + dr*dr
                    end do
                end do
                msd = msd / dble(nm)
                ! O MSD in Angstrom^2/ps
                msdO = 0.d0
                do i = 1,nm
                    ii = 3*(i-1) + 1
                    do k = 1,3
                        dr = r(k,ii,1)-r0(k,ii,1)
                        msdO = msdO + dr*dr
                    enddo
                enddo
                msdO = msdO / dble(nm)
                write(fh_MSD%unit(),'(3f12.6)')je*dtps,msd*toA*toA,msdO*toA*toA
            endif
        endif
#ifdef PARALLEL_BINDING
        endif
#endif
    enddo
#ifdef PARALLEL_BINDING
    if (myid.eq.0) then
#endif
    !closing of any handled file can be unconditional, because 
    !the routine throws no error if the file is closed
    worked = fh_vt%close()
    worked = fh_dipole_sys%close()
    worked = fh_dielec_conv%close()
    worked = fh_pressure%close()
    worked = fh_density%close()     
    if (print_traj) then
        do iter=lbound(vmd_traj_set,1),ubound(vmd_traj_set,1),1
            worked=vmd_traj_set(iter)%close()
        end do
        worked = fh_mol_dipole%close()
    endif
    if(print_beads)then
        do iter=lbound(vmd_bead_set,2),ubound(vmd_bead_set,2),1
            do ib = lbound(vmd_bead_set,1),ubound(vmd_bead_set,1),1
                worked=vmd_bead_set(ib,iter)%close()
            end do
        end do
    endif

    if (reftraj.eq.0) then
        tavang = tavang / dble(nrdf * nm * nb)
        tavoh = tavoh / dble(nrdf * nm * nb)

        tavdip = tavdip / dble(nrdf)
        tavdipx = tavdipx / dble(nrdf)
        tavdipy = tavdipy / dble(nrdf)
        tavdipz = tavdipz / dble(nrdf)
        tavdip2 = tavdip2 / dble(nrdf)

        tave  = tave / dble(nrdf * nm)
        tavee = tavee / dble(nrdf * nm)
        tq1ave = tq1ave / dble(nrdf * nm)
        tq2ave = tq2ave / dble(nrdf * nm)
        tq1aee = tq1aee / dble(nrdf * nm)
        tq2aee = tq2aee / dble(nrdf * nm)

        vave = vave / dble(ng * nm)
        v1ave = v1ave / dble(ng * nm)
        v2ave = v2ave / dble(ng * nm)
        v3ave = v3ave / dble(ng * nm)
        v1eeav = v1eeav / dble(nrdf * nm)
        v2eeav = v2eeav / dble(nrdf * nm)
        v3eeav = v3eeav / dble(nrdf * nm)
        veeav = veeav / dble(nrdf * nm)

        pav = pav / dble(nrdf)
        if (ens.eq.'NPT') then
            volav = volav / dble(nrdf)
            xav = xav / dble(nrdf)
            yav = yav / dble(nrdf)
            zav = zav / dble(nrdf)
            pvav = pvav / dble(nrdf * nm)
        end if

        rgoav = rgoav / dble(nrdf)
        rghav = rghav / dble(nrdf)
        rgcmav = rgcmav / dble(nrdf)

        ! Print average KE and V

        write(6,*)
        write(6,*)'* Energies '
        write(6,*)'----------------------------'
        if (nb.eq.1) then
            write(6,*)'<V> per molecule = ',toKjmol*vave,' KJ mol^-1'
            write(6,*)'<V>_inter = ', toKjmol*v1ave,' KJ mol^-1'
            write(6,*)'<V>_intra = ', toKjmol*v2ave,' KJ mol^-1'
            write(6,*)
            write(6,*)'<Virial KE> per molecule = ',toKjmol*tave, &
            ' KJ mol^-1'
            write(6,*)'<Virial KE> per molecule = ',toK*tave*(2.0/3.0),' K ' 
            write(6,*)'<Virial KE> per molecule = ',toKjmol*tave,' KJ mol^-1 '
            write(6,*)'<Virial KE>_inter = ', toKjmol*tq1ave,' KJ mol^-1'
            write(6,*)'<Virial KE>_intra = ', toKjmol*tq2ave,' KJ mol^-1'
        else
            write(6,*)'<V> per molecule = ',toKjmol*vave
            write(6,*)'<V>_ew = ', toKjmol*v1ave,' KJ mol^-1'
            write(6,*)'<V>_lj = ', toKjmol*v2ave,' KJ mol^-1'
            write(6,*)'<V>_inter = ', toKjmol*(v1ave+v2ave),' KJ mol^-1'
            write(6,*)'<V>_intra = ', toKjmol*v3ave,' KJ mol^-1'
            write(6,*)
            write(6,*)'<Virial KE> per molecule = ',toK*tave*(2.0/3.0),' K '
            write(6,*)'<Virial KE> per molecule = ',toKjmol*tave,' Kjmol-1 '
            write(6,*)'<Virial KE>_inter = ', toKjmol*tq1ave,' KJ mol^-1'
            write(6,*)'<Virial KE>_intra = ', toKjmol*tq2ave,' KJ mol^-1'
            if (itst(2).eq.1) then
                if ((nbdf1.gt.0).or.(nbdf2.gt.0)) then
                    write(6,*)
                    write(6,*)'* Energies using Exact Estimators '
                    write(6,*)'----------------------------'
                    write(6,*)'<V> per molecule = ',toKjmol* &
                    (v1eeav+v2eeav+v3eeav),' KJ mol^-1'
                    write(6,*)'<V>_ew = ', toKjmol*v1eeav,' KJ mol^-1'
                    write(6,*)'<V>_lj = ', toKjmol*v2eeav,' KJ mol^-1'
                    write(6,*)'<V>_inter = ', toKjmol*(v1eeav+v2eeav), &
                    ' KJ mol^-1'
                    write(6,*)'<V>_intra = ', toKjmol*v3eeav,' KJ mol^-1'
                    write(6,*)
                    write(6,*)'<Virial KE> per molecule = ', toKjmol*tavee, &
                    ' KJ mol^-1'
                    write(6,*)'<Virial KE>_inter = ', toKjmol*tq1aee, &
                    ' KJ mol^-1'
                    write(6,*)'<Virial KE>_intra = ', toKjmol*tq2aee, &
                    ' KJ mol^-1'
                endif
                if (rpmddft.eq.1) then
                    write(6,*)
                    write(6,*)'* Energies using Exact Estimators '
                    write(6,*)'----------------------------'
                    write(6,*)'<V> per molecule = ',toKjmol* veeav, &
                    ' KJ mol^-1'
                    write(6,*)'<Virial KE> per molecule = ', toKjmol*tavee, &
                    ' KJ mol^-1'
                endif
            endif
        endif
    
        ! Print Radius of Gyration
        if (nb.gt.1) then
            write(6,*)
            write(6,*)'* Radius of Gyration '
            write(6,*)'----------------------------'
            write(6,*)'RMS Rg Hydrogen = ', dsqrt(rghav),' bohr'
            write(6,*)'RMS Rg Oxygen   = ', dsqrt(rgoav),' bohr'
            write(6,*)'RMS Rg COM      = ', dsqrt(rgcmav),' bohr'
        endif
    
        ! Pressure, volume and enthalpy
    
        write(6,*)
        write(6,*) '* Volume,Pressure and Enthalpy'
        write(6,*) '-------------------------------'
        write(6,*) ' Average Pressure = ', pav*tobar, ' bar'
        if (ens.eq.'NPT') then
            denav = (dble(nm)*wmass*tokgcm3)/volav
            write(6,*) ' Average Volume   = ', volav,' bohr^3'
            write(6,*) ' Average Length x = ', xav ,' bohr'
            write(6,*) ' Average Length y = ', yav ,' bohr'
            write(6,*) ' Average Length z = ', zav ,' bohr'
            write(6,*) ' Average Density  = ', denav,' g cm^-3'
            write(6,*) ' Average PV       = ', tokjmol*pvav, ' KJ mol^-1'
        endif
        write(6,*)
    
        ! Dielectric constant
    
        vol = boxlxyz(1)*boxlxyz(2)*boxlxyz(3)
        fac = (4.d0*pi*beta)/(3.d0*vol)
        tavdipsq = tavdipx**2 + tavdipy**2 + tavdipz**2
        eps = 1.d0 + fac * (tavdip2-tavdipsq)
        write(6,*)' * Dielectric constant '
        write(6,*)'-------------------------------------------'
        write(6,*) ' <M_x>               =  ', dconv*tavdipx,' Debye'
        write(6,*) ' <M_y>               =  ', dconv*tavdipy, ' Debye'
        write(6,*) ' <M_z>               =  ', dconv*tavdipz, ' Debye'
        write(6,*) ' <M>^2               =  ', dconv2*tavdipsq,' Debye'
        write(6,*) ' <M^2>               =  ', dconv2*tavdip2, ' Debye'
        write(6,*) ' Dielectric constant =  ', eps
        write(6,*)
    
        ! Print molecular properties
    
        write(6,*)'* Molecular Properties: '
        write(6,*)'-------------------------------------------'
        write(6,*)' Average O-H bond length =  ',tavoh*ToA,' Angstroms '
        write(6,*)' Average bond angle      =  ',tavang,' degrees '
        write(6,*)' Average dipole          =  ',dconv*tavdip,' Debye '
        write(6,*)
    
        ! Print the RDF output
    
        write (6,*)'* Static properties calculation complete. '
        if (itst(1).eq.1) then
            call print_rdf(ihoo,ihoh,ihhh,na,nb,boxlxyz,nrdf,'')
        endif
    endif
#ifdef PARALLEL_BINDING
    endif
#endif
    deallocate (dist,ihhh,ihoo,ihoh,m_dipole)

    return
end subroutine md_static
