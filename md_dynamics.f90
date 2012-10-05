subroutine dynamics(nt,m,p,r,dvdr,dvdr2,na,nm,nb,boxlxyz,z,beta, &
dt,mass,irun,itcf,itst,pt,pb,print,iskip,ntherm,vacfac)
    use thermostat
    use intmod
    use avmod
    implicit none
  include 'globals.inc'
#ifdef PARALLEL_BINDING
  include 'mpif.h' !parallel
#endif
    ! ------------------------------------------------------------------
    ! Routine to calculate dynamical properties
    ! ------------------------------------------------------------------
    integer nt,m,na,nm,nb,irun,itcf(3),itst(3),pt,pb
    integer k,i,ib,print(3),iskip,ntherm,vacfac,ierr,myid,eefile
    real(8) r(3,na,nb),p(3,na,nb),dvdr(3,na,nb),dvdr2(3,na,nb)
    real(8) boxlxyz(3),z(na),mass(na),vir(3,3),vir_lf(3,3)
    real(8) beta,dt,dtfs,dtps,cscale,wt,tv,pe,dqq,dqx,dqy,dqz,ran2,thresh,ttaufs
    real(8) v,v1,v2,v3,dtv,dqqt,davx,davy,davz,dpe,sum
    real(8), allocatable :: ct(:),dct(:)
    real(8), allocatable :: dmo1(:,:),dmot1(:,:)
    real(8), allocatable :: rmtr(:),rmtv(:),dmtr(:),dmtv(:)
    real(8), allocatable :: irmtr(:,:,:),irmtv(:,:,:),idmtr(:,:,:),idmtv(:,:,:)
    real(8), allocatable :: dcvinter(:),dcvintra(:),dke(:)
    real(8), allocatable :: cvinter(:), cvintra(:),cke(:)
    integer(8), allocatable :: ihhh(:),ihoo(:),ihoh(:)
    character(128) file_name
    external ran2
    common /thinp/ ttaufs


    ! Define some local variable values
#ifdef PARALLEL_BINDING
	call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr)
!	write(*,*) "dynamics:",myid
#endif
    dtfs = dt/tofs
    dtps = dtfs*1.d-3
    cscale = (toA/(dtps/dt))**2
    wt = 1.d0 / dble(m)

    if (ttaufs.gt.0.d0 .and. (therm.eq.'AND' .or. therm.eq.'PRA')) then
        thresh = dtfs/ttaufs
    else
        thresh = 1.d0/dsqrt(dble(ntherm))   !for AND thermostat between trajectories
        thresh = max(thresh,0.01d0)         !same definition as in md_eq.f90
    end if

    ! Allocate relevant arrays

    allocate (ct(0:nt),dct(0:nt))
    allocate (rmtr(0:nt),rmtv(0:nt))
    allocate (dmtr(0:nt),dmtv(0:nt))
    allocate (irmtr(0:nt,4,2),irmtv(0:nt,4,2))
    allocate (idmtr(0:nt,4,2),idmtv(0:nt,4,2))
    allocate (dmo1(6,0:nt),dmot1(6,0:nt))
    allocate (cvinter(0:2*nt), dcvinter(0:2*nt))
    allocate (cvintra(0:2*nt), dcvintra(0:2*nt))
    allocate (cke(0:2*nt), dke(0:2*nt))
    allocate (ihhh(imaxbin),ihoo(imaxbin),ihoh(imaxbin))

    tv = 0.d0
    pe = 0.d0
    dqq = 0.d0
    dqx = 0.d0
    dqy = 0.d0
    dqz = 0.d0
    ct(:) = 0.d0
    rmtr(:) = 0.d0
    rmtv(:) = 0.d0
    irmtr(:,:,:) = 0.d0
    irmtv(:,:,:) = 0.d0
    dmo1(:,:) = 0.d0
    cvinter(:) = 0.d0
    cvintra(:) = 0.d0
    cke(:) = 0.d0

    ihoo(:) = 0
    ihoh(:) = 0
    ihhh(:) = 0

#ifdef PARALLEL_BINDING
if(myid.eq.0) then
#endif
    if (itcf(1).eq.1) then
        open(71,file='diff.out')
    endif

    if (itst(2).eq.1) then
        eefile=8765
        open(eefile,file='trajee.out')
    end if

    if (pt.gt.0) then
        if (print(1).eq.1) then
            open(36,file='vmd_traj.xyz')
        endif
        if (print(2).eq.1) then
            open(37,file='vmd_traj.frc')
        endif
        if (print(3).eq.1) then
            open(38,file='vmd_traj.vel')
        endif
    endif
    if (pb.gt.0) then
        if (print(1).eq.1) then
            do ib = 1,nb
                write(file_name,'(A,I0,A)') 'vmd_bead-',ib,'.xyz'
                open(300127+ib,file=file_name)
            enddo
        endif
        if (print(2).eq.1) then
            do ib = 1,nb
                write(file_name,'(A,I0,A)') 'vmd_bead-',ib,'.frc'
                open(400127+ib,file=file_name)
            enddo
        endif
        if (print(3).eq.1) then
            do ib = 1,nb
                write(file_name,'(A,I0,A)') 'vmd_bead-',ib,'.vel'
                open(500127+ib,file=file_name)
            enddo
        endif
    endif
#ifdef PARALLEL_BINDING
endif
#endif

    if (vacfac.ne.1) then
        call intinit(na,toA*boxlxyz,dtfs,(2*nt)/iskip+1)
        call avinit(m)
    end if

    do k = 1,m !Anzahl paralleler Trajektorien für Dynamik
#ifdef PARALLEL_BINDING
			if(myid.eq.0) then
#endif

			    write(eefile,*)
			    write(eefile,*)'Trajectory',k
			    write(eefile,*)'----------------------------'

        ! Find a new configuration
     
        if (therm.eq.'AND') then
            ! Resample momenta
            call sample(p,na,nb,mass,beta,irun,dt)
        endif

#ifdef PARALLEL_BINDING
endif
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
write(*,*) "myid in evolve:", myid
#endif

        !thermalize for ntherm steps if AND or PRA thermostat
        if ((therm.eq.'AND').or.(therm.eq.'PRA')) then
            do i = 1,ntherm
                call evolve(p,r,v,v1,v2,v3,dvdr,dvdr2,dt,mass,na,nb, &
                boxlxyz,z,beta,vir,vir_lf,irun,0)
#ifdef PARALLEL_BINDING
								if(myid.eq.0) then
#endif			
                !if ((therm.eq.'AND').or.(therm.eq.'PRA')) then
                !    if (ran2(irun,0.d0,1.d0) .lt. thresh) then
                !        call sample(p,na,nb,mass,beta,irun,dt)
                !    endif
                !endif
#ifdef PARALLEL_BINDING
								endif
		 						call MPI_bcast(p,SIZE(p),MPI_real8,0,MPI_COMM_WORLD,ierr)
#endif			
            enddo
        end if

        ! Run an NVT (NVE if therm == AND) trajectory of 2*nt time steps
        call trajectory(p,r,v,dvdr,dvdr2,dct,dtv,dmtr,dmtv,idmtr,idmtv, &
        dqqt,davx,davy,davz,dmot1,ihoo,ihoh,ihhh,nt,na,nb,boxlxyz,z, &
        beta,dt,mass,itcf,itst,nm,dpe,irun,dcvinter,dcvintra,&
        dke,pt,pb,print,iskip,vacfac,eefile)
#ifdef PARALLEL_BINDING
				if(myid.eq.0) then
#endif			
        write(6,*)'* Completed sampling: ',k,' of ', m
        call flush(6)
        tv = tv + wt*dtv
        pe = pe + wt*dpe
        dqq = dqq + dqqt
        dqx = dqx + davx
        dqy = dqy + davy
        dqz = dqz + davz

        ! Accumulate correlation functions

        ct(:) = ct(:) + dct(:)
        rmtr(:) = rmtr(:) + wt*dmtr(:)
        rmtv(:) = rmtv(:) + wt*dmtv(:)
        irmtr(:,:,:) = irmtr(:,:,:) + wt*idmtr(:,:,:)
        irmtv(:,:,:) = irmtv(:,:,:) + wt*idmtv(:,:,:)
        dmo1(:,:) = dmo1(:,:) + wt*dmot1(:,:)
        cvinter(:) = cvinter(:) + wt*dcvinter(:)
        cvintra(:) = cvintra(:) + wt*dcvintra(:)
        cke(:) = cke(:) + wt*dke(:)

        if (vacfac.ne.1) then
            write(6,*) 'Mean/Corr calculation and module reset...'
            call meancalc()
            write(6,*) 'Calculation done.'
            call avupdate()
            write(6,*) 'Average arrays updated.'
            call intreset()
            write(6,*) 'Module resetted.'
        end if

        ! Check convergence of the diffusion constant:

        if (itcf(1).eq.1) then
            call simp(nt+1,dtps,ct(0),sum)
            write(71,*)k,cscale*sum/(3.d0*dble(k))
            call flush(71)
        endif
#ifdef PARALLEL_BINDING
		endif
#endif	

    enddo

    if (vacfac.ne.1) then
        call avcalcm()
        call avintprint()
    end if
#ifdef PARALLEL_BINDING
				if(myid.eq.0) then
#endif	
    if (pt.gt.0) then
        if (print(1).eq.1) then
            close(unit=36)
        endif
        if (print(2).eq.1) then
            close(unit=37)
        endif
        if (print(3).eq.1) then
            close(unit=38)
        endif
    endif
    if (pb.gt.0) then
        if (print(1).eq.1) then
            do ib = 1,nb
                close(unit=300127+ib)
            enddo
        endif
        if (print(2).eq.1) then
            do ib = 1,nb
                close(unit=400127+ib)
            enddo
        endif
        if (print(3).eq.1) then
            do ib = 1,nb
                close(unit=500127+ib)
            enddo
        endif
    endif

    ! Print the total RDFs.
    if (itst(1).eq.1) then
      call print_rdf(ihoo,ihoh,ihhh,na,boxlxyz,nt,nb,(2*nt/iskip+1)*m)
    endif

    ! Print the (unit.potential) correlation functions.
    open(51,file='c1v.out',status='unknown')
    do i = 0, 2*nt
        write(51,'(4f14.8)')i*dtfs*1d-3,cvinter(i),cvintra(i),cke(i)
    enddo
    close(unit=51)

    ! Print the velocity auto-correlation function.

    if (itcf(1).eq.1) then
        ct(:) = ct(:) * wt
        call print_cvv(ct,nt,dtfs,dt)
    endif

    ! Print the dipole autocorrelation functions.

    if (itcf(2).eq.1) then
        dqq = dqq / dble(m)
        dqx = dqx / dble(m)
        dqy = dqy / dble(m)
        dqz = dqz / dble(m)
        call print_dipole(dtfs,nt,rmtr,rmtv,dqq,dqx,dqy,dqz, &
        boxlxyz,beta,dt,iskip)
        if (vacfac.ne.1) then
            call print_dipole_int(dtfs,nt,irmtr,irmtv,iskip)
        end if
    endif
   
   
    ! Print the orientational correlation functions.

    if (itcf(3).eq.1) then
        call print_orientation(dtfs,nt,dmo1,iskip)
    endif
    write(6,*)'-------------------------------------------'

    ! Print energies

    write(6,*)
    write(6,*)'* Calculation results  '
    write(6,*)
    write(6,'(1x,"Energies per atom:")' )
    tv = toKjmol*(tv/na)
    pe = toKjmol*(pe/na)
    write (6,65) tv
    write (6,66) pe
65  format(1x,'<KE> = ',f10.4,' Kj/mol')
66  format(1x,'<PE> = ',f10.4,' Kj/mol')

    deallocate (ct,dct,rmtr,irmtr,rmtv,irmtv,dmtr,idmtr,dmtv,idmtv,dmo1,dmot1)
#ifdef PARALLEL_BINDING
endif
#endif
    return
end subroutine dynamics

subroutine trajectory(p,r,v,dvdr,dvdr2,ct,dtv,dmtr,dmtv,idmtr,idmtv, &
dqq,davx,davy,davz,dmot1,ihoo,ihoh,ihhh,nt,na,nb,boxlxyz,z, &
beta,dt,mass,itcf,itst,nm,pe,irun,dcvinter, &
dcvintra,dke,pt,pb,print,iskip,vacfac,eefile)
    use intmod
    implicit none
  include 'globals.inc'
#ifdef PARALLEL_BINDING
	include 'mpif.h'
#endif
    ! -----------------------------------------------------------------
    ! Sample dynamic properties over NVE trajectory.
    ! -----------------------------------------------------------------
    integer na,nb,nt,itcf(3),itst(3),nm,irun,myid,ierr,eefile
    integer i,ib,j,jt,it,it10,jt10,k,pt,pb,print(3),iskip,vacfac,rpmddft,nbdf1,nbdf2,ibin,ii,jj
    real(8) dt,beta,alpha,alpha2,wm,wh,v,v1,v2,v3,w1,w2,w3,dx,dy,dz,rij,delr,boxmax
    real(8) em,tv,dipx,dipy,dipz,tbar,wt,wt1,vir(3,3),vir_lf(3,3)
    real(8) fac1,fac2,pe,dtq1,dtq2,ecut,tvxyz(3)
    real(8) dtv,dqq,davx,davy,davz,dip2,dipm,dotx,doty,dotz
    real(8) emi,vew,vlj,vint,rke
    real(8) tavee,tq1aee,tq2aee,v1eeav,v2eeav,v3eeav,tq1,tq2,tqe
    common /ew_param/ alpha,ecut,wm,wh
    common /RPMDDFT/ rpmddft
    common /beaddiabatic/ nbdf1,nbdf2

    ! Passed Arrays

    real(8) p(3,na,nb),r(3,na,nb),dvdr(3,na,nb),dvdr2(3,na,nb)
    real(8) boxlxyz(3),z(na),mass(na)
    real(8) ct(0:nt),dmot1(6,0:nt)
    real(8) dmtr(0:nt),dmtv(0:nt),idmtr(0:nt,4,2),idmtv(0:nt,4,2)
    real(8) dcvinter(0:2*nt), dcvintra(0:2*nt),dke(0:2*nt)
    integer(8) ihoo(imaxbin),ihoh(imaxbin),ihhh(imaxbin)

    ! Cvv, OCF and Dipole Working Arrays

    real(8), allocatable :: pc(:,:,:), pca(:,:), rca(:,:)
    real(8), allocatable :: ax1(:,:,:),ay1(:,:,:),az1(:,:,:)
    real(8), allocatable :: ax2(:,:,:),ay2(:,:,:),az2(:,:,:)
    real(8), allocatable :: dmxr(:),dmyr(:),dmzr(:)
    real(8), allocatable :: dmxv(:),dmyv(:),dmzv(:)
    real(8), allocatable :: idmxr(:,:,:),idmyr(:,:,:),idmzr(:,:,:)
    real(8), allocatable :: idmxv(:,:,:),idmyv(:,:,:),idmzv(:,:,:)
    real(8), allocatable :: dist(:,:),dvdre(:,:,:)

    ! Allocate arrays

    if (itcf(2).eq.1) then
        allocate (dmxr(0:2*nt),dmyr(0:2*nt),dmzr(0:2*nt))
        allocate (dmxv(0:2*nt),dmyv(0:2*nt),dmzv(0:2*nt))
        allocate (idmxr(0:2*nt,4,2),idmyr(0:2*nt,4,2),idmzr(0:2*nt,4,2))
        allocate (idmxv(0:2*nt,4,2),idmyv(0:2*nt,4,2),idmzv(0:2*nt,4,2))
    end if

    if (itcf(3).eq.1) then
        allocate (ax1(3,nm,0:(2*nt/iskip)))
        allocate (ay1(3,nm,0:(2*nt/iskip)))
        allocate (az1(3,nm,0:(2*nt/iskip)))
        allocate (ax2(6,nm,0:(2*nt/iskip)))
        allocate (ay2(6,nm,0:(2*nt/iskip)))
        allocate (az2(6,nm,0:(2*nt/iskip)))
    end if

    allocate (pc(3,nm,0:2*nt),pca(3,na),rca(3,na),dvdre(3,na,nb),dist(na,na))

    ! Define some useful constants and zero-out arrays

    alpha2 = 0.5d0 * (1.d0 - alpha)
    nm = na / 3
    em = mass(1) + mass(2) + mass(3)

    if (itst(2).eq.1) then

        v1eeav = 0.d0
        v2eeav = 0.d0
        v3eeav = 0.d0

        tavee  = 0.d0
        tq1aee = 0.d0
        tq2aee = 0.d0

    end if

    if (itcf(2).eq.1) then

        dmxr(:) = 0.d0
        dmyr(:) = 0.d0
        dmzr(:) = 0.d0
        dmxv(:) = 0.d0
        dmyv(:) = 0.d0
        dmzv(:) = 0.d0
  
        idmxr(:,:,:) = 0.d0
        idmyr(:,:,:) = 0.d0
        idmzr(:,:,:) = 0.d0
        idmxv(:,:,:) = 0.d0
        idmyv(:,:,:) = 0.d0
        idmzv(:,:,:) = 0.d0
    end if

    if (itcf(3).eq.1) then
        ax1(:,:,:) = 0.d0
        ay1(:,:,:) = 0.d0
        az1(:,:,:) = 0.d0
        ax2(:,:,:) = 0.d0
        ay2(:,:,:) = 0.d0
        az2(:,:,:) = 0.d0
    end if

    vir(:,:) = 0.d0
    vir_lf(:,:) = 0.d0

    ct(:) = 0.d0
    dmtr(:) = 0.d0
    dmtv(:) = 0.d0
    idmtr(:,:,:) = 0.d0
    idmtv(:,:,:) = 0.d0
    dmot1(:,:) = 0.d0
    dcvinter(:) = 0.d0
    dcvintra(:) = 0.d0
    dke(:) = 0.d0

    dqq  = 0.d0
    davx = 0.d0
    davy = 0.d0
    davz = 0.d0
    tv = 0.d0
    pe = 0.d0
    nm = na/3
    dist(:,:) = 0.d0

    if (3*nm .ne. na) stop 'trajectory : na is NOT 3*nm !'
#ifdef PARALLEL_BINDING
	call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr)
  call MPI_bcast(r,SIZE(r),MPI_real8,0,MPI_COMM_WORLD,ierr)
#endif
    ! Calculate initial forces....

    call full_forces(r,na,nb,v,vew,vlj,vint,vir,z,boxlxyz,dvdr,dvdr2)
    dcvinter(0) = vew + vlj ! inter
    dcvintra(0) = vint      ! intra

    ! Loop over 2 * nt time steps.
  
    do jt = 0,2*nt
     
        if (jt .gt. 0) then
        
            ! Evolve the system by dt
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
		 write(*,*) "myid in trajectory:", myid
#endif
           call evolve(p,r,v,v1,v2,v3,dvdr,dvdr2,dt,mass,na,nb, &
            boxlxyz,z,beta,vir,vir_lf,irun,0)
#ifdef PARALLEL_BINDING
					if(myid.eq.0) then
#endif
            if ((pt.gt.0).and.(mod(jt,pt).eq.0)) then
                if (print(1).eq.1) then
                    call print_vmd_full(r,nb,na,nm,boxlxyz,36)
                endif
                if (print(2).eq.1) then
                    call print_vmd_full_forces(dvdr,dvdr2,nb,na,boxlxyz,37)
                endif
                if (print(3).eq.1) then
                    call print_vmd_full_vels(p,mass,nb,na,boxlxyz,38)
                endif
            endif
            if ((pb.gt.0).and.(mod(jt,pb).eq.0)) then
                if (print(1).eq.1) then
                    do ib = 1,nb
                        call print_vmd_bead(r,nb,ib,na,nm,boxlxyz,300127+ib)
                    enddo
                endif
                if (print(2).eq.1) then
                    do ib = 1,nb
                        call print_vmd_bead_forces(dvdr,dvdr2,nb,ib,na,boxlxyz,400127+ib)
                    enddo
                endif
                if (print(3).eq.1) then
                    do ib = 1,nb
                        call print_vmd_bead_vels(p,mass,nb,ib,na,boxlxyz,500127+ib)
                    enddo
                endif
            endif

            call virial_ke(r,dvdr,dvdr2,dtv,tvxyz,dtq1,dtq2,beta,na,nb)
            call dipole (r,dipx,dipy,dipz,dip2,dipm,z,na,nb)
            davx = davx + dipx
            davy = davy + dipy
            davz = davz + dipz
            dqq = dqq + dip2
            tv = tv + dtv
            pe = pe + v

            ! Store PE contributions at time t.
            if (nb.gt.1) then
                dcvinter(jt) = v1 + v2 ! vew+vlj
                dcvintra(jt) = v3      ! vintra
            else
                dcvinter(jt) = v1 ! LF pot
                dcvintra(jt) = v2 ! HF pot
            endif
#ifdef PARALLEL_BINDING
					endif
#endif

        endif

#ifdef PARALLEL_BINDING
						if(myid.eq.0) then
#endif

        if (vacfac.ne.1) then
                if (mod(jt,iskip).eq.0) then
                    call center_atoms(r,rca,na,nb)
                    call intupdate(toA*rca)     !interface and hbond calculations
                end if
        end if

        ! KE
        rke = 0.d0
        do k = 1, nb
            do i = 1, na
                rke = rke + (p(1,i,k)**2+p(2,i,k)**2+p(3,i,k)**2)/(2.d0*mass(i))
            enddo
        enddo
        dke(jt) = rke

        ! RDF g_oo,g_oh,g_hh
        !------------------------
        if ((itst(1).eq.1) .and. (mod(jt,iskip).eq.0)) then
           do k = 1, nb
              do i = 2, na
                 do j = 1, i-1
                    dx = r(1,i,k) - r(1,j,k)
                    dy = r(2,i,k) - r(2,j,k)
                    dz = r(3,i,k) - r(3,j,k)
                    dx = dx-boxlxyz(1)*dble(nint(dx/boxlxyz(1)))
                    dy = dy-boxlxyz(2)*dble(nint(dy/boxlxyz(2)))
                    dz = dz-boxlxyz(3)*dble(nint(dz/boxlxyz(3)))
                    dist(i,j) = dsqrt(dx*dx+dy*dy+dz*dz)
                    dist(j,i) = dist(i,j)
                 enddo
              enddo

              do i = 1,na
                 dist(i,i) = 0.d0
              enddo

              ! calculate and bin O - O pair distances


              boxmax = max(boxlxyz(1),boxlxyz(2),boxlxyz(3))
              delr = dble(0.5d0*boxmax/dble(imaxbin))
              do i = 1, nm - 1
                 do j = i + 1, nm
                    ii = 3*i - 2
                    jj = 3*j - 2
                    rij = dist(ii,jj)
                    ibin = int(rij/delr) + 1
                    if (ibin.le.imaxbin) then
                       ihoo(ibin) = ihoo(ibin) + 2
                    endif
                 enddo
              enddo

              ! calculate and bin O - H pair distances

              do i = 1, nm
                 do j = 1, na, 3
                    ii = 3*i - 2
                    jj = j+1
                    rij = dist(ii,jj)
                    ibin = int( rij / delr) + 1
                    if (ibin.le.imaxbin) then
                       ihoh(ibin) = ihoh(ibin) + 1
                    endif
                    rij=dist(ii,jj+1)
                    ibin = int( rij / delr) + 1
                    if (ibin.le.imaxbin) then
                       ihoh(ibin) = ihoh(ibin) + 1
                    endif
                 enddo
              enddo

              ! calculate and bin H - H pair distances

              do i = 1, na-1
                 do j = i+1, na
                    if (mod(i-1,3).ne.0.and.mod(j-1,3).ne.0) then
                       rij = dist(i,j)
                       ibin = int(rij/delr) + 1
                       if (ibin.le.imaxbin) then
                          ihhh(ibin) = ihhh(ibin) + 2
                       endif
                    endif
                 enddo
              enddo
           enddo
        endif

        ! Exact Estimators for Energies
        ! -----------------------------

        if ((itst(2).eq.1) .and. (mod(jt,iskip).eq.0)) then   ! not suitable for AI-RPMD use
            if ((rpmddft.eq.0).and.((nbdf1.gt.0).or.(nbdf2.gt.0))) then

                ! Exact Estimators for Beadiabatic

                ! Ewald
                call forces(r,v1,dvdre,nb,na,boxlxyz,z,vir,2)
                call virial_ke_ee(r,dvdre,tv,tqe,beta,na,nb)

                ! LJ
                call forces(r,v2,dvdre,nb,na,boxlxyz,z,vir,3)
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
            endif
        endif


        ! Cvv - centroid momenta
        ! -----------------------

        if (itcf(1) .eq. 1) then
            call center_water(p,pc(1,1,jt),nm,nb)
        endif

        ! OCF
        ! ----

        if (itcf(3).eq.1.and.mod(jt,iskip).eq.0)then
            jt10 = jt / iskip
            call orientation(r,ax1,ay1,az1,ax2,ay2,az2, &
            jt10,na,nb,nm,nt,z,iskip)
        endif

        ! Dipole Spectrum
        ! ----------------

        ! Calculation of dipole and dipole derivatives using centroids.
        ! Note that using the centroids is NOT correct if the charges
        ! vary with molecular geometry!!!

        if (itcf(2) .eq. 1.and.mod(jt,iskip).eq.0) then
            call center_atoms(r,rca,na,nb)
            call center_atoms(p,pca,na,nb)

            ! Replace the oxygen coordinates with the M-site coordinates.

            do j = 1, na,3
                fac1 = (z(j+1)+alpha2*z(j))
                fac2 = (z(j+2)+alpha2*z(j))
                w1 = 1.d0 / mass(j)
                w2 = 1.d0 / mass(j+1)
                w3 = 1.d0 / mass(j+2)
                dmxr(jt) = dmxr(jt) + alpha*z(j)*rca(1,j)
                dmyr(jt) = dmyr(jt) + alpha*z(j)*rca(2,j)
                dmzr(jt) = dmzr(jt) + alpha*z(j)*rca(3,j)
                dmxr(jt) = dmxr(jt) + fac1*rca(1,j+1)
                dmyr(jt) = dmyr(jt) + fac1*rca(2,j+1)
                dmzr(jt) = dmzr(jt) + fac1*rca(3,j+1)
                dmxr(jt) = dmxr(jt) + fac2*rca(1,j+2)
                dmyr(jt) = dmyr(jt) + fac2*rca(2,j+2)
                dmzr(jt) = dmzr(jt) + fac2*rca(3,j+2)

                dmxv(jt) = dmxv(jt) + alpha*z(j)*pca(1,j)*w1
                dmyv(jt) = dmyv(jt) + alpha*z(j)*pca(2,j)*w1
                dmzv(jt) = dmzv(jt) + alpha*z(j)*pca(3,j)*w1
                dmxv(jt) = dmxv(jt) + fac1 * pca(1,j+1) * w2
                dmyv(jt) = dmyv(jt) + fac1 * pca(2,j+1) * w2
                dmzv(jt) = dmzv(jt) + fac1 * pca(3,j+1) * w2
                dmxv(jt) = dmxv(jt) + fac2 * pca(1,j+2) * w3
                dmyv(jt) = dmyv(jt) + fac2 * pca(2,j+2) * w3
                dmzv(jt) = dmzv(jt) + fac2 * pca(3,j+2) * w3
            enddo
        
            if (vacfac.ne.1) then
                do ib=1,2
                    do i=1,4
                        do j = 1, na,3

                            if (inint(j/3+1,i,ib,wntc)) then

                                fac1 = (z(j+1)+alpha2*z(j))
                                fac2 = (z(j+2)+alpha2*z(j))
                                w1 = 1.d0 / mass(j)
                                w2 = 1.d0 / mass(j+1)
                                w3 = 1.d0 / mass(j+2)
                                idmxr(jt,i,ib) = idmxr(jt,i,ib) + alpha*z(j)*rca(1,j)
                                idmyr(jt,i,ib) = idmyr(jt,i,ib) + alpha*z(j)*rca(2,j)
                                idmzr(jt,i,ib) = idmzr(jt,i,ib) + alpha*z(j)*rca(3,j)
                                idmxr(jt,i,ib) = idmxr(jt,i,ib) + fac1*rca(1,j+1)
                                idmyr(jt,i,ib) = idmyr(jt,i,ib) + fac1*rca(2,j+1)
                                idmzr(jt,i,ib) = idmzr(jt,i,ib) + fac1*rca(3,j+1)
                                idmxr(jt,i,ib) = idmxr(jt,i,ib) + fac2*rca(1,j+2)
                                idmyr(jt,i,ib) = idmyr(jt,i,ib) + fac2*rca(2,j+2)
                                idmzr(jt,i,ib) = idmzr(jt,i,ib) + fac2*rca(3,j+2)

                                idmxv(jt,i,ib) = idmxv(jt,i,ib) + alpha*z(j)*pca(1,j)*w1
                                idmyv(jt,i,ib) = idmyv(jt,i,ib) + alpha*z(j)*pca(2,j)*w1
                                idmzv(jt,i,ib) = idmzv(jt,i,ib) + alpha*z(j)*pca(3,j)*w1
                                idmxv(jt,i,ib) = idmxv(jt,i,ib) + fac1 * pca(1,j+1) * w2
                                idmyv(jt,i,ib) = idmyv(jt,i,ib) + fac1 * pca(2,j+1) * w2
                                idmzv(jt,i,ib) = idmzv(jt,i,ib) + fac1 * pca(3,j+1) * w2
                                idmxv(jt,i,ib) = idmxv(jt,i,ib) + fac2 * pca(1,j+2) * w3
                                idmyv(jt,i,ib) = idmyv(jt,i,ib) + fac2 * pca(2,j+2) * w3
                                idmzv(jt,i,ib) = idmzv(jt,i,ib) + fac2 * pca(3,j+2) * w3
                            end if
                        enddo
                    enddo
                enddo
            endif
        end if
#ifdef PARALLEL_BINDING
	endif
#endif
    enddo
#ifdef PARALLEL_BINDING
						if(myid.eq.0) then
#endif
    ! Average accumulators

    wt = 1.d0/dble(2*nt)
    tv = tv*wt
    pe = pe*wt
    dqq = dqq*wt
    davx = davx*wt
    davy = davy*wt
    davz = davz*wt

    wt=1.d0/dble((2*nt/iskip+1)*nm)
    tavee = wt*tavee
    tq1aee = wt*tq1aee
    tq2aee = wt*tq2aee
    v1eeav = wt*v1eeav
    v2eeav = wt*v2eeav
    v3eeav = wt*v3eeav

    if (itst(2).eq.1) then
        if ((rpmddft.eq.0).and.((nbdf1.gt.0).or.(nbdf2.gt.0))) then

            write(eefile,*)'<V> per molecule = ',toKjmol* &
            (v1eeav+v2eeav+v3eeav),' KJ mol^-1'
            write(eefile,*)'<V>_ew = ', toKjmol*v1eeav,' KJ mol^-1'
            write(eefile,*)'<V>_lj = ', toKjmol*v2eeav,' KJ mol^-1'
            write(eefile,*)'<V>_inter = ', toKjmol*(v1eeav+v2eeav), &
            ' KJ mol^-1'
            write(eefile,*)'<V>_intra = ', toKjmol*v3eeav,' KJ mol^-1'
            write(eefile,*)
            write(eefile,*)'<Virial KE> per molecule = ', toKjmol*tavee, &
            ' KJ mol^-1'
            write(eefile,*)'<Virial KE>_inter = ', toKjmol*tq1aee, &
            ' KJ mol^-1'
            write(eefile,*)'<Virial KE>_intra = ', toKjmol*tq2aee, &
            ' KJ mol^-1'
        endif
    endif

    ! Velocity auto-correlation function.
    ! -----------------------------------

    ! Convert momenta to velocities

    emi = 1.d0/em
    if (itcf(1).eq.1) then
        pc(:,:,:) = pc(:,:,:)*emi
    endif

    ! Calculate correlation function

    wt = 1.d0/((nm-1)*(nt+1))
    if (itcf(1).eq.1) then
        do jt = 0,nt
            do it = 0,nt
                do j = 1,nm
                    do i = 1,3
                        ct(jt) = ct(jt) + pc(i,j,it)*pc(i,j,it+jt)
                    enddo
                enddo
            enddo
            ct(jt) = wt*ct(jt)
        enddo
    endif
      
    ! Dipole auto-correlation functions.
    ! -----------------------------------

    wt1 = 1.d0 / dble((nt/iskip)+1)
    if (itcf(2).eq.1) then
        do jt = 0,nt, iskip
            do it = 0, nt, iskip
                dmtr(jt) = dmtr(jt) + (dmxr(it)*dmxr(it+jt)+ &
                dmyr(it)*dmyr(it+jt)+dmzr(it)*dmzr(it+jt))
                dmtv(jt) = dmtv(jt) + (dmxv(it)*dmxv(it+jt)+ &
                dmyv(it)*dmyv(it+jt)+dmzv(it)*dmzv(it+jt))
            enddo
            dmtr(jt) = wt1 * dmtr(jt)
            dmtv(jt) = wt1 * dmtv(jt)
        enddo
        if (vacfac.ne.1) then
            do ib=1,2
                do i=1,4
                    do jt = 0,nt,iskip
                        do it = 0,nt,iskip
                            idmtr(jt,i,ib) = idmtr(jt,i,ib) + (idmxr(it,i,ib)*idmxr(it+jt,i,ib)+ &
                            idmyr(it,i,ib)*idmyr(it+jt,i,ib)+idmzr(it,i,ib)*idmzr(it+jt,i,ib))
                            idmtv(jt,i,ib) = idmtv(jt,i,ib) + (idmxv(it,i,ib)*idmxv(it+jt,i,ib)+ &
                            idmyv(it,i,ib)*idmyv(it+jt,i,ib)+idmzv(it,i,ib)*idmzv(it+jt,i,ib))
                        enddo
                        idmtr(jt,i,ib) = wt1 * idmtr(jt,i,ib)
                        idmtv(jt,i,ib) = wt1 * idmtv(jt,i,ib)
                    enddo
                enddo
            enddo
        end if
    endif

    ! Orientational correlation functions.
    ! -------------------------------------
    ! First order Legendre polynomial in dmot(1-3,:).
    ! Second order Legendre polynomial in dmot(4-6,:).

    if (itcf(3) .eq. 1) then
        do jt = 0, nt, iskip
            jt10 = jt/iskip
            tbar = 0.d0
            do it = 0, nt, iskip
                it10 = it/iskip
                do j = 1, nm
                    dotx = 0.d0
                    doty = 0.d0
                    dotz = 0.d0
                    do i = 1,3
                        dotx = dotx + ax1(i,j,it10) * ax1(i,j,it10+jt10)
                        doty = doty + ay1(i,j,it10) * ay1(i,j,it10+jt10)
                        dotz = dotz + az1(i,j,it10) * az1(i,j,it10+jt10)
                    enddo
                    dmot1(1,jt) = dmot1(1,jt) + dotx
                    dmot1(2,jt) = dmot1(2,jt) + doty
                    dmot1(3,jt) = dmot1(3,jt) + dotz
              
                    dotx = 0.d0
                    doty = 0.d0
                    dotz = 0.d0
                    do i = 1,6
                        dotx = dotx + ax2(i,j,it10) * ax2(i,j,it10+jt10)
                        doty = doty + ay2(i,j,it10) * ay2(i,j,it10+jt10)
                        dotz = dotz + az2(i,j,it10) * az2(i,j,it10+jt10)
                    enddo
                    dmot1(4,jt) = dmot1(4,jt) + dotx
                    dmot1(5,jt) = dmot1(5,jt) + doty
                    dmot1(6,jt) = dmot1(6,jt) + dotz
                enddo
                tbar = tbar + 1.d0
            enddo
            do k = 1, 6
                dmot1(k,jt) = dmot1(k,jt) / (tbar * nm)
            enddo
        enddo
    endif

#ifdef PARALLEL_BINDING
	endif
#endif
    ! Deallocate arrays

    deallocate(pc,pca,rca,dvdre,dist)

    if (itcf(2).eq.1) then
        deallocate(dmxr,dmyr,dmzr,dmxv,dmyv,dmzv)
        deallocate(idmxr,idmyr,idmzr,idmxv,idmyv,idmzv)
    end if
    if (itcf(3).eq.1) then
        deallocate(ax1,ay1,az1,ax2,ay2,az2)
    end if

    return
end subroutine trajectory
