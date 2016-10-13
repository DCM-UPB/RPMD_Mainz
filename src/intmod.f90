module intmod
    implicit none
    include '2DPMF.inc'
    !!This module stores output data arrays concerning the instantaneous water interface (see Willard and Chandler)
    !!and hydrogen bonding in water. It contains routines to update those data arrays with
    !!trajectory data (on-the-fly from MD or data stream from file) and calculate means and
    !!correlation functions, when the trajectory is finished. A second storage module is suggested
    !!to average over several trajectories (e.g. trajavmod.f90 ).
    !!A input file named 'intinput' specifies the analysis configuration, while the initialization subroutine 'ininit'
    !!is called with the number of atoms, the box, and timestep and total number of 'intupdate' calls.
    real(8) au,pi
    !---input---
    integer cstep
    real(8) cgxi,dtpfs,gridspac,zspac,jcspac
    real(8) lprox(3),uprox(3),pmfdist(56),pmfcos(56)
    !--configuration--
    real(8) box(3),vbox,boxcenterz,cut,cutsq,cgxisq,dfac,expfac,exparg,maxrang,hbdist(3),hbcos(3)
    integer wna,wnm,totalsteps,wnt,wntc,xyzpts(3),ntpts,intpts,zpts,jcpts,yzpts,xypts,hbcorrlen
    real(8), allocatable :: grid(:,:), rmct(:,:,:)
    !--density arrays--
    real(8), allocatable :: meandfield(:),densnz(:),densnzmean(:)

    !--interface arrays--
    real(8), allocatable :: interface(:,:,:),meannvec(:,:,:),jcdist(:)
    real(8) meanint(2),rhobulk

    !--hbond arrays--
    real(8), allocatable :: intrabond(:,:),interbond(:,:),inintct(:,:,:),hbct(:,:),avbno(:,:)
    real(8), allocatable :: ddct(:,:),sdct(:,:),ndct(:,:),acct(:,:)
    logical, allocatable :: inint(:,:,:,:),hbonded(:,:,:,:)
    real(8)              :: hbgt1(4)

    !--correlation and post-loop arrays--
    real(8), allocatable :: intcorr(:),intmean(:,:),intsigma(:,:)
    real(8), allocatable :: inintcorr(:,:),bondcorr(:,:,:),restbcorr(:,:,:),brokrepcorr(:,:,:),bcrcorr(:,:,:)


    namelist/intinput/ cgxi,cstep,gridspac,zspac,jcspac,lprox,uprox



contains

    subroutine intinit (na,boxlxyz,dtfs,trajsteps)              !initialize configurational variables and data arrays
        implicit none

        integer na,trajsteps
        real(8) boxlxyz(3),dtfs
        integer iter,densi,densx,densy,densz

        wna=na
        wnm=wna/3
        box=boxlxyz
        dtpfs=dtfs
        totalsteps=trajsteps

        au=1.6605389d0        !atomic mass unit
        pi=2.d0*dacos(0.d0)

        !hbond definition

        hbdist(1)=3.1d0
        hbdist(2)=3.5d0
        hbdist(3)=4.0d0

        hbcos(1)=-Cos(160.d0/180.d0*pi)
        hbcos(2)=-Cos(140.d0/180.d0*pi)
        hbcos(3)=-Cos(130.d0/180.d0*pi)

        pmfdist(:)=pmfralph(1,:)
        pmfcos(:)=Cos(pmfralph(2,:)/180.d0*pi)

        !---data input---

        open(6234,file='intinput')
        read(6234,intinput)
        close (unit=6234)
        cgxisq=cgxi**2                      !for faster evaluation in critical loop later
        exparg=-0.5d0/cgxisq
        cut=4.d0*cgxi
        cutsq=cut**2

        !---initialization---
        if ((box(1)-int(box(1)).eq.0.d0).or.(box(2)-int(box(2)).eq.0.d0).or. &
        (box(3)-int(box(3)).eq.0.d0)) then
            write(6,*) 'Boxlengths which equal integer values result in divison by 0 due' &
            //'to the periodic grid coice. Stop.'
            stop
        end if
        vbox=box(1)*box(2)*box(3)
        boxcenterz=box(3)*0.5d0
        do iter=1,3
            xyzpts(iter)=int(box(iter)/gridspac)+1   !positions go from 0 to the int before box
        enddo
        yzpts=xyzpts(2)*xyzpts(3)
        xypts=xyzpts(1)*xyzpts(2)
        ntpts=xyzpts(1)*xyzpts(2)*xyzpts(3)     !number of total points per time step
        intpts=(totalsteps-1)/cstep+1
        maxrang=0.5d0*sqrt(box(1)**2+box(2)**2+box(3)**2) !greatest possible box distance to surface
                                                                      !with minimum image convention applied
        zpts=int(2.d0*maxrang/zspac)+1
        jcpts=int(2.d0/jcspac)+1
        hbcorrlen=int(20000.d0/dtpfs)   !calculate hbond correlation for 20 ps
        if (hbcorrlen.gt.totalsteps/2) then
            hbcorrlen=totalsteps/2
        end if
        dfac=18.d0*au/(box(1)*box(2)*zspac)  !normalizes dprofiles to g/cm^3
        expfac=(18.d0*au)*(2.d0*pi*cgxisq)**(-1.5d0)  !for coarse graining gaussian


        write(6,*) 'Total Time Steps:',totalsteps             !number of time steps in traj
        write(6,*) 'Density Calculation Steps:',intpts        !number of density and interface evaluation steps
        write(6,*) 'HB-Corr. Length in time steps:',hbcorrlen !Number of hb-correlation time steps
        write(6,*)
                                                                                                 !array allocation
        allocate (rmct(3,wnm,intpts),grid(3,ntpts),meannvec(3,xypts,2), &
        jcdist(jcpts*jcpts*zpts),densnz(zpts),densnzmean(zpts), &
        meandfield(ntpts), &
        intcorr(intpts),interface(xypts,2,intpts),hbct(4,4),avbno(4,4), &
        intmean(xypts,2),intsigma(xypts,2),hbonded(4,3,wnm,totalsteps),intrabond(4,3),interbond(4,3), &
        inint(wnm,4,2,intpts),inintct(4,2,intpts),inintcorr(4,hbcorrlen),bondcorr(4,4,hbcorrlen), &
        restbcorr(4,4,hbcorrlen),brokrepcorr(4,4,hbcorrlen),bcrcorr(4,4,hbcorrlen), &
        ddct(4,4),sdct(4,4),ndct(4,4),acct(4,4))

        !sum up arrays go zero

        call intreset()

        densi=0
        do densx=1,xyzpts(1)                    !grid initialization
            do densy=1,xyzpts(2)
                do densz=1,xyzpts(3)
                    densi=densi+1
                    grid(1,densi)=gridspac*dble(densx-1)     !densx=1 represents x=0 and so on
                    grid(2,densi)=gridspac*dble(densy-1)
                    grid(3,densi)=gridspac*dble(densz-1)
                enddo
            enddo
        enddo

    end subroutine

    subroutine intreset()
        implicit none

        wnt=0
        wntc=0
        meannvec(:,:,:)=0.d0
        jcdist(:)=0.d0
        meandfield(:)=0.d0
        densnz(:)=0.d0
        densnzmean(:)=0.d0

        intrabond(:,:)=0.d0
        interbond(:,:)=0.d0
        hbonded(:,:,:,:)=.false.
        inint(:,:,:,:)=.false.
        inintct(:,:,:)=0.d0
        hbct(:,:)=0.d0
        avbno(:,:)=0.d0
        ddct(:,:)=0.d0
        sdct(:,:)=0.d0
        ndct(:,:)=0.d0
        acct(:,:)=0.d0
        hbgt1(:)=0.d0

        intmean(:,:)=0.d0
        intsigma(:,:)=0.d0
        meanint(:)=0.d0
        rhobulk=0.d0

        intcorr(:)=0.d0
        inintcorr(:,:)=0.d0
        bondcorr(:,:,:)=0.d0
        restbcorr(:,:,:)=0.d0
        brokrepcorr(:,:,:)=0.d0
        bcrcorr(:,:,:)=0.d0
    end subroutine

    subroutine intupdate(ra)            !prepare coordinate data, if cstep update dfield and int data, always update hb data
        implicit none
        real(8) ra(3,wna)

        real(8) rct(3,wnm),rot(3,wnm),bondvec(3,2,wnm),instdfield(ntpts)

        integer im,ia,ib,iv, int12,iter,inti,densx,densy,densz,densi,densi2,modz,modz2,maxdensz
        real(8) distvec(3),dist,rhomax,hbulk,hvec(3,2),boxsh,boxsh2,cosangl(2)
        real(8) hprox(xypts),intnvec(3,xypts),proximity(wnm)
        integer cthb(4,2)
        logical foundhb,hbflag(4)


        wnt=wnt+1       !counts internal calculation step number

        rct(:,:)=0.d0

        hbulk=0.d0
        instdfield(:)=0.d0

        ia=0
        do im=1,wnm          !preparation of coordinate data
            ia=ia+1
            rot(:,im)=ra(:,ia)                      !oxygen data
            rct(:,im)=rct(:,im)+16.d0*rot(:,im)

            do ib=1,2                               !hydrogen data
                ia=ia+1
                rct(:,im)=rct(:,im)+ra(:,ia)

                bondvec(:,ib,im)=ra(:,ia)-rot(:,im)     !normalized bond vectors
                bondvec(:,ib,im)=bondvec(:,ib,im)/sqrt(bondvec(1,ib,im)**2+bondvec(2,ib,im)**2+bondvec(3,ib,im)**2)
            enddo

            rct(:,im)=rct(:,im)/18.d0
        enddo

        do iv=1,3       !apply period boundary conditions to COMs and O
            rct(iv,:)=rct(iv,:)-box(iv)*anint(rct(iv,:)/box(iv)-0.5d0)
            rot(iv,:)=rot(iv,:)-box(iv)*anint(rot(iv,:)/box(iv)-0.5d0)
        enddo

        if (mod(wnt-1,cstep).eq.0) then

            wntc=wntc+1 !density calculation step number

            rmct(:,:,wntc)=rct(:,:)  !saves COM trajectory for certain mean interface calculations

                                                !density field calculation routine (consumes much time)
            do densi=1,ntpts                    !approx 100 million loop steps for 11*11*11 molecules and VF 3
                do im=1,wnm
                    distvec(:)=rct(:,im)-grid(:,densi)
                    distvec(:)=distvec(:)-box(:)*anint(distvec(:)/box(:))       !minimum image convention

                    if (abs(distvec(3)).lt.cut) then    !faster for large z-direction boxes (as vacuum boxes are)
                        dist=dot_product(distvec,distvec)
                        if (dist .lt. cutsq) then       !cutoff 4 cgxi -> 16 cgxisq

                            !coarse graining with gaussian (2*pi*cgxi^2)^(-3/2) * Exp[-((r-rm)/cgxi)^2 / 2]
                            instdfield(densi)=instdfield(densi)+exp(exparg*dist)

                        end if
                    end if
                enddo
            enddo

            instdfield(:)=expfac*instdfield(:)

            !interface calculation routine

            hbulk=0.d0
            iv=0

            do inti=1,xypts         !calculate the density in bulk for next step
                maxdensz=MAXLOC(instdfield((inti-1)*xyzpts(3)+1:inti*xyzpts(3)),1)
                rhomax=0.9d0*MAXVAL(instdfield((inti-1)*xyzpts(3)+1:inti*xyzpts(3)))
                                                                  ! find average density in
                densi=1                                           ! the 90% of max zone
                densi2=1 !to avoid compiler warnings
                do densz=maxdensz-xyzpts(3)/2,maxdensz-1
                    modz=modulo(densz,xyzpts(3))+1
                    densi=(inti-1)*xyzpts(3)+modz
                    if (instdfield(densi) .gt. rhomax) then
                        exit
                    end if
                enddo

                do densz=-(maxdensz+xyzpts(3)/2),-(maxdensz+1)
                    modz2=modulo(-densz,xyzpts(3))+1
                    densi2=(inti-1)*xyzpts(3)+modz2
                    if (instdfield(densi2) .gt. rhomax) then
                        exit
                    end if
                enddo

                if (densi2.gt.densi) then                           !sum bulk densities
                    hbulk=hbulk+sum(instdfield(densi:densi2))
                    iv=iv+densi2-densi+1
                else if (densi2.lt.densi) then
                    hbulk=hbulk+sum(instdfield(densi2:densi))
                    iv=iv+densi-densi2+1
                else
                    write(6,*) 'Warning: densi2 and densi in bulk density calculation are equal!'
                end if
            enddo

            hbulk=hbulk/dble(iv)
            rhobulk=rhobulk+hbulk
            hbulk=0.5d0*hbulk           !half the bulk density

            if (maxval(instdfield).gt.hbulk) then
                do inti=1,xypts

                    !linear interpolation between the z points including the interface
                    !y(z):= [(y2-y1)/(z2-z1)] * (z-z1) + y1 && y(z0)==hbulk  -->  z0=z1+(z1-z2)/(y1-y2)*(hbulk-y1)

                    maxdensz=MAXLOC(instdfield((inti-1)*xyzpts(3)+1:inti*xyzpts(3)),1)

                    do densz=maxdensz-xyzpts(3)/2,maxdensz-1
                        modz=modulo(densz-1,xyzpts(3))+1
                        modz2=modulo(densz,xyzpts(3))+1

                        densi=(inti-1)*xyzpts(3)+modz                       !interface 1 ("left")
                        densi2=(inti-1)*xyzpts(3)+modz2

                        if (instdfield(densi2) .gt. hbulk) then
                            if (modz.gt.maxdensz) then
                                boxsh=grid(3,densi)-box(3)
                            else
                                boxsh=grid(3,densi)
                            end if
                            if (modz2.gt.maxdensz) then
                                boxsh2=grid(3,densi2)-box(3)
                            else
                                boxsh2=grid(3,densi2)
                            end if

                            !if you replace gridspac with z2-z1 you get the interpolation formula from above
                            interface(inti,1,wntc)=boxsh+(boxsh2-boxsh)*(hbulk-instdfield(densi)) / &
                            (instdfield(densi2)-instdfield(densi))
                            intmean(inti,1)=intmean(inti,1) + interface(inti,1,wntc)
                            exit
                        end if
                    enddo

                    do densz=-(maxdensz+xyzpts(3)/2),-(maxdensz+1)
                        modz=modulo(-densz-1,xyzpts(3))+1
                        modz2=modulo(-densz,xyzpts(3))+1

                        densi=(inti-1)*xyzpts(3)+modz                       !interface 2 ("right"),inversed direction
                        densi2=(inti-1)*xyzpts(3)+modz2

                        if (instdfield(densi) .gt. hbulk) then
                            if (modz.lt.maxdensz) then
                                boxsh=grid(3,densi)+box(3)
                            else
                                boxsh=grid(3,densi)
                            end if
                            if (modz2.lt.maxdensz) then
                                boxsh2=grid(3,densi2)+box(3)
                            else
                                boxsh2=grid(3,densi2)
                            end if

                            interface(inti,2,wntc)=boxsh+(boxsh2-boxsh)*(hbulk-instdfield(densi)) / &
                            (instdfield(densi2)-instdfield(densi))
                            intmean(inti,2)=intmean(inti,2) + interface(inti,2,wntc)
                            exit
                        end if
                    enddo
                enddo

                                !The zero components of this plane vectors are known
                hvec(2,1)=0.d0
                hvec(1,2)=0.d0

                !calculate interfacial properties

                do int12=1,2
                                    !calculate interface normal vector field
                    do inti=1,xypts
                        !For surface normal vector calculation the interface is taken in every point (x0,y0) seperately
                        !as a regular(!) plain parametric surface P:(x,y)->(x,y,mx*x + my*y [+c]),
                        !where the slopes mx and my are calculated by interface height values at points
                        !(x1,y0)&(x2,y0) and (x0,y1)&(x0,y2),where x1,x2,y1,y2 are the next to x0,y0.
                        !This is an approximation to the tangential plane in (x0,y0)
                        !and therefore suitable for normal vector calculation.
                        !For a regular surface one can obtain the surface normal vector through: Fx= D(P,x),Fy=D(P,y),
                        !N=(Fx x Fy)/Norm(Fx x Fy)

                        densi=modulo(inti-xyzpts(2)-1,xypts)+1
                        densi2=modulo(inti+xyzpts(2)-1,xypts)+1                                 !Fx
                        dist=grid(1,densi2*xyzpts(3))-grid(1,densi*xyzpts(3))
                        dist=dist-box(1)*anint(dist/box(1))     !needed for border case
                        hvec(1,1)=dist
                        hvec(3,1)=(interface(densi2,int12,wntc)-interface(densi,int12,wntc))/dist
                        densi=((inti-1)/xyzpts(2))*xyzpts(2)+modulo(inti-2,xyzpts(2))+1
                        densi2=((inti-1)/xyzpts(2))*xyzpts(2)+modulo(inti,xyzpts(2))+1          !Fy
                        dist=grid(2,densi2*xyzpts(3))-grid(2,densi*xyzpts(3))
                        dist=dist-box(2)*anint(dist/box(2))
                        hvec(2,2)=dist
                        hvec(3,2)=(interface(densi2,int12,wntc)-interface(densi,int12,wntc))/dist

                        intnvec(:,inti)=crossp(hvec(:,mod(int12,2)+1),hvec(:,mod(int12-1,2)+1))

                        intnvec(:,inti)=intnvec(:,inti)/sqrt(intnvec(1,inti)**2 + &
                        intnvec(2,inti)**2 + intnvec(3,inti)**2)

                    enddo

                    do im=1,wnm    !now calculate properties depending on normal vectors

                        do inti=1,xypts             !find closest surface point
                            distvec(1)=grid(1,inti*xyzpts(3))-rct(1,im)
                            distvec(2)=grid(2,inti*xyzpts(3))-rct(2,im)
                            distvec(3)=interface(inti,int12,wntc)-rct(3,im)
                            distvec(:)=distvec(:)-box(:)*anint(distvec(:)/box(:)) !min img conv
                            hprox(inti)=dot_product(distvec,distvec)
                        enddo

                        inti=MINLOC(hprox,1)        !here

                        distvec(1)=grid(1,inti*xyzpts(3))-rct(1,im)
                        distvec(2)=grid(2,inti*xyzpts(3))-rct(2,im)
                        distvec(3)=interface(inti,int12,wntc)-rct(3,im)
                        distvec(:)=distvec(:)-box(:)*anint(distvec(:)/box(:))

                        proximity(im)=distvec(1)*intnvec(1,inti)+distvec(2)*intnvec(2,inti)+distvec(3)*intnvec(3,inti)

                        densz=int((proximity(im)+maxrang)/zspac)+1
                        densnz(densz)=densnz(densz) + 1.d0

                        !joint conditional distribution calculation

                        densi=int((bondvec(1,1,im)*intnvec(1,inti)+bondvec(2,1,im)*intnvec(2,inti)+ &     ! (a.b)/|a|/|b|=Cos(angle)
                        bondvec(3,1,im)*intnvec(3,inti)+1.d0)/jcspac)*jcpts + &
                        int((bondvec(1,2,im)*intnvec(1,inti)+bondvec(2,2,im)*intnvec(2,inti)+ &
                        bondvec(3,2,im)*intnvec(3,inti)+1.d0)/jcspac)

                        densi2=int((proximity(im)+maxrang)/zspac)*jcpts*jcpts + 1
                        densi=densi+densi2

                        jcdist(densi)=jcdist(densi) + 1.d0

                        if (proximity(im).gt.lprox(1) .and. proximity(im).lt.uprox(1)) then
                            inint(im,1,int12,wntc)=.true.
                            inintct(1,int12,wntc)=inintct(1,int12,wntc)+1.d0
                        end if

                        if (proximity(im).gt.lprox(2) .and. proximity(im).lt.uprox(2)) then
                            inint(im,2,int12,wntc)=.true.
                            inintct(2,int12,wntc)=inintct(2,int12,wntc)+1.d0
                        end if

                        if (proximity(im).gt.lprox(3) .and. proximity(im).lt.uprox(3)) then
                            inint(im,3,int12,wntc)=.true.
                            inintct(3,int12,wntc)=inintct(3,int12,wntc)+1.d0
                        end if

                    enddo

                enddo

                inint(:,4,1,wntc) = (.not.inint(:,1,1,wntc)) .and. (.not.inint(:,2,1,wntc)) .and. &
                (.not.inint(:,3,1,wntc))
                inint(:,4,2,wntc) = (.not.inint(:,1,2,wntc)) .and. (.not.inint(:,2,2,wntc)) .and. &
                (.not.inint(:,3,2,wntc))
                inint(:,4,1,wntc) = inint(:,4,1,wntc) .and. inint(:,4,2,wntc)
                inint(:,4,2,wntc) = inint(:,4,1,wntc)
                inintct(4,1,wntc)=wnm-sum(inintct(1:3,1,wntc))-sum(inintct(1:3,2,wntc))
                inintct(4,2,wntc)=inintct(4,1,wntc)
            else
                write(6,*) 'Warning: maxval(instdfield).lt.hbulk !'
            end if

            meandfield(:)=meandfield(:) + instdfield(:)
            instdfield(:)=0.d0

            if (((totalsteps-1)/cstep+1)/10.gt.0) then
                if (wnt.ne.0 .and. mod((wnt-1)/cstep,((totalsteps-1)/cstep+1)/10).eq.0) then
                    write(6,*) 'IntStep',wntc,'of',totalsteps/cstep+1,'completed.'
                end if
            end if

            !---INTERFACE CALCULATION STEP COMPLETED---
        end if

        do im=1,wnm                                            !hydrogen bonding based on fixed definitions (1-3) or 2d-pmf (4)
            cthb(:,:)=0
            do iter=1,wnm
                if (iter.ne.im) then

                    distvec(:)=rot(:,iter)-rot(:,im)
                    distvec(:)=distvec(:)-box(:)*anint(distvec(:)/box(:))
                    if (abs(distvec(3)).lt.hbdist(3)) then    !faster for large z slabs (as waterslabs with surface often are)
                        dist=sqrt(dot_product(distvec,distvec))
                        if (dist.lt.hbdist(3)) then                     !hbdist(3) is weaker than all pmf distances !!
                            cosangl(1)=(bondvec(1,1,im)*distvec(1)+bondvec(2,1,im)*distvec(2)+ &     ! (a.b)/|a|/|b|=Cos(angle)
                            bondvec(3,1,im)*distvec(3))/dist
                            cosangl(2)=(bondvec(1,2,im)*distvec(1)+bondvec(2,2,im)*distvec(2)+ &
                            bondvec(3,2,im)*distvec(3))/dist

                                          !check fixed definition hbonds
                            do densi=1,2
                                hbflag(:)=.false. !says whether the current im,iter,densi,densx config is a bond

                                if (cosangl(densi).gt.hbcos(3)) then    !hbcos(3) is weaker than all pmf angles !!
                                    cthb(3,densi)=cthb(3,densi)+1
                                    hbflag(3)=.true.
                                    hbonded(3,densi,im,wnt)=.true.
                                    hbonded(3,3,iter,wnt)=.true.        !3 is the index for accepting bond
                                    if (dist.lt.hbdist(2)) then
                                        if (cosangl(densi).gt.hbcos(2)) then
                                            cthb(2,densi)=cthb(2,densi)+1
                                            hbflag(2)=.true.
                                            hbonded(2,densi,im,wnt)=.true.
                                            hbonded(2,3,iter,wnt)=.true.
                                            if (dist.lt.hbdist(1)) then
                                                if (cosangl(densi).gt.hbcos(1)) then
                                                    cthb(1,densi)=cthb(1,densi)+1
                                                    hbflag(1)=.true.
                                                    hbonded(1,densi,im,wnt)=.true.
                                                    hbonded(1,3,iter,wnt)=.true.
                                                end if
                                            end if
                                        end if
                                    end if

                                        !check pmf hbonds

                                    densi2=int((dist-pmfdist(1))/pmfdelt)+1

                                    if (densi2.lt.pmflen) then !!.lt.
                                        if (densi2.lt.1) then

                                            if (cosangl(densi).gt.pmfcos(1)) then
                                                cthb(4,densi)=cthb(4,densi)+1
                                                hbflag(4)=.true.
                                                hbonded(4,densi,im,wnt)=.true.
                                                hbonded(4,3,iter,wnt)=.true.
                                            end if

                                        else
                                            boxsh=(dist-pmfdist(densi2))/pmfdelt    !is positive and between 0 and 1
                                                                                    !due to densi2 calculation
                                            if (cosangl(densi).gt.((1.d0-boxsh)*pmfcos(densi2)+ &
                                            boxsh*pmfcos(densi2+1))) then     !lin. interpolation
                                                                                !of PMF definition
                                                cthb(4,densi)=cthb(4,densi)+1
                                                hbflag(4)=.true.
                                                hbonded(4,densi,im,wnt)=.true.
                                                hbonded(4,3,iter,wnt)=.true.

                                            end if
                                        end if
                                    end if

                                        !check directionality of hbonds in terms of layers

                                    do densx=1,4
                                        if(hbflag(densx)) then
                                            foundhb=.false.
                                            do int12=1,2
                                                do inti=1,4
                                                    if (inint(im,inti,int12,wntc)) then
                                                        do densy=1,2
                                                            do densz=1,4
                                                                if(inint(iter,densz,densy,wntc)) then


                                                                    if (inti.ne.4) then
                                                                        if (inti.eq.densz) then
                                                                            intrabond(densx,inti)= &
                                                                            intrabond(densx,inti)+1.d0
                                                                        else
                                                                            interbond(densx,inti)= &
                                                                            interbond(densx,inti)+1.d0
                                                                        end if
                                                                    end if
                                                                    avbno(densx,inti)= &
                                                                    avbno(densx,inti)+1.d0
                                                                    avbno(densx,densz)= &
                                                                    avbno(densx,densz)+1.d0


                                                                    foundhb=.true.

                                                                    exit
                                                                end if
                                                            enddo
                                                            if(foundhb) then    !needed because layer 4 is bulk
                                                                exit            !and therefore int12=1,2 are same
                                                            end if
                                                        enddo
                                                    end if
                                                enddo
                                                if(foundhb) then    !needed because layer 4 is bulk
                                                    exit            !and therefore int12=1,2 are same
                                                end if              !-> no double counting
                                            enddo

                                        end if
                                    enddo



                                end if
                            enddo
                        end if
                    end if
                end if
            enddo

            !check how often an oh donates more than 1 hbond per definition

            do densi=1,2
                !write(6,*) cthb(1,densi),cthb(2,densi),cthb(3,densi)
                do densx=1,4
                    if (cthb(densx,densi).gt.1) then
                        hbgt1(densx)=hbgt1(densx)+1
                    end if
                enddo
            enddo
        enddo

    end subroutine




    subroutine  meancalc ()
        implicit none

        integer im, iter,int12,inti,densx,densy,densz,densi,densi2
        real(8) distvec(3),dist,hvec(3,2),sumint(2)
        real(8) hproxm(xypts),proxm(wnm),ffthelp(intpts),intcorrh(intpts)
        logical restcorr(4,2)

        if (wnt.ne.totalsteps) then
            write(6,*) 'Warning: wnt!=totalsteps!'
        end if
        if (wntc.ne.intpts) then
            write(6,*) 'Warning: wntc!=intpts!'
        end if

        rhobulk=rhobulk/intpts !get back to full density
        meandfield(:)=meandfield(:)/intpts      !only evaluated for wnt%ctstep=0
        densnz(:)=0.5d0*dfac*densnz(:)/intpts
        write(6,*) rhobulk
        hvec(2,1)=0.d0
        hvec(1,2)=0.d0

        do densz=1,zpts
            do inti=(densz-1)*jcpts*jcpts+1,densz*jcpts*jcpts
                if (densnz(densz).ne.0.d0) then
                    jcdist(inti)=jcdist(inti)/intpts/densnz(densz)
                else
                    jcdist(inti)=0.d0
                end if
            enddo
        enddo

        do int12=1,2

            intmean(:,int12)=intmean(:,int12)/intpts
            meanint(int12)=sum(intmean(:,int12))/xypts

            do inti=1,xypts
                do iter=1,intpts
                    intsigma(inti,int12)=intsigma(inti,int12)+(interface(inti,int12,iter)-intmean(inti,int12))**2
                enddo
                intsigma(inti,int12)=sqrt(intsigma(inti,int12)/(intpts-1))

                !mean interface normal vector field

                densi=modulo(inti-xyzpts(2)-1,xypts)+1
                densi2=modulo(inti+xyzpts(2)-1,xypts)+1                                 !Fx
                dist=grid(1,densi2*xyzpts(3))-grid(1,densi*xyzpts(3))
                dist=dist-box(1)*anint(dist/box(1))     !needed for border case
                hvec(1,1)=dist
                hvec(3,1)=(intmean(densi2,int12)-intmean(densi,int12))/dist

                densi=((inti-1)/xyzpts(2))*xyzpts(2)+modulo(inti-2,xyzpts(2))+1
                densi2=((inti-1)/xyzpts(2))*xyzpts(2)+modulo(inti,xyzpts(2))+1          !Fy
                dist=grid(2,densi2*xyzpts(3))-grid(2,densi*xyzpts(3))
                dist=dist-box(2)*anint(dist/box(2))
                hvec(2,2)=dist
                hvec(3,2)=(intmean(densi2,int12)-intmean(densi,int12))/dist

                meannvec(:,inti,int12)=crossp(hvec(:,mod(int12,2)+1),hvec(:,mod(int12-1,2)+1))
                meannvec(:,inti,int12)=meannvec(:,inti,int12)/sqrt(meannvec(1,inti,int12)**2 + &
                meannvec(2,inti,int12)**2 + meannvec(3,inti,int12)**2)

            enddo
        enddo

        do iter=1,intpts
            do int12=1,2
                do im=1,wnm    !now calculate properties depending on normal vectors

                    do inti=1,xypts             !find closest surface point
                        distvec(1)=grid(1,inti*xyzpts(3))-rmct(1,im,iter)
                        distvec(2)=grid(2,inti*xyzpts(3))-rmct(2,im,iter)
                        distvec(3)=intmean(inti,int12)-rmct(3,im,iter)
                        distvec(:)=distvec(:)-box(:)*anint(distvec(:)/box(:)) !min img conv
                        hproxm(inti)=dot_product(distvec,distvec)
                    enddo

                    inti=MINLOC(hproxm,1)        !here

                    distvec(1)=grid(1,inti*xyzpts(3))-rmct(1,im,iter)
                    distvec(2)=grid(2,inti*xyzpts(3))-rmct(2,im,iter)
                    distvec(3)=intmean(inti,int12)-rmct(3,im,iter)
                    distvec(:)=distvec(:)-box(:)*anint(distvec(:)/box(:))
                    proxm(im)=distvec(1)*meannvec(1,inti,int12)+distvec(2)*meannvec(2,inti,int12)+distvec(3)*meannvec(3,inti,int12)

                    densz=int((proxm(im)+maxrang)/zspac)+1
                    densnzmean(densz)=densnzmean(densz) + 1.d0
                enddo
            enddo

        enddo

        densnzmean(:)=0.5d0*dfac*densnzmean(:)/intpts

        do int12=1,2
            do inti=1,xypts
                ffthelp(:)=interface(inti,int12,:)-intmean(inti,int12)
                call FastAutoCorr(ffthelp,intcorrh,intpts)
                intcorr(:)=intcorr(:)+intcorrh(:)
            enddo
        enddo
        intcorr(:)=intcorr(:)/intcorr(1)

        do densy=1,totalsteps-hbcorrlen                     !calculate hbond correlation functions
            do int12=1,2
                do iter=1,4
                    do im=1,wnm

                        if (inint(im,iter,int12,(densy-1)/cstep+1)) then

                            do densz=1,hbcorrlen
                                if (inint(im,iter,int12,(densy+densz-2)/cstep+1)) then
                                    inintcorr(iter,densz)=inintcorr(iter,densz)+1.d0
                                end if
                            enddo
                            restcorr(:,:)=.true.

                            do densi=1,2

                                do densx=1,4
                                    if (hbonded(densx,densi,im,densy)) then
                                        do densz=1,hbcorrlen
                                            if (hbonded(densx,densi,im,densy+densz-1)) then
                                                bondcorr(densx,iter,densz)=bondcorr(densx,iter,densz)+1.d0
                                            else
                                                restcorr(densx,densi)=.false.
                                            end if

                                            if (restcorr(densx,densi)) then
                                                restbcorr(densx,iter,densz)=restbcorr(densx,iter,densz)+1.d0
                                            end if

                                            if (.not.hbonded(densx,mod(densi,2)+1,im, &
                                            densy) .and. (.not.hbonded(densx,densi,im,densy+densz-1) .and. &
                                            hbonded(densx,mod(densi,2)+1,im,densy+densz-1))) then

                                                bcrcorr(densx,iter,densz)=bcrcorr(densx,iter,densz)+1.d0

                                            end if
                                        enddo
                                    else if (densy.gt.1) then
                                        if (hbonded(densx,densi,im,densy-1)) then
                                            do densz=1,hbcorrlen
                                                if (hbonded(densx,densi,im,densy+densz-1)) then
                                                    brokrepcorr(densx,iter,densz)=brokrepcorr(densx,iter,densz)+1.d0
                                                    exit
                                                end if
                                            enddo
                                        end if
                                    end if
                                enddo

                            enddo

                        end if
                    enddo
                enddo
            enddo
        enddo

        do densy=1,totalsteps
            do im=1,wnm
                do int12=1,2
                    do iter=1,4
                        do densx=1,4
                            do densi=1,2
                                if (hbonded(densx,densi,im,densy).and.inint(im,iter,int12,(densy-1)/cstep+1)) then
                                    hbct(densx,iter)=hbct(densx,iter)+1.d0          !count total number of bonded oh
                                end if
                            enddo
                            if (inint(im,iter,int12,(densy-1)/cstep+1)) then                    !count donor configurations

                                if (.not.hbonded(densx,1,im,densy) .and. .not.hbonded(densx,2,im,densy)) then
                                    ndct(densx,iter)=ndct(densx,iter)+1.d0
                                else if (XOR(hbonded(densx,1,im,densy),hbonded(densx,2,im,densy))) then
                                    sdct(densx,iter)=sdct(densx,iter)+1.d0
                                else if (hbonded(densx,1,im,densy) .and. hbonded(densx,2,im,densy)) then
                                    ddct(densx,iter)=ddct(densx,iter)+1.d0
                                end if

                                if (hbonded(densx,3,im,densy)) then
                                    acct(densx,iter)=acct(densx,iter)+1.d0
                                end if

                            end if
                        enddo
                    enddo
                enddo
            enddo
        enddo

        do iter=1,4                     !normalize counters
            sumint(:)=0.d0
            do densy=1,totalsteps
                sumint(:)=sumint(:)+inintct(iter,:,(densy-1)/cstep+1)
            enddo
            sumint(1)=sumint(1)+sumint(2)
            sumint(2)=sum(inintct(iter,1,:))
            sumint(2)=sumint(2)+sum(inintct(iter,2,:))
            inintcorr(iter,:)=inintcorr(iter,:)/inintcorr(iter,1)
            do densx=1,4

                bondcorr(densx,iter,:)=2.d0*bondcorr(densx,iter,:)/hbct(densx,iter)
                restbcorr(densx,iter,:)=2.d0*restbcorr(densx,iter,:)/hbct(densx,iter)
                bcrcorr(densx,iter,:)=2.d0*bcrcorr(densx,iter,:)/sdct(densx,iter)
                brokrepcorr(densx,iter,:)=2.d0*bondcorr(densx,iter,:)/hbct(densx,iter)
                avbno(densx,iter)=avbno(densx,iter)/sumint(1)
                if (iter.eq.4) then
                    avbno(densx,iter)=2.d0*avbno(densx,iter)
                else
                    interbond(densx,iter)=interbond(densx,iter)/hbct(densx,iter)
                    intrabond(densx,iter)=intrabond(densx,iter)/hbct(densx,iter)
                end if
                hbct(densx,iter)=0.5d0*hbct(densx,iter)/sumint(1)
                ndct(densx,iter)=ndct(densx,iter)/sumint(1)
                sdct(densx,iter)=sdct(densx,iter)/sumint(1)
                ddct(densx,iter)=ddct(densx,iter)/sumint(1)
                acct(densx,iter)=acct(densx,iter)/sumint(1)


            enddo
        enddo

    end subroutine

    subroutine FastAutoCorr(input,output,n)
    implicit none
    integer n
    double precision    :: input(n), output(n)

    double complex :: help(n/2+1)
    integer(8)          :: plan,plani

    output(:)=0.d0
    help(:)=(0.d0,0.d0)

    call dfftw_plan_dft_r2c_1d(plan,n,input,help,'FFTW_ESTIMATE')

    call dfftw_execute(plan)

    help=help*CONJG(help)

    call dfftw_plan_dft_c2r_1d(plani,n,help,output,'FFTW_ESTIMATE')

    call dfftw_execute(plani)

    call dfftw_destroy_plan(plan)

    call dfftw_destroy_plan(plani)

    return

end subroutine FastAutoCorr


    function crossp(vec1,vec2) result (cpvec)
        implicit none
        real(8) vec1(3),vec2(3),cpvec(3)

        cpvec(1)=vec1(2)*vec2(3)-vec1(3)*vec2(2)
        cpvec(2)=vec1(3)*vec2(1)-vec1(1)*vec2(3)
        cpvec(3)=vec1(1)*vec2(2)-vec1(2)*vec2(1)

    end function

end module intmod
