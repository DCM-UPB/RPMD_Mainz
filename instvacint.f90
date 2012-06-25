subroutine instvacint(ra_bohr,inint,na,boxlxyz_bohr)
    implicit none
    include 'globals.inc'

    !---passed variables--- Units in Bohr !!
    integer na
    real(8) boxlxyz_bohr(3)
    real(8) ra_bohr(3,na)
    logical inint(na/3,3,2)

    !--- Angström working arrays---
    real(8) boxlxyz(3),ra(3,na)

    !---    hard coded parameters   ---
    real(8) cgxi,gridspac,zspac,au,pi,hbulk
    real(8) lprox1,uprox1,lprox2,uprox2

    !---    general variables   ---
    integer ia,im,densx,densy,densz,densi,densi2,inti,iter,int12    !iterators
    integer xyzpts(3),ntpts,zpts,yzpts,xypts,maxdensz
    real(8), allocatable :: rct(:,:),grid(:,:)
    real(8) vbox,boxcenterz,watercomz

    !-- density field vars  --
    real(8), allocatable :: instdfield(:)
    real(8) distvec(3),dist,cut,cutsq,cgxisq,exparg

    !-- interface vars  --
    real(8), allocatable :: interface(:,:),intnvec(:,:)
    real(8) boxsh,boxsh2,maxrang,hvec1(3),hvec2(3)
    integer modz,modz2
    real(8), allocatable :: bondvec1(:,:),bondvec2(:,:),proximity(:),hprox(:)

    boxlxyz(:)=toA*boxlxyz_bohr(:)
    ra(:,:)=toA*ra_bohr(:,:)

    cgxi=2.4d0
    gridspac=1.d0
    zspac=0.1d0
    au=1.6605389d0
    pi=3.14159d0
    hbulk=0.5d0*0.997d0

    lprox1=-2.0d0
    uprox1=0.5d0
    lprox2=0.5d0
    uprox2=3.0d0

    cgxisq=cgxi**2                      !for faster evaluation in critical loop later
    exparg=-0.5d0/cgxisq
    cut=4.d0*cgxi
    cutsq=cut**2

    vbox=boxlxyz(1)*boxlxyz(2)*boxlxyz(3)
    boxcenterz=boxlxyz(3)*0.5d0
    do iter=1,3
        xyzpts(iter)=int(boxlxyz(iter)/gridspac)+1   !positions go from 0 to the int before boxlxyz
    enddo
    yzpts=xyzpts(2)*xyzpts(3)
    xypts=xyzpts(1)*xyzpts(2)
    ntpts=xyzpts(1)*xyzpts(2)*xyzpts(3)     !number of total points per time step
    maxrang=0.5d0*sqrt(boxlxyz(1)**2+boxlxyz(2)**2+boxlxyz(3)**2) !greatest possible box distance to surface
                                                                  !with minimum image convention applied
    zpts=int(2.d0*maxrang/zspac)+1 !zpts is misunderstandable.. it represents the range of proximity values

    allocate (rct(3,na/3),grid(3,ntpts),intnvec(3,xypts), &
    bondvec1(3,na/3),bondvec2(3,na/3),proximity(na/3),hprox(xypts), &
    instdfield(ntpts),interface(xypts,2))

    instdfield(:)=0.d0


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

    im=0
    do ia=1,na,3
        im=im+1
        rct(:,im)=(16.d0*ra(:,ia)+ra(:,ia+1)+ra(:,ia+2))/18.d0
        !bondvec1(:,im)=ra(:,ia+1)-ra(:,ia)
        !bondvec2(:,im)=ra(:,ia+2)-ra(:,ia)     !not in use for the moment

        rct(1,im)=rct(1,im)-boxlxyz(1)*anint(rct(1,im)/boxlxyz(1)-0.5d0)
        rct(2,im)=rct(2,im)-boxlxyz(2)*anint(rct(2,im)/boxlxyz(2)-0.5d0)
        !bondvec1(:,im)=bondvec1(:,im)/sqrt(bondvec1(1,im)**2+bondvec1(2,im)**2+bondvec1(3,im)**2)
        !bondvec2(:,im)=bondvec2(:,im)/sqrt(bondvec2(1,im)**2+bondvec2(2,im)**2+bondvec2(3,im)**2)
    enddo
    watercomz=sum(rct(3,:))/na*3.d0
    rct(3,:)=rct(3,:)-boxlxyz(3)*anint((rct(3,:)-watercomz)/boxlxyz(3)) !shift molecules in a box around the COM (in z)
    rct(3,:)=rct(3,:)+boxcenterz-watercomz                              !move the created box to the center of the simulation box

                                        !density field calculation routine (consumes much time)
    do densi=1,ntpts                    !approx 100 million loop steps for 11*11*11 molecules and VF 3
        do im=1,na/3
            distvec(:)=rct(:,im)-grid(:,densi)
            distvec(:)=distvec(:)-boxlxyz(:)*anint(distvec(:)/boxlxyz(:))       !minimum image convention

            if (abs(distvec(3)).lt.cut) then    !faster for large z-direction boxes (as vacuum boxes are)
                dist=dot_product(distvec,distvec)
                if (dist .lt. cutsq) then       !cutoff 4 cgxi -> 16 cgxisq

                    !coarse graining with gaussian (2*pi*cgxi^2)^(-3/2) * Exp[-((r-rm)/cgxi)^2 / 2]
                    instdfield(densi)=instdfield(densi)+exp(exparg*dist)

                end if
            end if
        enddo
    enddo

    instdfield(:)= na*6.d0*au/vbox / sum(instdfield)*ntpts * instdfield(:)  !normalization and unit transformation to g/cm^3

                     !interface calculation routine
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
                        boxsh=grid(3,densi)-boxlxyz(3)
                    else
                        boxsh=grid(3,densi)
                    end if
                    if (modz2.gt.maxdensz) then
                        boxsh2=grid(3,densi2)-boxlxyz(3)
                    else
                        boxsh2=grid(3,densi2)
                    end if
                    !if you replace gridspac with z2-z1 you get the interpolation formula from above

                    interface(inti,1)=boxsh+(boxsh2-boxsh)*(hbulk-instdfield(densi)) / &
                    (instdfield(densi2)-instdfield(densi))

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
                        boxsh=grid(3,densi)+boxlxyz(3)
                    else
                        boxsh=grid(3,densi)
                    end if
                    if (modz2.lt.maxdensz) then
                        boxsh2=grid(3,densi2)+boxlxyz(3)
                    else
                        boxsh2=grid(3,densi2)
                    end if

                    interface(inti,2)=boxsh+(boxsh2-boxsh)*(hbulk-instdfield(densi)) / &
                    (instdfield(densi2)-instdfield(densi))

                    exit
                end if
            enddo
        enddo

        hvec1(1)=1.d0
        hvec1(2)=0.d0
        hvec2(1)=0.d0
        hvec2(2)=1.d0

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
                dist=dist-boxlxyz(1)*anint(dist/boxlxyz(1))     !needed for border case
                hvec1(3)=(interface(densi2,int12)-interface(densi,int12))/dist

                densi=((inti-1)/xyzpts(2))*xyzpts(2)+modulo(inti-2,xyzpts(2))+1
                densi2=((inti-1)/xyzpts(2))*xyzpts(2)+modulo(inti,xyzpts(2))+1          !Fy
                dist=grid(2,densi2*xyzpts(3))-grid(2,densi*xyzpts(3))
                dist=dist-boxlxyz(2)*anint(dist/boxlxyz(2))
                hvec2(3)=(interface(densi2,int12)-interface(densi,int12))/dist

                intnvec(1,inti)=-hvec1(3)       !thats what is left over from the hvec1 x hvec2 if all 0 and 1 are erased
                intnvec(2,inti)=-hvec2(3)
                intnvec(3,inti)=1.d0
                intnvec(:,inti)=intnvec(:,inti)/sqrt(intnvec(1,inti)**2 + &
                intnvec(2,inti)**2 + intnvec(3,inti)**2)*(-1)**int12

            enddo

            do im=1,na/3    !now calculate properties depending on normal vectors

                do inti=1,xypts             !find closest surface point
                    distvec(1)=grid(1,inti*xyzpts(3))-rct(1,im)
                    distvec(2)=grid(2,inti*xyzpts(3))-rct(2,im)
                    distvec(3)=interface(inti,int12)-rct(3,im)
                    distvec(:)=distvec(:)-boxlxyz(:)*anint(distvec(:)/boxlxyz(:)) !min img conv
                    hprox(inti)=dot_product(distvec,distvec)
                enddo

                inti=MINLOC(hprox,1)        !here

                distvec(1)=grid(1,inti*xyzpts(3))-rct(1,im)
                distvec(2)=grid(2,inti*xyzpts(3))-rct(2,im)
                distvec(3)=interface(inti,int12)-rct(3,im)
                distvec(:)=distvec(:)-boxlxyz(:)*anint(distvec(:)/boxlxyz(:))

                proximity(im)=distvec(1)*intnvec(1,inti)+distvec(2)*intnvec(2,inti)+distvec(3)*intnvec(3,inti)

                if (proximity(im).gt.lprox1 .and. proximity(im).lt.uprox1) then
                    inint(im,1,int12)=.true.
                else
                    inint(im,1,int12)=.false.
                end if

                if (proximity(im).gt.lprox2 .and. proximity(im).lt.uprox2) then
                    inint(im,2,int12)=.true.
                else
                    inint(im,2,int12)=.false.
                end if
                
                inint(im,3,int12)=inint(im,1,int12).or.inint(im,2,int12)
            enddo
        enddo
    end if
	
    deallocate (rct,grid,intnvec,bondvec1,bondvec2,proximity,hprox,instdfield,interface)
    return
end subroutine instvacint
