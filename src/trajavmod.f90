module avmod
    use intmod
    implicit none

    integer ntraj,itraj

    !--density arrays--
    real(8), allocatable :: avmeandfield(:),avdensnz(:),avdensnzmean(:)

    !--interface arrays--
    real(8), allocatable :: avmeannvec(:,:,:),avjcdist(:)
    real(8) avmeanint(2),avrhobulk

    !--hbond arrays--
    real(8), allocatable :: avintrabond(:,:),avinterbond(:,:),avinintct(:,:,:),avhbct(:,:),avavbno(:,:)
    real(8), allocatable :: avddct(:,:),avsdct(:,:),avndct(:,:),avacct(:,:)
    real(8)              :: avhbgt1(4)

    !--correlation and post-loop arrays--
    real(8), allocatable :: avintcorr(:),avintmean(:,:),avintsigma(:,:)
    real(8), allocatable :: avinintcorr(:,:),avbondcorr(:,:,:),avrestbcorr(:,:,:)
    real(8), allocatable :: avbrokrepcorr(:,:,:),avbcrcorr(:,:,:)

contains

    subroutine avinit(trajno)
        implicit none

        integer trajno
        ntraj=trajno
        itraj=0

        allocate (avmeannvec(3,xypts,2), &
        avjcdist(jcpts*jcpts*zpts),avdensnz(zpts),avdensnzmean(zpts), &
        avmeandfield(ntpts), &
        avintcorr(intpts),avhbct(4,4),avavbno(4,4), &
        avintmean(xypts,2),avintsigma(xypts,2),avintrabond(4,3),avinterbond(4,3), &
        avinintct(4,2,intpts),avinintcorr(4,hbcorrlen),avbondcorr(4,4,hbcorrlen), &
        avrestbcorr(4,4,hbcorrlen),avbrokrepcorr(4,4,hbcorrlen),avbcrcorr(4,4,hbcorrlen), &
        avddct(4,4),avsdct(4,4),avndct(4,4),avacct(4,4))

        avmeannvec(:,:,:)=0.d0
        avjcdist(:)=0.d0
        avmeandfield(:)=0.d0
        avdensnz(:)=0.d0
        avdensnzmean(:)=0.d0

        avintrabond(:,:)=0.d0
        avinterbond(:,:)=0.d0
        avinintct(:,:,:)=0.d0
        avhbct(:,:)=0.d0
        avavbno(:,:)=0.d0
        avddct(:,:)=0.d0
        avsdct(:,:)=0.d0
        avndct(:,:)=0.d0
        avacct(:,:)=0.d0
        avhbgt1(:)=0.d0

        avintmean(:,:)=0.d0
        avintsigma(:,:)=0.d0
        avmeanint(:)=0.d0
        avrhobulk=0.d0

        avintcorr(:)=0.d0
        avinintcorr(:,:)=0.d0
        avbondcorr(:,:,:)=0.d0
        avrestbcorr(:,:,:)=0.d0
        avbrokrepcorr(:,:,:)=0.d0
        avbcrcorr(:,:,:)=0.d0
    end subroutine

    subroutine avupdate()
        implicit none

        itraj=itraj+1
        avmeannvec(:,:,:)=avmeannvec(:,:,:)+meannvec(:,:,:)
        avjcdist(:)=avjcdist(:)+jcdist(:)
        avmeandfield(:)=avmeandfield(:)+meandfield(:)
        avdensnz(:)=avdensnz(:)+densnz(:)
        avdensnzmean(:)=avdensnzmean(:)+densnzmean(:)

        avintrabond(:,:)=avintrabond(:,:)+intrabond(:,:)
        avinterbond(:,:)=avinterbond(:,:)+interbond(:,:)
        avinintct(:,:,:)=avinintct(:,:,:)+inintct(:,:,:)
        avhbct(:,:)=avhbct(:,:)+hbct(:,:)
        avavbno(:,:)=avavbno(:,:)+avbno(:,:)
        avddct(:,:)=avddct(:,:)+ddct(:,:)
        avsdct(:,:)=avsdct(:,:)+sdct(:,:)
        avndct(:,:)=avndct(:,:)+ndct(:,:)
        avacct(:,:)=avacct(:,:)+acct(:,:)
        avhbgt1(:)=avhbgt1(:)+hbgt1(:)

        avintmean(:,:)=avintmean(:,:)+intmean(:,:)
        avintsigma(:,:)=avintsigma(:,:)+intsigma(:,:)
        avmeanint(:)=avmeanint(:)+meanint(:)
        avrhobulk=avrhobulk+rhobulk

        avintcorr(:)=avintcorr(:)+intcorr(:)
        avinintcorr(:,:)=avinintcorr(:,:)+inintcorr(:,:)
        avbondcorr(:,:,:)=avbondcorr(:,:,:)+bondcorr(:,:,:)
        avrestbcorr(:,:,:)=avrestbcorr(:,:,:)+restbcorr(:,:,:)
        avbrokrepcorr(:,:,:)=avbrokrepcorr(:,:,:)+brokrepcorr(:,:,:)
        avbcrcorr(:,:,:)=avbcrcorr(:,:,:)+bcrcorr(:,:,:)
    end subroutine

    subroutine avcalcm()
        implicit none
        if (itraj.ne.ntraj) then
            write(6,*) 'Warning: itraj!=ntraj'
        end if
        avmeannvec(:,:,:)=avmeannvec(:,:,:)/ntraj
        avjcdist(:)=avjcdist(:)/ntraj
        avmeandfield(:)=avmeandfield(:)/ntraj
        avdensnz(:)=avdensnz(:)/ntraj
        avdensnzmean(:)=avdensnzmean(:)/ntraj

        avintrabond(:,:)=avintrabond(:,:)/ntraj
        avinterbond(:,:)=avinterbond(:,:)/ntraj
        avinintct(:,:,:)=avinintct(:,:,:)/ntraj
        avhbct(:,:)=avhbct(:,:)/ntraj
        avavbno(:,:)=avavbno(:,:)/ntraj
        avddct(:,:)=avddct(:,:)/ntraj
        avsdct(:,:)=avsdct(:,:)/ntraj
        avndct(:,:)=avndct(:,:)/ntraj
        avacct(:,:)=avacct(:,:)/ntraj
        avhbgt1(:)=avhbgt1(:)/ntraj

        avintmean(:,:)=avintmean(:,:)/ntraj
        avintsigma(:,:)=avintsigma(:,:)/ntraj
        avmeanint(:)=avmeanint(:)/ntraj
        avrhobulk=avrhobulk/ntraj

        avintcorr(:)=avintcorr(:)/ntraj
        avinintcorr(:,:)=avinintcorr(:,:)/ntraj
        avbondcorr(:,:,:)=avbondcorr(:,:,:)/ntraj
        avrestbcorr(:,:,:)=avrestbcorr(:,:,:)/ntraj
        avbrokrepcorr(:,:,:)=avbrokrepcorr(:,:,:)/ntraj
        avbcrcorr(:,:,:)=avbcrcorr(:,:,:)/ntraj
    end subroutine

    subroutine avintprint ()
        implicit none
        integer densi,densx,densy,densz,iter,int12,inti
        real(8) sumint(2)


        open(1001,file='meandfield.txt')
        open(1005,file='intcorr.txt')
        open(1007,file='intmeansig1.txt')
        open(1008,file='intmeansig2.txt')
        open(1009,file='densint.txt')
        open(1011,file='densmean.txt')
        open(1013,file='jcdist.txt')
        open(1015,file='bondcorr1.txt')
        open(1016,file='bondcorr2.txt')
        open(1017,file='bondcorr3.txt')
        open(1018,file='bondcorr4.txt')
        open(2015,file='bcrcorr1.txt')
        open(2016,file='bcrcorr2.txt')
        open(2017,file='bcrcorr3.txt')
        open(2018,file='bcrcorr4.txt')
        open(3015,file='restbcorr1.txt')
        open(3016,file='restbcorr2.txt')
        open(3017,file='restbcorr3.txt')
        open(3018,file='restbcorr4.txt')
        open(4015,file='brokrepcorr1.txt')
        open(4016,file='brokrepcorr2.txt')
        open(4017,file='brokrepcorr3.txt')
        open(4018,file='brokrepcorr4.txt')
        open(5015,file='inintcorr1.txt')
        open(5016,file='inintcorr2.txt')
        open(5017,file='inintcorr3.txt')
        open(5018,file='inintcorr4.txt')
        open(1020,file='interbond.txt')
        open(1021,file='intrabond.txt')
        open(1022,file='totalbond.txt')
        open(1023,file='donor.txt')
        open(1025,file='avbno.txt')

        do densi=1,ntpts
            write(1001,*) grid(1,densi),grid(2,densi),grid(3,densi),avmeandfield(densi)
        enddo

        do densi=1,intpts
        write(1005,*) (densi-1)*cstep*dtpfs,avintcorr(densi)
        enddo

        do iter=1,4
            do densi=1,hbcorrlen
                write(1014+iter,*) (densi-1)*dtpfs*cstep,avbondcorr(:,iter,densi)
                write(2014+iter,*) (densi-1)*dtpfs*cstep,avbcrcorr(:,iter,densi)
                write(3014+iter,*) (densi-1)*dtpfs*cstep,avrestbcorr(:,iter,densi)
                write(4014+iter,*) (densi-1)*dtpfs*cstep,avbrokrepcorr(:,iter,densi)
                write(5014+iter,*) (densi-1)*dtpfs*cstep,avinintcorr(iter,densi)
            enddo
        enddo

        do iter=1,3
            sumint(1)=sum(avinintct(iter,1,:))
            sumint(2)=sum(avinintct(iter,2,:))
            sumint(1)=0.5d0*(sumint(1)+sumint(2))/intpts
            write(1020,*) iter,sumint(1),avinterbond(:,iter)
            write(1021,*) iter,sumint(1),avintrabond(:,iter)
        enddo

        do iter=1,4
            write(1022,*) iter,avhbct(:,iter)
            write(1025,*) iter,avavbno(:,iter)
        enddo



        do iter=1,4
            write(1023,*) 'Layer',iter
            do densx=1,4
                write(1023,*) densx,avndct(densx,iter),avsdct(densx,iter),avddct(densx,iter),avacct(densx,iter)
            enddo
        enddo


        do int12=1,2
            do inti=1,xypts
                write(1006+int12,*) grid(1,inti*xyzpts(3)),grid(2,inti*xyzpts(3)),avintmean(inti,int12),&
                avintsigma(inti,int12)
            enddo
        enddo
        do densz=1,zpts
            write(1009,*) (densz-0.5d0)*zspac-maxrang,avdensnz(densz)
        enddo
        do densz=1,zpts
            write(1011,*) (densz-0.5d0)*zspac-maxrang,avdensnzmean(densz)
        enddo
        do densz=1,zpts
            if ((densz-0.5d0)*zspac-maxrang.gt.-2.d0 .and. &
            (densz-0.5d0)*zspac-maxrang.lt.abs(meanint(2)-meanint(1))/2.d0) then
                write(1013,*) 'Dist',(densz-0.5d0)*zspac-maxrang,jcpts*jcpts
                do densx=1,jcpts
                    if ((densx-0.5d0)*jcspac-1.d0 .lt. 1.d0) then
                        do densy=1,jcpts
                            if ((densy-0.5d0)*jcspac-1.d0 .lt. 1.d0) then
                                write(1013,*) (densx-0.5d0)*jcspac-1.d0,(densy-0.5d0)*jcspac-1.d0, &
                                avjcdist(densy+(densx-1)*jcpts+(densz-1)*jcpts*jcpts)
                            else
                                write(1013,*) (densx-0.5d0)*jcspac-1.d0,0.5d0*(densy-1)*jcspac, &
                                avjcdist(densy+(densx-1)*jcpts+(densz-1)*jcpts*jcpts)
                            end if
                        enddo
                    else
                        do densy=1,jcpts
                            if ((densy-0.5d0)*jcspac-1.d0 .lt. 1.d0) then
                                write(1013,*) 0.5d0*(densx-1)*jcspac,(densy-0.5d0)*jcspac-1.d0, &
                                avjcdist(densy+(densx-1)*jcpts+(densz-1)*jcpts*jcpts)
                            else
                                write(1013,*) 0.5d0*(densx-1)*jcspac,0.5d0*(densy-1)*jcspac, &
                                avjcdist(densy+(densx-1)*jcpts+(densz-1)*jcpts*jcpts)
                            end if
                        enddo
                    end if
                enddo
            end if
        enddo

    end subroutine

end module avmod
