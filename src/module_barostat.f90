module barostat
  implicit none
  integer nacc,nmove
  real(8) taub,patm,dv
  character(len=3) baro
contains
 
  subroutine beren_driver(vir,tv,tvxyz,dt,r,boxlxyz,na,nb)
    implicit none
    ! ------------------------------------------------------------------
    ! Berendesen Driver
    ! ------------------------------------------------------------------
    integer na,nb
    real(8) r(3,na,nb),vir(3,3),tvxyz(3),boxlxyz(3)
    real(8) dt,tv
    logical iamcub
    common /symmetry/ iamcub
    
    if (iamcub) then
       call beren_baro(vir,tv,dt,r,boxlxyz,na,nb)
    else
       call beren_baro_nc(vir,tvxyz,dt,r,boxlxyz,na,nb)
    endif
    
    return
  end subroutine beren_driver
  
  subroutine beren_baro_nc(vir,tvxyz,dt,r,boxlxyz,na,nb)
    implicit none
    ! ------------------------------------------------------------------
    ! Berendsen Barostat - Non-Isotropic Version
    !
    ! compi is the isothermal compressibility of water in atomic units
    ! ------------------------------------------------------------------
    integer na,nm,nb,i,j,k,ii,jj,kk,nacc,nmove,rpmddft
    real(8) r(3,na,nb),boxlxyz(3),scale(3),pres(3),patm,dt
    real(8) vir(3,3),tvxyz(3),compi,bstat,dv,ptail,vol
    character(len=3) baro
    real(8) oo_eps,oo_sig,oo_gam,rcut
    real(8), allocatable :: rcm(:,:)
    common /oo_param/ oo_eps,oo_sig,oo_gam,rcut
    common /RPMDDFT/ rpmddft
    ! Current pressure in x, y and z directions
    
    if (oo_gam.eq.0.d0.and.rpmddft.eq.0) then !!! LJ correction nicht bei RPMDDFT
       call pres_lj_tail(ptail,boxlxyz,na)
    else
       ptail = 0.d0
    endif
    
    vol = boxlxyz(1)*boxlxyz(2)*boxlxyz(3)
    do j = 1,3
       pres(j) = (-vir(j,j) + 2.d0*tvxyz(j))/vol
       pres(j) = pres(j) + ptail
    enddo
    
    ! Non-Isotropic Berendesen Barostat
    
    nm = na/3
    
    compi = 13531.d0
    bstat = compi/(3.d0*taub)
    
    scale(:) = abs(1.d0 - (dt*bstat)*(patm - pres(:)))
    boxlxyz(:) = scale(:) * boxlxyz(:)
    
    ! Scale water COMs
    
    allocate (rcm(3,nm))
    call center_water(r,rcm,nm,nb)
    do i = 1, nm
       ii = 3*i-2
       jj = 3*i-1
       kk = 3*i
       do k = 1, nb
          do j = 1,3
             r(j,ii,k) = r(j,ii,k) + rcm(j,i)*(scale(j)-1.d0)
             r(j,jj,k) = r(j,jj,k) + rcm(j,i)*(scale(j)-1.d0)
             r(j,kk,k) = r(j,kk,k) + rcm(j,i)*(scale(j)-1.d0)
          enddo
       enddo
    enddo
    deallocate (rcm)
    
    return
  end subroutine beren_baro_nc
  
  subroutine beren_baro(vir,tv,dt,r,boxlxyz,na,nb)
    implicit none
    ! ------------------------------------------------------------------
    ! Berendsen Barostat - Isotropic Version
    !
    ! compi is the isothermal compressibility of water in atomic units
    ! ------------------------------------------------------------------
    integer na,nm,nb,i,k,ii,jj,kk
    real(8) r(3,na,nb),boxlxyz(3),scale,dt
    real(8) vir(3,3),tv,compi,bstat,pres
    real(8), allocatable :: rcm(:,:)
    
    ! Current pressure
    
    call pressure(pres,vir,tv,na,boxlxyz)
    
    ! Isotropic Berendesen Barostat
    
    nm = na/3
    compi = 13531.d0
    bstat = compi/(3.d0*taub)
    scale = abs(1.d0 - (dt*bstat)*(patm - pres))
    boxlxyz(:) = scale * boxlxyz(:)

    ! Scale water COMs
    
    allocate (rcm(3,nm))
    call center_water(r,rcm,nm,nb)  ! Wenn Fehler in NPT dann hier !!!! da nicht wirklich center of water sonder mit Massen gewichtet!
    do i = 1, nm 
       ii = 3*i-2
       jj = 3*i-1
       kk = 3*i
       do k = 1, nb
          r(1:3,ii,k) = r(1:3,ii,k) + rcm(1:3,i)*(scale-1.d0)
          r(1:3,jj,k) = r(1:3,jj,k) + rcm(1:3,i)*(scale-1.d0)
          r(1:3,kk,k) = r(1:3,kk,k) + rcm(1:3,i)*(scale-1.d0)
       enddo
    enddo
    deallocate (rcm)
    
    return
  end subroutine beren_baro

  subroutine mc_baro(r,dvdr,dvdr2,vir,eo,z,beta,boxlxyz, &
                     na,nb,irun)
    implicit none
    ! ------------------------------------------------------------------
    ! Monte Carlo stochastic volume move
    ! 
    ! dv      is the maximum move size in ln(volume)
    ! naccept is the number of accepted moves
    ! nmove   is the number of moves attempted
    ! ------------------------------------------------------------------
    integer na,nm,nb,irun,i,k,ii,jj,kk
    real(8) r(3,na,nb),dvdr(3,na,nb),dvdr2(3,na,nb),vir(3,3)
    real(8) boxlxyz(3),z(na),boxlxyzn(3),virn(3,3),vew,vlj,vint
    real(8) vo,vn,vnlog,eo,en,beta,test,scale,ran2
    real(8) thresh
    logical accept
    external ran2
    real(8), allocatable :: rcm(:,:),rn(:,:,:),dvdrn(:,:,:)
    real(8), allocatable :: dvdr2n(:,:,:)

    nm = na/3
    thresh = 1.d0/taub

    if (ran2(irun,0.d0,1.d0) .lt. thresh) then

      call full_forces(r,na,nb,eo,vew,vlj,vint,vir,z,boxlxyz, &   !dvdr und dvdr2 wird neu berechnet
                       dvdr,dvdr2)

       nmove = nmove + 1
       vo = boxlxyz(1)*boxlxyz(2)*boxlxyz(3)

       ! Pick a random change in ln(volume)
       ! vnlog = dlog(vo) + dv*ran2(irun,-1.d0,1.d0)
       ! vn = dexp(vnlog)

       ! Pick a random change in volume

       vn = vo + dv*ran2(irun,-1.d0,1.d0)

       scale = (vn/vo)**(1.d0/3.d0)
       boxlxyzn(:) = scale * boxlxyz(:)

       ! Scale water COMs

       allocate (rcm(3,na/3))
       allocate (rn(3,na,nb),dvdrn(3,na,nb),dvdr2n(3,na,nb))

       call center_water(r,rcm,nm,nb)
       do i = 1, nm
          ii = 3*i-2
          jj = 3*i-1
          kk = 3*i
          do k = 1, nb
             rn(1:3,ii,k) = r(1:3,ii,k) + rcm(1:3,i)*(scale-1.d0)
             rn(1:3,jj,k) = r(1:3,jj,k) + rcm(1:3,i)*(scale-1.d0)
             rn(1:3,kk,k) = r(1:3,kk,k) + rcm(1:3,i)*(scale-1.d0)
          enddo
       enddo

       ! Calculate new forces and energy

       call full_forces(rn,na,nb,en,vew,vlj,vint,virn,z,boxlxyzn, &
                        dvdrn,dvdr2n)                                 

       ! Acceptance test

       test = (en-eo) + patm*(vn-vo) - dble(nm+1)*dlog(vn/vo)/beta
       test = (en-eo) + patm*(vn-vo) - dble(nm)*dlog(vn/vo)/beta
       test = dexp(-beta*test)
       accept = test.gt.ran2(irun,0.d0,1.d0)

       if (accept) then

          ! Replace energies, forces, etc. with new values

          nacc = nacc + 1
          eo = en
          boxlxyz(:) = boxlxyzn(:)
          r(:,:,:) = rn(:,:,:)
          dvdr(:,:,:) = dvdrn(:,:,:)
          dvdr2(:,:,:) = dvdr2n(:,:,:)
          vir(:,:) = virn(:,:)
       endif
       deallocate (rcm)
       deallocate (rn,dvdrn,dvdr2n)
    endif

    return

  end subroutine mc_baro
end module barostat
