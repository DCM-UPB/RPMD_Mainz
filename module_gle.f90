module gle
  use thermostat 
  implicit none
  real(8), allocatable, save ::  gS(:,:), gT(:,:), gp(:,:,:,:), ngp(:,:,:,:)
  integer ns
contains
  
  subroutine therm_pro(p,mass,beta,ctau,dt,na,nb,irun)
    implicit none
    include 'globals.inc'
    ! ------------------------------------------------------------------
    ! "Critical damping" white noise Langevin for internal ring modes
    ! ------------------------------------------------------------------
    integer na,nb,irun,k,j,i,ic
    real(8), intent(inout):: p(3,na,nb)
    real(8) mass(na),gamma,dt,beta, ctau
    real(8) c1,c2,ck(3),c3,gaussian
    real(8) pi,pibyn,twown,betan,radtwo,wk
    real(8) mm, pcom(3)
    external gaussian
    
    betan=beta/nb
    pi = dacos(-1.d0)
    twown = 2.d0/(betan*hbar)
    pibyn = pi/nb
    radtwo=dsqrt(2.0d0)
    
    do j = 1, 3
       ck(j)=dsqrt(mass(j))
    enddo
    
    do k = 1,nb
       wk=twown*dsin((k-1)*pibyn)
       if (k.eq.1) then
          gamma=2.*pi/ctau   !thermostat on the centroid is given from options
       else
          gamma=2.*wk       !thermostat on necklace mode, "optimal damping" choice
       endif
       c1 = dexp (-gamma*dt)
       c2 = dsqrt((1.d0 - c1*c1)*0.5d0/betan)
       if(k.eq.1 .or. ((mod(nb,2).eq.0) .and. k.eq.(nb/2+1))) c2=radtwo*c2
       
       do j = 1, na
          ic = mod(j,3)
          if (ic.eq.0) then
             ic = 3
          endif
          c3 = c2 * ck(ic)
          do i = 1, 3
             p(i,j,k) = c1 * p(i,j,k) + c3*gaussian(irun,1.d0)
          enddo
       enddo
    enddo
    
    !remove COM from the centroid mode
    pcom=0.d0
    mm=0.d0
    do j = 1, na
       do i = 1, 3
          pcom(i)=pcom(i)+p(i,j,1)
       enddo
       mm=mm+mass(j)
    enddo
    pcom=pcom*(1.0d0/mm)
    do j = 1, na
       mm=mass(j)
       do i = 1, 3
          p(i,j,1)=p(i,j,1)-pcom(i)*mm
       enddo
    enddo
    
    return
  end subroutine therm_pro
  
    subroutine therm_pro_global(p,mass,beta,ctau,dt,na,nb,irun)
    implicit none
    include 'globals.inc'
    ! ------------------------------------------------------------------
    ! "Critical damping" white noise Langevin for internal ring modes
    ! and global centroid thermostat
    ! ------------------------------------------------------------------
    integer na,nb,irun,k,j,i,ic
    real(8), intent(inout):: p(3,na,nb)
    real(8) mass(na),gamma,dt,beta, ctau
    real(8) c1,c2,ck(3),c3,gaussian
    real(8) pi,pibyn,twown,betan,radtwo,wk
    real(8) mm, pcom(3)
    external gaussian
    
    betan=beta/nb
    pi = dacos(-1.d0)
    twown = 2.d0/(betan*hbar)
    pibyn = pi/nb
    radtwo=dsqrt(2.0d0)
    
    do j = 1, 3
       ck(j)=dsqrt(mass(j))
    enddo
    
    ! Critical damping white noise for non centroid modes
    
    do k = 2,nb
       wk=twown*dsin((k-1)*pibyn)
       gamma=2.*wk       !thermostat on necklace mode, "optimal damping" choice
       c1 = dexp (-gamma*dt)
       c2 = dsqrt((1.d0 - c1*c1)*0.5d0/betan)
       if(((mod(nb,2).eq.0) .and. k.eq.(nb/2+1))) c2=radtwo*c2
       
       do j = 1, na
          ic = mod(j,3)
          if (ic.eq.0) then
             ic = 3
          endif
          c3 = c2 * ck(ic)
          do i = 1, 3
             p(i,j,k) = c1 * p(i,j,k) + c3*gaussian(irun,1.d0)
          enddo
       enddo
    enddo
    
    ! Centroid mode thermostatted by global Bussi-Parrinello thermostat
    
    call parinello_therm(p(1,1,1),mass,ctau,na,1,dt,irun,betan)
    
    !remove COM from the centroid mode
    pcom=0.d0
    mm=0.d0
    do j = 1, na
       do i = 1, 3
          pcom(i)=pcom(i)+p(i,j,1)
       enddo
       mm=mm+mass(j)
    enddo
    pcom=pcom*(1.0d0/mm)
    do j = 1, na
       mm=mass(j)
       do i = 1, 3
          p(i,j,1)=p(i,j,1)-pcom(i)*mm
       enddo
    enddo
    
    return
  end subroutine therm_pro_global
  
  
  subroutine therm_gle_init(ttau,na,nb,dt,irun,beta)
    implicit none
    real(8), intent(in)  :: dt, beta, ttau
    integer, intent(in) :: na, nb, irun
    real(8), allocatable :: gA(:,:), gC(:,:), gr(:)
    integer i, j, k, h, cns, ios
    real(8), external :: gaussian     !this is in the pH2 RPMD code, replace it with appropiate normal deviate sampler
    
    !reads in matrices
    !reads A (in units of the maximum frequency present)
    open(121,file='GLE-A',status='OLD',iostat=ios)
    read(121,*) ns
    write(6,*)  "GLE INIT CALLED, ns=",ns," na=",na," nb=",nb
    !allocate everything we need
    allocate(gA(ns+1,ns+1))
    allocate(gC(ns+1,ns+1))
    allocate(gS(ns+1,ns+1))
    allocate(gT(ns+1,ns+1))
    allocate(gp(3,na,nb,ns+1))   
    allocate(ngp(3,na,nb,ns+1))   
    allocate(gr(ns+1))
    
    if (ios.ne.0) write(0,*) "COULD NOT FIND GLE-A!"
    do i=1,ns+1
       read(121,*) gA(i,:)
    enddo
    close(121)
    gA=gA*(2.*3.1415927/ttau)
    !reads C (in K)
    open(121,file='GLE-C',status='OLD',iostat=ios)
    if (ios.ne.0) then            
       gC=0.d0
       do i=1,ns+1
          gC(i,i)=(nb/beta)
       enddo
    else           
       read(121,*) cns
       if (cns.ne.ns) write(0,*) "NS MISMATCH BETWEEN GLE-A and GLE-C!!"
       do i=1,ns+1
          read(121,*) gC(i,:)
       enddo
       gC=gC*(nb*3.166829d-6) !goes in energy units
    end if
    
    ! the deterministic part of the propagator is obtained in a second
    call matrix_exp(-dt*gA, ns+1,15,15,gT)
    
    ! the stochastic part is just as easy. we use gA as a temporary array
    gA=gC-matmul(gT,matmul(gC,transpose(gT)))
    call cholesky(gA, gS, ns+1)
    
    ! then, we must initialize the auxiliary vectors. we keep general - as we might be 
    ! using non-diagonal C to break detailed balance - and we use cholesky decomposition
    ! of C. again, since one would like to initialize correctly the velocities in 
    ! case of generic C, we use an extra slot for gp for the physical momentum, as we 
    ! could then use it to initialize the momentum in the calling code
    gA=gC   
    call cholesky(gA, gC, ns+1)
    
    do h=1,3
       do k=1,na
          do j=1,nb
             do i=1,ns+1
                gr(i)=gaussian(irun,1.0d0)
             enddo
             gp(h,k,j,:)=matmul(gC,gr)
          enddo
       enddo
    end do
    
    write(6,*) "GLE INITIALIZATION DONE"
    deallocate(gA)
    deallocate(gC)
    deallocate(gr)
    ! debug: prints out the propagator matrices
    ! open(121,file='CMP-T',status='NEW',iostat=ios)
    ! write(121,*) gT
    ! close(121)
    ! open(121,file='CMP-S',status='NEW',iostat=ios)
    ! write(121,*) gS
    ! close(121)
  end subroutine therm_gle_init

  subroutine therm_gle(p,dheat,mass,na,nb,irun)
    implicit none
    integer,intent(in)  :: na, nb,irun
    real(8), intent(in) :: mass(na)
    real(8), intent(inout) :: p(3,na,nb)
    real(8), intent(out) :: dheat
    real(8), external :: gaussian
    integer i, j, k, h, n
    real(8) mfac
    n=3*na*nb
    dheat=0.d0
    ! write(6,*) "GLE PROPAGATOR HAS BEEN CALLED"
    do k=1,na
       mfac=1.0/dsqrt(mass(k))
       do h=1,3
          do j=1,nb  ! go to mass-scaled coordinates
             gp(h,k,j,1)=p(h,k,j)*mfac
          enddo
       enddo
    end do
    !we pretend that gp is a (3 na nb)x(ns+1) matrix, which should be fine....
    call dgemm('n','t',n,ns+1,ns+1,1.0d0,gp,n,gT,ns+1,0.0d0,ngp,n)
    
    !now, must compute random part. 
    !first, fill up gp of random n
    do i=1,ns+1
       do j=1,nb
          do k=1,na
             do h=1,3
                gp(h,k,j,i)=gaussian(irun,1.0d0)
             end do
          end do
       end do
    end do
    
    call dgemm('n','t',n,ns+1,ns+1,1.0d0,gp,n,gS,ns+1,1.0d0,ngp,n)
    gp=ngp
    
    !clean out COM component including ALL THE ADDITIONAL DOF, in mass-scaled coordinates
    !this is not strictly necessary, but ensures that the "pseudo-normal-modes" are consistent
    !between physical and fictitious DOF, so that the thermostat can "feel" the appropriate
    !dynamics.
    !NB: for the same reason, if working with rigid bodies, all the fictitious momenta should
    !kept orthogonal to the constraints. - THIS IS NOT IMPLEMENTED RIGHT NOW!!!
    call pshift_gle(gp,mass,na,nb,ns)
    
    do k=1,na
       mfac=dsqrt(mass(k))
       do h=1,3
          do j=1,nb  ! back to conventional momenta
             p(h,k,j)=gp(h,k,j,1)*mfac
          enddo
       enddo
    end do
  end subroutine therm_gle

  subroutine pshift_gle(gp,mass,na,nb,ns)
    integer,intent(in)  :: na, nb, ns
    real(8), intent(in) :: mass(na)
    real(8), intent(inout) :: gp(3,na,nb,ns+1)
    integer i,j,k,h
    real(8) mfac, mm
    real(8) pcom(3)
    
    mm=0.d0
    do i=1,na
       mm=mm+mass(i)
    enddo
    mm=1./sqrt(nb*mm)
    
    do i=1,ns+1
       pcom=0.d0
       do j=1,na
          mfac=dsqrt(mass(j))
          do k=1,nb
             do h=1,3
                pcom(h)=pcom(h)+gp(h,j,k,i)*mfac
             enddo
          enddo
       enddo
       pcom=pcom*mm*mm
       
       do j=1,na
          mfac=dsqrt(mass(j))
          do k=1,nb
             do h=1,3
                gp(h,j,k,i)=gp(h,j,k,i)-pcom(h)*mfac
             enddo
          enddo
       enddo
    enddo
  end subroutine pshift_gle

  ! matrix exponential by scale & square      
  subroutine matrix_exp(M, n, j, k, EM)
    integer, intent(in)  :: n, j, k
    real(8), intent(in)   :: M(n,n)
    real(8), intent(out)   :: EM(n,n)
    
    real(8) :: tc(j+1), SM(n,n)
    integer p, i
    tc(1)=1.d0
    do i=1,j
       tc(i+1)=tc(i)/dble(i)
    enddo
    
    !scale
    SM=M*(1./2.**k)
    EM=0.d0
    do i=1,n
       EM(i,i)=tc(j+1)
    enddo
    
    !taylor exp of scaled matrix
    do p=j,1,-1
       EM=matmul(SM,EM);
       do i=1,n
          EM(i,i)=EM(i,i)+tc(p)
       enddo
    enddo
    
    !square
    do p=1,k
       EM=matmul(EM,EM)
    enddo
  end subroutine matrix_exp
  
  ! brute-force cholesky decomposition
  subroutine cholesky(SST, S, n)
    integer, intent(in)  :: n
    real(8), intent(in)   :: SST(n,n)
    real(8), intent(out)   :: S(n,n)
    integer i,j,k
    S=0.d0
    S(1,1)=sqrt(SST(1,1))
    do i=1,n
       do j=1,i-1
          S(i,j)=SST(i,j);
          do k=1,j-1
             S(i,j)=S(i,j)-S(i,k)*S(j,k)
          enddo
          S(i,j)=S(i,j)/S(j,j)
       enddo
       S(i,i)=SST(i,i)
       do k=1,i-1
          S(i,i)=S(i,i)-S(i,k)**2
       end do
       S(i,i)=sqrt(S(i,i))
    enddo
  end subroutine cholesky
  
end module gle
