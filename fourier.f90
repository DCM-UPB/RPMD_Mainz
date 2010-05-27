
!**********************************************************************!
!                                                                      !
! realft    - FFT of m real arrays (if mode = 1) or complex Hermitian  !
!             arrays in real storage (if mode = -1), using -lfftw3.    !
!             Support multiple plans- up to max_plan                   !
!                                                                      !
!                                                                      !
! Arguments: data(m,n) --> array to transform                          !
!            m         --> array size 1                                !
!            n         --> array size 2                                !
!            mode      --> +1 or -1, defining forwards or backwards    !
!                          Fourier transformation.                     !
!                                                                      !
!**********************************************************************!
Subroutine realft (data,m,n,mode)
  implicit none
  integer :: m, n, mode, k, j, nplan, iplan
  integer, parameter :: nmax = 256
  integer, parameter :: max_plan = 10
  integer*8 :: plana(max_plan),planb(max_plan)
  integer :: np(max_plan)
  double precision :: data(m,n)
  double precision :: copy(nmax,max_plan)
  double precision :: scale_f(max_plan)
  data np /max_plan * 0/
  data nplan /0/
  save copy,scale_f,plana,planb,np,nplan

  iplan = 0
  do j = 1,nplan
     if (n.eq.np(j)) then
        iplan = j
     endif
  enddo

  if (iplan .eq. 0) then
     nplan = nplan + 1
     iplan = nplan
     if (n .gt. nmax) stop 'realft 1'
     if (iplan .gt. max_plan) stop 'realft 2'
     scale_f(iplan) = dsqrt(1.d0/dble(n))
     call dfftw_plan_r2r_1d (plana(iplan),n,copy(1,iplan),copy(1,iplan),0,64)
     call dfftw_plan_r2r_1d (planb(iplan),n,copy(1,iplan),copy(1,iplan),1,64)
     np(iplan) = n
  endif
  do k = 1,m
     do j = 1,n
        copy(j,iplan) = data(k,j)
     enddo
     if (mode .eq. 1) then
        call dfftw_execute (plana(iplan))
     else if (mode .eq. -1) then
        call dfftw_execute (planb(iplan))
     else 
        stop 'realft 3'
     endif
     do j = 1,n
        data(k,j) = scale_f(iplan)*copy(j,iplan)
     enddo
  enddo
  return
end Subroutine realft