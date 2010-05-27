function gaussian(idum,sigma)
  implicit none
  ! ------------------------------------------------------------------
  ! Function returns a random number from Gaussian distribution with
  ! mean zero and standard deviation sigma.
  !  ------------------------------------------------------------------
  integer iset, idum
  real(8) gset, sigma,v1,v2,fac,gaussian,rsq,ran2
  save iset,gset
  data iset/0/
  external ran2
  if (iset.eq.0) then
1    v1=2.d0*ran2(idum,0.d0,1.d0)-1.d0
     v2=2.d0*ran2(idum,0.d0,1.d0)-1.d0
     rsq=v1**2+v2**2  
     if(rsq.ge.1.d0.or.rsq.eq.0.d0)goto 1 
     fac=sqrt(-2.d0*log(rsq)/rsq) 
     gset=v1*fac
     gaussian=sigma * v2 * fac
     iset=1
  else 
     gaussian=sigma * gset
     iset=0 
  endif
  return
end function gaussian

function ran2(idum,rmin,rmax)
  implicit none
  ! ------------------------------------------------------------------
  ! Long period random number generator of L'Ecuyer with Bays-Durham shuffle
  ! and added safeguards. Returns a uniform random deviate between
  ! 0 .0 and 1.0 (exclusive of the endpoint values).
  ! Call with idum a negative integer to initialize; thereafter, do not
  ! alter idum between successive deviates in a sequence.
  ! rnmx should approximate the largest floating value that is less than 1.
  ! ------------------------------------------------------------------
  integer im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
  parameter ( im1=2147483563,im2=2147483399,imm1=im1-1,   &
       ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211, & 
       ir2=3791,ntab=32,ndiv=1+imm1/ntab )
  real(8) am, eps, rnmx
  parameter ( am=1./im1, eps=1.2e-7,rnmx=1.-eps )
  real(8) ran2,rmin,rmax
  integer iv(ntab),j,k,iy,idum,idum2
  save iv,iy,idum2
  data idum2/123456789/, iv/ntab*0/, iy/0/
    
  if(idum.le.0) then 
     idum=max(-idum,1) 
     idum2=idum
     do j=ntab+8,1,-1 
        k=idum/iq1
        idum=ia1*(idum-k*iq1)-k*ir1
        if (idum.lt.0) idum=idum+im1
        if (j.le.ntab) iv(j)=idum
     enddo
     iy=iv(1)
  endif
  k=idum/iq1 
  idum=ia1*(idum-k*iq1)-k*ir1 
  if(idum.lt.0) idum=idum+im1
  k=idum2/iq2
  idum2=ia2*(idum2-k*iq2)-k*ir2 
  if(idum2.lt.0) idum2=idum2+im2
  j=1+iy/ndiv 
  iy=iv(j)-idum2 
  iv(j)=idum
  if(iy.lt.1)iy=iy+imm1
  ran2=((rmax-rmin)*(min(am*iy,rnmx)))+rmin      
  return
end function ran2

