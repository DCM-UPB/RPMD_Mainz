subroutine pressure(pres,vir,tv,na,boxlxyz)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! Pressure calculation including tail correction
	! 01.02.2012 If RPMD-DFT then get pressure from CP2K and convert it 
	! 					 from Bar to internal units
  ! ------------------------------------------------------------------
  integer na, rpmddft,ierr
  real(8) vir(3,3),boxlxyz(3),tv,vol,pres,ptail,w
  real(8) oo_eps,oo_sig,oo_gam,rcut
  common /oo_param/ oo_eps,oo_sig,oo_gam,rcut
	common /RPMDDFT/ rpmddft

	if(rpmddft.eq.0) then
  	w = vir(1,1) + vir(2,2) + vir(3,3)
  	w = -w/3.d0
  	vol = boxlxyz(1)*boxlxyz(2)*boxlxyz(3)
  	pres = (2.d0*tv/3.d0 + w) / vol
  
  	if (oo_gam.eq.0.d0) then
     	call pres_lj_tail(ptail,boxlxyz,na)
  	else
     	ptail = 0.d0
  	endif
  
  	pres = pres + ptail
  else
			w = (vir(1,1) + vir(2,2) + vir(3,3)) 
  		w = -w/3.d0
  		vol = boxlxyz(1)*boxlxyz(2)*boxlxyz(3)
  		pres = (2.d0*tv/3.d0 + w) / vol
	endif

  return
end subroutine pressure

subroutine pres_lj_tail(ptail,boxlxyz,na)
  implicit none
  ! ------------------------------------------------------------------
  ! LJ tail correction to pressure
  ! ------------------------------------------------------------------
  integer na,nm
  real(8) ptail,boxlxyz(3),vol,pi,rho,prefac
  real(8) oo_eps,oo_sig,oo_gam,rcut
  common /oo_param/ oo_eps,oo_sig,oo_gam,rcut
  
  nm = na/3
  pi = dacos(-1.d0)
  vol = boxlxyz(1)*boxlxyz(2)*boxlxyz(3)
  rho = dble(nm) / vol
  prefac = (16.d0*pi*(rho**2)*oo_eps*(oo_sig**3)) / (3.d0)
  ptail  = prefac*( (2.d0/3.d0)*(oo_sig/rcut)**9 &
       - (oo_sig/rcut)**3)
  
  return
end subroutine pres_lj_tail
