subroutine evolve(p,r,v,v1,v2,v3,dvdr,dvdr2,dt,mass,na,nb, &
                  boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
  implicit none
  ! ------------------------------------------------------------------
  ! Driver routine for evolution
  ! ------------------------------------------------------------------
  integer na,nb,irun,nbaro,rpmddft,nbdf3,rctdk
  real(8) p(3,na,nb),r(3,na,nb),dvdr(3,na,nb),dvdr2(3,na,nb)
  real(8) mass(na),z(na),boxlxyz(3),vir(3,3),vir_lf(3,3)
  real(8) v,v1,v2,v3,dt,beta
  logical iamrigid
  common /structure/ iamrigid
  common /RPMDDFT/ rpmddft,nbdf3,rctdk

  if (nb.eq.1) then
    if (rpmddft.eq.1) then
      call evolve_cl_RPMDDFT(p,r,v,v1,v2,dvdr,dvdr2,dt,mass,na,nb, &
                            boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
    else
      if (iamrigid) then
         call evolve_rig_cl(p,r,v,v1,v2,dvdr,dvdr2,dt,mass,na,nb, &
                            boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
      else
         call evolve_cl(p,r,v,v1,v2,dvdr,dvdr2,dt,mass,na,nb, &
                       boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
      endif
    endif
  else
    if (rpmddft.eq.1) then
      if (rctdk.eq.1) then
        call evolve_cl_pi_RPMDDFT(p,r,v,v1,v2,v3,dvdr,dvdr2,dt,mass,na,nb, &
                           boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
      else      
        if (nbdf3.eq.0) then
					if (iamrigid) then
          call evolve_rigid_pi_RPMDDFT(p,r,v,v1,v2,dvdr,dvdr2,dt,mass,na,nb, &
                           boxlxyz,z,beta,vir,vir_lf,irun,nbaro)					
					else
          call evolve_pi_RPMDDFT(p,r,v,v1,v2,v3,dvdr,dvdr2,dt,mass,na,nb, &
                           boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
					endif
        else
          call evolve_pi_rc_RPMDDFT(p,r,v,v1,v2,v3,dvdr,dvdr2,dt,mass,na,nb, &
                           boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
        endif
      endif
    else
     if (iamrigid) then
        call evolve_rig_pi(p,r,v,v1,v2,dvdr,dvdr2,dt,mass,na,nb, &
                           boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
     else
        call evolve_pi(p,r,v,v1,v2,v3,dvdr,dvdr2,dt,mass,na,nb, &
                       boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
     endif
    endif
  endif  

  return
end subroutine evolve
