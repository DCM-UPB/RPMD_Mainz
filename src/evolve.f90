! ------------------------------------------------------------------
!<summary>
!   Driver routine for evolution
!</summary>
!<param name="p"> the momentum of all particles and beads</param>
!<param name="r"> the coordinates of all particles and beads</param>
!<param name="mass"> the masses of all atoms </param>
!<param name="na">the number of atoms</param>
!<param name="nb">the number of beads (per atom?)</param>
!<param name="boxlxyz">The length of the box?</param>
!<param name=beta"> the inverse temperature of the system</param>
! ------------------------------------------------------------------
subroutine evolve(p,r,v,v1,v2,v3,dvdr,dvdr2,dt,mass,na,nb, &
                  boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
    implicit none
    integer,intent(in)    :: na,nb
    integer irun,nbaro,rpmddft,nbdf3,rctdk
    real(8),intent(inout) :: p(3,na,nb),r(3,na,nb)
    real(8) dvdr(3,na,nb),dvdr2(3,na,nb)
    real(8) mass(na),z(na),boxlxyz(3),vir(3,3),vir_lf(3,3)
    real(8) v,v1,v2,v3,dt,beta
    logical iamrigid
    common /structure/ iamrigid
    common /RPMDDFT/ rpmddft,nbdf3,rctdk

    if (nb.eq.1) then
        if (rpmddft.eq.1) then
            if (iamrigid) then
                call evolve_rigid_cl_RPMDDFT(p,r,v,v1,v2,dvdr,dvdr2,dt,mass,&
                na,nb,boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
            else
                call evolve_cl_RPMDDFT(p,r,v,v1,v2,dvdr,dvdr2,dt,mass,na,nb, &
                boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
            end if
        else
            if (iamrigid) then
                call evolve_rig_cl(p,r,v,v1,v2,dvdr,dvdr2,dt,mass,na,nb, &
                boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
            else
                call evolve_cl(p,r,v,v1,v2,dvdr,dvdr2,dt,mass,na,nb, &
                boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
            end if
        end if
    else
        if (rpmddft.eq.1) then
            if (rctdk.eq.1) then
                call evolve_cl_pi_RPMDDFT(p,r,v,v1,v2,v3,dvdr,dvdr2,dt,mass,&
                na,nb,boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
            else      
                if (nbdf3.eq.0) then
                    if (iamrigid) then
                        call evolve_rigid_pi_RPMDDFT(p,r,v,v1,v2,dvdr,dvdr2,&
                        dt,mass,na,nb, boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
                    else
                        call evolve_pi_RPMDDFT(p,r,v,v1,v2,v3,dvdr,dvdr2,dt,&
                        mass,na,nb,boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
                    end if
                else
                    call evolve_pi_rc_RPMDDFT(p,r,v,v1,v2,v3,dvdr,dvdr2,dt,&
                    mass,na,nb,boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
                end if
            end if
        else
            if (iamrigid) then
                call evolve_rig_pi(p,r,v,v1,v2,dvdr,dvdr2,dt,mass,na,nb, &
                boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
            else
                call evolve_pi(p,r,v,v1,v2,v3,dvdr,dvdr2,dt,mass,na,nb, &
                boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
            end if
        end if
    end if  
    return
end subroutine evolve
