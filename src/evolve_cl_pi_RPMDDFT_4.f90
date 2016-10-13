subroutine evolve_cl_pi_RPMDDFT_4(p,r,v,vew,vlj,vint,dvdr,dvdr2,dt,mass,na,nb, &
                     boxlxyz,z,beta,vir,vir_lf,irun,nbaro)
  use thermostat
  use barostat
  use gle
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! RPMD / RPMDFT evolution
  ! 
  ! ------------------------------------------------------------------
  integer na,nb,irun,k,j,i,nbdf1,nbdf2,im,ic,mts,nbaro,rpmddft,nbdf3
  real(8) p(3,na,nb),r(3,na,nb),dvdr(3,na,nb),dvdr2(3,na,nb),dvdre(3,na,nb),dvdrl(3,na,nb)
  real(8) mass(na),z(na),boxlxyz(3)
  real(8) vir(3,3),vir_lj(3,3),vir_ew(3,3),tvxyz(3),tvxyz_itr(3)
  real(8) vir_itr(3,3),vir_ewc(3,3),vir_lf(3,3),vir_hf(3,3),virMM(3,3),virCP2K(3,3),virMMh(3,3),virMMl(3,3)
  real(8) pprime,halfdtsmall,halfdt,om,gaussian
  real(8) dt,v,beta,dtsmall,vew,vlj,vint,sig,ve,vCP2K,vMM,vMMh,vMMl
  real(8) tv,tv_itr,tq1,tq2,c1,c2,c3
  real(8) dheat, comx, comy, comz, mm !!GLE
  character(len=4) type
  external gaussian
  common /path_i/ om,type
  common /beaddiabatic/ nbdf1,nbdf2
  common /multiple_ts/ mts
  common /correct/ sig
  common /RPMDDFT/ rpmddft,nbdf3
  real(8) rnm(3,na,nb),rb(3,na,nbdf3),dvdrCP2K(3,na,nb),dvdrMM(3,na,nb),dvdrb(3,na,nbdf3),dvdrMMl(3,na,nb),dvdrMMh(3,na,nb)


  ! Zero constants and arrays
  dvdrl(:,:,:) = 0.d0
  dvdre(:,:,:) = 0.d0
  dvdrMM(:,:,:) = 0.d0
  dvdrMMh(:,:,:) = 0.d0
  dvdrMMl(:,:,:) = 0.d0
  dvdrCP2K(:,:,:) = 0.d0
  dvdrb(:,:,:) = 0.d0
  rnm(:,:,:) = 0.d0
  rb(:,:,:) = 0.d0
  vCP2K = 0.d0
  vMM = 0.d0
  vMMh = 0.d0
  vMMl = 0.d0
  
  vint = 0.d0
  vlj = 0.d0
  vir_hf(:,:) = 0.d0
  vir_lj(:,:) = 0.d0

  vir_ew(:,:) = 0.d0
  vir_hf(:,:) = 0.d0
  virCP2K(:,:) = 0.d0
  virMM(:,:) = 0.d0
	virMMh(:,:) = 0.d0
	virMMl(:,:) = 0.d0
  tv = 0.d0
  tvxyz(:) = 0.d0

  halfdt = 0.5d0*dt

  if (type.eq.'RPMD') then
     if (therm.eq.'PRG') then
        call parinello_therm(p,mass,ttau,na,nb,halfdt,irun,beta)
     else if (therm.eq.'PRL') then
        call parinello_therm_loc(p,mass,ttau,na,nb, &
                                 dt,irun,beta)
     else if (therm.eq.'GLE') then     !! GLE
         call therm_gle(p,dheat,mass,na,nb,irun)
     endif
  endif

  p(:,:,:) = p(:,:,:) - halfdt*(dvdr(:,:,:)+dvdr2(:,:,:)) !hier sind alle Kräfte drin!

  !  get new coordinates
  call freerp_rpmd(p,r,dt,mass,na,nb,beta)

  ! Form Bead Positions and evaluate high-frequency forces
  call forces(r,vint,dvdr2,nb,na,boxlxyz,z,vir_hf,4)
  call virial_ke(r,dvdr,dvdr2,tv,tvxyz, &
                  tq1,tq2,beta,na,nb,mass)

  p(:,:,:) = p(:,:,:) - halfdt*dvdr2(:,:,:)
     
   vir(:,:) = vir_lf(:,:) + vir_hf(:,:)   
  ! Barostat
  ! (note:// COMs scaled therefore do not need to recalculate
  !  intramolecular forces as they remain the same)

  if (nbaro.eq.1) then
     if (baro.eq.'BER') then
        call beren_driver(vir,tv,tvxyz,dt,r,boxlxyz,na,nb)
     else if (baro.eq.'MCI') then
        call mc_baro(r,dvdr,dvdr2,vir,v,z,beta,boxlxyz,na,nb,irun)
     else 
        write(6,*) ' ** Invalid barostat **'
        stop
     endif
  endif

  ! Evaluation of the low frequency forces
  ! --------------------------------------------------------------
	call forces(r,vlj,dvdr,nb,na,boxlxyz,z,vir_lj,1)

 
  ! Evaluate CP2K and MM, nb=nbdf3 force:
  if (nbdf3.eq.0) then 
    write(6,*) '* nbdf3 has to be nonzero if cl_pi is used'
    stop
  else

  ! Form contracted RP coordinates
    ! Convert the positions to the normal mode representation
    rnm(:,:,:) = r(:,:,:)
    call realft (rnm,3*na,nb,+1)

    call ring_contract(rnm,rb,na,nb,nbdf3)

  ! Evaluate forces on the bead coordinates rb

    call forces(rb,vCP2K,dvdrb,nbdf3,na,boxlxyz,z,virCP2K,9)
    call force_contract(dvdrCP2K,dvdrb,na,nb,nbdf3)
    dvdrb(:,:,:) = 0.d0


    !call forces(rb,vMM,dvdrb,nbdf3,na,boxlxyz,z,virMM,0) !PRoblem, das bei full forces hier RPMDDFT_forece aufgerufen wird, da rpmddft =1
		call forces(rb,vMMl,dvdrb,nbdf3,na,boxlxyz,z,virMMl,1)
		call force_contract(dvdrMMl,dvdrb,na,nb,nbdf3)
    dvdrb(:,:,:) = 0.d0
		call forces(rb,vMMh,dvdrb,nbdf3,na,boxlxyz,z,virMMh,4)
		call force_contract(dvdrMMh,dvdrb,na,nb,nbdf3)
		dvdrMM(:,:,:) = dvdrMMl(:,:,:)+dvdrMMh(:,:,:)
		virMM(:,:) = virMMl(:,:) + virMMh(:,:)
		vMM = vMMl + vMMh
  endif

  ! Convert dvdrCP2K and dvdrMM back to bead representation
  call realft(dvdrCP2K,3*na,nb,-1)
  call realft(dvdrMM,3*na,nb,-1)


	!write(*,*) "dvdr:", dvdr(:,1,3)
	!write(*,*) "dvdrCP2K:", dvdrCP2K(:,1,3)
	!write(*,*) "dvdrMM:", dvdrMM(:,1,3)
	!write(*,*) "dvdrCP2K-dvdrMM", dvdrCP2K(:,1,3)-dvdrMM(:,1,3)
	!write(*,*) "dvdr+dvdrCP2K-dvdrMM", dvdr(:,1,3)+dvdrCP2K(:,1,3)-dvdrMM(:,1,3)
	dvdr(:,:,:) = dvdr(:,:,:)+dvdrCP2K(:,:,:)-dvdrMM(:,:,:)

  ! Evolve the momenta under the low+CP2K+MM forces

  p(:,:,:) = p(:,:,:) - halfdt*dvdr(:,:,:)


  
  if (type.eq.'RPMD') then
     if (therm.eq.'PRG') then
        call parinello_therm(p,mass,ttau,na,nb,halfdt,irun,beta)
     else if (therm.eq.'PRL') then
        call parinello_therm_loc(p,mass,ttau,na,nb, &
                                       dt,irun,beta)
     else if (therm.eq.'GLE') then     !! GLE
         call therm_gle(p,dheat,mass,na,nb,irun)
     endif
  endif

  ! Potential Energy

  v = vlj + vint - vMM + vCP2K ! ??????????????????

  ! Virial

  vir_lf(:,:) = vir_lj(:,:) + virCP2K(:,:) - virMM(:,:)    !????
  vir(:,:) =  vir_lf(:,:) + vir_hf(:,:)  


  !compute COM and outputs to make sure that we removed COM velocity correctly....
  comx=0.d0
  comy=0.d0
  comz=0.d0
  mm=0.0d0
  do j = 1,na
    do k = 1,nb
      comx=comx+mass(j)*r(1,j,k)
      comy=comy+mass(j)*r(2,j,k)
      comz=comz+mass(j)*r(3,j,k)
    enddo
    mm=mm+mass(j)
  enddo
  comx=comx/mm
  comy=comy/mm
  comz=comz/mm
!  write(*,*) "COM: ",comx,comy,comz
  return
end subroutine evolve_cl_pi_RPMDDFT_4
