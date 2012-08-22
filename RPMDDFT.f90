#ifdef CP2K_BINDING
subroutine RPMDDFT_force(r,dvdr,na,nb,v,vir,boxlxyz,bead)
  use f_env
  implicit none
  include 'globals.inc'
#ifdef PARALLEL_BINDING
  include 'mpif.h' !parallel
#endif
  ! ------------------------------------------------------------------
  ! RPMD-DFT force modul
  ! ------------------------------------------------------------------

  integer na,nb,bead,ierr,k,j,jj,kk,rctdk,rpmddft,nbdf3,myid

  real(8) r(3,na),dvdr(3,na),vir(3,3),rtest(3,na),boxlxyz(3),cell(3,3)
  real(8) v,dx1,dy1,dz1
  character(len=3) ens

  common /ensemble/ ens

  common /RPMDDFT/ rpmddft,nbdf3,rctdk

  v = 0.d0
  vir(:,:) = 0.d0
  dvdr(:,:) = 0.d0
  cell(:,:) = 0.d0
#ifdef PARALLEL_BINDING
	call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr)
!write(*,*) "myid in RPMDDFTFORCE:", myid
	call MPI_bcast(r,SIZE(r),MPI_real8,0,MPI_COMM_WORLD,ierr)
	call MPI_bcast(boxlxyz,SIZE(boxlxyz),MPI_real8,0,MPI_COMM_WORLD,ierr)
#endif
  ! If Barostat is used set new cell
  if(ens.eq."NPT") then  
    cell(1,1) = boxlxyz(1)
    cell(2,2) = boxlxyz(2)
    cell(3,3) = boxlxyz(3)
    call cp_set_cell(f_env_id(bead),cell,ierr) !funktioniert und wurde getestet
    if (ierr.ne.0) STOP "set_cell"
  endif
  ! If rctdk ....
  !if(rctdk.eq.1) then
  !  cell(1,1) = boxlxyz(1)
  !  cell(2,2) = boxlxyz(2)
  !  cell(3,3) = boxlxyz(3)
  !  call cp_set_cell(f_env_id(bead),cell,ierr) !funktioniert und wurde getestet
  !  if (ierr.ne.0) STOP "set_cell"
  !endif

  ! Set Positions in CP2K
  call cp_set_pos(f_env_id(bead),r,SIZE(r),ierr)
  if (ierr.ne.0) STOP "set_pos"
!  write(*,*) "bead = ", bead, "f_env_id =", f_env_id(bead)

  ! Calculate new energy and force
  call cp_calc_energy_force(f_env_id(bead), r, SIZE(r), v, dvdr, SIZE(dvdr), ierr)
  if (ierr.ne.0) STOP "calc_energy_force"
#ifdef PARALLEL_BINDING
if(myid.eq.0) then
#endif
  ! Get force and energy
  call cp_get_force(f_env_id(bead),dvdr,SIZE(dvdr),ierr)
  if (ierr.ne.0) STOP "get_force"
  call cp_get_energy(f_env_id(bead),v,ierr)
  if (ierr.ne.0) STOP "get_energy"
  ! Get virial
  call cp_get_virial(f_env_id(bead),vir,SIZE(vir)*na,ierr)  ! gibt vir !!!!!!![vorher] in Bar zurück
  if (ierr.ne.0) STOP "get_virial"

  ! transform vir, because force has to be transformed to dvdr
  vir(:,:) = -vir(:,:) 
  
  ! transform force to dvdr
  dvdr(:,:) = -dvdr(:,:)
#ifdef PARALLEL_BINDING
endif
	call MPI_bcast(dvdr,SIZE(dvdr),MPI_real8,0,MPI_COMM_WORLD,ierr)
	call MPI_bcast(vir,SIZE(vir),MPI_real8,0,MPI_COMM_WORLD,ierr)
#endif
  return 
end subroutine RPMDDFT_force
#endif
