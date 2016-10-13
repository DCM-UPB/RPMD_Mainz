subroutine epsr_run(r,boxlxyz)
  implicit none
  include "globals.inc"
  ! ------------------------------------------------------------------
  ! Calculate EPSR correction
  ! ------------------------------------------------------------------
  real(8) r(3,ina),boxlxyz(3)
  character(len=255) :: cwd,cmd
  integer epsr_update
  logical epsr
  real(8) pos(5000),potOO(5000),frcOO(5000),potOH(5000),frcOH(5000),potHH(5000),frcHH(5000)
  common /EPSR/ epsr, epsr_update, pos, potOO, frcOO, potOH, frcOH, potHH, frcHH

  !write(6,*) "Calculating EPSR correction"

  ! Save current box in subfolder
  open(123456,file='EPSRrun/vmd_current.xyz',STATUS='replace')
  call print_vmd_bead(r,inb,1,ina,inm,boxlxyz,123456)
  close(123456)

  ! Run epsr
  call chdir("EPSRrun")
  call getcwd(cwd)
  ! TODO Add global bool to /EPSR/ to we can turn on/off /dev/null on the fly
  write(cmd, '(A,A,A)') 'readxyz ', trim(cwd), '/ readxyz vmd_current.xyz 0.65 > /dev/null'
  call system(cmd)
  write(cmd, '(A,A,A)') 'epsr ', trim(cwd), '/ epsr vmd_current.EPSR.inp > /dev/null'
  call system(cmd)
  call chdir("..")
end subroutine


subroutine epsr_read()
  implicit none
  include "globals.inc"
  real(8) diff
  character(len=255) :: line
  integer i, io_err
  integer epsr_update
  logical epsr, file_exists
  real(8) pos(5000),potOO(5000),frcOO(5000),potOH(5000),frcOH(5000),potHH(5000),frcHH(5000)
  common /EPSR/ epsr, epsr_update, pos, potOO, frcOO, potOH, frcOH, potHH, frcHH
  ! Get results
  potOO(:) = 0.0d0
  potOH(:) = 0.0d0
  potHH(:) = 0.0d0
  frcOO(:) = 0.0d0
  frcOH(:) = 0.0d0
  frcHH(:) = 0.0d0


#ifdef EPSR_STARTUP_READ
  inquire(FILE="vmd_current.EPSR.p01", EXIST=file_exists)
  if (file_exists.eqv..true.) then
    write(6,*) "Reading EPSR correction file now"
    open(123456, file="vmd_current.EPSR.p01", action="read")
  else
    write(6,*) "EPSR correction file vmd_current.EPSR.p01 not found!"
    write(6,*) "Setting everything to zero..."
    do i=1,5000
        pos(i) = i
    enddo
    return
  endif
#else
  open(123456, file="EPSRrun/vmd_current.EPSR.p01", action="read")
#endif
  read(123456,*) line
  do i=1,5000
    read(123456,*,iostat=io_err) pos(i),potOO(i),frcOO(i),potOH(i),frcOH(i),potHH(i),frcHH(i)
    if(io_err.ne.0) then
      exit
    endif
  enddo
  do i=2,5000
    pos(i-1) = pos(i)
    potOO(i-1) = potOO(i)
    potOH(i-1) = potOH(i)
    potHH(i-1) = potHH(i)
    frcOO(i-1) = frcOO(i)
    frcOH(i-1) = frcOH(i)
    frcHH(i-1) = frcHH(i)
  enddo
  close(123456)
  potOO(:) = potOO(:)/toKjmol
  potOH(:) = potOH(:)/toKjmol
  potHH(:) = potHH(:)/toKjmol
  pos(:) = pos(:)/toA

  ! calculate force as it is not (yet) implemented in EPSR
  do i = 2,4999
    diff = (pos(i+1)-pos(i-1))
    if (diff.gt.1.0d-5.or.diff.lt.-1.0d-5) then
      frcOO(i) = -(potOO(i+1)-potOO(i-1))/diff
      frcOH(i) = -(potOH(i+1)-potOH(i-1))/diff
      frcHH(i) = -(potHH(i+1)-potHH(i-1))/diff
    endif
  enddo

  open(123456, file="internal_epsr_potential.dat", action="write")
  do i=1,5000
    write(123456,*) pos(i),potOO(i),frcOO(i),potOH(i),frcOH(i),potHH(i),frcHH(i)
  enddo
  close(123456)

end subroutine epsr_read

subroutine epsr_driver(r,dvdr,v,vir,list,point,na,boxlxyz,njump)
  implicit none
  include 'globals.inc'
  ! ------------------------------------------------------------------
  ! EPSR Driver
  ! ------------------------------------------------------------------
  integer na,point(na+3),list(maxnab*na),njump
  real(8) r(3,na),dvdr(3,na),vir(3,3),boxlxyz(3)
  real(8) v,oo_eps,oo_sig,oo_gam,rcut,boxmax
  logical epsr
  integer epsr_update
  real(8) pos(5000),potOO(5000),frcOO(5000),potOH(5000),frcOH(5000),potHH(5000),frcHH(5000)
  common /EPSR/ epsr, epsr_update, pos, potOO, frcOO, potOH, frcOH, potHH, frcHH
  common /oo_param/ oo_eps,oo_sig,oo_gam,rcut

#ifdef EPSR_STARTUP_READ
  ! Only read once on startup
  if (istep.eq.0) then
    call epsr_read()
  else if (mod(istep,epsr_update).eq.0) then
    call epsr_run(r,boxlxyz)
  endif
#else
  if (mod(istep,epsr_update).eq.0) then
    call epsr_run(r,boxlxyz)
    call epsr_read()
  endif
#endif

  dvdr(:,:) = 0.0d0
  v = 0.0d0
  vir(:,:) = 0.0d0

  ! Cubic

  !if (list(1).ne.0) then
  !   call epsr_list(r,dvdr,v,vir,na,boxlxyz,list,point,njump)
  !else
  !   boxmax = maxval(boxlxyz)
  !   if (3.d0*rcut .gt. boxmax) then
        call epsr_basic(r,dvdr,v,vir,na,boxlxyz,njump)
  !   else
  !
  !      ! Use linked cell list
  !
  !      call epsr_cell(r,v,vir,dvdr,na,boxlxyz,njump)
  !   endif
  !endif
  
end subroutine epsr_driver


subroutine epsr_basic(r,dvdr,v,vir,na,boxlxyz,njump)
  implicit none
  ! ----------------------------------------------------------------
  ! EPSR - General Version
  ! ----------------------------------------------------------------
  integer na,i,j,njump, bin
  real(8) r(3,na),dvdr(3,na),vir(3,3),boxlxyz(3),v
  real(8) onboxx,onboxy,onboxz
  real(8) boxx,boxy,boxz
  real(8) drsq,dfx,dfy,dfz
  real(8) dx,dy,dz,sq
  logical epsr
  integer epsr_update, which_H
  real(8) pos(5000),potOO(5000),frcOO(5000),potOH(5000),frcOH(5000),potHH(5000),frcHH(5000)
  common /EPSR/ epsr, epsr_update, pos, potOO, frcOO, potOH, frcOH, potHH, frcHH

  boxx = boxlxyz(1)
  boxy = boxlxyz(2)
  boxz = boxlxyz(3)
  onboxx = 1.d0/boxx
  onboxy = 1.d0/boxy
  onboxz = 1.d0/boxz

  v = 0.d0
  vir(:,:) = 0.d0
  dvdr(:,:) = 0.d0

  do j = 1+njump,na,njump
     do i = 1,j-njump,njump
        ! O--O
        dx = r(1,i)-r(1,j)
        dy = r(2,i)-r(2,j)
        dz = r(3,i)-r(3,j)
        dx = dx - boxx*dble(nint(onboxx*dx))
        dy = dy - boxy*dble(nint(onboxy*dy))
        dz = dz - boxz*dble(nint(onboxz*dz))
        drsq = dx*dx + dy*dy + dz*dz
        sq = sqrt(drsq)
        bin = int(sq/(pos(5)-pos(4)))
        if (bin .lt. 5000 .and. bin.gt.0) then
           v = v + potOO(bin)
           dfx = frcOO(bin) * dx*dx/drsq
           dfy = frcOO(bin) * dy*dy/drsq
           dfz = frcOO(bin) * dz*dz/drsq
           !write(6,*) bin, frcOO(bin), dfx,dfy,dfz
           !write(6,*) dx, dy, dz, sq
           dvdr(1,i) = dvdr(1,i) - dfx
           dvdr(2,i) = dvdr(2,i) - dfy
           dvdr(3,i) = dvdr(3,i) - dfz
           dvdr(1,j) = dvdr(1,j) + dfx
           dvdr(2,j) = dvdr(2,j) + dfy
           dvdr(3,j) = dvdr(3,j) + dfz
           vir(1,1) = vir(1,1) - dx * dfx
           vir(2,2) = vir(2,2) - dy * dfy
           vir(3,3) = vir(3,3) - dz * dfz
           vir(1,2) = vir(1,2) - dx * dfy
           vir(1,3) = vir(1,3) - dx * dfz
           vir(2,1) = vir(2,1) - dy * dfx
           vir(2,3) = vir(2,3) - dy * dfz
           vir(3,1) = vir(3,1) - dz * dfx
           vir(3,2) = vir(3,2) - dz * dfy
        endif

        ! O--H
        do which_H=1,2
          dx = r(1,i)-r(1,j+which_H)
          dy = r(2,i)-r(2,j+which_H)
          dz = r(3,i)-r(3,j+which_H)
          dx = dx - boxx*dble(nint(onboxx*dx))
          dy = dy - boxy*dble(nint(onboxy*dy))
          dz = dz - boxz*dble(nint(onboxz*dz))
          drsq = dx*dx + dy*dy + dz*dz
          sq = sqrt(drsq)
          bin = int(sq/(pos(5)-pos(4)))
          if (bin .lt. 5000 .and. bin.gt.0) then
             v = v + potOH(bin)
             dfx = frcOH(bin) * dx*dx/drsq
             dfy = frcOH(bin) * dy*dy/drsq
             dfz = frcOH(bin) * dz*dz/drsq
             !write(6,*) bin, frcOO(bin), dfx,dfy,dfz
             !write(6,*) dx, dy, dz, sq
             dvdr(1,i) = dvdr(1,i) - dfx
             dvdr(2,i) = dvdr(2,i) - dfy
             dvdr(3,i) = dvdr(3,i) - dfz
             dvdr(1,j+which_H) = dvdr(1,j+which_H) + dfx
             dvdr(2,j+which_H) = dvdr(2,j+which_H) + dfy
             dvdr(3,j+which_H) = dvdr(3,j+which_H) + dfz
             vir(1,1) = vir(1,1) - dx * dfx
             vir(2,2) = vir(2,2) - dy * dfy
             vir(3,3) = vir(3,3) - dz * dfz
             vir(1,2) = vir(1,2) - dx * dfy
             vir(1,3) = vir(1,3) - dx * dfz
             vir(2,1) = vir(2,1) - dy * dfx
             vir(2,3) = vir(2,3) - dy * dfz
             vir(3,1) = vir(3,1) - dz * dfx
             vir(3,2) = vir(3,2) - dz * dfy
          endif

          dx = r(1,i+which_H)-r(1,j)
          dy = r(2,i+which_H)-r(2,j)
          dz = r(3,i+which_H)-r(3,j)
          dx = dx - boxx*dble(nint(onboxx*dx))
          dy = dy - boxy*dble(nint(onboxy*dy))
          dz = dz - boxz*dble(nint(onboxz*dz))
          drsq = dx*dx + dy*dy + dz*dz
          sq = sqrt(drsq)
          bin = int(sq/(pos(5)-pos(4)))
          if (bin .lt. 5000 .and. bin.gt.0) then
             v = v + potOH(bin)
             dfx = frcOH(bin) * dx*dx/drsq
             dfy = frcOH(bin) * dy*dy/drsq
             dfz = frcOH(bin) * dz*dz/drsq
             !write(6,*) bin, frcOH(bin), dfx,dfy,dfz
             !write(6,*) dx, dy, dz, sq
             dvdr(1,i+which_H) = dvdr(1,i+which_H) - dfx
             dvdr(2,i+which_H) = dvdr(2,i+which_H) - dfy
             dvdr(3,i+which_H) = dvdr(3,i+which_H) - dfz
             dvdr(1,j) = dvdr(1,j) + dfx
             dvdr(2,j) = dvdr(2,j) + dfy
             dvdr(3,j) = dvdr(3,j) + dfz
             vir(1,1) = vir(1,1) - dx * dfx
             vir(2,2) = vir(2,2) - dy * dfy
             vir(3,3) = vir(3,3) - dz * dfz
             vir(1,2) = vir(1,2) - dx * dfy
             vir(1,3) = vir(1,3) - dx * dfz
             vir(2,1) = vir(2,1) - dy * dfx
             vir(2,3) = vir(2,3) - dy * dfz
             vir(3,1) = vir(3,1) - dz * dfx
             vir(3,2) = vir(3,2) - dz * dfy
          endif
        enddo ! O--H

!!!        ! H--H on same molecule
!!!        dx = r(1,i+1)-r(1,i+2)
!!!        dy = r(2,i+1)-r(2,i+2)
!!!        dz = r(3,i+1)-r(3,i+2)
!!!        dx = dx - boxx*dble(nint(onboxx*dx))
!!!        dy = dy - boxy*dble(nint(onboxy*dy))
!!!        dz = dz - boxz*dble(nint(onboxz*dz))
!!!        drsq = dx*dx + dy*dy + dz*dz
!!!        sq = sqrt(drsq)
!!!        bin = int(sq/(pos(5)-pos(4)))
!!!        if (bin .lt. 5000 .and. bin.gt.0) then
!!!           v = v + potOO(bin)
!!!           dfx = frcHH(bin) * dx*dx/drsq
!!!           dfy = frcHH(bin) * dy*dy/drsq
!!!           dfz = frcHH(bin) * dz*dz/drsq
!!!           !write(6,*) bin, frcHH(bin), dfx,dfy,dfz
!!!           !write(6,*) dx, dy, dz, sq
!!!           dvdr(1,i+1) = dvdr(1,i+1) - dfx
!!!           dvdr(2,i+1) = dvdr(2,i+1) - dfy
!!!           dvdr(3,i+1) = dvdr(3,i+1) - dfz
!!!           dvdr(1,i+2) = dvdr(1,i+2) + dfx
!!!           dvdr(2,i+2) = dvdr(2,i+2) + dfy
!!!           dvdr(3,i+2) = dvdr(3,i+2) + dfz
!!!           vir(1,1) = vir(1,1) - dx * dfx
!!!           vir(2,2) = vir(2,2) - dy * dfy
!!!           vir(3,3) = vir(3,3) - dz * dfz
!!!           vir(1,2) = vir(1,2) - dx * dfy
!!!           vir(1,3) = vir(1,3) - dx * dfz
!!!           vir(2,1) = vir(2,1) - dy * dfx
!!!           vir(2,3) = vir(2,3) - dy * dfz
!!!           vir(3,1) = vir(3,1) - dz * dfx
!!!           vir(3,2) = vir(3,2) - dz * dfy
!!!        endif
!!!        dx = r(1,j+1)-r(1,j+2)
!!!        dy = r(2,j+1)-r(2,j+2)
!!!        dz = r(3,j+1)-r(3,j+2)
!!!        dx = dx - boxx*dble(nint(onboxx*dx))
!!!        dy = dy - boxy*dble(nint(onboxy*dy))
!!!        dz = dz - boxz*dble(nint(onboxz*dz))
!!!        drsq = dx*dx + dy*dy + dz*dz
!!!        sq = sqrt(drsq)
!!!        bin = int(sq/(pos(5)-pos(4)))
!!!        if (bin .lt. 5000 .and. bin.gt.0) then
!!!           v = v + potOO(bin)
!!!           dfx = frcHH(bin) * dx*dx/drsq
!!!           dfy = frcHH(bin) * dy*dy/drsq
!!!           dfz = frcHH(bin) * dz*dz/drsq
!!!           !write(6,*) bin, frcHH(bin), dfx,dfy,dfz
!!!           !write(6,*) dx, dy, dz, sq
!!!           dvdr(1,j+1) = dvdr(1,j+1) - dfx
!!!           dvdr(2,j+1) = dvdr(2,j+1) - dfy
!!!           dvdr(3,j+1) = dvdr(3,j+1) - dfz
!!!           dvdr(1,j+2) = dvdr(1,j+2) + dfx
!!!           dvdr(2,j+2) = dvdr(2,j+2) + dfy
!!!           dvdr(3,j+2) = dvdr(3,j+2) + dfz
!!!           vir(1,1) = vir(1,1) - dx * dfx
!!!           vir(2,2) = vir(2,2) - dy * dfy
!!!           vir(3,3) = vir(3,3) - dz * dfz
!!!           vir(1,2) = vir(1,2) - dx * dfy
!!!           vir(1,3) = vir(1,3) - dx * dfz
!!!           vir(2,1) = vir(2,1) - dy * dfx
!!!           vir(2,3) = vir(2,3) - dy * dfz
!!!           vir(3,1) = vir(3,1) - dz * dfx
!!!           vir(3,2) = vir(3,2) - dz * dfy
!!!        endif

        ! H--H on different molecule
        do which_H=1,2
          dx = r(1,i+1)-r(1,j+which_H)
          dy = r(2,i+1)-r(2,j+which_H)
          dz = r(3,i+1)-r(3,j+which_H)
          dx = dx - boxx*dble(nint(onboxx*dx))
          dy = dy - boxy*dble(nint(onboxy*dy))
          dz = dz - boxz*dble(nint(onboxz*dz))
          drsq = dx*dx + dy*dy + dz*dz
          sq = sqrt(drsq)
          bin = int(sq/(pos(5)-pos(4)))
          if (bin .lt. 5000 .and. bin.gt.0) then
             v = v + potHH(bin)
             dfx = frcHH(bin) * dx*dx/drsq
             dfy = frcHH(bin) * dy*dy/drsq
             dfz = frcHH(bin) * dz*dz/drsq
             !write(6,*) bin, frcHH(bin), dfx,dfy,dfz
             !write(6,*) dx, dy, dz, sq
             dvdr(1,i+1) = dvdr(1,i+1) - dfx
             dvdr(2,i+1) = dvdr(2,i+1) - dfy
             dvdr(3,i+1) = dvdr(3,i+1) - dfz
             dvdr(1,j+which_H) = dvdr(1,j+which_H) + dfx
             dvdr(2,j+which_H) = dvdr(2,j+which_H) + dfy
             dvdr(3,j+which_H) = dvdr(3,j+which_H) + dfz
             vir(1,1) = vir(1,1) - dx * dfx
             vir(2,2) = vir(2,2) - dy * dfy
             vir(3,3) = vir(3,3) - dz * dfz
             vir(1,2) = vir(1,2) - dx * dfy
             vir(1,3) = vir(1,3) - dx * dfz
             vir(2,1) = vir(2,1) - dy * dfx
             vir(2,3) = vir(2,3) - dy * dfz
             vir(3,1) = vir(3,1) - dz * dfx
             vir(3,2) = vir(3,2) - dz * dfy
          endif

          dx = r(1,i+2)-r(1,j+which_H)
          dy = r(2,i+2)-r(2,j+which_H)
          dz = r(3,i+2)-r(3,j+which_H)
          dx = dx - boxx*dble(nint(onboxx*dx))
          dy = dy - boxy*dble(nint(onboxy*dy))
          dz = dz - boxz*dble(nint(onboxz*dz))
          drsq = dx*dx + dy*dy + dz*dz
          sq = sqrt(drsq)
          bin = int(sq/(pos(5)-pos(4)))
          if (bin .lt. 5000 .and. bin.gt.0) then
             v = v + potHH(bin)
             dfx = frcHH(bin) * dx*dx/drsq
             dfy = frcHH(bin) * dy*dy/drsq
             dfz = frcHH(bin) * dz*dz/drsq
             !write(6,*) bin, frcHH(bin), dfx,dfy,dfz
             !write(6,*) dx, dy, dz, sq
             dvdr(1,i+2) = dvdr(1,i+2) - dfx
             dvdr(2,i+2) = dvdr(2,i+2) - dfy
             dvdr(3,i+2) = dvdr(3,i+2) - dfz
             dvdr(1,j+which_H) = dvdr(1,j+which_H) + dfx
             dvdr(2,j+which_H) = dvdr(2,j+which_H) + dfy
             dvdr(3,j+which_H) = dvdr(3,j+which_H) + dfz
             vir(1,1) = vir(1,1) - dx * dfx
             vir(2,2) = vir(2,2) - dy * dfy
             vir(3,3) = vir(3,3) - dz * dfz
             vir(1,2) = vir(1,2) - dx * dfy
             vir(1,3) = vir(1,3) - dx * dfz
             vir(2,1) = vir(2,1) - dy * dfx
             vir(2,3) = vir(2,3) - dy * dfz
             vir(3,1) = vir(3,1) - dz * dfx
             vir(3,2) = vir(3,2) - dz * dfy
          endif
        enddo ! O--H
     enddo
  enddo

  return
end subroutine epsr_basic

subroutine epsr_list(r,dvdr,v,vir,na,boxlxyz,list,point,njump)
  implicit none
  include 'globals.inc'
  ! ----------------------------------------------------------------
  ! EPSR - Neighbour list version
  ! ----------------------------------------------------------------
  integer na,i,j,point(na+3),list(maxnab*na),ibeg,iend,inab,njump
  real(8) r(3,na),dvdr(3,na),vir(3,3),boxlxyz(3),v
  real(8) oo_eps,oo_sig,oo_gam,rcut
  real(8) sigsq,rcutsq,vij
  real(8) drsq,onr2,fij,dfx,dfy,dfz,sr2,sr6,wij
  real(8) dx,dy,dz,vscale,dscale
  real(8) rxj,ryj,rzj,dvdxj,dvdyj,dvdzj
  real(8) onboxx,onboxy,onboxz,boxx,boxy,boxz
  common /oo_param/ oo_eps,oo_sig,oo_gam,rcut
  logical epsr
  integer epsr_update
  real(8) pos(5000),potOO(5000),frcOO(5000),potOH(5000),frcOH(5000),potHH(5000),frcHH(5000)
  common /EPSR/ epsr, epsr_update, pos, potOO, frcOO, potOH, frcOH, potHH, frcHH

  ! clear arrays

  v = 0.d0
  vir(:,:) = 0.d0
  dvdr(:,:) = 0.d0

  rcutsq = rcut*rcut
  boxx = boxlxyz(1)
  boxy = boxlxyz(2)
  boxz = boxlxyz(3)
  onboxx = 1.d0/boxx
  onboxy = 1.d0/boxy
  onboxz = 1.d0/boxz

  ! use the neighbour list to calculate energy and forces

  sigsq = oo_sig*oo_sig
  do j = 1,na,njump
     ibeg = point(j)
     iend = point(j+njump) - 1
     if (ibeg.le.iend) then
        rxj = r(1,j)
        ryj = r(2,j)
        rzj = r(3,j)
        dvdxj = dvdr(1,j)
        dvdyj = dvdr(2,j)
        dvdzj = dvdr(3,j)
        do inab = ibeg,iend
           i = list(inab)
           dx = r(1,i) - rxj
           dy = r(2,i) - ryj
           dz = r(3,i) - rzj
           dx = dx - boxx*dble(nint(onboxx*dx))
           dy = dy - boxy*dble(nint(onboxy*dy))
           dz = dz - boxz*dble(nint(onboxz*dz))
           drsq = dx*dx + dy*dy + dz*dz
           if (drsq .lt. rcutsq) then
              onr2 = 1.d0/drsq
              sr2 = sigsq * onr2
              sr6 = sr2 * sr2 * sr2
              vij = sr6 * (sr6-1.d0)
              v = v + vij
              wij = sr6 * (sr6-0.5d0)
              fij = wij * onr2
              dfx = fij * dx
              dfy = fij * dy
              dfz = fij * dz
              dvdr(1,i) = dvdr(1,i) - dfx
              dvdr(2,i) = dvdr(2,i) - dfy
              dvdr(3,i) = dvdr(3,i) - dfz
              dvdxj = dvdxj + dfx
              dvdyj = dvdyj + dfy
              dvdzj = dvdzj + dfz
              vir(1,1) = vir(1,1) - dx * dfx
              vir(2,2) = vir(2,2) - dy * dfy
              vir(3,3) = vir(3,3) - dz * dfz
              vir(1,2) = vir(1,2) - dx * dfy
              vir(1,3) = vir(1,3) - dx * dfz
              vir(2,1) = vir(2,1) - dy * dfx
              vir(2,3) = vir(2,3) - dy * dfz
              vir(3,1) = vir(3,1) - dz * dfx
              vir(3,2) = vir(3,2) - dz * dfy
           endif
        enddo
        dvdr(1,j) = dvdxj
        dvdr(2,j) = dvdyj
        dvdr(3,j) = dvdzj
     endif
  enddo

  vscale = 4.d0*oo_eps
  dscale = 48.d0*oo_eps

  v = vscale * v
  vir(:,:) = dscale * vir(:,:)
  do i = 1,na,njump
     dvdr(1,i) = dscale*dvdr(1,i)
     dvdr(2,i) = dscale*dvdr(2,i)
     dvdr(3,i) = dscale*dvdr(3,i)
  enddo

  return
end subroutine epsr_list

subroutine epsr_cell(r,v,vir,dvdr,na,boxlxyz,njump)
  implicit none
  ! ------------------------------------------------------------------
  ! EPSR - Linked Cell Version
  ! ------------------------------------------------------------------
  integer na,ix,iy,iz,k,i,j,incx,incy,incz,njump
  integer jx,jy,jz,nn,ncellxyz(3),mcellxyz(3)
  integer list(na)
  real(8) r(3,na),dvdr(3,na),vir(3,3),boxlxyz(3)
  real(8) oo_eps,oo_sig,oo_gam,rcut,rcutsq,vij,sigsq,v
  real(8) drsq,onr2,fij,dfx,dfy,dfz,sr2,sr6,wij
  real(8) dx,dy,dz,vscale,dscale
  real(8) boxx,boxy,boxz,onboxx,onboxy,onboxz
  integer, allocatable :: link(:,:,:),ink(:,:)
  integer, allocatable :: mapx(:,:),mapy(:,:),mapz(:,:)
  common /oo_param/ oo_eps,oo_sig,oo_gam,rcut
  logical epsr
  integer epsr_update
  real(8) pos(5000),potOO(5000),frcOO(5000),potOH(5000),frcOH(5000),potHH(5000),frcHH(5000)
  common /EPSR/ epsr, epsr_update, pos, potOO, frcOO, potOH, frcOH, potHH, frcHH

  ! generate linked cell list

  call cell_setup(boxlxyz,rcut,mcellxyz,ncellxyz,nn)

  allocate (link(ncellxyz(1),ncellxyz(2),ncellxyz(3)))
  allocate (ink(3,nn))
  allocate (mapx(ncellxyz(1),-mcellxyz(1):mcellxyz(1)))
  allocate (mapy(ncellxyz(2),-mcellxyz(2):mcellxyz(2)))
  allocate (mapz(ncellxyz(3),-mcellxyz(3):mcellxyz(3)))

  call cell_list(r,list,link,ink,mapx,mapy,mapz,nn, &
                 na,boxlxyz,ncellxyz,mcellxyz,njump)

  ! and evaluate the potential energy and forces

  boxx = boxlxyz(1)
  boxy = boxlxyz(2)
  boxz = boxlxyz(3)

  onboxx = 1.d0/boxx
  onboxy = 1.d0/boxy
  onboxz = 1.d0/boxz

  v = 0.d0
  vir(:,:) = 0.d0
  dvdr(:,:) = 0.d0

  rcutsq = rcut*rcut
  sigsq = oo_sig*oo_sig

  do iz = 1,ncellxyz(3)
     do iy = 1,ncellxyz(2)
        do ix = 1,ncellxyz(1)
           
           ! including contributions from within the same cell
           
           i = link(ix,iy,iz)
           do while (i .gt. 0)
              j = list(i)
              do while (j .gt. 0)
                 dx = r(1,i)-r(1,j)
                 dy = r(2,i)-r(2,j)
                 dz = r(3,i)-r(3,j)
                 dx = dx-boxx*dble(nint(onboxx*dx))
                 dy = dy-boxy*dble(nint(onboxy*dy))
                 dz = dz-boxz*dble(nint(onboxz*dz))
                 drsq = dx*dx + dy*dy + dz*dz
                 if (drsq .lt. rcutsq) then
                    onr2 = 1.d0/drsq
                    sr2 = sigsq*onr2
                    sr6 = sr2*sr2*sr2
                    vij = sr6*(sr6-1.d0)
                    v = v + vij
                    wij = sr6*(sr6-0.5d0)
                    fij = wij*onr2
                    dfx = fij*dx
                    dfy = fij*dy
                    dfz = fij*dz
                    dvdr(1,i) = dvdr(1,i) - dfx
                    dvdr(2,i) = dvdr(2,i) - dfy
                    dvdr(3,i) = dvdr(3,i) - dfz
                    dvdr(1,j) = dvdr(1,j) + dfx
                    dvdr(2,j) = dvdr(2,j) + dfy
                    dvdr(3,j) = dvdr(3,j) + dfz
                    vir(1,1) = vir(1,1) - dx * dfx
                    vir(2,2) = vir(2,2) - dy * dfy
                    vir(3,3) = vir(3,3) - dz * dfz
                    vir(1,2) = vir(1,2) - dx * dfy
                    vir(1,3) = vir(1,3) - dx * dfz
                    vir(2,1) = vir(2,1) - dy * dfx
                    vir(2,3) = vir(2,3) - dy * dfz
                    vir(3,1) = vir(3,1) - dz * dfx
                    vir(3,2) = vir(3,2) - dz * dfy
                 endif
                 j = list(j)
              enddo
              i = list(i)
           enddo
           
           ! neigbouring cells

           do k = 1,nn
              incx = ink(1,k)
              incy = ink(2,k)
              incz = ink(3,k)
              jx = mapx(ix,incx)
              jy = mapy(iy,incy)
              jz = mapz(iz,incz)

              i = link(ix,iy,iz)
              do while (i .gt. 0)
                 j = link(jx,jy,jz)
                 do while (j .gt. 0)
                    dx = r(1,i)-r(1,j)
                    dy = r(2,i)-r(2,j)
                    dz = r(3,i)-r(3,j)
                    dx = dx-boxx*dble(nint(onboxx*dx))
                    dy = dy-boxy*dble(nint(onboxy*dy))
                    dz = dz-boxz*dble(nint(onboxz*dz))
                    drsq = dx*dx+dy*dy+dz*dz
                    if (drsq .lt. rcutsq) then
                       onr2 = 1.d0/drsq
                       sr2 = sigsq*onr2
                       sr6 = sr2*sr2*sr2
                       vij = sr6*(sr6-1.d0)
                       v = v + vij
                       wij = sr6*(sr6-0.5d0)
                       fij = wij*onr2
                       dfx = fij*dx
                       dfy = fij*dy
                       dfz = fij*dz
                       dvdr(1,i) = dvdr(1,i) - dfx
                       dvdr(2,i) = dvdr(2,i) - dfy
                       dvdr(3,i) = dvdr(3,i) - dfz
                       dvdr(1,j) = dvdr(1,j) + dfx
                       dvdr(2,j) = dvdr(2,j) + dfy
                       dvdr(3,j) = dvdr(3,j) + dfz
                       vir(1,1) = vir(1,1) - dx * dfx
                       vir(2,2) = vir(2,2) - dy * dfy
                       vir(3,3) = vir(3,3) - dz * dfz
                       vir(1,2) = vir(1,2) - dx * dfy
                       vir(1,3) = vir(1,3) - dx * dfz
                       vir(2,1) = vir(2,1) - dy * dfx
                       vir(2,3) = vir(2,3) - dy * dfz
                       vir(3,1) = vir(3,1) - dz * dfx
                       vir(3,2) = vir(3,2) - dz * dfy
                    endif
                    j = list(j)
                 enddo
                 i = list(i)
              enddo
           enddo
        enddo
     enddo
  enddo
  
  ! common factors

  vscale = 4.d0*oo_eps
  dscale = 48.d0*oo_eps
  v = vscale*v
  vir(:,:) = dscale*vir(:,:)
  do i = 1,na,njump
     dvdr(1,i) = dscale*dvdr(1,i)
     dvdr(2,i) = dscale*dvdr(2,i)
     dvdr(3,i) = dscale*dvdr(3,i)
  enddo

  deallocate (link,ink,mapx,mapy,mapz)

  return
end subroutine epsr_cell
