!Subroutine for conversion of rpmd variables to e3b variables
   subroutine scon(na,nb,nm,r,boxlxyz,pos,box)
      implicit none
      include 'globals.inc'      

      integer i,j,k,l,n
      integer na,nb,nm
      real(8) r(3,na,nb),boxlxyz(3),pos(na/3,4,3),box(3),m(3,na/3,nb),alpha,alpha2
      common /ew_param/ alpha

      !Conversion boxlxyz (bohr) to box (nanometer)
       box(1) = boxlxyz(1) * ToA / 10
       box(2) = boxlxyz(2) * ToA / 10
       box(3) = boxlxyz(3) * ToA / 10
      

 
    !Calculation of M-position
     alpha2 = 0.5d0 * (1.d0 - alpha)
     n = 0
       do i = 1, na, 3
         n=n+1
        do j = 1, 3
           m(j,n,nb) = alpha * r(j,i,nb) + alpha2*(r(j,i+1,nb)+r(j,i+2,nb))
        End do
      End do

     !Generate Array pos(molecule,site,xyz)
      do k = 1, nm        
        do l=1,3
             pos(k,l,1) = r(1,k*3-3+l,nb)* ToA / 10
             pos(k,l,2) = r(2,k*3-3+l,nb)* ToA / 10
             pos(k,l,3) = r(3,k*3-3+l,nb)* ToA / 10
        End do
       pos(k,4,1) = m(1,k,nb)* ToA / 10
       pos(k,4,2) = m(2,k,nb)* ToA / 10
       pos(k,4,3) = m(3,k,nb)* ToA / 10
     End do
      
    End

!**********************************************************************
!     Example subroutine to demonstrate how the two-body correction and 
!     three-body terms are calculated for the E3B model.  
!     E3B uses TIP4P as a reference potential.  We used
!     Gromacs to calculate the potential, forces and virial for TIP4P.
!     A (paralellized) version of this code was then added to Gromacs (after
!     do_force in md.c for instance).  This subroutine calculates the potential
!     energy, forces and virial for the two-body correction and three-body
!     terms which are added to what Gromacs calculates for TIP4P.  
!     (See JCP, 134, 184501)
!
!     Written by C. J. Tainter - 6/27/2012
!
!     nm is number of molecules
!     pos is input coordinates in nanometers - indexed by (molecule,site,xyz)
!     sites are: O=1, H1=2, H2=3, M=4
!     f3B is output forces in kJ/(mol*nm) - indexed same as pos
!     vir is output virial in kJ/mol 
!     box is input box size in nanometers 
!     PE is output potential energy in kJ/mol
!     (I did not divide PE by nm, Gromacs does this with the -nmol option of g_energy) 
!
!     Subroutine split into two loops
!
!     First calculates two-body correction and terms for three-body loop
!     Terms for three-body loop are as follows
!     pair: specifies indices of two molecules 
!     es: exponential*switching function for all O-H combinations of given pair
!     des: derivative of es with respect to r
!     vec: vector distances 
!     total: keeps track of the sum of three body terms for each site 
!
!     Second loop calculates three-body PE/forces/virial.
!     We avoid a triple loop by using the total array and subtracting
!     appropriate term(s) depending on which type of three-body term
!     is being calculated.  
!**********************************************************************
      subroutine tainter_e3b(nm,pos,f3B,vir,box,PE)
      implicit none
      include 'globals.inc'
      
      integer i,j,k,kk,l,ia,ib,ja,jb
      integer nm,fac,ntmp
      parameter(fac=40) !estimate #3B terms (based on density) to guess array size
      real(8) c3B(4),exp3B,exp2B,rs,rf,rmax,rs2,rf2,sw_fac,sw_num
      real(8) pos(nm,4,3),f3B(nm,4,3),vir(3,3),box(3),PE
      real(8) invBox(3),total(nm,3)
      real(8) es(fac*nm,3,3),des(fac*nm,3,3),vec(fac*nm,3,3,3)
      real(8) r,tmp,ftmp,dis,pos1(3),pos2(3),xtmp(3)
      integer pair(fac*nm,2)

!     E3B parameters
      c3B(1)=1745.7           ! Ea  (kJ/mol)
      c3B(2)=-4565.0          ! Eb  (kJ/mol)
      c3B(3)=7606.8           ! Ec  (kJ/mol)
      c3B(4)=2.349*1000000.0  ! E2  (kJ/mol)
      exp2B=4.872*10.0        ! k2  (nm^-1)
      exp3B=1.907*10.0        ! k3  (nm^-1)
      rs=0.5
      rs2=rs*rs
      rf=0.52
      rf2=rf*rf
      rmax=rf+0.09572
!     sw_fac and sw_num are numbers needed for the switching function
      sw_fac=1.0/((rf-rs)**3.0D0)
      sw_num=rf-3.0*rs
      PE=0.0

!     Initialize arrays
      Do i=1,3
        invBox(i)=1.0/box(i)
        Do j=1,3
           vir(i,j)=0.0
        End Do
      End Do
      Do i=1,nm
        Do j=1,3
          total(i,j)=0.0
          Do k=1,4
            f3B(i,k,j)=0.0
          End Do
        End Do
      End Do
      Do i=1,fac*nm
        Do j=1,2
          pair(i,j)=0
        End Do
        Do j=1,3
          Do k=1,3
            es(i,j,k)=0.0
            des(i,j,k)=0.0
            Do l=1,3
              vec(i,j,k,l)=0.0
            End Do
          End Do
        End Do
      End Do

      ntmp=1
!**********************************************************************
!     First Loop
!     First calculates Two-Body Correction
!     Second, because three-body ineractions (intermolecular O-H distances)
!     are zero after 0.52 nm, we know there will not be a contribution to 
!     three-body potential if O-O distance is greater than rmax=0.52+0.09572
!**********************************************************************
      Do ia=1,nm-1
        Do ib=ia+1,nm

!         Two-Body Correction
!         Save xyz coords of oxygens on molecules ia/ib to pos1/pos2
          Do k=1,3
            pos1(k)=pos(ia,1,k)
            pos2(k)=pos(ib,1,k)
          End Do
!         Calculate nearest image distance between the two oxygens
!         Vector is saved in xtmp
!         Update potential energy, forces, and virial
          r=dis(pos1,pos2,box,invBox,xtmp)
          tmp=c3B(4)*exp(-exp2B*r)
          PE=PE+tmp
          Do k=1,3
            ftmp=xtmp(k)*tmp*exp2B/r
            call virial(vir,xtmp,ftmp,k)
            f3B(ia,1,k)=f3B(ia,1,k)+ftmp
            f3B(ib,1,k)=f3B(ib,1,k)-ftmp
          End Do
          
!         Now, if O-O distance is less than rmax, calculate terms
!         needed for 3B interactions.  ntmp specifies this ia,ib pair.
!         Note: some of these still may be zero. 
          If (r.le.rmax)Then
            If (ntmp.gt.fac*nm)Then
              print*,'# interactions > array size, increase fac'
              stop
            End If 
            pair(ntmp,1)=ia
            pair(ntmp,2)=ib
  
            Do jb=2,3  ! Hydrogens on molecule ib (oxygen on ia)
!             Save xyz coords and calculate distance
              Do k=1,3
                pos1(k)=pos(ia,1,k)
                pos2(k)=pos(ib,jb,k)
              End Do
              r=dis(pos2,pos1,box,invBox,xtmp)
!             Calculate 3B terms if r<rf (else they will remain zero)
              If (r.le.rf) Then
                es(ntmp,1,jb)=exp(-exp3B*r)
                des(ntmp,1,jb)=exp3B*es(ntmp,1,jb)/r
!               If 0.5<r<0.52, we need to do additional stuff with switching function
                If (r.gt.rs) Then
                  tmp=sw_fac*(rf-r)*(rf-r)*(sw_num+2.0*r)
                  des(ntmp,1,jb)=des(ntmp,1,jb)*tmp
                  des(ntmp,1,jb)=des(ntmp,1,jb)-es(ntmp,1,jb)* &
                                6.0*sw_fac*(rf-r)*(rs-r)/r
                  es(ntmp,1,jb)=es(ntmp,1,jb)*tmp
                End If ! r.gt.rs 
!               Save vector information in vec
                vec(ntmp,1,jb,1)=xtmp(1)
                vec(ntmp,1,jb,2)=xtmp(2)
                vec(ntmp,1,jb,3)=xtmp(3)
!               Keep track of totals
                total(ia,1)=total(ia,1)+es(ntmp,1,jb)
                total(ib,jb)=total(ib,jb)+es(ntmp,1,jb)
              End If ! r.le.rf
            End Do ! jb
            Do ja=2,3  ! Hydrogens on molecule ia (oxygen on ib)
!             Save xyz coords and calculate distance
              Do k=1,3
                 pos1(k)=pos(ia,ja,k)
                 pos2(k)=pos(ib,1,k)
              End Do
              r=dis(pos1,pos2,box,invBox,xtmp)   
!             Calculate 3B terms if r<rf (else they will remain zero)
              If (r.le.rf) Then
                es(ntmp,ja,1)=exp(-exp3B*r)
                des(ntmp,ja,1)=exp3B*es(ntmp,ja,1)/r
!               If 0.5<r<0.52, we need to do additional stuff with switching function
                If (r.gt.rs) Then
                  tmp=sw_fac*(rf-r)*(rf-r)*(sw_num+2.0*r)
                  des(ntmp,ja,1)=des(ntmp,ja,1)*tmp
                  des(ntmp,ja,1)=des(ntmp,ja,1)-es(ntmp,ja,1)* &
                                6.0*sw_fac*(rf-r)*(rs-r)/r
                  es(ntmp,ja,1)=es(ntmp,ja,1)*tmp
                End If ! r.gt.rs
!               Save vector information in vec
                vec(ntmp,ja,1,1)=xtmp(1)
                vec(ntmp,ja,1,2)=xtmp(2)
                vec(ntmp,ja,1,3)=xtmp(3)
!               Keep track of totals
                total(ib,1)=total(ib,1)+es(ntmp,ja,1)
                total(ia,ja)=total(ia,ja)+es(ntmp,ja,1)
              End If !r.le.rf
            End Do !ja
            ntmp=ntmp+1
          End If ! r.le.rmax         
        End Do !ib
      End Do !ia

!**********************************************************************
!     Second Loop
!     Here we use the information saved in the fist loop (es,des,total,vec)
!     to calculate PE/forces/virial for the 3B terms. 
!**********************************************************************
      Do i=1,ntmp-1  
        ia=pair(i,1)
        ib=pair(i,2)
        Do k=2,3
          kk=mod(k-1,2)+2 ! k=2(3), kk=3(2)
          Do j=1,3
            pos1(j)=vec(i,k,1,j)
            pos2(j)=vec(i,1,k,j)
          End Do
!         Type A
          tmp=c3B(1)*(total(ia,kk)-es(i,kk,1))
          PE=PE+es(i,k,1)*tmp/2.0
          Do j=1,3
            ftmp=tmp*des(i,k,1)*pos1(j)
            call virial(vir,pos1,ftmp,j)
            f3B(ib,1,j)=f3B(ib,1,j)-ftmp
            f3B(ia,k,j)=f3B(ia,k,j)+ftmp
          End Do
          tmp=c3B(1)*(total(ib,kk)-es(i,1,kk))
          PE=PE+es(i,1,k)*tmp/2.0
          Do j=1,3
            ftmp=tmp*des(i,1,k)*pos2(j)
            call virial(vir,pos2,ftmp,j)
            f3B(ia,1,j)=f3B(ia,1,j)-ftmp
            f3B(ib,k,j)=f3B(ib,k,j)+ftmp
          End Do
!         Type B
          tmp=c3B(2)*(total(ia,1)-es(i,1,k)-es(i,1,kk))
          PE=PE+es(i,k,1)*tmp/2.0
          Do j=1,3
            ftmp=tmp*des(i,k,1)*pos1(j)
            call virial(vir,pos1,ftmp,j)
            f3B(ib,1,j)=f3B(ib,1,j)-ftmp
            f3B(ia,k,j)=f3B(ia,k,j)+ftmp
          End Do
          tmp=c3B(2)*(total(ib,1)-es(i,k,1)-es(i,kk,1))
          PE=PE+es(i,1,k)*tmp/2.0
          Do j=1,3
            ftmp=tmp*des(i,1,k)*pos2(j)
            call virial(vir,pos2,ftmp,j)
            f3B(ia,1,j)=f3B(ia,1,j)-ftmp
            f3B(ib,k,j)=f3B(ib,k,j)+ftmp
          End Do
          tmp=c3B(2)*(total(ib,k)-es(i,1,k)+total(ib,kk)-es(i,1,kk))
          PE=PE+es(i,k,1)*tmp/2.0
          Do j=1,3
            ftmp=tmp*des(i,k,1)*pos1(j)
            call virial(vir,pos1,ftmp,j)
            f3B(ib,1,j)=f3B(ib,1,j)-ftmp
            f3B(ia,k,j)=f3B(ia,k,j)+ftmp
          End Do
          tmp=c3B(2)*(total(ia,k)-es(i,k,1)+total(ia,kk)-es(i,kk,1))
          PE=PE+es(i,1,k)*tmp/2.0
          Do j=1,3
            ftmp=tmp*des(i,1,k)*pos2(j)
            call virial(vir,pos2,ftmp,j)
            f3B(ia,1,j)=f3B(ia,1,j)-ftmp
            f3B(ib,k,j)=f3B(ib,k,j)+ftmp
          End Do
!         Typce C
          tmp=c3B(3)*(total(ib,1)-es(i,k,1)-es(i,kk,1))
          PE=PE+es(i,k,1)*tmp/2.0
          Do j=1,3
            ftmp=tmp*des(i,k,1)*pos1(j)
            call virial(vir,pos1,ftmp,j)
            f3B(ib,1,j)=f3B(ib,1,j)-ftmp
            f3B(ia,k,j)=f3B(ia,k,j)+ftmp
          End Do
          tmp=c3B(3)*(total(ia,1)-es(i,1,k)-es(i,1,kk))
          PE=PE+es(i,1,k)*tmp/2.0
          Do j=1,3
            ftmp=tmp*des(i,1,k)*pos2(j)
            call virial(vir,pos2,ftmp,j)
            f3B(ia,1,j)=f3B(ia,1,j)-ftmp
            f3B(ib,k,j)=f3B(ib,k,j)+ftmp
          End Do
        End Do  ! k
      End Do  ! i
     
      End
!**********************************************************************
!     function to calculate nearest image distances
!**********************************************************************
      function dis(pos1,pos2,box,invBox,vec)
      implicit none
      include 'globals.inc'
  
      real(8) pos1(3),pos2(3),vec(3),box(3),invBox(3)
      real(8) dx,dy,dz,ddx,ddy,ddz,dis
      integer tx,ty,tz

      
      dx=pos2(1)-pos1(1)
      dy=pos2(2)-pos1(2)
      dz=pos2(3)-pos1(3)
    
      tx=int(dx*invBox(1)+1.5)-1
      ty=int(dy*invBox(2)+1.5)-1
      tz=int(dz*invBox(3)+1.5)-1

      ddx=tx*box(1)-dx
      ddy=ty*box(2)-dy
      ddz=tz*box(3)-dz

      vec(1)=ddx
      vec(2)=ddy
      vec(3)=ddz

      dis=sqrt(ddx*ddx+ddy*ddy+ddz*ddz)
      End
!**********************************************************************
!     subroutine to add to virial 
!**********************************************************************
      subroutine virial(vir,xtmp,ftmp,k)
      implicit none
      include 'globals.inc'
      integer k
      real(8) vir(3,3),xtmp(3),ftmp

      vir(1,k)=vir(1,k)-0.5D0*xtmp(1)*ftmp
      vir(2,k)=vir(2,k)-0.5D0*xtmp(2)*ftmp
      vir(3,k)=vir(3,k)-0.5D0*xtmp(3)*ftmp

      End


subroutine econ(nm,f3B,vir_e3b,PE,dvdr_e3b)
 implicit none
  include 'globals.inc'
  integer nm,h,i
  real(8) f3B(nm,4,3),vir_e3b(3,3),PE,dvdr_e3b(3,nm*3) 

  !Conversion of potential Energy PE
  PE=PE/toKjmol
  
  !Conversion of Virial
  do i=1,3
        vir_e3b(i,1) = vir_e3b(i,1)/toKjmol
        vir_e3b(i,2) = vir_e3b(i,2)/toKjmol
        vir_e3b(i,3) = vir_e3b(i,3)/toKjmol
  end do
      
  !Conversion Forces
    do  h=1,nm
      do i=1,3
        dvdr_e3b(1,i+3*h-3)=f3B(h,i,1)/toKjmol*toA/10
        dvdr_e3b(2,i+3*h-3)=f3B(h,i,2)/toKjmol*toA/10
        dvdr_e3b(3,i+3*h-3)=f3B(h,i,3)/toKjmol*toA/10
      end do
    end do
  End

