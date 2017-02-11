module RDFMod
use useful, only : lrtrim
implicit none
private
character(*),parameter :: moduleN='RDFMod'
type                                    :: RDF
    real(8),Dimension(:,:),allocatable  :: gofr 
end type
public :: RDF,calc_rdf
contains
    !-------------------------------------------------------------------
    !<summary>
    !   calculate the radial distribution function for an n bead
    !   system between the atoms center_type and partner_type
    !</summary>
    !<param name="geo">the positions of all beads of all atoms</param>
    !<param name="element"> the element of all atoms </param>
    !<param name="center_type">what element shall be at the center</param>
    !<param name="partner_type">what element shall be the partner</param>
    !<param name="rdf_lbound">the lower value of the rdf_array</param>
    !<param name="rdf_ubound">the upper value of the rdf_array</param>
    !<param name="rdf">the array with the histogram setting</param>
    !-------------------------------------------------------------------
    subroutine calc_rdf(geo,element,center_type,partner_type,box,rdf_lbound,&
    rdf_ubound,rdf)
        implicit none
        real(8),Dimension(:,:,:),intent(in) :: geo !Dimension,atom,bead
        character(LEN=2),Dimension(:),intent(in) :: element
        character(LEN=2),intent(in) :: center_type, partner_type
        real(8),Dimension(:),intent(in) :: box
        real(8),intent(in) :: rdf_lbound,rdf_ubound
        real(8),Dimension(:),intent(inout) :: rdf
        integer :: iBead, iCenter,iPartner
        real(8),Dimension(size(box)) :: dist_vec,inv_box
        real(8) :: distance
        integer :: bin
        !this might be a good spot to parallelize over beads
        !make nBead copies of rdf, iterate and then sum up later
        inv_box(:)=1/box(:)
        do iBead=lbound(geo,3),ubound(geo,3),1
            do iCenter=lbound(geo,2),ubound(geo,2),1
                if( lrtrim(element(iCenter))/=lrtrim(center_type) )cycle
                do iPartner=lbound(geo,2),ubound(geo,2),1
                    if(iCenter.eq.iPartner)cycle
                    if( lrtrim(element(iPartner)).ne.lrtrim(partner_type) )cycle
                    dist_vec=geo(:,iBead,iPartner)-geo(:,iBead,iCenter)
                    dist_vec = dist_vec*dble(nint(dist_vec*inv_box))
                    distance=sqrt(sum(dist_vec**2))
                    bin=inBin(distance,rdf_lbound,rdf_ubound,size(rdf,1))
                    rdf(bin)=rdf(bin)+1
                end do
            end do
        end do
    end subroutine
    !-------------------------------------------------------------------
    !<summary>
    !   determine the correct bin to place an elemet in a histogram
    !</summary>
    !<param name="bin_value">
    !    the value for which the bin shall be deterimine
    !</param>
    !<param name="lbound">the lowest histiogam value</param>
    !<param name="ubound">the highest histogram value</param>
    !<param naem="size">the size of the historgram>/param>
    !<returns> the number of the bin</returns>
    !-------------------------------------------------------------------
    elemental function inBin(bin_value,lbound,ubound,size)
        implicit none 
        integer                 :: inBin
        real(8),intent(in)      :: bin_value,lbound,ubound
        integer,intent(in)      :: size
        inBin=floor((bin_value-lbound)/(ubound-lbound))*size
        inBin=inBin+1!becuase fortran arrays start at 
    end function
end module
