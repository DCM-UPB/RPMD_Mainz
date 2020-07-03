module useful
implicit none
private  
character(*), parameter :: moduleN="useful"
integer :: unit_search_start=2136
public :: get_Unit,toChar,lrtrim
interface toChar
    module procedure toChar_R
    module procedure toChar_int
end interface
contains
    !-----------------------------------------------------------------------
    !<summary>
    !returns the ioUnit of an open file or supplies a new ioUnit if the file
    !is not opended yet
    !</summary>
    !<param name="fName">the file name</param>
    !<param name="ioUnit">returns the ioUnit number of the file</param>
    !<param name="ioIsOpen">returns wheter the file is open</param>
    !-----------------------------------------------------------------------
    subroutine get_Unit(fName,ioUnit,ioIsOpen)
        implicit none
        character(LEN=*),intent(in)                    :: fName
        integer,intent(out)                            :: ioUnit
        logical,intent(out),optional                   :: ioIsOpen
        logical                                        :: loc_ioIsOpen

        loc_ioIsOpen=.FALSE.
        inquire(FILE=fName,OPENED=loc_ioIsOpen,NUMBER=ioUnit)
        !-----------------------------------------------
        !if file is not open, search for a ioUnit not
        !already in use.
        if(.NOT.loc_ioIsOpen)then
            loc_ioIsOpen=.TRUE.
            ioUnit=unit_search_start
            do while(loc_ioIsOpen)
                inquire(UNIT=ioUnit,OPENED=loc_ioIsOpen)
                if(loc_ioIsOpen)ioUnit=ioUnit+1
            end do
            unit_search_start=ioUnit
        end if
        if(present(ioIsOpen))ioIsOpen=loc_ioIsOpen
    end subroutine
    pure function toChar_int(an_int)
        implicit none
        integer,intent(in)                          :: an_int
        character(LEN=getSize_int(an_int))          :: toChar_int
        character(LEN=32)                           :: format_string
        write(format_string,*)"(I",len(toChar_int),")"
        write(toChar_int,FMT=format_string)an_int
        toChar_int=trim(adjustl(toChar_int))
    end function
    !-------------------------------------------------------------------
    !<summary>
    !   retun a number as a character of appropriate length
    !</summary>
    !-------------------------------------------------------------------
    pure function toChar_R(a_real)
        implicit none
        double precision,intent(in)                 :: a_real
        character(LEN=getSize_real_form(a_real))    :: toChar_r
        character(LEN=32)                           :: format_string
        write(format_string,*)"(F",len(toChar_r),".6)"
        write(toChar_R,FMT=format_string)a_real
        toChar_R=trim(adjustl(toChar_r))
    end function
    !-------------------------------------------------------------------
      !<summary>
      !    removes spaces from left and right side of a string
      !</summary>
      !<param name="str">the string to be trimmed</param>
      !<returns>a string withoug preceding or following spaces</returns>
    !-------------------------------------------------------------------
    pure function lrtrim(str)
        implicit none
        character(*),intent(in)                     :: str
        character(LEN=:),allocatable                :: lrtrim
        allocate(character(LEN=len_trim(adjustl(str))) :: lrtrim)
        if(len(lrtrim).gt.0)then
            write(lrtrim,FMT='(A'//toChar(len(lrtrim))//')')trim(adjustl(str))
        end if
    end function
    !-----------------------------------------------------------------------
    !<summary>
    !    Funktionen zum Festlegen der Stringlänge einer Zahl
    !</summary>
    !<param name="an_int">a integer number</param>
    !<retruns> the number of digits the integer needs for print</returns>
    !-----------------------------------------------------------------------
    pure function getSize_int(an_int)result(sizeof)
        implicit none
        integer,intent(in)                             :: an_int
        integer                                        :: calc_int
        integer                                        :: sizeof

        sizeof=1
        calc_int=an_int
        do while (abs((calc_int/10)).gt.0)
            calc_int=calc_int/10
            sizeof=sizeof+1
        end do
        if(an_int.lt.0)then
            sizeof=sizeof+1
        end if
    end function
    !-----------------------------------------------------------------
    !<summary>
    !   Returns actual length of a real of a given format spec
    !   Exploits knowlegde of the how an IEEE 753 double precision
    !   is organized
    !-----------------------------------------------------------------
    pure function getSize_real_form(a_real,passed_formats)result(sizeof)
        implicit none
        character(*),parameter :: FName='getSize_real_form'
        character(*),parameter :: CName=moduleN//":"//FName
        double precision,intent(in)                     :: a_real
        character(LEN=*),intent(in),optional            :: passed_formats
        integer                                         :: sizeof
        integer(KIND=8)                                 :: InpBitPat
        integer(KIND=8)                                 :: MantBitPat
        integer(KIND=8)                                 :: ExpBitPat
        integer                                         :: ExpSize
        integer                                         :: MantSize
        integer                                         :: SgnSize
        character(LEN=:),allocatable                    :: formats
        character                                       :: CFormats
        logical                                         :: Exponents

        InpBitPat  = transfer(a_real,InpBitPat)    !transfer bit pattern to int(8)
        ExpBitPat  = ibits(InpBitPat,52,11)        !get the exponent bits
        MantBitPat = ibits(InpBitPat,0,52)         !get the mantissa bits
        sgnSize    = merge(1,0,btest(InpBitPat,63))!check if the sign bit is set
        !-----------------------------------------------
        !get custom length of mantissa
        if(present(passed_formats))then
            allocate(character(LEN=5) :: formats)
            formats=adjustl(passed_formats)
        else
            allocate(character(LEN=5) :: formats)
            formats='F12.6'
        end if
        read(formats,FMT='(A1,2x,1x,I1)')CFormats,MantSize
        Exponents=(CFormats=="E")                  !Test if Exponents exist
        !-----------------------------------------------
        !Exponents: 4 chars for '0.E+', Float: 1 char for '.'
        MantSize = merge(MantSize+4,MantSize+1,Exponents)
        select case (ExpBitPat)
            case(int(z'000', KIND=8))!handle denormal numbers
                ExpSize=merge(3,1,Exponents)
            case(int(z'7ff', KIND=8))!represents infinity or Nan depending on Mantissa
                sizeof=merge(8+sgnSize,3,MantBitPat.eq.int(z'00000', KIND=8))
                return
            case default
                if(Exponents)then
                    ExpSize=2   !total length is 4 but 2 already counted
                else
                  !------------------------------------
                  !determine how many digit(10**x),x>=0)
                  !a floating point has. Then add +1 b/c
                  !10**0 also occupies a digit.
                  ExpSize=floor(log(abs(a_real))/log(1.d1))
                  ExpSize=merge(ExpSize,0,ExpSize.gt.0)+1
             end if
        end select
        sizeof=MantSize+ExpSize+sgnSize
    end function
end module

