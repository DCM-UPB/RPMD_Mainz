module FileMod
use, intrinsic :: iso_fortran_env, only : error_unit
use useful, only : obtain_unit=>get_unit
implicit none
character(*),parameter :: moduleN='FileMod'
type                                    :: FileHandle
    private
        character(LEN=64)               :: my_name
        integer                         :: my_unit
        logical                         :: opened
    contains
        procedure,pass(this),public     :: open
        procedure,pass(this),public     :: close
        procedure,pass(this),public     :: unit=>get_unit
        procedure,pass(this),public     :: isOpen
end type    
interface FileHandle
    module procedure fh_cnstr1
end interface
contains
    !------------------------------------------------------------------
    !<summary>
    !   Constructor for type file handle
    !</summary>
    !<param name="fName">name of the file connected to the handler</param>
    !<return> a new file handle object</return>
    !------------------------------------------------------------------
    function fh_cnstr1(fName)
        implicit none
        character(*),parameter :: routineN='fh_cnstr1'
        character(*),parameter :: routineP=moduleN//":"//routineN
        type(FileHandle)                :: fh_cnstr1
        character(LEN=*),intent(in)     :: fName
        fh_cnstr1%my_name=fName
        fh_cnstr1%my_unit=-1
        fh_cnstr1%opened=.FALSE.
    end function
    !------------------------------------------------------------------
    !<summary>
    !   open the file connected to the handler 
    !</summary>
    !<param name="this"> the binding type</param>
    !<param name="status">referring to the preexistence of file</param>
    !<param name="action">what is to be done with the file</param>
    !<param name="access">set r/w permissions for file</param>
    !<param name="iostat">return the IO status number </param>
    !<param name="iomsg"> return the IO message associated with IO stat </param>
    !<param name="position">set the file pointer </param>
    !<returns>whether the file was successfully opened</returns>
    !------------------------------------------------------------------
    logical function open(this,status,action,access,iostat,iomsg,position)
        implicit none    
        character(*),parameter :: routineN='open'
        character(*),parameter :: routineP=moduleN//":"//routineN
        class(FileHandle),intent(inout) :: this
        character(*),intent(in),optional:: status, action,access,position
        character(*),intent(out),optional:: iomsg
        integer,intent(out),optional    :: iostat
        character(LEN=:),allocatable    :: loc_status,loc_action,loc_access
        character(LEN=:),allocatable    :: loc_iomsg, loc_position
        integer                         :: loc_iostat
        integer                         :: ioUnit
        logical                         :: ioIsOpen

        allocate(character(LEN=25) :: loc_status,loc_action,loc_access)
        allocate(character(LEN=25) :: loc_iomsg,loc_position)
        loc_status='UNKNOWN'
        loc_action='READWRITE'
        loc_access='SEQUENTIAL'
        loc_position='ASIS'
        if(present(status))loc_status=status
        if(present(action))loc_action=action
        if(present(access))loc_access=access
        if(present(position))loc_position=position

        call obtain_unit(this%my_name,ioUnit,ioIsOpen)
        if(ioIsOpen)then
            open=.FALSE.
            return 
        end if    
        open(ioUnit,FILE=this%my_name,STATUS=loc_status,ACTION=loc_action,&
        ACCESS=loc_access,POSITION=loc_position,IOSTAT=loc_iostat,&
        IOMSG=loc_iomsg)
        if(present(iostat))iostat=loc_iostat
        if(present(iomsg))iomsg=loc_iomsg
        !------------------------------o
        !the following passage attempts to imitate statndard 
        !fortran behavior. Hence, stop program if iostat is not
        !specified, else return the error but continue
        if( loc_iostat.ne.0 .AND. .not. present(iostat) )then
            write(error_unit,*)loc_iomsg
            stop
        else if(loc_iostat.ne.0)then
            open=.FALSE.
        else
            this%my_unit=ioUnit
            this%opened=.TRUE.
            open=.TRUE.
        end if
    end function
    !------------------------------------------------------------------
    !<summary>
    !   close an open file
    !</summary>
    !<param name=this>the binding type</param>
    !<param name="status">what shall happen with file upon closing</param>
    !<return>whether the closing succedes</return>
    !------------------------------------------------------------------
    logical function close(this,status)
        implicit none
        character(*),parameter :: routineN='close'
        character(*),parameter :: routineP=moduleN//":"//routineN
        class(FileHandle),intent(inout)         :: this
        character(LEN=*),optional,intent(in)    :: status
        character(LEN=:),allocatable            :: loc_status
        allocate(character(LEN=25) :: loc_status)
        loc_status='KEEP'
        if(present(status))loc_status=status
        if (this%opened)then
            close(this%my_unit,STATUS=loc_status)
            close=.TRUE. 
            this%opened=.FALSE.
            this%my_unit=-1
        else
            close=.FALSE.
        end if
    end function
    !------------------------------------------------------------------
    !<summary>
    !   returns a positive number if the file is open and a negativ 
    !   number if it is not.
    !</summary>
    !<param name="this"> the binding type</param>
    !<return>the unit number of the associated file</return>
    !------------------------------------------------------------------
    integer pure function get_unit(this)
        implicit none
        character(*),parameter :: routineN='get_unit'
        character(*),parameter :: routineP=moduleN//":"//routineN
        class(FileHandle),intent(in)    :: this
        get_unit=this%my_unit
    end function
    !------------------------------------------------------------------
    !<summary>
    !   get procedure for the opened property
    !</summary>
    !<param name="this"> the binding type</param>
    !<return>the value of the type's opened field</return>
    !------------------------------------------------------------------
    logical pure function isOpen(this)
        implicit none
        character(*),parameter :: routineN='isOpen'
        character(*),parameter :: routineP=moduleN//":"//routineN
        class(FileHandle),intent(in)    :: this
        isOpen=this%opened
    end function
end module
