! module parse_ini
!
! Purpose: contains routines to parse INI files
!
! General usage:
!   Open an INI file to read with openIniFile
!   Read values with getIniString, getIniInt, getIniReal, and getIniBool
!   Close the INI file with closeIniFile
!   
!   The INI file format should be as follows:
!     [Section1]
!     ; Comment
!     keyword1=data1
!     keyword2=data2
!     [Section2]
!     keyword3=data3
!     etc...
!  
!   All comments, blank lines, etc. are ignored.
!
! Author: Jason Meiring
! Date: 2005/10/08
      module parse_ini
      
      implicit none
      
!     Scope Settings
! **** note that the 'private' at the beginning of this routine forces
! **** you to EXPLICITLY define variables that you want to be global
! **** and avoids unexpected references to variables within this module
      private
      public :: openIniFile, closeIniFile
      public :: getIniString, getIniInt, getIniReal, getIniPetscTruth
      
      ! Module constants *** NOTE *** THESE ARE PRIVATE TO THIS MODULE
      integer, parameter :: maxlen = 512  ! maximum line length in a file
      
      ! Module variables
      character(len=maxlen), allocatable, dimension(:) :: iniLines  ! holds file data in memory
      
      contains
      
      ! subroutine openIniFile
      !
      ! Purpose: read an INI file containing input parameters into memory
      !   This function must be called before any others. Only one INI file
      !   can be read at a time. 
      ! 
      ! Input:
      !   iniFilename: name of INI file to read
      !
      ! Output:
      !   stat (optional): true if successful read, otherwise false
      subroutine openIniFile(iniFilename, stat)
      
        implicit none
        
        ! Input/output parameters
        character(len=*), intent(in) :: iniFilename
        logical, intent(out), optional :: stat
        
        ! Local function variables
        integer :: istat                ! error status of file I/O operations
        integer :: astat                ! error status of memory allocation  
        integer :: numlines             ! number of lines in the INI file
        integer :: i                    ! loop counter
        character(len=maxlen) :: buffer ! read buffer
        
        ! Only one INI file can be opened at a time. Check to make sure
        ! no other file is open before proceeding
        if (allocated(iniLines)) then
          write(*,*) 'File I/O error!'
          write(*,*) '  openIniFile(): File already open'
          write(*,*) ''
          stop
        end if
        
        ! Try to open the file
        open(unit=27, file=iniFilename, iostat=istat, form='formatted',  &
              access='sequential', status='old', action='read')
      
        ! Error checking
        if(istat /= 0) then
          if(present(stat)) then
            stat = .false.
            return
          else
            write(*,*) 'File I/O error!'
            write(*,*) '  openIniFile(): Unable to open ',               &
              trim(iniFilename), ': istat=', istat
            write(*,*) ''
            stop
          end if
        end if
        
        ! Read through the file to determine how much memory to allocate
        numlines = 0
        
        do
          read(unit=27, fmt='(A)', iostat=istat) buffer    
          if (istat /= 0) exit   ! exit at EOF
          
          ! Increase the line counter
          numlines = numlines + 1
             
        end do
        
        if (numlines == 0) then
          if(present(stat)) then
            stat = .false.
            return
          else
            write(*,*) 'File I/O error!'
            write(*,'(1X,3A)') '  openIniFile(): File ',   &
              trim(iniFilename), ' is empty'
            write(*,*) ''
            close(unit=27)
            stop
          end if
        end if
        
        ! Try to allocate memory to hold file 
        allocate(iniLines(1:numlines), stat=astat)
        
        ! Error checking
        if (astat /= 0) then
          write(*,*) 'Memory error!'
          write(*,*) '  openIniFile(): Unable to allocate enough memory'
          write(*,*) ''
          close(unit=27)
          stop
        end if
        
        ! Now, read the file into memory
        rewind(unit=27)  ! rewind the file
        
        ! Read in the data, chopping off any space on the right and left
        do i=1,numlines
          read(27, '(A)') buffer
          iniLines(i) = adjustl(trim(buffer))
        end do
        
        ! All finished  
        close(unit=27)
        
        if(present(stat)) stat = .true.   ! successful open
        
      end subroutine openIniFile
      
      
      
      ! subroutine closeIniFile()
      !
      ! Purpose: deallocates memory used for a previously open INI file
      !
      ! Inputs:
      !   none
      ! 
      ! Output:
      !   stat (optional): true if deallocation was successful, otherwise false
      subroutine closeIniFile(stat)
      
        implicit none
        
        ! Output parameter
        logical, intent(out), optional :: stat
      
        ! Local function variable
        integer :: astat           ! status of array deallocation
        
        ! Check if the array was ever allocated. If not, just return true and exit
        ! Note: doesn't seem to work with Intel compilers
        if (.not. allocated(iniLines)) then
          if(present(stat)) stat = .true.
          return  
        end if
        
        deallocate(iniLines, stat=astat)
        
        ! Error checking
        if (astat /= 0) then
          write(*,*) 'Memory error!'
          write(*,*) '  closeIniFile(): Unable to deallocate memory'
          write(*,*) ''
          stop
        end if
        
        if(present(stat)) stat = .true.
        
      end subroutine closeIniFile
      
      
      
      ! subroutine findIniValue()
      !
      ! Purpose: search the INI file for 'keyword' in 'section'
      !   If found, returns the value associated with keyword.
      !   This subroutine is intended to be used only by the 
      !   getIniXXX subroutines in this module
      !
      ! Inputs:
      !   section: [section] in which to look for keyword
      !   keyword: name of string, given as keyword=string
      !
      ! Output:
      !   value: keyword value (returns '' if not found)
      subroutine findIniValue(section, keyword, value)
        
        use string_utils
        implicit none
        
        ! Input/output parameters
        character(len=*), intent(in) :: section
        character(len=*), intent(in) :: keyword
        character(len=*), intent(out) :: value
          
        ! Local function variables
        character(len=maxlen) :: line                 ! current line
        integer :: i                                  ! loop counter
        integer :: l                                  ! line length
        integer :: s                                  ! substring location  
        integer :: curLine                            ! current line number
        logical :: sectionFound                       ! true if given section is found
        logical :: keywordFound                       ! true if given keyword is found
          
        ! Check that the keyword is a non-empty string
        if(trim(keyword) == '') then
          value = ''
          return
        end if
        
        ! First find the section
        sectionFound = .false.
        curLine = 1
           
        do i=1,size(iniLines)
          line = iniLines(i)
          l = len_trim(line)
      
          ! Skip blank or short lines
          if(l < 3) cycle
       
          ! Skip comment lines
          if(line(1:1) == ';' .or. line(1:1) == '*') cycle
          
          ! Check if the line contains [ ]
          if(line(1:1) == '[' .and. line(l:l) == ']') then
            if(strUpCase(line(2:(l-1))) == strUpCase(section)) then
              ! Found the specified section
              sectionFound = .true.
              curLine = i
              exit
            end if
          end if
        end do
          
        if(.not. sectionFound) then
          value = ''
          return
        end if
        
        ! Next, search for the keyword
        keywordFound = .false.
         
        do i=curLine+1,size(iniLines)
          line = iniLines(i)
          l = len_trim(line)
          
          ! Skip blank or short lines
          if(l < 3) cycle
        
          ! Skip comment lines
          if(line(1:1) == ';' .or. line(1:1) == '*') cycle
          
          ! Make sure that we don't hit a new section. If we do, then quit searching
          if(line(1:1) == '[' .and. line(l:l) == ']') then
            exit
          end if
          
          ! Search for the keyword
          s = index(line, '=')
          
          if(s > 0) then
            ! An '=' sign was found. Now see if the keyword exists
            ! to the left of the '=' sign.      
           if(trim(adjustl(strUpCase(line(1:s-1)))) == strUpCase(keyword))then
              ! The keyword was found, so quit searching
              keywordFound = .true.
              curLine = i
              exit
            end if  
          end if
        end do
        
        ! Return an empty string if the keyword was not found
        if(.not. keywordFound) then
          value = ''
          return
        end if
        
        ! Get the keyword value and remove extra spaces
        l = len_trim(iniLines(curLine))
        value = adjustl(trim(iniLines(curLine)(s+1:l)))
        
      end subroutine findIniValue
      
      
      
      ! subroutine getIniString()
      !
      ! Purpose: Find string data for the given keyword in [section].
      !
      ! Inputs:
      !   section: [section] in which to look for keyword
      !   keyword: name of string, given as keyword=string
      !   default: value to use if keyword is not found
      !
      ! Outputs:
      !   value: keyword value (returns '' if error)
      !   stat (optional): true if successful, otherwise false
      subroutine getIniString(section, keyword, value, default, stat)
      
        implicit none
       
        ! Input/output parameters
        character(len=*), intent(in) :: section
        character(len=*), intent(in) :: keyword
        character(len=*), intent(out) :: value
        character(len=*), intent(in) :: default
        logical, intent(out), optional :: stat
        
        ! Local function variables
        character(len=maxlen) :: valStr      ! string representation of value
      
        ! Verify that an INI file is open
        if (.not. allocated(iniLines)) then
          if(present(stat)) then
            stat = .false.
            value = ''
            return
          else
            write(*,*) 'File I/O error!'
            write(*,*) '  getIniString(): INI file not open'
            stop
          end if
        end if
           
        call findIniValue(section, keyword, valStr)
      
        if (trim(valStr) == '') then
          value = default
        else
          value = valStr
        end if
        
        ! Successful search
        if(present(stat)) stat = .true.
      
      end subroutine getIniString
      
      
      
      ! subroutine getIniInt()
      !
      ! Purpose: Find integer data for the given keyword in [section].
      !
      ! Inputs:
      !   section: [section] in which to look for keyword
      !   keyword: name of integer, given as keyword=integer
      !   default: value to use if keyword is not found
      !
      ! Outputs:
      !   value: keyword value (returns 0 if error)
      !   stat (optional): true if successful, otherwise false
      subroutine getIniInt(section, keyword, value, default, stat)
      
        implicit none
       
        ! Input/output parameters
        character(len=*), intent(in) :: section
        character(len=*), intent(in) :: keyword
        integer, intent(out) :: value
        integer, intent(in) :: default
        logical, intent(out), optional :: stat
      
        ! Local function variables
        character(len=maxlen) :: valStr        ! string representation of value
        
        ! Verify that an INI file is open
        if (.not. allocated(iniLines)) then
          write(*,*) 'File I/O warning!'
          write(*,*) '  getIniString(): INI file not open'
          write(*,*) '  returning default'
          value = default
          return
        end if
        
        call findIniValue(section, keyword, valStr)
      
        if (valStr == '') then
          value = default
        else
          ! Convert string to integer
          read(valStr, *) value 
        end if
        
        ! Successful search
        if(present(stat)) stat = .true.
        
      end subroutine getIniInt
      
      
      
      ! subroutine getIniReal()
      !
      ! Purpose: Find floating point data for the given keyword in [section].
      !
      ! Inputs:
      !   section: [section] in which to look for keyword
      !   keyword: name of real value, given as keyword=real
      !   default: value to use if keyword is not found
      !
      ! Outputs:
      !   value: keyword value (returns 0.0 if error)
      !   stat (optional): true if successful, otherwise false
      subroutine getIniReal(section, keyword, value, default, stat)
      
        implicit none
       
        ! Input/output parameters
        character(len=*), intent(in) :: section
        character(len=*), intent(in) :: keyword
        double precision, intent(out) :: value
        double precision, intent(in) :: default
        logical, intent(out), optional:: stat
      
        ! Local function variables
        character(len=maxlen) :: valStr        ! string representation of value
        
        ! Verify that an INI file is open
        if (.not. allocated(iniLines)) then
          write(*,*) 'File I/O warning!'
          write(*,*) '  getIniString(): INI file not open'
          write(*,*) '  returning default'
          value = default
          return
        end if
         
        call findIniValue(section, keyword, valStr)
      
        if (valStr == '') then
          value = default
        else
          ! Convert string to real
          read(valStr, *) value 
        end if
        
        ! Successful search
        if(present(stat)) stat = .true.
        
      end subroutine getIniReal
      
      
      
!      ! subroutine getIniBool()
!      !
!      ! Purpose: Find boolean data for the given keyword in [section].
!      !
!      ! Inputs:
!      !   section: [section] in which to look for keyword
!      !   keyword: name of boolean value, given as keyword=boolean
!      !   default: value to use if keyword is not found
!      !
!      ! Outputs:
!      !   value: keyword value (returns .false. if error)
!      !   stat (optional): true if successful, otherwise false
!      subroutine getIniBool(section, keyword, value, default, stat)
!      
!        implicit none
!       
!        ! Input/output parameters
!        character(len=*), intent(in) :: section
!        character(len=*), intent(in) :: keyword
!        logical, intent(out) :: value
!        logical, intent(in) :: default
!        logical, intent(out), optional :: stat
!      
!        ! Local function variables
!        character(len=maxlen) :: valStr        ! string representation of value
!        
!        ! Verify that an INI file is open
!        if (.not. allocated(iniLines)) then
!          if(present(stat)) then
!            stat = .false.
!            value = ''
!            return
!          else
!            write(*,*) 'File I/O error!'
!            write(*,*) '  getIniString(): INI file not open'
!            stop
!          end if
!        end if
!         
!        call findIniValue(section, keyword, valStr)
!      
!        if (valStr == '') then
!          value = default
!        else
!          ! add periods on both sides
!          valStr = '.'//trim(valStr)//'.'
!          ! Convert string to logical
!  strange error here          read(valStr, *) value 
!        end if
!        
!        ! Successful search
!        if(present(stat)) stat = .true.
!        
!      end subroutine getIniBool
      
      ! subroutine getIniPetscTruth()
      !
      ! Purpose: Find PetscTruth data for the given keyword in [section].
      !
      ! Inputs:
      !   section: [section] in which to look for keyword
      !   keyword: name of PetscTruth value, given as keyword=PetscTruth
      !   default: value to use if keyword is not found
      !
      ! Outputs:
      !   value: keyword value (returns .false. if error)
      !   stat (optional): true if successful, otherwise false
      subroutine getIniPetscTruth(section, keyword, value)
      
        implicit none
#include "finclude/petscdef.h"
       
        ! Input/output parameters
        character(len=*), intent(in) :: section
        character(len=*), intent(in) :: keyword
        PetscTruth, intent(out) :: value
        PetscTruth, external :: stringcomp
      
        ! Local function variables
        character(len=maxlen) :: valStr        ! string representation of value
        
        ! Verify that an INI file is open
        if (.not. allocated(iniLines)) then
            write(*,*) 'File I/O error!'
            write(*,*) '  getIniString(): INI file not open'
            stop
        end if
         
        call findIniValue(section, keyword, valStr)
      
        value = stringcomp(trim(valStr)//char(0),"true"//char(0))

      end subroutine getIniPetscTruth
      
      end module parse_ini
