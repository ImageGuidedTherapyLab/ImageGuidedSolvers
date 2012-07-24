! module string_utils
!
! Purpose: contains methods for converting strings to uppercase and lowercase
!   in a portable fashion. Adapted from: 
!     Figure 3.5B, pg 80, "Upgrading to Fortran 90", by Cooper Redwine,
!     1995 Springer-Verlag, New York. 
!
! Author: Cooper Redwine, adapted by Jason Meiring
! Date: 2005/10/07
      module string_utils
      
      implicit none
      
      ! Scope settings
      private
      public :: strUpCase, strLowCase
      
      ! Module constants
      character(len=*), private, parameter ::  &
                        lower_case = 'abcdefghijklmnopqrstuvwxyz'
      character(len=*), private, parameter ::  &
                        upper_case = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      
      contains
      
      ! function strLowCase()
      !
      ! Purpose: convert a string to all lowercase
      ! 
      ! Input:
      !   inputString: string to convert to lowercase
      !
      ! Output:
      !   outputString: lowercase version of inputString
      function strUpCase (inputString) result (outputString)
        
        implicit none
        
        ! Input parameters
        character(len=*), intent(in) :: inputString
           
        ! Output parameters
        character(len(inputString)) :: outputString
           
        ! Local function variables
        integer :: i, n
      
        ! Copy input string
        outputString = inputString
           
        ! Loop over string elements
        do i = 1, len(outputString )
             
          ! Find location of letter in lower case constant string
          n = index(lower_case, outputString(i:i))
             
          ! If current substring is a lower case letter, make it upper case
          if (n /= 0) outputString(i:i) = upper_case(n:n)
        end do
           
      end function strUpCase
      
      
      
      ! function strLowCase()
      !
      ! Purpose: convert a string to all lowercase
      ! 
      ! Input:
      !   inputString: string to convert to lowercase
      !
      ! Output:
      !   outputString: lowercase version of inputString
      function strLowCase (inputString) result (outputString)
      
        implicit none
      
        ! Input parameters
        character(len=*), intent(in) :: inputString
        
        ! Output parameters
        character(len(inputString)) :: outputString
        
        ! Local function variables
        integer :: i, n
      
        ! Copy input string
        outputString = inputString
        
        ! Loop over string elements
        do i = 1, len(outputString )
        
          ! Find location of letter in upper case constant string
          n = index(upper_case, outputString(i:i))
             
          ! If current substring is an upper case letter, make it lower case
          if (n /= 0) outputString(i:i) = lower_case(n:n)
        
        end do
        
      end function strLowCase
      
      end module string_utils 
