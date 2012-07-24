#include "finclude/petscdef.h"
! ensure same default values between C and Fortran
#include "parser_defaults.h" 
! ensure same default values between C and Fortran
#include "mri_params.h"
! ensure same default values between C and Fortran
#include "variable_map.h"
!----------------------------------------------------------------------
!
!   module  - pennes_model     (latest revision: Aug 08)
!
!   purpose - store constitutive data for Pennes model with MRTI for
!             the inverse problem
!
!----------------------------------------------------------------------
      module pennes_model
      implicit none
! the save attribute specifies that upon use of the module all values 
! of the variables remain as they were upon the previous instance of
! the module
      save
      ! element-wise error estimates and adjoint gradient contributions
      integer, parameter ::  NUMHPELEMVAR = 4
      PetscScalar,allocatable :: HPELMFLD(:,:)
      ! element-wise blood perfusion and thermal conductivity
      ! the 0th entry in W_0 and K_0 are used to store the part of 
      ! the field that is constant amongst elements (if it exists)
      ! if there is no constant portion of the feild the first entry
      ! just gives the first elements W_0 and K_0
      PetscScalar, allocatable , dimension(:) :: W_0 , K_0 
      PetscInt  :: ID_W_0 , ID_K_0 
!     Module variables
      PetscScalar :: RHO,C_P                 ! tissue density and specific heat
      PetscScalar :: C_BLOOD, U_A            ! specific heat of blood & arterial temp
      PetscScalar :: COEFF_COOL, U_INF       ! coeff of cooling & ambient temp for cauchy bc
      PetscScalar :: G_FLUX                  ! heat flux for neumann bc
      PetscScalar :: K_1,K_2,K_3             ! k(u,x) = K_0(x) + K_1 * 2/pi* atan(K_2(u-K_3))
      PetscScalar :: W_N,W_I ,W_D            ! w(u,x) = W_0(x) + (W_N+W_D)/2 +
      PetscScalar :: W_2,W_NI,W_ID           !  2/pi * [ (W_I-W_N)/2*atan(W_2(u-W_NI)) 
                                             !          -(W_I-W_D)/2*atan(W_2(u-W_ID)) ]
      PetscScalar :: X_0, Y_0, Z_0           ! position of laser
      PetscScalar :: U_INIT,U_PROBE          ! variable for initial conditions
      PetscScalar :: MU_A,MU_S,ANFACT,MU_TR,MU_EFF     ! laser parameter coeff 
      PetscScalar :: S_0,S_1,S_2,S_3         ! s(u) = S_0(x) + S_1 * 2/pi* atan(S_2(u-S_3))
      PetscScalar :: K_0_LB,K_1_LB,K_2_LB,K_3_LB       ! bounds on parameters
      PetscScalar :: K_0_UB,K_1_UB,K_2_UB,K_3_UB       ! bounds on parameters
      PetscScalar :: W_0_LB,W_0_UB                     ! bounds on parameters
      PetscScalar :: W_N_LB,W_I_LB,W_D_LB,W_2_LB       ! bounds on parameters
      PetscScalar :: W_N_UB,W_I_UB,W_D_UB,W_2_UB       ! bounds on parameters
      PetscScalar :: W_NID_LB, W_NID_MD, W_NID_UB      ! bounds on parameters
      PetscScalar :: X_0_LB, Y_0_LB, Z_0_LB, POW_LB    ! bounds on parameters 
      PetscScalar :: X_0_UB, Y_0_UB, Z_0_UB, POW_UB    ! bounds on parameters
      PetscScalar :: MU_A_LB,MU_A_UB,MU_S_LB,MU_S_UB   ! bounds on parameters
      ! perfusion field and thermal conductivity field parameters
      PetscTruth  :: K_0_FIELD  ! TRUE ==> let K_0 vary element wise
      PetscTruth  :: W_0_FIELD  ! TRUE ==> let W_0 vary element wise
      PetscInt :: NEXACTVERIFNUMBER ! test suite id

      ! Arrhenius damage model coefficients
      PetscScalar :: Arr_R,Arr_A,Arr_Ea

      ! Two State damage model coefficients
      PetscScalar :: TS_h,TS_alpha,TS_beta

      type :: POWER
      ! power as a function of time
          PetscScalar :: time
          PetscScalar :: power
      end type POWER
      type(POWER) , allocatable , dimension(:):: POW
      ! the data structure is setup as follows and ASSUMES equispace time 
      !     distances. ie POW(i)%time - POW(i-1)%time = IDEAL_DT   \forall i
      !  
      ! 
      !               t = 0      POW(0)%power SHOULD NEVER BE USED POW(0)%power 
      !                 |        is set to 1.0d0 for plotting purposes
      !                 |
      !                 |
      !           POW(1)%power      the power between t = 0 and 
      !                 |            t = POW(1)%time is POW(1)%power
      !                 |
      ! POW(1)%time ---------
      !                 |
      !                 |           the power between t = POW(1)%time and 
      !           POW(2)%power       t = POW(2)%time is POW(2)%power
      !                 |
      !                 |
      ! POW(2)%time ---------
      !                 |
      !                 |           the power between t = POW(2)%time and 
      !           POW(3)%power       t = POW(3)%time is POW(3)%power
      !                 |
      !                 |
      ! POW(3)%time ---------
      !           .
      !           .
      !           .
      !           .
      ! create a buffer to write the avs field file
      PetscInt :: IDEAL_NZERO,IDEAL_NTIME  ! time slice of current optimization
      PetscInt :: ISTEPS_PER_IDEAL ! fem steps per ideal step
      PetscScalar ::  IDEAL_DT  ! time step of ideal data

      PetscScalar:: TAU  ! final time stored for verif problems only

      ! characteristic function data for goals algorithm
      PetscScalar :: CHAR_X_0, CHAR_Y_0, CHAR_Z_0, CHAR_RAD

      ! ideal function data temp_control, dam_control, hsp_control
      PetscScalar :: IDEAL_X_0, IDEAL_Y_0, IDEAL_Z_0, IDEAL_RAD,  &
                     IDEAL__IN, IDEAL_OUT 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      end module pennes_model
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      ! debugging utility f90 not working well with gdb
      ! from debugger call debug_pennes_model()
      subroutine debug_pennes_model() 
        use pennes_model
        write(*,*) "POW               =" 
        do i=0,size(POW)-1
          write(*,'(e12.6,x,e12.6)') POW(i)%time,POW(i)%power
        end do
        write(*,*) "K_0               =", K_0
        write(*,*) "W_0               =", W_0 
        write(*,*) "ID_W_0            =", ID_W_0 
        write(*,*) "ID_K_0            =", ID_K_0 
        write(*,*) "RHO               =", RHO
        write(*,*) "C_P               =", C_P                 
        write(*,*) "C_BLOOD           =", C_BLOOD
        write(*,*) "U_A               =", U_A            
        write(*,*) "COEFF_COOL        =", COEFF_COOL
        write(*,*) "U_INF             =", U_INF       
        write(*,*) "G_FLUX            =", G_FLUX               
        write(*,*) "K_1               =", K_1
        write(*,*) "K_2               =", K_2
        write(*,*) "K_3               =", K_3             
        write(*,*) "W_N               =", W_N
        write(*,*) "W_I               =", W_I 
        write(*,*) "W_D               =", W_D            
        write(*,*) "W_2               =", W_2
        write(*,*) "W_NI              =", W_NI
        write(*,*) "W_ID              =", W_ID           
        write(*,*) "X_0               =", X_0
        write(*,*) "Y_0               =", Y_0
        write(*,*) "Z_0               =", Z_0           
        write(*,*) "U_INIT            =", U_INIT 
        write(*,*) "MU_A              =", MU_A
        write(*,*) "MU_S              =", MU_S
        write(*,*) "ANFACT            =", ANFACT
        write(*,*) "MU_TR             =", MU_TR
        write(*,*) "MU_EFF            =", MU_EFF    
        write(*,*) "K_0_LB            =", K_0_LB
        write(*,*) "K_1_LB            =", K_1_LB
        write(*,*) "K_2_LB            =", K_2_LB
        write(*,*) "K_3_LB            =", K_3_LB      
        write(*,*) "K_0_UB            =", K_0_UB
        write(*,*) "K_1_UB            =", K_1_UB
        write(*,*) "K_2_UB            =", K_2_UB
        write(*,*) "K_3_UB            =", K_3_UB      
        write(*,*) "W_0_LB            =", W_0_LB
        write(*,*) "W_0_UB            =", W_0_UB
        write(*,*) "W_N_LB            =", W_N_LB
        write(*,*) "W_I_LB            =", W_I_LB
        write(*,*) "W_D_LB            =", W_D_LB
        write(*,*) "W_2_LB            =", W_2_LB      
        write(*,*) "W_N_UB            =", W_N_UB
        write(*,*) "W_I_UB            =", W_I_UB
        write(*,*) "W_D_UB            =", W_D_UB
        write(*,*) "W_2_UB            =", W_2_UB      
        write(*,*) "W_NID_LB          =", W_NID_LB
        write(*,*) "W_NID_MD          =", W_NID_MD
        write(*,*) "W_NID_UB          =", W_NID_UB     
        write(*,*) "X_0_LB            =", X_0_LB
        write(*,*) "Y_0_LB            =", Y_0_LB
        write(*,*) "Z_0_LB            =", Z_0_LB
        write(*,*) "POW_LB            =", POW_LB   
        write(*,*) "X_0_UB            =", X_0_UB
        write(*,*) "Y_0_UB            =", Y_0_UB
        write(*,*) "Z_0_UB            =", Z_0_UB
        write(*,*) "POW_UB            =", POW_UB   
        write(*,*) "MU_A_LB           =", MU_A_LB
        write(*,*) "MU_A_UB           =", MU_A_UB
        write(*,*) "MU_S_LB           =", MU_S_LB
        write(*,*) "MU_S_UB           =", MU_S_UB  
        write(*,*) "K_0_FIELD         =", K_0_FIELD  
        write(*,*) "W_0_FIELD         =", W_0_FIELD  
        write(*,*) "NEXACTVERIFNUMBER =", NEXACTVERIFNUMBER 
        write(*,*) "Arr_R             =", Arr_R
        write(*,*) "Arr_A             =", Arr_A
        write(*,*) "Arr_Ea            =", Arr_Ea
        write(*,*) "TS_h              =", TS_h
        write(*,*) "TS_alpha          =", TS_alpha
        write(*,*) "TS_beta           =", TS_beta
        write(*,*) "IDEAL_NZERO       =", IDEAL_NZERO
        write(*,*) "IDEAL_NTIME       =", IDEAL_NTIME  
        write(*,*) "ISTEPS_PER_IDEAL  =", ISTEPS_PER_IDEAL 
        write(*,*) "IDEAL_DT          =", IDEAL_DT  
        write(*,*) "TAU               =", TAU  
        call FLUSH()
      end subroutine debug_pennes_model
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine forwardexact(Xpoint,Time,Texact,Gradexact)
        use pennes_model
        use global_params, only:  L,R,PI
        implicit none
        PetscScalar, intent(in)  :: Xpoint(3),Time
        PetscScalar, intent(out) :: Texact,Gradexact(3)
        PetscScalar :: k_4,T_infinity,q0,usln_half,t,  &
          theta,gamma,R3(3,3),R2(3,3),graddum(3),xdum(3),radius,xcoord
        PetscInt :: iii,jjj,kkk

            radius=sqrt(Xpoint(1)**2+Xpoint(2)**2+Xpoint(3)**2)
            xcoord= Xpoint(1)
            t=Time
            q0=POW(1)%power
            !DEFAULT : exact solution not known return return ZERO
            Texact=0.0d0
            select case(NEXACTVERIFNUMBER) 
            case(1)
                Texact=(L-xcoord)**2/2.0d0+COEFF_COOL*xcoord+K_0(0)
                gradexact(1)= -L+xcoord+COEFF_COOL
                gradexact(2)=0.0d0
                gradexact(3)=0.0d0
            case(2)
                Texact=q0*(L-xcoord)**2*t+K_0(0)
                gradexact(1)=-2.0d0 * q0*t*(L-xcoord)
                gradexact(2)=0.0d0
                gradexact(3)=0.0d0
            case(22)
                theta=-PI/6.0d0;
                gamma=-PI/4.0d0;
                R3(1,1)=cos(theta)
                R3(1,2)=-sin(theta)
                R3(1,3)=0.0d0
                R3(2,1)=sin(theta)
                R3(2,2)=cos(theta)
                R3(2,3)=0.0d0
                R3(3,1)=0.0d0
                R3(3,2)=0.0d0
                R3(3,3)=1.0d0
                R2(1,1)=cos(gamma)
                R2(1,2)=0.0d0
                R2(1,3)=-sin(gamma)
                R2(2,1)=0.0d0
                R2(2,2)=1.0d0
                R2(2,3)=0.0d0
                R2(3,1)=sin(gamma)
                R2(3,2)=0.0d0
                R2(3,3)=cos(gamma)
                xdum=0.0d0
                gradexact=0.0d0
                do iii=1,3
                   do jjj=1,3
                       do kkk=1,3
                 xdum(iii)=xdum(iii)+R3(iii,kkk)*R2(kkk,jjj)*Xpoint(jjj)
                       enddo
                   enddo
                enddo
                Texact=U_INF+q0*L/COEFF_COOL + q0*L/K_0(0)*xdum(1)-        &
                     q0/2.0d0/K_0(0)*xdum(1)*xdum(1)
                  graddum(1)=q0*L/K_0(0)-q0/K_0(0)*xdum(1)
                  graddum(2)=0.0d0
                  graddum(3)=0.0d0
                do iii=1,3
                   do jjj=1,3
                       do kkk=1,3
                          gradexact(iii)=gradexact(iii)+                   &
                                    R2(kkk,iii)*R3(jjj,kkk)*graddum(jjj)
                       enddo
                   enddo
                enddo
            case(3)
                T_infinity=U_INF
                k_4=(q0*R-K_0(0)+q0*q0*R/COEFF_COOL)/T_infinity
                Texact=q0/k_4*radius-K_0(0)/k_4;
            case(4)
                T_infinity=U_INF
                k_4=(q0*q0*R**3/2.0d0/COEFF_COOL-K_0(0)+q0*R**2/2.0d0)/T_infinity
                Texact=q0/k_4*radius**2/2.0d0-K_0(0)/k_4
            case(5)
                Texact= q0/2.0d0 *Xpoint(1)*Xpoint(1) - q0*L*Xpoint(1) + U_INF
            case(6)
                Texact= radius*t/K_2+K_3
            case(7)
                Texact= (L-Xpoint(1))**2*exp(-t)+K_3
            case(8)
                Texact= (L-Xpoint(1))**2*exp(-t)+K_3
            case(9)
!                Texact= (L-xcoord)**2*exp(-t)/K_0(0)*q0*W_0
!                Texact= (L-xcoord)**2*t/K_0(0)*q0*W_0
            end select
      end subroutine forwardexact
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine adjointexact(Xpoint,Idstep,Dualexact,Gradient,Deltat)
        use pennes_model
        use global_params , ONLY: L,R
        implicit none
        PetscScalar, intent(in)  :: Xpoint(3), Deltat
        PetscInt , intent(in)  :: Idstep
        PetscScalar, intent(out) :: Dualexact,Gradient(3)
        PetscScalar :: t,xcoord,q0,D,omega

        xcoord = Xpoint(1)
        q0=POW(1)%power
            !DEFAULT : exact solution not known return return ZERO
            Dualexact= 0.0d0
            Gradient= 0.0d0
            select case(NEXACTVERIFNUMBER) 
            case(1)
               ! should be g = k dudx(L) = 0 neuman condition at x = L
               ! should be h = (1-exp(L/k)) to satisfy cauchy bc at x=0
               !     k dudx(0) = h u(0)
               dualexact = exp(xcoord/K_0(0))- 1.0d0/K_0(0)*exp(L/K_0(0))*xcoord
            case(2)
               if(Idstep.ne.0) then
                 t=(Idstep-1) * Deltat
               else
                 t=Idstep * Deltat
               endif
               dualexact = q0*(L-xcoord)**2*(TAU-t)
            case(8)
               if(Idstep.ne.0) then
                 t=(Idstep + (Idstep-1) ) * 0.5d0 * Deltat
               else
                 t=Idstep * Deltat
               endif
               dualexact = ((xpoint(1)-L)**3*COEFF_COOL/L**2/          &
                       (3.0d0*K_0(0)+3.0d0*K_1*atan(K_2/exp(t)*L**2)   &
                        +COEFF_COOL*L)+1.0d0)/COEFF_COOL*(TAU-t)
            case(9)
!               D  = RHO*C_P/K_0(0)**2*DELTAT*q0*(W_0-30.0d0)
!               omega=sqrt(2.0d0/K_0(0)*(rho*c_p/Deltat+W_0*c_blood/2.0d0))
!               dualexact = 
!     .     exp(-omega*L)*exp(omega*xcoord)+exp(omega*L)*exp(-omega*xcoor
!     #d)+(D*L**2+2*D/omega**2)/omega**2-2*D*L/omega**2*xcoord+D/omega**2
!     #*xcoord**2
            end select
      end subroutine adjointexact
! 
! the next set of routines are various function evaluations at the guass
! points to evaluate the Function , Jacobian, and Gradient of the QOI
! 
!--------------------------------------------------------------------------- 
      subroutine nonlinpennes(Deltat,Usln,Usln_prev,Adiff,Creac)
        use global_params, only: TWODPI
        use pennes_model , only: K_0,K_1,K_2,K_3,W_0,W_N,W_I,W_D,W_2, &
                                 W_NI,W_ID,ID_W_0,ID_K_0,             &
                                 RHO,C_P,U_A,C_BLOOD
        implicit none
        PetscScalar,intent(in)  :: Deltat,Usln,Usln_prev
        PetscScalar,intent(out) :: Adiff,Creac
        PetscScalar :: usln_half
        ! compute default coefficients for nonlinear pennes w/ iso laser
        usln_half  = (Usln+Usln_prev)*0.5d0
        !diffusion coefficient
        Adiff         = K_0(ID_K_0)+K_1*TWODPI*atan(K_2*(usln_half-K_3))
        !reaction TERM
        Creac=RHO*C_P*(usln-usln_prev)/deltat                             &
               +C_BLOOD*( W_0(ID_W_0) + (W_N+W_D)*0.5d0 + TWODPI *        &
                        ( (W_I-W_N)*0.5d0*atan(W_2*(usln_half-W_NI))      &
                         -(W_I-W_D)*0.5d0*atan(W_2*(usln_half-W_ID)) ) )  &
                                                      *(usln_half - U_A)
      end subroutine nonlinpennes
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine getverifformfcn(Xpoint,Idstep,Deltat,Soln,Nsize,     &
                                                  Adiff,Creac,Source)
        use global_params , only: L,R
        use pennes_model
        implicit none
        PetscScalar,intent(in) :: Xpoint(3),Deltat,Soln(0:Nsize-1)
        PetscInt,intent(in)    :: Idstep,Nsize
        PetscScalar,intent(out):: Adiff,Creac,Source
        PetscInt :: idpow
        PetscScalar :: t_infinity,k_4,radius,q0,t,usln_half,xcoord,dist, &
                       usln,usln_prev

        call nonlinpennesisolaser(Xpoint,Idstep,Deltat,Soln,Nsize,     &
                                                  Adiff,Creac,Source)
        ! compute pennes coefficient
        Usln     = Soln(0)
        Usln_prev= Soln(1)

        ! modify for verification suite
        xcoord=Xpoint(1)
        radius=sqrt(Xpoint(1)**2+Xpoint(2)**2+Xpoint(3)**2)
        t=0.5d0*( Idstep+Idstep-1 )  * Deltat
        usln_half=(Usln+Usln_prev)/2.0d0
        q0 = POW(1)%power
        select case(NEXACTVERIFNUMBER) 
        case(1)
            Adiff = K_0(0)
            Creac=0.0d0
            Source=-K_0(0)
        case(2)
            Adiff = K_0(0)
            Creac = RHO*C_P*(Usln-Usln_prev)/DELTAT
            Source= RHO*C_P*q0*(L-xcoord)**2-2.0d0*K_0(0)*q0*t
        case(3)
            t_infinity=U_INF
            k_4=(q0*R-K_0(0)+q0*q0*R/COEFF_COOL)/t_infinity
            Source=-q0*q0*3.0d0/k_4
            Adiff=K_0(0)+k_4*usln_half
            Creac=0.0d0
        case(4)
            t_infinity=U_INF
            k_4=(q0*q0*R**3/2.0d0/COEFF_COOL-K_0(0)+q0*R**2/2.0d0)/t_infinity
            Source=-5.0d0/2.0d0*q0*q0/k_4*                          &
                           (Xpoint(1)**2+Xpoint(2)**2+Xpoint(3)**2)
            Adiff=K_0(0)+k_4*usln_half
            Creac=0.0d0
        case(5)
            t_infinity=U_INF
            k_4=-K_0(0)/t_infinity
            Source=-k_4*(q0*Xpoint(1)-q0*L)**2-(k_4*(q0/2.0d0*          &
                     Xpoint(1)**2-q0*L*Xpoint(1)+t_infinity)+K_0(0))*q0
            Adiff=K_0(0)+k_4*usln_half
            Creac=0.0d0
!        case(6)
!            Source= rho*c_p/K_2*radius-1.0d0/radius**2*
!     .            (2.0d0*radius*(K_0(0)+K_1*atan(radius*t))/K_2*t+
!     .             radius**2*K_1/K_2*t**2/(1.0d0+radius**2*t**2))+
!     .             (W_0(0)+W_1*atan(W_2*(1.0d0/K_2*radius*t+K_3-W_3)))*
!     .                          c_blood*(1.0d0/K_2*radius*t+K_3-U_A)
!     .             + c_blood*U_A*(W_0(0)+W_1*atan(W_2*(usln_half-W_3)) )
!     .             + 2.0d0*rho*c_p/DELTAT*usln_prev
!        case(7)
!            Source=
!     .-rho*c_p*(L-xcoord)**2*exp(-t)-4*K_1*K_2*(L-xcoord)**2*exp(-t
!     #)**2/(1+K_2**2*(L-xcoord)**4*exp(-t)**2)-(2*K_0(0)+2*K_1*atan(K_2*
!     #(L-xcoord)**2*exp(-t)))*exp(-t)+(W_0(0)+W_1*atan(W_2*((L-xcoord)**
!     #2*exp(-t)+K_3-W_3)))*c_blood*((L-xcoord)**2*exp(-t)+K_3-U_A)
!        case(8)
!            Source= -rho*c_p*(L-Xpoint(1))**4*exp(-t)-16*K_1*K_2*(L-
!     #Xpoint(1))**6*exp(-t)**2/(1+K_2**2*(L-Xpoint(1))**8*exp(-t)**2)-(1
!     #2*K_0(0)+12*K_1*atan(K_2*(L-Xpoint(1))**4*exp(-t)))*
!     #(L-Xpoint(1))**2*exp(-t)+(W_0(0)+W_1*atan(W_2*((L-Xpoint(1))**4*
!     #exp(-t)+K_3-W_3)))*c_blood*((L-Xpoint(1))**4*exp(-t)+K_3 -U_A)
        case(9)
             Source= q0*(L-xcoord)
!            Source=-rho*c_p*(L-xcoord)**2*exp(-t)/K_0(0)*q0*W_0(0)-2*exp(-t)*
!     .        q0*W_0(0)+W_0(0)*c_blood*((L-xcoord)**2*exp(-t)/K_0(0)*q0*W_0(0)-U_A)
!            Source=rho*c_p*(L-xcoord)**2/K_0(0)*q0*W_0(0)-2*t*q0*W_0(0)+
!     .              W_0(0)*c_blood*((L-xcoord)**2*t/K_0(0)*q0*W_0(0)-U_A)
        end select
      end subroutine getverifformfcn
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine nonlinpennesisolaser(Xpoint,Idstep,Deltat,Soln,Nsize,     &
                                                  Adiff,Creac,Source)
        use global_params, only: PI
        use pennes_model , only: X_0,Y_0,Z_0,MU_A,MU_TR,MU_EFF,MU_S,  &
                                 POW,ANFACT,ISTEPS_PER_IDEAL
        implicit none
        PetscInt,intent(in)    :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Deltat,Soln(0:Nsize-1)
        PetscScalar,intent(out):: Adiff,Creac,Source
        PetscScalar :: usln_half,q0,dist, Usln,   Usln_prev
        PetscInt :: idpow

        ! compute pennes coefficient
        Usln     = Soln(0)
        Usln_prev= Soln(1)
       
        !DEC$ ATTRIBUTES INLINE :: nonlinpennes
        call nonlinpennes(Deltat,Usln,Usln_prev,Adiff,Creac)

        ! NOTE the index "idpow" is for a FORTRAN data structure
        ! INDEXING BEGINS WITH ONE!!!!!!!!!!!!!!!!!
        if( modulo(Idstep,ISTEPS_PER_IDEAL) .ne. 0) then
           idpow= Idstep / ISTEPS_PER_IDEAL  + 1
        else 
           idpow= Idstep / ISTEPS_PER_IDEAL 
        endif
        q0 = POW( idpow )%power
        !SOURCE TERM
        dist=sqrt(  (xpoint(1)-X_0)**2+ (xpoint(2)-Y_0)**2+ (xpoint(3)-Z_0)**2)
        MU_TR=MU_A+MU_S*(1.0d0-ANFACT)
        MU_EFF=sqrt(3.0d0*MU_A*MU_TR)
        Source=0.75d0*q0*MU_A*MU_TR*exp(-MU_EFF*dist)/PI/dist
      end subroutine nonlinpennesisolaser
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine nonlinpennesmonte(Xpoint,Idstep,Deltat,Soln,Nsize,     &
                                                  Adiff,Creac,Source)
        use global_params, only: PI
        use pennes_model , only: X_0,Y_0,Z_0,MU_A,MU_TR,MU_EFF,MU_S,  &
                                 POW,ANFACT,ISTEPS_PER_IDEAL
        implicit none
        PetscInt,intent(in)    :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Deltat,Soln(0:Nsize-1)
        PetscScalar,intent(out):: Adiff,Creac,Source
        PetscScalar :: laserorigin(0:2),q0,usln,usln_prev
        PetscInt :: idpow

        ! compute pennes coefficient
        Usln     = Soln(0)
        Usln_prev= Soln(1)
        ! inline
        call nonlinpennes(Deltat,Usln,Usln_prev,Adiff,Creac)

        ! NOTE the index "idpow" is for a FORTRAN data structure
        ! INDEXING BEGINS WITH ONE!!!!!!!!!!!!!!!!!
        if( modulo(Idstep,ISTEPS_PER_IDEAL) .ne. 0) then
           idpow= Idstep / ISTEPS_PER_IDEAL  + 1
        else 
           idpow= Idstep / ISTEPS_PER_IDEAL 
        endif
        q0 = POW( idpow )%power
        laserorigin(0)=X_0;
        laserorigin(1)=Y_0;
        laserorigin(2)=Z_0;
        !SOURCE TERM
        !Source= q0 * mcheat_nopt(Xpoint,laserorigin)
        Source= q0 
      end subroutine nonlinpennesmonte
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine penaltytemperature(Xpoint,Idstep,Deltat,Soln,Nsize,     &
                                                  Adiff,Creac,Source)
        use global_params, only: DirichletPenalty
        use pennes_model , only: U_PROBE
        implicit none
        PetscInt,intent(in)    :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Deltat,Soln(0:Nsize-1)
        PetscScalar,intent(out):: Adiff,Creac,Source
        PetscScalar :: Usln, Usln_prev, usln_half

        ! penalty term for temperature of the form
        !   \int_{\Omega_P} Penalty ( u - U_PROBE ) \phi dx = 0 
        !                                                   \forall \phi
        Usln     = Soln(0)
        Usln_prev= Soln(1)
        usln_half  = (Usln+Usln_prev)*0.5d0

        ! inline
        Adiff  =   0.0
        Creac  = usln_half * DirichletPenalty
        Source = U_PROBE   * DirichletPenalty 
      end subroutine penaltytemperature
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine penaltyvoltage(Xpoint,Idstep,Deltat,Soln,Nsize,     &
                                                  Adiff,Creac,Source)
        use global_params, only: DirichletPenalty
        use pennes_model , only: ISTEPS_PER_IDEAL  ,POW
        implicit none
        PetscInt,intent(in)    :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Deltat,Soln(0:Nsize-1)
        PetscScalar,intent(out):: Adiff,Creac,Source
        PetscScalar :: voltage,Volt_cur
        PetscInt    :: idpow
        ! penalty term for voltage (z)  of the form
        ! 
        ! \int_{\Omega_P} Penalty ( Creac * z - voltage ) \psi  dx = 0
        !                                                      \forall \psi

        Volt_cur = Soln(5)

        Adiff  =   0.0
        Creac  =   Volt_cur * DirichletPenalty
        ! NOTE the index "idpow" is for a FORTRAN data structure
        ! INDEXING BEGINS WITH ONE!!!!!!!!!!!!!!!!!
        if( modulo(Idstep,ISTEPS_PER_IDEAL) .ne. 0) then
           idpow= Idstep / ISTEPS_PER_IDEAL  + 1
        else 
           idpow= Idstep / ISTEPS_PER_IDEAL 
        endif
        voltage = POW( idpow )%power
        Source  = voltage * DirichletPenalty 
      end subroutine penaltyvoltage
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine nonlinpennesrf(Xpoint,Idstep,Deltat,Soln,Nsize,     &
                                                  Adiff,Creac,Source)
        use global_params, only: PI
        use pennes_model , only: S_0, S_1
        implicit none
        PetscInt,intent(in)    :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Deltat,Soln(0:Nsize-1)
        PetscScalar,intent(out):: Adiff,Creac,Source
        PetscScalar :: gradz(3),usln_half,q0,dist,Usln,Usln_prev
        PetscInt :: idpow
        ! compute pennes coefficient
        Usln     = Soln(0)
        Usln_prev= Soln(1)
        gradz = Soln(7:9)
        ! inline
        call nonlinpennes(Deltat,Usln,Usln_prev,Adiff,Creac)
        ! compute default coefficients for nonlinear pennes w/ iso laser
        usln_half  = (Usln+Usln_prev)*0.5d0
        Source=(S_0 + S_1*usln_half)*(gradz(1)*gradz(1)+ &
                                      gradz(2)*gradz(2)+ &
                                      gradz(3)*gradz(3))
      end subroutine nonlinpennesrf
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine pennesidealmonte(Xpoint,Idstep,Deltat,Soln,Nsize,     &
                                                  Adiff,Creac,Source)
        use global_params, only: PI
        use pennes_model , only: X_0,Y_0,Z_0,MU_A,MU_TR,MU_EFF,MU_S,  &
                                 POW,ANFACT,ISTEPS_PER_IDEAL
        implicit none
        PetscInt,intent(in)    :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Deltat,Soln(0:Nsize-1)
        PetscScalar,intent(out):: Adiff,Creac,Source
        PetscScalar :: laserorigin(0:2),q0,usln,usln_prev
        PetscInt :: idpow

        ! compute pennes coefficient
        Usln     = Soln(0)
        Usln_prev= Soln(1)
        ! inline
        call nonlinpennes(Deltat,Usln,Usln_prev,Adiff,Creac)

        ! NOTE the index "idpow" is for a FORTRAN data structure
        ! INDEXING BEGINS WITH ONE!!!!!!!!!!!!!!!!!
        if( modulo(Idstep,ISTEPS_PER_IDEAL) .ne. 0) then
           idpow= Idstep / ISTEPS_PER_IDEAL  + 1
        else 
           idpow= Idstep / ISTEPS_PER_IDEAL 
        endif
        q0 = POW( idpow )%power
        laserorigin(0)=X_0;
        laserorigin(1)=Y_0;
        laserorigin(2)=Z_0;
        !SOURCE TERM
        !Source= q0 * mcheat_ideal(Xpoint,laserorigin)
        Source= q0 
      end subroutine pennesidealmonte
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine sourceemrf(Xpoint,Idstep,Deltat,Soln,Nsize,     &
                                                  Adiff,Creac,Source)
        use global_params, only: PI
        use pennes_model , only: S_0, S_1
        implicit none
        PetscInt,intent(in)    :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Deltat,Soln(0:Nsize-1)
        PetscScalar,intent(out):: Adiff,Creac,Source
        PetscScalar :: usln_half,usln,usln_prev
        PetscInt :: idpow
        ! retrieve solution values
        Usln     = Soln(0)
        Usln_prev= Soln(1)
        ! compute default coefficients for nonlinear pennes w/ iso laser
        usln_half  = (Usln+Usln_prev)*0.5d0
        Adiff =S_0 + S_1*usln_half
        Creac =0.0
        Source=0.0
      end subroutine sourceemrf
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine sourceisolaser(Xpoint,Idstep,Deltat,Soln,Nsize,     &
                                                  Adiff,Creac,Source)
        use global_params, only: PI
        use pennes_model , only: X_0,Y_0,Z_0,MU_A,MU_TR,MU_EFF,MU_S,  &
                                 POW,ANFACT,ISTEPS_PER_IDEAL
        implicit none
        PetscInt,intent(in)    :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Deltat,Soln(0:Nsize-1)
        PetscScalar,intent(out):: Adiff,Creac,Source
        PetscScalar :: Z_cur,q0,dist
        PetscInt :: idpow
        ! retrieve solution values
        Z_cur     = Soln(5)
        ! use this uncoupled variable to plot the laser source term
        Adiff=0.0
        Creac= Z_cur
        ! NOTE the index "idpow" is for a FORTRAN data structure
        ! INDEXING BEGINS WITH ONE!!!!!!!!!!!!!!!!!
        if( modulo(Idstep,ISTEPS_PER_IDEAL) .ne. 0) then
           idpow= Idstep / ISTEPS_PER_IDEAL  + 1
        else 
           idpow= Idstep / ISTEPS_PER_IDEAL 
        endif
        q0 = POW( idpow )%power
        !SOURCE TERM
        dist=sqrt(  (xpoint(1)-X_0)**2+ (xpoint(2)-Y_0)**2+ (xpoint(3)-Z_0)**2)
        MU_TR=MU_A+MU_S*(1.0d0-ANFACT)
        MU_EFF=sqrt(3.0d0*MU_A*MU_TR)
        Source=0.75d0*q0*MU_A*MU_TR*exp(-MU_EFF*dist)/PI/dist
      end subroutine sourceisolaser
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine heattransbc(Time,Newtoncoeff,T_infinity,Neumanflux,Temp,Penalty)
        use global_params, only: DirichletPenalty
        use pennes_model , only:  COEFF_COOL,U_INF,G_FLUX,U_INIT
        implicit none
        PetscScalar,intent(in)::Time
        PetscScalar,intent(out)::Newtoncoeff,T_infinity,Neumanflux,Temp,Penalty

        Newtoncoeff=COEFF_COOL
        T_infinity =U_INF
        Neumanflux =G_FLUX
        Temp       =U_INIT
        Penalty    =0.0
      end subroutine heattransbc
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine rfvoltagebc(Time,Newtoncoeff,V_infinity,V_flux,Volt,Penalty)
        use global_params, only: DirichletPenalty
        use pennes_model , only: ISTEPS_PER_IDEAL,POW,IDEAL_DT
        implicit none
        PetscScalar,intent(in)::Time
        PetscScalar,intent(out)::Newtoncoeff,V_infinity,V_flux,Volt,Penalty
        PetscInt :: idpow,idstep

        Newtoncoeff=DirichletPenalty
        V_infinity =0.0
        V_flux     =0.0
        idstep = Time / IDEAL_DT * ISTEPS_PER_IDEAL  
        ! NOTE the index "idpow" is for a FORTRAN data structure
        ! INDEXING BEGINS WITH ONE!!!!!!!!!!!!!!!!!
        if( modulo(Idstep,ISTEPS_PER_IDEAL) .ne. 0) then
           idpow= Idstep / ISTEPS_PER_IDEAL  + 1
        else 
           idpow= Idstep / ISTEPS_PER_IDEAL 
        endif
        Volt = POW( idpow )%power
        Penalty    =DirichletPenalty
      end subroutine rfvoltagebc
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine getverifbc(Time,Newtoncoeff,T_infinity,Neumanflux,Penalty)
        use global_params, only:  L
        use pennes_model , only:  K_0,K_1,K_2,K_3,COEFF_COOL,U_INF,    &
                                   G_FLUX,POW,NEXACTVERIFNUMBER
        implicit none
        PetscScalar,intent(in)::Time
        PetscScalar,intent(out)::Newtoncoeff,T_infinity,Neumanflux,Penalty
        PetscScalar :: t,q0

        !initialize with data input from control.ini file
        call heattransbc(Time,Newtoncoeff,T_infinity,Neumanflux,Penalty)

        !modify bc data for specific verification problem
          !         t=0.5d0*( gettime(Idstep)+gettime(Idstep-1) )
          t=Time
          q0=POW(1)%power
          select case(NEXACTVERIFNUMBER) 
          case(1)
              ! cauchy bc at x=0
              T_infinity= L**2/2.0d0+K_0(0)/COEFF_COOL*L
              ! neumann bc at x=L
              Neumanflux=-K_0(0)*COEFF_COOL
          case(2)
              T_infinity=q0*L**2*t+K_0(0)+2*K_0(0)*q0*L*t/COEFF_COOL
          case(6)
              Newtoncoeff= -(K_0(0)+K_1*atan(t))*t/(t+K_3*K_2-U_INF*K_2)
          case(7)
              Newtoncoeff=  &
      -2*(K_0(0)+K_1*atan(K_2*exp(-t)*L**2))*L*exp(-t)/(exp(-t)*L**2+K_3-U_INF)
          case(8)
              T_infinity= exp(-t)*L**2+K_3+(2.0d0*K_0(0)+                  &
               2.0d0*K_1*atan(K_2*exp(-t)*L**2))*L*exp(-t)/Newtoncoeff
          case(9)
!              T_infinity=L**2*exp(-t)/K_0(0)*q0*W_0(0)+2*L*exp(-t)*q0*
!     .                                                 W_0(0)/COEFF_COOL
!              T_infinity=t/K_0(0)*q0*W_0(0)*L**2+2*L*t*q0*W_0(0)/COEFF_COOL
          end select
      end subroutine getverifbc
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      ! inline
      subroutine nonlinpennesjac(Deltat,Usln,Usln_prev, & 
                                        Adiff,Bconv,Creac)
        use global_params, only: TWODPI
        use pennes_model , only:  K_0,K_1,K_2,K_3,W_0,W_N,W_I,W_D,W_2,  &
                                  W_NI,W_ID,ID_W_0,ID_K_0,              &
                                  RHO,C_P,C_BLOOD,U_A
        implicit none
        PetscScalar, intent(in)  :: Deltat,Usln,Usln_prev
        PetscScalar, intent(out) :: Adiff,Bconv,Creac
        PetscScalar :: usln_half,dwdu,dkdu
        usln_half = (Usln+Usln_prev)*0.5d0
        dwdu=TWODPI*                                                &
             ( (W_I/2-W_N/2)*W_2/(1+W_2**2*(usln_half-W_NI)**2) -   &
               (W_I/2-W_D/2)*W_2/(1+W_2**2*(usln_half-W_ID)**2)   )
        dkdu=TWODPI*K_1*K_2/(1.0d0+(K_2*(usln_half-K_3))**2)
        !elliptic coefficients 
        !diffusion TERM 
        Adiff= ( K_0(ID_K_0) + K_1*TWODPI*atan(K_2*(usln_half-K_3)) )*0.5d0
        !convective TERM
        Bconv         = dkdu*0.5d0
        !reaction TERM 
        Creac=RHO*C_P/deltat +(                                         &
             W_0(ID_W_0) + (W_N+W_D)*0.5d0 + TWODPI *                   &
                    ( (W_I-W_N)*0.5d0*atan(W_2*(usln_half-W_NI))        &
                     -(W_I-W_D)*0.5d0*atan(W_2*(usln_half-W_ID)) ) )    &
                                      *C_BLOOD*0.5d0                    &
                     +dwdu*C_BLOOD*(usln_half-U_A)*0.5d0
      end subroutine nonlinpennesjac
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine jacpennesrf(Deltat,Usln,Usln_prev,Gradz, &
                             Adiff,Bconv,Creac,Sigma,DsigmaDu)
        use pennes_model , only : S_0, S_1
        implicit none
        PetscScalar, intent(in)  :: Deltat,Usln,Usln_prev,Gradz(3)
        PetscScalar, intent(out) :: Adiff,Bconv,Creac,Sigma,DsigmaDu
        PetscScalar :: usln_half
        ! get solution
        usln_half = (Usln+Usln_prev)*0.5d0
        ! inline
        call nonlinpennesjac(Deltat,Usln,Usln_prev,Adiff,Bconv,Creac)
        Creac = Creac + S_1 * (Gradz(1) * Gradz(1)+ &
                               Gradz(2) * Gradz(2)+ &
                               Gradz(3) * Gradz(3)) 
        ! get solution
        Sigma    = (S_0 + S_1 * usln_half)
        DsigmaDu =  S_1 
      end subroutine jacpennesrf
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine getverifjac(Deltat,Usln,Usln_prev, & 
                                    Adiff,Bconv,Creac)
        use global_params, only:  L,R
        use pennes_model
        implicit none
        PetscScalar, intent(in)  :: Deltat,Usln,Usln_prev
        PetscScalar, intent(out) :: Adiff,Bconv,Creac
        PetscScalar :: k_4,T_infinity,q0,usln_half

        usln_half = (Usln+Usln_prev)*0.5d0
        !initialize default values
        call nonlinpennesjac(Deltat,Usln,Usln_prev,Adiff,Bconv,Creac)

        !modify for verification problems
        select case(NEXACTVERIFNUMBER) 
        case(1)
             Adiff=K_0(0)
             Bconv=0.0d0
             Creac=0.0d0
        case(2)
             Adiff=K_0(0)/2.0d0
             Bconv=0.0d0
             Creac=RHO* C_P/Deltat
        case(3)
             q0 = POW(1)%power
             T_infinity=U_INF
             k_4=(q0*R-K_0(0)+q0*q0*R/COEFF_COOL)/T_infinity
             Adiff=(K_0(0)+k_4*usln_half)/2.0d0
             Bconv = k_4/2.d0
             Creac=0.0d0
        case(4)
             q0 = POW(1)%power
             T_infinity=U_INF
             k_4=(q0*q0*R**3/2.0d0/COEFF_COOL-K_0(0)+q0*R**2/2.0d0)/T_infinity
             Adiff=(K_0(0)+k_4*usln_half)/2.0d0
             Bconv = k_4/2.d0
             Creac=0.0d0
        case(5)
             T_infinity=U_INF
             k_4=-K_0(0)/T_infinity
             Adiff=(K_0(0)+k_4*usln_half)/2.0d0
             Bconv = k_4/2.d0
             Creac=0.0d0
        end select
      end subroutine getverifjac
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine getverifdualjac(Deltat,Usln,Usln_prev, & 
                                        Adiff,Bconv,Creac)
        use global_params, only:  L,R
        use pennes_model 
        implicit none
        PetscScalar, intent(in)  :: Deltat,Usln,Usln_prev
        PetscScalar, intent(out) :: Adiff,Bconv,Creac
        PetscScalar :: k_4,T_infinity,q0,usln_half

        !initialize default values
        call nonlinpennesjac(Deltat,Usln,Usln_prev,Adiff,Bconv,Creac)

        !modify for verification problems
        usln_half = (Usln+Usln_prev)/2.0d0
        select case(NEXACTVERIFNUMBER) 
        case(1)
             Adiff=K_0(0)
             Bconv=0.0d0
             !mass term was subtracted for special construction of adjoint !matrix
             !Creac=-RHO*C_P/Deltat 
             !NOT DONE ANYMORE
             Creac=0.0d0
        case(2)
             Adiff=K_0(0)/2.0d0
             Bconv=0.0d0
             !Creac=0.0d0
             Creac=RHO*C_P/Deltat 
        case(3)
             q0 = POW(1)%power
             T_infinity=U_INF
             k_4=(q0*R-K_0(0)+q0*q0*R/COEFF_COOL)/T_infinity
             Adiff=(K_0(0)+k_4*usln_half)/2.0d0
             Bconv = k_4/2.d0
             Creac=0.0d0
        case(4)
             q0 = POW(1)%power
             T_infinity=U_INF
             k_4=(q0*q0*R**3/2.0d0/COEFF_COOL-K_0(0)+q0*R**2/2.0d0)/T_infinity
             Adiff=(K_0(0)+k_4*usln_half)/2.0d0
             Bconv = k_4/2.d0
             Creac=0.0d0
        case(5)
             T_infinity=U_INF
             k_4=-K_0(0)/T_infinity
             Adiff=(K_0(0)+k_4*usln_half)/2.0d0
             Bconv = k_4/2.d0
             Creac=0.0d0
        end select
      end subroutine getverifdualjac
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine getverifbcdual(Time,Newtoncoeff,T_infty,Neumanflux)
        use global_params, only:  L,R
        use pennes_model 
        implicit none
        PetscScalar , intent(in) :: Time
        PetscScalar,intent(out) :: Newtoncoeff,T_infty,Neumanflux
        PetscScalar :: t,q0

        !initialize with data input from control.ini file
        call heattransbc(Time,Newtoncoeff,T_infty,Neumanflux)

          !t=0.5d0*( gettime(Idstep)+gettime(Idstep-1) )
          t=Time

          q0=POW(1)%power
          select case(NEXACTVERIFNUMBER) 
          case(1)
              T_infty= L**2/4.0d0+K_0(0)/COEFF_COOL*L/2.0d0
              Neumanflux=-K_0(0)*COEFF_COOL/2.0d0
          case(2)
              Newtoncoeff = -2*K_0(0)/L
          case(6)
              Newtoncoeff= -(K_0(0)+K_1*atan(t))*t/(t+K_3*K_2-U_INF*K_2)
          case(7)
              Newtoncoeff=   &
      -2*(K_0(0)+K_1*atan(K_2*exp(-t)*L**2))*L*exp(-t)/(exp(-t)*L**2+K_3-U_INF)
          case(8)
              T_infty= exp(-t)*L**2+K_3+(2.0d0*K_0(0)+                 &
               2.0d0*K_1*atan(K_2*exp(-t)*L**2))*L*exp(-t)/Newtoncoeff
          end select
      end subroutine getverifbcdual
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !Error Estimate
      function error_estimate(Xpoint,Idstep,Deltat,Soln,Nsize) 
        implicit none
        PetscInt   ,intent(in) :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Soln(0:Nsize-1),Deltat  
        PetscScalar :: Grad_u_kmhalf(3),Grad_v_k__now(3),&
                       Usln,Usln_prev,V_k__now
        PetscScalar :: Adiff,Creac,Source
        PetscScalar :: error_estimate

        ! compute default values
        call nonlinpennesisolaser(Xpoint,Idstep,Deltat,Soln,Nsize,     &
                                                  Adiff,Creac,Source)
        ! map solution
        Usln          = Soln(0)
        Usln_prev     = Soln(1)
        Grad_u_kmhalf = Soln(2:4)
        V_k__now      = Soln(5)
        Grad_v_k__now = Soln(6:8)

        ! NOTE the minus sign is b/c the pde part of the gradient is subtracted
        error_estimate=-(Adiff *( Grad_u_kmhalf(1)*Grad_v_k__now(1) +      &
                                 Grad_u_kmhalf(2)*Grad_v_k__now(2) +      &
                                 Grad_u_kmhalf(3)*Grad_v_k__now(3) ) +    &
                        (Creac - Source )* V_k__now)

      end function error_estimate
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !Error Estimate call of verification problem
      function error_estimate_verif(Xpoint,Idstep,Deltat,Soln,Nsize) 
        implicit none
        PetscInt   ,intent(in) :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Soln(0:Nsize-1),Deltat  
        PetscScalar :: Grad_u_kmhalf(3),Grad_v_k__now(3),&
                       Usln,Usln_prev,V_k__now
        PetscScalar :: Adiff,Creac,Source
        PetscScalar :: error_estimate_verif

        ! compute default values
        call getverifformfcn(Xpoint,Idstep,Deltat,Soln,Nsize,     &
                                                  Adiff,Creac,Source)
        ! map solution
        Usln          = Soln(0)
        Usln_prev     = Soln(1)
        Grad_u_kmhalf = Soln(2:4)
        V_k__now      = Soln(5)
        Grad_v_k__now = Soln(6:8)

        ! NOTE the minus sign is b/c the pde part of the gradient is subtracted
        error_estimate_verif=                                             &
                      -(Adiff *( Grad_u_kmhalf(1)*Grad_v_k__now(1) +      &
                                 Grad_u_kmhalf(2)*Grad_v_k__now(2) +      &
                                 Grad_u_kmhalf(3)*Grad_v_k__now(3) ) +    &
                       (Creac - Source )* V_k__now)

      end function error_estimate_verif
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !Gradient Evaluation K_0
      function dkdk0 (Xpoint,Idstep,Deltat,Soln,Nsize)
        use pennes_model , only:  K_0,K_1,K_2,K_3,ID_W_0,ID_K_0
        implicit none
        PetscInt   ,intent(in) :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Soln(0:Nsize-1),Deltat
        PetscScalar :: Grad_u_kmhalf(3),Grad_v_k__now(3),&
                       Usln,Usln_prev,V_k__now
        PetscScalar :: usln_half,Adiff
        PetscScalar :: dkdk0
        ! map solution
        Usln          = Soln(0)
        Usln_prev     = Soln(1)
        Grad_u_kmhalf = Soln(2:4)
        V_k__now      = Soln(5)
        Grad_v_k__now = Soln(6:8)

        usln_half  = (Usln+Usln_prev)/2.0d0
        Adiff = 1.0d0
        dkdk0 = Adiff *( Grad_u_kmhalf(1)*Grad_v_k__now(1) +      &
                         Grad_u_kmhalf(2)*Grad_v_k__now(2) +      &
                         Grad_u_kmhalf(3)*Grad_v_k__now(3) ) 
      end function dkdk0
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !Gradient Evaluation K_1
      function dkdk1 (Xpoint,Idstep,Deltat,Soln,Nsize)
        use global_params, only: TWODPI
        use pennes_model , only:  K_0,K_1,K_2,K_3,ID_W_0,ID_K_0
        implicit none
        PetscInt   ,intent(in) :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Soln(0:Nsize-1),Deltat
        PetscScalar :: Grad_u_kmhalf(3),Grad_v_k__now(3),&
                       Usln,Usln_prev,V_k__now
        PetscScalar :: usln_half,Adiff
        PetscScalar :: dkdk1
        ! map solution
        Usln          = Soln(0)
        Usln_prev     = Soln(1)
        Grad_u_kmhalf = Soln(2:4)
        V_k__now      = Soln(5)
        Grad_v_k__now = Soln(6:8)

        usln_half  = (Usln+Usln_prev)/2.0d0
        Adiff = TWODPI*atan(K_2*(usln_half-K_3))
        dkdk1 = Adiff *( Grad_u_kmhalf(1)*Grad_v_k__now(1) +      &
                         Grad_u_kmhalf(2)*Grad_v_k__now(2) +      &
                         Grad_u_kmhalf(3)*Grad_v_k__now(3) ) 
      end function dkdk1
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !Gradient Evaluation K_2
      function dkdk2 (Xpoint,Idstep,Deltat,Soln,Nsize)
        use global_params, only: TWODPI
        use pennes_model , only:  K_0,K_1,K_2,K_3,ID_W_0,ID_K_0
        implicit none
        PetscInt   ,intent(in) :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Soln(0:Nsize-1),Deltat
        PetscScalar :: Grad_u_kmhalf(3),Grad_v_k__now(3),&
                       Usln,Usln_prev,V_k__now
        PetscScalar :: usln_half,Adiff
        PetscScalar :: dkdk2
        ! map solution
        Usln          = Soln(0)
        Usln_prev     = Soln(1)
        Grad_u_kmhalf = Soln(2:4)
        V_k__now      = Soln(5)
        Grad_v_k__now = Soln(6:8)

        usln_half  = (Usln+Usln_prev)/2.0d0
        Adiff = K_1*TWODPI*(usln_half-K_3)/ &
                          (1.0d0+K_2**2*(usln_half-K_3)**2)
        dkdk2 = Adiff *( Grad_u_kmhalf(1)*Grad_v_k__now(1) +      &
                         Grad_u_kmhalf(2)*Grad_v_k__now(2) +      &
                         Grad_u_kmhalf(3)*Grad_v_k__now(3) ) 
      end function dkdk2
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !Gradient Evaluation K_3
      function dkdk3 (Xpoint,Idstep,Deltat,Soln,Nsize)
        use global_params, only: TWODPI
        use pennes_model , only:  K_0,K_1,K_2,K_3,ID_W_0,ID_K_0
        implicit none
        PetscInt   ,intent(in) :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Soln(0:Nsize-1),Deltat
        PetscScalar :: Grad_u_kmhalf(3),Grad_v_k__now(3),&
                       Usln,Usln_prev,V_k__now
        PetscScalar :: usln_half,Adiff
        PetscScalar :: dkdk3 
        ! map solution
        Usln          = Soln(0)
        Usln_prev     = Soln(1)
        Grad_u_kmhalf = Soln(2:4)
        V_k__now      = Soln(5)
        Grad_v_k__now = Soln(6:8)

        usln_half  = (Usln+Usln_prev)/2.0d0
        Adiff = -TWODPI*K_1*K_2/(1.0d0+K_2**2*(usln_half-K_3)**2)
        dkdk3 = Adiff *( Grad_u_kmhalf(1)*Grad_v_k__now(1) +      &
                         Grad_u_kmhalf(2)*Grad_v_k__now(2) +      &
                         Grad_u_kmhalf(3)*Grad_v_k__now(3) ) 
      end function dkdk3
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !Gradient Evaluation W_0
      function dwdw0(Xpoint,Idstep,Deltat,Soln,Nsize)
        use pennes_model , only:  W_0,W_N,W_I,W_D,W_2,W_NI,W_ID, &
                                  ID_W_0,C_BLOOD,U_A 
        implicit none
        PetscInt   ,intent(in) :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Soln(0:Nsize-1),Deltat
        PetscScalar :: Grad_u_kmhalf(3),Grad_v_k__now(3),&
                       Usln,Usln_prev,V_k__now
        PetscScalar :: usln_half,Creac
        PetscScalar :: dwdw0 

        ! map solution
        Usln          = Soln(0)
        Usln_prev     = Soln(1)
        Grad_u_kmhalf = Soln(2:4)
        V_k__now      = Soln(5)
        Grad_v_k__now = Soln(6:8)

        usln_half  = (Usln+Usln_prev)/2.0d0
        Creac =  1.0d0*C_BLOOD*(usln_half-U_A)
        dwdw0 =  Creac * V_k__now
      end function dwdw0
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !Gradient Evaluation W_N
      function dwdwn (Xpoint,Idstep,Deltat,Soln,Nsize)
        use global_params, only: TWODPI
        use pennes_model , only:  W_0,W_N,W_I,W_D,W_2,W_NI,W_ID, &
                                  ID_W_0,C_BLOOD,U_A
        implicit none
        PetscInt   ,intent(in) :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Soln(0:Nsize-1),Deltat
        PetscScalar :: Grad_u_kmhalf(3),Grad_v_k__now(3),&
                       Usln,Usln_prev,V_k__now
        PetscScalar :: usln_half,Creac
        PetscScalar :: dwdwn
        ! map solution
        Usln          = Soln(0)
        Usln_prev     = Soln(1)
        Grad_u_kmhalf = Soln(2:4)
        V_k__now      = Soln(5)
        Grad_v_k__now = Soln(6:8)

        usln_half  = (Usln+Usln_prev)/2.0d0
        Creac = (1.D0-TWODPI*atan(W_2*(usln_half-W_NI)) ) *0.5d0 &
                                             *C_BLOOD*(usln_half-U_A)
        dwdwn =  Creac * V_k__now
      end function dwdwn
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !Gradient Evaluation W_D
      function dwdwd (Xpoint,Idstep,Deltat,Soln,Nsize)
        use global_params, only: TWODPI
        use pennes_model , only:  W_0,W_N,W_I,W_D,W_2,W_NI,W_ID, &
                                  ID_W_0,C_BLOOD,U_A
        implicit none
        PetscInt   ,intent(in) :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Soln(0:Nsize-1),Deltat
        PetscScalar :: Grad_u_kmhalf(3),Grad_v_k__now(3),&
                       Usln,Usln_prev,V_k__now
        PetscScalar :: usln_half,Creac
        PetscScalar :: dwdwd 
        ! map solution
        Usln          = Soln(0)
        Usln_prev     = Soln(1)
        Grad_u_kmhalf = Soln(2:4)
        V_k__now      = Soln(5)
        Grad_v_k__now = Soln(6:8)

        usln_half  = (Usln+Usln_prev)/2.0d0
        Creac = (1.D0+TWODPI*atan(W_2*(usln_half-W_ID)) ) *0.5d0 &
                                           *C_BLOOD*(usln_half-U_A)
        dwdwd =  Creac * V_k__now
      end function dwdwd
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !Gradient Evaluation W_I
      function dwdwi (Xpoint,Idstep,Deltat,Soln,Nsize)
        use global_params, only: TWODPI
        use pennes_model , only:  W_0,W_N,W_I,W_D,W_2,W_NI,W_ID, &
                                  ID_W_0,C_BLOOD,U_A
        implicit none
        PetscInt   ,intent(in) :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Soln(0:Nsize-1),Deltat
        PetscScalar :: Grad_u_kmhalf(3),Grad_v_k__now(3),&
                       Usln,Usln_prev,V_k__now
        PetscScalar :: usln_half,Creac
        PetscScalar :: dwdwi
        ! map solution
        Usln          = Soln(0)
        Usln_prev     = Soln(1)
        Grad_u_kmhalf = Soln(2:4)
        V_k__now      = Soln(5)
        Grad_v_k__now = Soln(6:8)

        usln_half  = (Usln+Usln_prev)/2.0d0
        Creac = 0.5d0*TWODPI*                                              &
              ( atan(W_2*(usln_half-W_NI))-atan(W_2*(usln_half-W_ID)) )     &
                                               *C_BLOOD*(usln_half-U_A)
        dwdwi =  Creac * V_k__now
      end function dwdwi
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !Gradient Evaluation W_2
      function dwdw2 (Xpoint,Idstep,Deltat,Soln,Nsize)
        use global_params, only: TWODPI
        use pennes_model , only:  W_0,W_N,W_I,W_D,W_2,W_NI,W_ID, &
                                  ID_W_0,C_BLOOD,U_A
        implicit none
        PetscInt   ,intent(in) :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Soln(0:Nsize-1),Deltat
        PetscScalar :: Grad_u_kmhalf(3),Grad_v_k__now(3),&
                       Usln,Usln_prev,V_k__now
        PetscScalar :: usln_half,Creac
        PetscScalar :: dwdw2
        ! map solution
        Usln          = Soln(0)
        Usln_prev     = Soln(1)
        Grad_u_kmhalf = Soln(2:4)
        V_k__now      = Soln(5)
        Grad_v_k__now = Soln(6:8)

        usln_half  = (Usln+Usln_prev)/2.0d0
        Creac = TWODPI*(                                                   &
        (W_I/2-W_N/2)*(usln_half-W_NI)/(1+W_2**2*(usln_half-W_NI)**2)       &
       -(W_I/2-W_D/2)*(usln_half-W_ID)/(1+W_2**2*(usln_half-W_ID)**2) )     &
                                               *C_BLOOD*(usln_half-U_A)
        dwdw2 =  Creac * V_k__now
      end function dwdw2
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !Gradient Evaluation W_NI
      function dwdwni (Xpoint,Idstep,Deltat,Soln,Nsize)
        use global_params, only: TWODPI
        use pennes_model , only:  W_0,W_N,W_I,W_D,W_2,W_NI,W_ID, &
                                  ID_W_0,C_BLOOD,U_A
        implicit none
        PetscInt   ,intent(in) :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Soln(0:Nsize-1),Deltat
        PetscScalar :: Grad_u_kmhalf(3),Grad_v_k__now(3),&
                       Usln,Usln_prev,V_k__now
        PetscScalar :: usln_half,Creac
        PetscScalar :: dwdwni
        ! map solution
        Usln          = Soln(0)
        Usln_prev     = Soln(1)
        Grad_u_kmhalf = Soln(2:4)
        V_k__now      = Soln(5)
        Grad_v_k__now = Soln(6:8)

        usln_half  = (Usln+Usln_prev)/2.0d0
        Creac = -TWODPI*(W_I/2-W_N/2)*W_2/(1+W_2**2*(usln_half-W_NI)**2) &
                                           *C_BLOOD*(usln_half-U_A)
        dwdwni =  Creac * V_k__now
      end function dwdwni
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !Gradient Evaluation W_ID
      function dwdwid (Xpoint,Idstep,Deltat,Soln,Nsize)
        use global_params, only: TWODPI
        use pennes_model , only:  W_0,W_N,W_I,W_D,W_2,W_NI,W_ID, &
                                  ID_W_0,C_BLOOD,U_A
        implicit none
        PetscInt   ,intent(in) :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Soln(0:Nsize-1),Deltat
        PetscScalar :: Grad_u_kmhalf(3),Grad_v_k__now(3),&
                       Usln,Usln_prev,V_k__now
        PetscScalar :: usln_half,Creac
        PetscScalar :: dwdwid
        ! map solution
        Usln          = Soln(0)
        Usln_prev     = Soln(1)
        Grad_u_kmhalf = Soln(2:4)
        V_k__now      = Soln(5)
        Grad_v_k__now = Soln(6:8)

        usln_half  = (Usln+Usln_prev)/2.0d0
        Creac = -TWODPI*(W_D/2-W_I/2)*W_2/(1+W_2**2*(usln_half-W_ID)**2) &
                                           *C_BLOOD*(usln_half-U_A)
        dwdwid =  Creac * V_k__now
      end function dwdwid
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !Gradient Evaluation X_0 (isotropic source)
      function dqlaserdx (Xpoint,Idstep,Deltat,Soln,Nsize)
        use global_params , only: PI
        use pennes_model , only:  X_0,Y_0,Z_0,MU_A,MU_S,MU_TR,MU_EFF, &
                                  POW,ANFACT,ID_W_0,ID_K_0,ISTEPS_PER_IDEAL
        implicit none
        PetscInt   ,intent(in) :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Soln(0:Nsize-1),Deltat
        PetscScalar :: Grad_u_kmhalf(3),Grad_v_k__now(3),&
                       Usln,Usln_prev,V_k__now
        PetscScalar :: usln_half,Source
        PetscScalar :: dqlaserdx
        PetscScalar :: dist,q0
        PetscInt :: idpow
        ! map solution
        Usln          = Soln(0)
        Usln_prev     = Soln(1)
        Grad_u_kmhalf = Soln(2:4)
        V_k__now      = Soln(5)
        Grad_v_k__now = Soln(6:8)

            dist=sqrt((xpoint(1)-X_0)**2+       &
                      (xpoint(2)-Y_0)**2+        &
                      (xpoint(3)-Z_0)**2)
        ! NOTE the index "idpow" is for a FORTRAN data structure
        ! INDEXING BEGINS WITH ONE!!!!!!!!!!!!!!!!!
        if( modulo(Idstep,ISTEPS_PER_IDEAL) .ne. 0) then
           idpow= Idstep / ISTEPS_PER_IDEAL  + 1
        else 
           idpow= Idstep / ISTEPS_PER_IDEAL 
        endif
        q0 = POW( idpow )%power
        Source= 0.75d0*q0*MU_A*MU_TR*(MU_EFF * dist+1.0d0) & 
                      *exp(-MU_EFF*dist)/PI/dist**3*(xpoint(1)-X_0)
        dqlaserdx = - Source * V_k__now
      end function dqlaserdx
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !Gradient Evaluation Y_0 (isotropic source)
      function dqlaserdy (Xpoint,Idstep,Deltat,Soln,Nsize)
        use global_params , only: PI
        use pennes_model , only:  X_0,Y_0,Z_0,MU_A,MU_S,MU_TR,MU_EFF, &
                                  POW,ANFACT,ID_W_0,ID_K_0,ISTEPS_PER_IDEAL
        implicit none
        PetscInt   ,intent(in) :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Soln(0:Nsize-1),Deltat
        PetscScalar :: Grad_u_kmhalf(3),Grad_v_k__now(3),&
                       Usln,Usln_prev,V_k__now
        PetscScalar :: usln_half,Source
        PetscScalar :: dqlaserdy
        PetscScalar :: dist,q0
        PetscInt :: idpow
        ! map solution
        Usln          = Soln(0)
        Usln_prev     = Soln(1)
        Grad_u_kmhalf = Soln(2:4)
        V_k__now      = Soln(5)
        Grad_v_k__now = Soln(6:8)

            dist=sqrt((xpoint(1)-X_0)**2+       &
                      (xpoint(2)-Y_0)**2+        &
                      (xpoint(3)-Z_0)**2)
        ! NOTE the index "idpow" is for a FORTRAN data structure
        ! INDEXING BEGINS WITH ONE!!!!!!!!!!!!!!!!!
        if( modulo(Idstep,ISTEPS_PER_IDEAL) .ne. 0) then
           idpow= Idstep / ISTEPS_PER_IDEAL  + 1
        else 
           idpow= Idstep / ISTEPS_PER_IDEAL 
        endif
        q0 = POW( idpow )%power
        Source= 0.75d0*q0*MU_A*MU_TR*(MU_EFF * dist+1.0d0) &
                    *exp(-MU_EFF*dist)/PI/dist**3*(xpoint(2)-Y_0)
        dqlaserdy = - Source * V_k__now
      end function dqlaserdy
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !Gradient Evaluation Z_0 (isotropic source)
      function dqlaserdz (Xpoint,Idstep,Deltat,Soln,Nsize)
        use global_params , only: PI
        use pennes_model , only:  X_0,Y_0,Z_0,MU_A,MU_S,MU_TR,MU_EFF, &
                                  POW,ANFACT,ID_W_0,ID_K_0,ISTEPS_PER_IDEAL
        implicit none
        PetscInt   ,intent(in) :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Soln(0:Nsize-1),Deltat
        PetscScalar :: Grad_u_kmhalf(3),Grad_v_k__now(3),&
                       Usln,Usln_prev,V_k__now
        PetscScalar :: usln_half,Source
        PetscScalar :: dqlaserdz
        PetscScalar :: dist,q0
        PetscInt :: idpow
        ! map solution
        Usln          = Soln(0)
        Usln_prev     = Soln(1)
        Grad_u_kmhalf = Soln(2:4)
        V_k__now      = Soln(5)
        Grad_v_k__now = Soln(6:8)

            dist=sqrt((xpoint(1)-X_0)**2+  &
                      (xpoint(2)-Y_0)**2+   &
                      (xpoint(3)-Z_0)**2)
        ! NOTE the index "idpow" is for a FORTRAN data structure
        ! INDEXING BEGINS WITH ONE!!!!!!!!!!!!!!!!!
        if( modulo(Idstep,ISTEPS_PER_IDEAL) .ne. 0) then
           idpow= Idstep / ISTEPS_PER_IDEAL  + 1
        else 
           idpow= Idstep / ISTEPS_PER_IDEAL 
        endif
        q0 = POW( idpow )%power
        Source= 0.75d0*q0*MU_A*MU_TR*(MU_EFF * dist+1.0d0) &
                   *exp(-MU_EFF*dist)/PI/dist**3*(xpoint(3)-Z_0)
        dqlaserdz = - Source * V_k__now
      end function dqlaserdz
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !Gradient Evaluation POW (isotropic source)
      function dqlaserdpow (Xpoint,Idstep,Deltat,Soln,Nsize)
        use global_params , only: PI
        use pennes_model , only:  X_0,Y_0,Z_0,MU_A,MU_S,MU_TR,MU_EFF, &
                                  POW,ANFACT,ID_W_0,ID_K_0
        implicit none
        PetscInt   ,intent(in) :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Soln(0:Nsize-1),Deltat
        PetscScalar :: Grad_u_kmhalf(3),Grad_v_k__now(3),&
                       Usln,Usln_prev,V_k__now
        PetscScalar :: usln_half,Source
        PetscScalar :: dqlaserdpow
        PetscScalar :: dist,dqdP
        ! map solution
        Usln          = Soln(0)
        Usln_prev     = Soln(1)
        Grad_u_kmhalf = Soln(2:4)
        V_k__now      = Soln(5)
        Grad_v_k__now = Soln(6:8)

            dist=sqrt((xpoint(1)-X_0)**2+  &
                      (xpoint(2)-Y_0)**2+   &
                      (xpoint(3)-Z_0)**2)
        dqdP=1.0d0 
        Source= 0.75d0*dqdP*MU_A*MU_TR*exp(-MU_EFF*dist)/PI/dist
        dqlaserdpow = - Source * V_k__now
      end function dqlaserdpow
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !Gradient Evaluation MU_A
      function dqlaserdmu_a (Xpoint,Idstep,Deltat,Soln,Nsize)
        use global_params , only: PI
        use pennes_model , only:  X_0,Y_0,Z_0,MU_A,MU_S,MU_TR,MU_EFF, &
                                  POW,ANFACT,ID_W_0,ID_K_0,ISTEPS_PER_IDEAL
        implicit none
        PetscInt   ,intent(in) :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Soln(0:Nsize-1),Deltat
        PetscScalar :: Grad_u_kmhalf(3),Grad_v_k__now(3),&
                       Usln,Usln_prev,V_k__now
        PetscScalar :: usln_half,Source
        PetscScalar :: dqlaserdmu_a
        PetscScalar :: dist,q0
        PetscInt :: idpow
        ! map solution
        Usln          = Soln(0)
        Usln_prev     = Soln(1)
        Grad_u_kmhalf = Soln(2:4)
        V_k__now      = Soln(5)
        Grad_v_k__now = Soln(6:8)

          dist=sqrt((xpoint(1)-X_0)**2+   &
                    (xpoint(2)-Y_0)**2+    &
                    (xpoint(3)-Z_0)**2)
        ! NOTE the index "idpow" is for a FORTRAN data structure
        ! INDEXING BEGINS WITH ONE!!!!!!!!!!!!!!!!!
        if( modulo(Idstep,ISTEPS_PER_IDEAL) .ne. 0) then
           idpow= Idstep / ISTEPS_PER_IDEAL  + 1
        else 
           idpow= Idstep / ISTEPS_PER_IDEAL 
        endif
        q0 = POW( idpow )%power
        MU_TR=MU_A+MU_S*(1.0d0-ANFACT)
        MU_EFF=sqrt(3.0d0*MU_A*MU_TR)
        Source= 0.75d0*q0*exp(-MU_EFF*dist)/PI/dist* &
                   (MU_TR+MU_A-MU_TR*MU_A*dist*1.5d0*(MU_A+MU_TR)/MU_EFF)
        dqlaserdmu_a = - Source * V_k__now
      end function dqlaserdmu_a
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !Gradient Evaluation MU_S
      function dqlaserdmu_s (Xpoint,Idstep,Deltat,Soln,Nsize)
        use global_params , only: PI
        use pennes_model , only:  X_0,Y_0,Z_0,MU_A,MU_S,MU_TR,MU_EFF, &
                                  POW,ANFACT,ID_W_0,ID_K_0,ISTEPS_PER_IDEAL
        implicit none
        PetscInt   ,intent(in) :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Soln(0:Nsize-1),Deltat
        PetscScalar :: Grad_u_kmhalf(3),Grad_v_k__now(3),&
                       Usln,Usln_prev,V_k__now
        PetscScalar :: usln_half,Source
        PetscScalar :: dqlaserdmu_s
        PetscScalar :: dist,q0
        PetscInt :: idpow
        ! map solution
        Usln          = Soln(0)
        Usln_prev     = Soln(1)
        Grad_u_kmhalf = Soln(2:4)
        V_k__now      = Soln(5)
        Grad_v_k__now = Soln(6:8)

          dist=sqrt((xpoint(1)-X_0)**2+   &
                    (xpoint(2)-Y_0)**2+    &
                    (xpoint(3)-Z_0)**2)
        ! NOTE the index "idpow" is for a FORTRAN data structure
        ! INDEXING BEGINS WITH ONE!!!!!!!!!!!!!!!!!
        if( modulo(Idstep,ISTEPS_PER_IDEAL) .ne. 0) then
           idpow= Idstep / ISTEPS_PER_IDEAL  + 1
        else 
           idpow= Idstep / ISTEPS_PER_IDEAL 
        endif
        q0 = POW( idpow )%power
        MU_TR=MU_A+MU_S*(1.0d0-ANFACT)
        MU_EFF=sqrt(3.0d0*MU_A*MU_TR)
        Source= 0.75d0*q0*exp(-MU_EFF*dist)/PI/dist* &
                  MU_A*(1.0d0-ANFACT)*(1.0d0-1.5d0*MU_TR*dist*MU_A/MU_EFF)
        dqlaserdmu_s = - Source * V_k__now
      end function dqlaserdmu_s
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !Gradient Evaluation X_0 (monte carlo source)
      function dqmontedx (Xpoint,Idstep,Deltat,Soln,Nsize)
        use global_params , only: PI
        use pennes_model , only:  X_0,Y_0,Z_0,MU_A,MU_S,MU_TR,MU_EFF, &
                                  POW,ANFACT,ID_W_0,ID_K_0,ISTEPS_PER_IDEAL
        implicit none
        PetscInt   ,intent(in) :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Soln(0:Nsize-1),Deltat
        PetscScalar :: Grad_u_kmhalf(3),Grad_v_k__now(3),&
                       Usln,Usln_prev,V_k__now
        PetscScalar :: q0,laserorigin(0:2),Source
        PetscScalar :: dqmontedx
        PetscInt :: idpow,izero
        PetscScalar , external :: dmcheat_noptdvector 
        ! map solution
        Usln          = Soln(0)
        Usln_prev     = Soln(1)
        Grad_u_kmhalf = Soln(2:4)
        V_k__now      = Soln(5)
        Grad_v_k__now = Soln(6:8)

        laserorigin(0)=X_0;
        laserorigin(1)=Y_0;
        laserorigin(2)=Z_0;
        ! NOTE the index "idpow" is for a FORTRAN data structure
        ! INDEXING BEGINS WITH ONE!!!!!!!!!!!!!!!!!
        if( modulo(Idstep,ISTEPS_PER_IDEAL) .ne. 0) then
           idpow= Idstep / ISTEPS_PER_IDEAL  + 1
        else 
           idpow= Idstep / ISTEPS_PER_IDEAL 
        endif
        q0 = POW( idpow )%power
        izero=0
        !Source= q0*dmcheat_noptdvector(Xpoint,laserorigin,izero)
        dqmontedx = - Source * V_k__now
      end function dqmontedx
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !Gradient Evaluation Y_0 (monte carlo source)
      function dqmontedy (Xpoint,Idstep,Deltat,Soln,Nsize)
        use global_params , only: PI
        use pennes_model , only:  X_0,Y_0,Z_0,MU_A,MU_S,MU_TR,MU_EFF, &
                                  POW,ANFACT,ID_W_0,ID_K_0,ISTEPS_PER_IDEAL
        implicit none
        PetscInt   ,intent(in) :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Soln(0:Nsize-1),Deltat
        PetscScalar :: Grad_u_kmhalf(3),Grad_v_k__now(3),&
                       Usln,Usln_prev,V_k__now
        PetscScalar :: q0,laserorigin(0:2),Source
        PetscScalar :: dqmontedy
        PetscInt :: idpow,ione
        PetscScalar , external :: dmcheat_noptdvector 
        ! map solution
        Usln          = Soln(0)
        Usln_prev     = Soln(1)
        Grad_u_kmhalf = Soln(2:4)
        V_k__now      = Soln(5)
        Grad_v_k__now = Soln(6:8)

        laserorigin(0)=X_0;
        laserorigin(1)=Y_0;
        laserorigin(2)=Z_0;
        ! NOTE the index "idpow" is for a FORTRAN data structure
        ! INDEXING BEGINS WITH ONE!!!!!!!!!!!!!!!!!
        if( modulo(Idstep,ISTEPS_PER_IDEAL) .ne. 0) then
           idpow= Idstep / ISTEPS_PER_IDEAL  + 1
        else 
           idpow= Idstep / ISTEPS_PER_IDEAL 
        endif
        q0 = POW( idpow )%power
        ione = 1
        !Source= q0*dmcheat_noptdvector(Xpoint,laserorigin,ione)
        dqmontedy = - Source * V_k__now
      end function dqmontedy
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !Gradient Evaluation Z_0 (monte carlo source)
      function dqmontedz (Xpoint,Idstep,Deltat,Soln,Nsize)
        use global_params , only: PI
        use pennes_model , only:  X_0,Y_0,Z_0,MU_A,MU_S,MU_TR,MU_EFF, &
                                  POW,ANFACT,ID_W_0,ID_K_0,ISTEPS_PER_IDEAL
        implicit none
        PetscInt   ,intent(in) :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Soln(0:Nsize-1),Deltat
        PetscScalar :: Grad_u_kmhalf(3),Grad_v_k__now(3),&
                       Usln,Usln_prev,V_k__now
        PetscScalar :: q0,laserorigin(0:2),Source
        PetscScalar :: dqmontedz
        PetscInt :: idpow,itwo
        PetscScalar , external :: dmcheat_noptdvector 
        ! map solution
        Usln          = Soln(0)
        Usln_prev     = Soln(1)
        Grad_u_kmhalf = Soln(2:4)
        V_k__now      = Soln(5)
        Grad_v_k__now = Soln(6:8)

        laserorigin(0)=X_0;
        laserorigin(1)=Y_0;
        laserorigin(2)=Z_0;
        ! NOTE the index "idpow" is for a FORTRAN data structure
        ! INDEXING BEGINS WITH ONE!!!!!!!!!!!!!!!!!
        if( modulo(Idstep,ISTEPS_PER_IDEAL) .ne. 0) then
           idpow= Idstep / ISTEPS_PER_IDEAL  + 1
        else 
           idpow= Idstep / ISTEPS_PER_IDEAL 
        endif
        q0 = POW( idpow )%power
        itwo = 2 
        !Source= q0*dmcheat_noptdvector(Xpoint,laserorigin,itwo)
        dqmontedz = - Source * V_k__now
      end function dqmontedz
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !Gradient Evaluation POW (monte carlo source)
      function dqmontedpow (Xpoint,Idstep,Deltat,Soln,Nsize)
        use global_params , only: PI
        use pennes_model , only:  X_0,Y_0,Z_0,MU_A,MU_S,MU_TR,MU_EFF, &
                                  POW,ANFACT,ID_W_0,ID_K_0
        implicit none
        PetscInt   ,intent(in) :: Idstep,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Soln(0:Nsize-1),Deltat
        PetscScalar :: Grad_u_kmhalf(3),Grad_v_k__now(3),&
                       Usln,Usln_prev,V_k__now
        PetscScalar :: laserorigin(0:2) , dqdP,Source
        PetscScalar :: dqmontedpow
        PetscScalar , external :: mcheat_nopt 
        ! map solution
        Usln          = Soln(0)
        Usln_prev     = Soln(1)
        Grad_u_kmhalf = Soln(2:4)
        V_k__now      = Soln(5)
        Grad_v_k__now = Soln(6:8)

        laserorigin(0)=X_0;
        laserorigin(1)=Y_0;
        laserorigin(2)=Z_0;
        dqdP=1.0d0 
        !Source= dqdP * mcheat_nopt(Xpoint,laserorigin)
        dqmontedpow = - Source * V_k__now
      end function dqmontedpow
!----------------------------------------------------------------------
!
!   subroutine  - write_power/read_power     (latest revision: Jun 06)
!
!   purpose     -  solve information is communicated between processor
!                  groups by writing the power field out to a file.
!                  the next group to read it in must read in the same format
!                  that was written out. read_power reads in the same format
!                  that write_power writes
!      
!----------------------------------------------------------------------
      subroutine write_power(FileID)  
        use pennes_model , only: POW
        implicit none
        PetscInt, intent(in) :: FileID
        character(len=MAXLEN):: powerfile 
        character(len=16)  :: charfileID
        PetscInt :: i
        !nplot must be different than in write_visualase_file so that 
        !    two threads do try to open a file with the same unit number.
        integer :: nplot=27  !DO NOT MAKE THIS INTO A PARAMETER WILL CRASH CODE
        integer :: ierr

        write(charfileID,*) FileID
        powerfile='files/power'//trim(adjustl(charfileID))//'.dat'
        ! Try to open the file
        open(unit=nplot, file=trim(powerfile),                 &
             iostat=ierr,form='formatted',access='sequential', &
             status='unknown',action='write')
        ! Error checking
        if(ierr /= 0) then
            write(*,*) "File I/O error! iostat= ",ierr
            write(*,*) 'write_power: Unable to open '//trim(powerfile)
            call abort()
        end if
        ! write out power data structure
        do i=0,size(POW)-1
          write(nplot,'(e12.6,x,e12.6)') POW(i)%time,POW(i)%power
        end do
        close(unit=nplot)
      end subroutine write_power
! 
      subroutine read_power(FileID,Nzero)
        use pennes_model , only: POW
        implicit none
#include "cinout.blk"
        PetscInt, intent(in) :: FileID, &
                                Nzero ! Nzero .ne. 0 is used w/ data_server 
                                      !  to prevent power history overwrite
        character(len=MAXLEN):: powerfile 
        character(len=16)  :: charfileID
        PetscInt :: i
        !nplot must be different than in write_visualase_file so that 
        !    two threads do try to open a file with the same unit number.
        integer :: nplot=27  !DO NOT MAKE THIS INTO A PARAMETER WILL CRASH CODE
        integer :: ierr
        PetscScalar :: dum ! dummy variable

        write(charfileID,*) FileID
        powerfile='files/power'//trim(adjustl(charfileID))//'.dat'
        ! Try to open the file
        open(unit=nplot, file=trim(powerfile),                      &
             iostat=ierr,form='formatted', access='sequential',     &
             status='old',action='read')
        ! Error checking
        if(ierr /= 0) then
            write(*,*) 'File I/O error! iostat=',ierr
            write(*,*) ' read_power: Unable to open '//trim(powerfile)
            call abort()
        end if
        ! read in power data structure
        ! Nzero .ne. 0 is used w/ data_server to prevent power history overwrite
        do i=0,Nzero-1 !must read past unneeded data
          read(nplot,'(e12.6,x,e12.6)') dum,dum
        end do
        do i=Nzero,size(POW)-1
          read(nplot,'(e12.6,x,e12.6)') POW(i)%time,POW(i)%power
        end do
        close(unit=nplot)
        !echo parameters to a file
        do i=0,size(POW)-1
          write(NOUT,'(a,a,I4,a,e12.6,x,e12.6)') trim(powerfile),      &
                     ': POWER(',i,' ,:)=',POW(i)%time,POW(i)%power
        end do
      end subroutine read_power
!----------------------------------------------------------------------
      ! write a file for vis of the power v.s. time profile
      subroutine vis_power(Idfile)
        use plot_params    
        use global_params , only: GroupID,PROFILEID
        use pennes_model   , only: POW
        implicit none
#include "finclude/petsc.h"
        PetscInt,intent(in) :: Idfile
        real,allocatable :: powbuffer(:) 
        character(len=MAXLEN):: avspowerfile,avsheader
        character(len=16)    :: plotsuffix,powlen,grpid
        integer  :: amode,state(MPI_STATUS_SIZE),err,j,ndumplot
        !nplot must be different than in read_power and write_power so that 
        !    two threads do try to open a file with the same unit number.
        integer  :: nplot=69  !DO NOT MAKE THIS INTO A PARAMETER WILL CRASH CODE
        write(grpid,*) GroupID
        write(plotsuffix,*) Idfile ! id of plot
        ndumplot  =NFLDPLOT-len(trim(adjustl(plotsuffix))) 
        plotsuffix=CZEROS(1:ndumplot  )//trim(adjustl(plotsuffix))
        avspowerfile= trim(COMPWRITEMRI)//"/"//trim(PROFILEID)//   &
                      "power_qoi_"//trim(adjustl(grpid))//BYTEID//  &
                      trim(plotsuffix)//".fld"

        !append the plot file
        amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY) 
        call MPI_File_open(PETSC_COMM_SELF,trim(avspowerfile),amode,MPI_INFO_NULL,nplot,err)
        ! copy to buffer 
        allocate(powbuffer(0:SIZE(POW)-1),stat=err)
        do j =  0,SIZE(POW)-1
           powbuffer(j)=POW(j)%power  ! implicit type conversion
        enddo
        !if(BYTESWAP .eqv. PETSC_TRUE) call swapbyteorderf90(powbuffer,SIZE(POW))
        write(powlen,*) SIZE(POW) ! array bounds info
        ! char(10) = '\n'
        avsheader="# AVS field file "                        //char(10)    &
               // "ndim=1   # no. of dimensions in the field"//char(10)    &
               // "nspace=1 # no. of phys. coords per point" //char(10)    &
               // "dim1="//trim(adjustl(powlen))//" #pow dim"//char(10)    &
               // "veclen=1    # no of comps at each point"  //char(10)    &
               // "label= Power # data labels"               //char(10)    &
               // "unit = Watts # data units "               //char(10)    &
               // "data=float   # native format binary"      //char(10)    &
               // "field=uniform # field type"               //char(10)    &
               //char(12)//char(12) ! char(12) = '\f'
        call MPI_File_write(nplot,avsheader,len(trim(avsheader)), &
                                                MPI_BYTE,state,err)
        call MPI_File_write(nplot,powbuffer,SIZE(POW),MPI_REAL,   &
                                                             state,err)
        deallocate(powbuffer,stat=err)
        !@@@@@ close the plot file
        call MPI_File_close(nplot,err) 
      end subroutine vis_power
!----------------------------------------------------------------------
!   subroutine  -  write_visualase_file(Idpow)  (latest revision: Nov 07)
!
!   purpose     -  write out: 
!                   1) file to control the visualase
!                   2) file to log the power output and time instances written
!                   3) file to vis the power as a function of time in AVS
!
!   arguments   -  Idpow: time instance
!
!----------------------------------------------------------------------
      subroutine write_visualase_file(Idpow,MaxTemp)  
        use global_params , only : VISUALASE_FILE,VISUALASE_MAX_POW,VISUALASE_MAX_TEMP
        use pennes_model  , only : POW
        use plot_params 
        implicit none
#include "finclude/petsc.h"
#include "cinout.blk"
        PetscInt    , intent(in) :: Idpow
        PetscScalar , intent(in) :: MaxTemp
        !nplot must be different than in read_power and write_power so that 
        !    two threads do try to open a file with the same unit number.
        integer  :: nplot=69  !DO NOT MAKE THIS INTO A PARAMETER WILL CRASH CODE
        integer :: ierr
        PetscScalar :: percentage,override_power, original_power
        PetscTruth  :: override
        real :: time

        !save the orignal power as a log and to recover if necessary
        original_power = POW(Idpow)%power

        ! fail safe. shut down when max temp exceeded
        if(MaxTemp > VISUALASE_MAX_TEMP) then 
           percentage = 0.0d0
           POW(Idpow)%power=0.0d0
           write(*,*)           "MaxTemp detected for Idpow =",Idpow
           write(NOUT,*) "vlase: MaxTemp detected for Idpow =",Idpow
        endif

        ! check if user wants to manually override the power
        call dynamic_laser_control(override,override_power);

        ! if all fails manually add/subtract a wattage to the power profile
        !  *note* in case of bad filtering and power being cut off
        !         can override with override_power=ZERO
        !         to recover the original power
        if(override .eqv. PETSC_TRUE) POW(Idpow)%power = original_power + override_power 

        ! check if the power is out of bounds
        if( POW(Idpow)%power > VISUALASE_MAX_POW) then
           percentage = 1.0d0
           POW(Idpow)%power=VISUALASE_MAX_POW
        else if( POW(Idpow)%power <       0.0d0      ) then
           percentage = 0.0d0
           POW(Idpow)%power=0.0d0
        else 
           percentage = POW(Idpow)%power/VISUALASE_MAX_POW
        endif


        percentage = 100.0d0 * percentage

        ! Try to open the file
        open(unit=nplot, file='../'//VISUALASE_FILE,              &
             iostat=ierr,form='formatted',access='sequential',    &
             status='unknown',action='write')
        ! Error checking
        if(ierr /= 0) then
            write(*,*) "File I/O error! iostat= ",ierr
            write(*,*) 'write_visualase: File I/O error!'
            call abort()
        end if

        ! write out power percentage
        write(nplot,'(f5.1)') percentage

        ! close file
        close(unit=nplot)
        ! keep a log of the data written to the visualase 
        ! what times the data written out, max temperature
        ! original power, and the override power
        call cpu_time(time)
        write(NOUT,'(a,f5.1,a,I8,a,f5.1,a,e12.5,a,I4,a,e12.5,a,e12.5)') &
        "vlase: t=",time , " Idpow=",Idpow,                             &
        " percentage=",percentage," MaxTemp=",MaxTemp,                  &
        " override=",override,                                          &
        " original_power=",original_power,                              &
        " override_power=",override_power

        ! write power vis file
        call vis_power(Idpow)

      end subroutine write_visualase_file
!!!----------------------------------------------------------------------
!!!   subroutine  -  alloc_hpelemfield/pack_hpelemfield/gatherdealloc_hpelemfield
!!!                                              (latest revision: Aug 06)
!!!
!!!   purpose     -  solve information is communicated between processor
!!!                  groups by writing the w_0 and k_0 fields out to a file.
!!!                  gather_hpelemfield gathers all processor field info onto
!!!                  rank 0 to write to a file 
!!!                  scatter_hpelemfield scatters field info read in from a file
!!!                  to each processor
!!!
!!!   arguments   -  Nsize: Array size
!!!                  IndivprocHPelm: array who's ith entry contains the 
!!!                                  # of elements on rank i
!!!                  HPelm_prev:     array who's ith entry contains the 
!!!                                  sum # of elements on rank 0 thru i-1
!!!
!!!----------------------------------------------------------------------
!!      subroutine alloc_hpelemfield(NelemHP)
!!        use pennes_model , only: HPELMFLD,NUMHPELEMVAR
!!        use data_structure
!!        implicit none
!!        PetscInt , intent(in) :: NelemHP
!!        PetscInt :: astat
!!
!!        allocate(HPELMFLD(NUMHPELEMVAR,0:NelemHP-1) , stat = astat)
!!        if(astat .ne. 0 ) then
!!           write(*,*) "error allocating field buffer"
!!           call abort
!!        endif
!!      end subroutine alloc_hpelemfield
!!!
!!      subroutine pack_hpelemfield(Idbuff,Idelem,IdMdle)
!!        use data_structure
!!        use pennes_model , only: HPELMFLD, K_0 , W_0
!!        implicit none
!!#include "finclude/petsc.h"
!!        PetscInt, intent(in) :: Idbuff,Idelem,IdMdle
!!        HPELMFLD(1,Idbuff) = W_0(ELEMB(Idelem)%map_w_0)
!!        HPELMFLD(2,Idbuff) = K_0(ELEMB(Idelem)%map_k_0)
!!        HPELMFLD(3,Idbuff) = ELEMB(Idelem)%error
!!        HPELMFLD(4,Idbuff) = IdMdle  ! pack the rank for visualization
!!      end subroutine pack_hpelemfield
!!
!!      subroutine gatherdealloc_hpelemfield(ControlComm,
!!     .                               Optimize_W_0,Optimize_K_0,
!!     .                               Nsize,IndivprocHPelm,HPelem_prev)
!!        use data_structure
!!        use pennes_model , only: HPELMFLD,NUMHPELEMVAR, K_0 , W_0
!!        implicit none
!!#include "finclude/petsc.h"
!!        PetscMPIInt, intent(in) :: ControlComm
!!        PetscInt, intent(in) :: Nsize !# of tasks in the computational group
!!        ! field optimizations
!!        PetscTruth, intent(in) :: Optimize_W_0,Optimize_K_0
!!        ! communication is between computational group AND control task
!!        !  NOTE that size Nsize is the size of the computational group
!!        !    the # of tasks in ControlComm = Nsize + 1 
!!        PetscInt, intent(in) :: IndivprocHPelm(0:Nsize), 
!!     .                          HPelem_prev(0:Nsize) 
!!        PetscInt :: mdle ,iii , jjj, idelem , icontrolrank
!!        PetscScalar :: dum
!!        PetscErrorCode :: ierr
!!        
!!        call MPI_Comm_rank(ControlComm,icontrolrank,ierr);
!!
!!        !gather the fields from all processors
!!        if(Control_Task) then
!!          ! NRELEB = 0 on rank 0
!!          call MPI_Gatherv(dum,
!!     .                     IndivprocHPelm(icontrolrank)*NUMHPELEMVAR,
!!     .                     MPIU_SCALAR,HPELMFLD,
!!     .                     IndivprocHPelm*NUMHPELEMVAR,
!!     .                     HPelem_prev*NUMHPELEMVAR,
!!     .                     MPIU_SCALAR,0,ControlComm,ierr)
!!          !to prevent overwriting only store if the group optimized the field
!!          if(Optimize_W_0)then
!!             mdle=0
!!             do iii = 1, NRELEB
!!                call nelconb(mdle,mdle)
!!                idelem = -NODEB(mdle)%father
!!                W_0(ELEMB(idelem)%map_w_0) =  HPELMFLD(1,iii-1)  !C -indexing!!!
!!             enddo
!!          endif
!!          !to prevent overwriting only store if the group optimized the field
!!          if(Optimize_K_0)then
!!             mdle=0
!!             do iii = 1, NRELEB
!!                call nelconb(mdle,mdle)
!!                idelem = -NODEB(mdle)%father
!!                K_0(ELEMB(idelem)%map_k_0) =  HPELMFLD(2,iii-1)  !C -indexing!!!
!!             enddo
!!          endif
!!          mdle=0
!!          do iii = 1, NRELEB
!!             call nelconb(mdle,mdle)
!!             idelem = -NODEB(mdle)%father
!!             ELEMB(idelem)%error        =     HPELMFLD(3,iii-1)  !C -indexing!!!
!!          enddo
!!          !retrieve mdle #'s owned by group to visualize domain decomposition
!!          do iii = 1, Nsize
!!            do jjj = 0 , IndivprocHPelm(iii)-1 !C -indexing!!!
!!              NODEB(int(HPELMFLD(4,HPelem_prev(iii)+jjj)))%irankown=iii
!!            enddo
!!          enddo
!!        else
!!          call MPI_Gatherv(HPELMFLD,
!!     .                     IndivprocHPelm(icontrolrank)*NUMHPELEMVAR,
!!     .                     MPIU_SCALAR,dum,
!!     .                     IndivprocHPelm*NUMHPELEMVAR,
!!     .                     HPelem_prev*NUMHPELEMVAR,
!!     .                     MPIU_SCALAR,0,ControlComm,ierr)
!!        endif
!!        ! deallocate temp storage
!!        deallocate(HPELMFLD, stat = ierr)
!!        if(ierr .ne. 0 ) then
!!           write(*,*) "error allocating field buffer"
!!           call abort
!!        endif
!!
!!      end subroutine gatherdealloc_hpelemfield
!----------------------------------------------------------------------
!
!   subroutine  - write_w_0field/read_w_0field  (latest revision: Jun 06)
!
!   purpose     -  solve information is communicated between processor
!                  groups by writing the w_0 field out to a file.
!                  the next group to read it in must read in the same format
!                  that was written out. read_w_0field reads in the same
!                  format that write_w_0field writes
!   arguments   -  FileID:  file name identification
!      
!----------------------------------------------------------------------
      subroutine write_w_0field(FileID)
        use pennes_model , only: W_0
        implicit none
        PetscInt, intent(in) :: FileID
        character(len=MAXLEN):: w_0fieldfile 
        character(len=16)  :: charfileID
        PetscInt:: i
        integer :: istat

        write(charfileID,*) FileID
        w_0fieldfile='files/w_0field'//trim(adjustl(charfileID))//'.dat'

        ! Try to open the file
        open(unit=27, file=trim(w_0fieldfile), &
             iostat=istat,form='formatted',access='sequential', &
                                   status='unknown',action='write')
        ! Error checking
        if(istat /= 0) then
            write(*,*) 'File I/O error!'
            write(*,*) 'write_w_0: Unable to open',trim(w_0fieldfile)
            write(*,*) ''
            call abort
        end if
        ! write out element field
        do i=0,size(W_0)-1 
          write(27,'(e12.6,x,e12.6)') W_0(i)
        end do
        close(unit=27)
      end subroutine write_w_0field
! 
      subroutine read_w_0field(FileID)
        use pennes_model , only: W_0
        implicit none
#include "cinout.blk"
        PetscInt, intent(in) :: FileID
        character(len=MAXLEN):: w_0fieldfile 
        character(len=16)    :: charfileID
        PetscInt:: i
        integer :: istat

        write(charfileID,*) FileID
        w_0fieldfile='files/w_0field'//trim(adjustl(charfileID))//'.dat'
        ! Try to open the file
        open(unit=27, file=trim(w_0fieldfile), &
             iostat=istat,form='formatted', access='sequential', &
                                        status='old',action='read')
        ! Error checking
        if(istat /= 0) then
            write(*,*) 'File I/O error!'
            write(*,*) 'read_w_0: Unable to open',trim(w_0fieldfile)
            write(*,*) ''
            call abort
        end if
        ! read in element field
        do i=0,size(W_0)-1 
          read(27,'(e12.6,x,e12.6)') W_0(i)
        end do
        close(unit=27)
        ! echo parameters read in
        write(NOUT,*) 'read_w_0field(',FileID,'):W_0=',W_0(0:size(W_0)-1)
      end subroutine read_w_0field
!----------------------------------------------------------------------
!
!   subroutine  - write_k_0field/read_k_0field  (latest revision: Jun 06)
!
!   purpose     -  solve information is communicated between processor
!                  groups by writing the k_0 field out to a file.
!                  the next group to read it in must read in the same format
!                  that was written out. read_k_0field reads in the same
!                  format that write_k_0field writes
!   arguments   -  GroupID: computational group Identification
!                  FileID:  file name identification
!      
!----------------------------------------------------------------------
      subroutine write_k_0field(FileID)
        use pennes_model , only: K_0
        implicit none
        PetscInt, intent(in) :: FileID
        character(len=MAXLEN):: k_0fieldfile 
        character(len=16)  ::   charfileID
        PetscInt:: i
        integer :: istat

        write(charfileID,*) FileID
        k_0fieldfile='files/k_0field'//trim(adjustl(charfileID))//'.dat'
        ! Try to open the file
        open(unit=27, file=trim(k_0fieldfile), &
             iostat=istat,form='formatted',access='sequential', &
             status='unknown',action='write')
        ! Error checking
        if(istat /= 0) then
            write(*,*) 'File I/O error!'
            write(*,*) 'write_k_0: Unable to open',trim(k_0fieldfile)
            write(*,*) ''
            call abort
        end if
        ! write out element field
        do i=0,size(K_0)-1 
          write(27,'(e12.6,x,e12.6)') K_0(i)
        end do
        close(unit=27)
      end subroutine write_k_0field
! 
      subroutine read_k_0field(FileID)
        use pennes_model , only: K_0
        implicit none
#include "cinout.blk"
        PetscInt, intent(in) :: FileID
        character(len=MAXLEN):: k_0fieldfile 
        character(len=16)  :: charfileID
        PetscInt:: i
        integer :: istat

        write(charfileID,*) FileID
        k_0fieldfile='files/k_0field'//trim(adjustl(charfileID))//'.dat'
        ! Try to open the file
        open(unit=27, file=trim(k_0fieldfile), &
             iostat=istat,form='formatted', access='sequential', &
                                        status='old',action='read')
        ! Error checking
        if(istat /= 0) then
            write(*,*) 'File I/O error!'
            write(*,*) 'read_k_0: Unable to open',trim(k_0fieldfile)
            write(*,*) ''
            call abort
        end if
        ! read in element field
        do i=0,size(K_0)-1 
          read(27,'(e12.6,x,e12.6)') K_0(i)
        end do
        close(unit=27)
        write(NOUT,*) 'read_k_0field(',FileID,'):K_0=',K_0(0:size(K_0)-1)
      end subroutine read_k_0field
!----------------------------------------------------------------------
!
!    the following are various small routines that pass consitutive data
!      from hp3d to C++ for the parameter optimization of W_0
!      
!----------------------------------------------------------------------
      !get a parameter buffer from the w_0 array
      subroutine getparam_w_0(Parambuff,Id)
        use pennes_model , ONLY: W_0
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = w_0(Id) 
      end subroutine getparam_w_0
      !get w_0_lb
      subroutine getparam_w_0_lb(Parambuff,Id)
        use pennes_model , ONLY: W_0_LB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = W_0_LB
      end subroutine getparam_w_0_lb
      !get w_0_ub
      subroutine getparam_w_0_ub(Parambuff,Id)
        use pennes_model , ONLY: W_0_UB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = W_0_UB
      end subroutine getparam_w_0_ub
      !put a parameter buffer into the w_0 array
      subroutine putparam_w_0(Parambuff,Id)
        use pennes_model , ONLY: W_0
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(in)  :: Parambuff
        w_0(Id) = Parambuff
      end subroutine putparam_w_0
      !get field variation state for w_0 
      function get_w_0_field()
        use pennes_model , ONLY: W_0_FIELD
        implicit none
        PetscTruth :: get_w_0_field
        get_w_0_field = W_0_FIELD
      end function get_w_0_field 
!----------------------------------------------------------------------
!
!    the following are various small routines that pass consitutive data
!      from hp3d to C++ for the parameter optimization of K_0
!      
!----------------------------------------------------------------------
      !get a parameter buffer from the k_0 array
      subroutine getparam_k_0(Parambuff,Id)
        use pennes_model , ONLY: K_0
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = k_0(Id) 
      end subroutine getparam_k_0
      !get k_0_lb 
      subroutine getparam_k_0_lb(Parambuff,Id)
        use pennes_model , ONLY: K_0_LB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = K_0_LB
      end subroutine getparam_k_0_lb
      !get k_0_ub
      subroutine getparam_k_0_ub(Parambuff,Id)
        use pennes_model , ONLY: K_0_UB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = K_0_UB
      end subroutine getparam_k_0_ub
      !put a parameter buffer into the k_0 array
      subroutine putparam_k_0(Parambuff,Id)
        use pennes_model , ONLY: K_0
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(in)  :: Parambuff
        k_0(Id) = Parambuff
      end subroutine putparam_k_0
      !get field variation state for k_0 
      function get_k_0_field()
        use pennes_model , ONLY: K_0_FIELD
        implicit none
        PetscTruth :: get_k_0_field
        get_k_0_field = K_0_FIELD
      end function get_k_0_field 
!----------------------------------------------------------------------
!
!    the following are various small routines that pass consitutive data
!      from hp3d to C++ for the parameter optimization of POW
!      
!----------------------------------------------------------------------
      !map for the POW array is built from the size
      function get_pow_size()
        use pennes_model , only: IDEAL_NZERO,IDEAL_NTIME
        implicit none
        PetscInt :: get_pow_size
        get_pow_size = IDEAL_NTIME - IDEAL_NZERO 
      end function get_pow_size
      !get a parameter buffer from the pow array
      subroutine getparam_pow(Parambuff,Id)
        use pennes_model , ONLY: POW,IDEAL_NZERO
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        PetscInt :: jlo
        jlo=IDEAL_NZERO+1
        Parambuff = POW( jlo + Id)%power 
      end subroutine getparam_pow
      !get POW_LB
      subroutine getparam_pow_lb(Parambuff,Id)
        use pennes_model , ONLY: POW_LB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = POW_LB
      end subroutine getparam_pow_lb
      !get POW_UB
      subroutine getparam_pow_ub(Parambuff,Id)
        use pennes_model , ONLY: POW_UB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = POW_UB
      end subroutine getparam_pow_ub
      !put a parameter buffer into the pow array
      subroutine putparam_pow(Parambuff,Id)
        use pennes_model , ONLY: POW , IDEAL_NZERO
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(in)  :: Parambuff
        PetscInt :: jlo
        jlo=IDEAL_NZERO+1
        POW( jlo + Id)%power = Parambuff
      end subroutine putparam_pow
!----------------------------------------------------------------------
!
!    the following are various small routines that pass consitutive data
!      from hp3d to C++ for the parameter optimization of K_1
!      
!----------------------------------------------------------------------
      !get a parameter buffer from the k_1
      subroutine getparam_k_1(Parambuff,Id)
        use pennes_model , ONLY: K_1
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = k_1
      end subroutine getparam_k_1
      !get k_1_LB
      subroutine getparam_k_1_lb(Parambuff,Id)
        use pennes_model , ONLY: K_1_LB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = K_1_LB
      end subroutine getparam_k_1_lb
      !get k_1_UB
      subroutine getparam_k_1_ub(Parambuff,Id)
        use pennes_model , ONLY: K_1_UB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = K_1_UB
      end subroutine getparam_k_1_ub
      !put a parameter buffer into the k_1 array
      subroutine putparam_k_1(Parambuff,Id)
        use pennes_model , ONLY: K_1
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(in)  :: Parambuff
        k_1 = Parambuff
      end subroutine putparam_k_1
!----------------------------------------------------------------------
!
!    the following are various small routines that pass consitutive data
!      from hp3d to C++ for the parameter optimization of K_2
!      
!----------------------------------------------------------------------
      !get a parameter buffer from the k_2
      subroutine getparam_k_2(Parambuff,Id)
        use pennes_model , ONLY: K_2
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = k_2
      end subroutine getparam_k_2
      !get k_2_LB
      subroutine getparam_k_2_lb(Parambuff,Id)
        use pennes_model , ONLY: K_2_LB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = K_2_LB
      end subroutine getparam_k_2_lb
      !get k_2_UB
      subroutine getparam_k_2_ub(Parambuff,Id)
        use pennes_model , ONLY: K_2_UB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = K_2_UB
      end subroutine getparam_k_2_ub
      !put a parameter buffer into the k_2 array
      subroutine putparam_k_2(Parambuff,Id)
        use pennes_model , ONLY: K_2
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(in)  :: Parambuff
        k_2 = Parambuff
      end subroutine putparam_k_2
!----------------------------------------------------------------------
!
!    the following are various small routines that pass consitutive data
!      from hp3d to C++ for the parameter optimization of K_3
!      
!----------------------------------------------------------------------
      !get a parameter buffer from the k_3
      subroutine getparam_k_3(Parambuff,Id)
        use pennes_model , ONLY: K_3
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = k_3
      end subroutine getparam_k_3
      !get k_3_LB
      subroutine getparam_k_3_lb(Parambuff,Id)
        use pennes_model , ONLY: K_3_LB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = k_3_LB
      end subroutine getparam_k_3_lb
      !get a parameter buffer from the k_3_UB
      subroutine getparam_k_3_ub(Parambuff,Id)
        use pennes_model , ONLY: K_3_UB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = k_3_UB
      end subroutine getparam_k_3_ub
      !put a parameter buffer into the k_3 array
      subroutine putparam_k_3(Parambuff,Id)
        use pennes_model , ONLY: K_3
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(in)  :: Parambuff
        k_3 = Parambuff
      end subroutine putparam_k_3
!----------------------------------------------------------------------
!
!    the following are various small routines that pass consitutive data
!      from hp3d to C++ for the parameter optimization of W_N
!      
!----------------------------------------------------------------------
      !get a parameter buffer from the w_n
      subroutine getparam_w_n(Parambuff,Id)
        use pennes_model , ONLY: W_N
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = w_n
      end subroutine getparam_w_n
      !get w_n_lb
      subroutine getparam_w_n_lb(Parambuff,Id)
        use pennes_model , ONLY: W_N_LB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = W_N_LB
      end subroutine getparam_w_n_lb
      !get w_n_ub
      subroutine getparam_w_n_ub(Parambuff,Id)
        use pennes_model , ONLY: W_N_UB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = W_N_UB
      end subroutine getparam_w_n_ub
      !put a parameter buffer into the w_n array
      subroutine putparam_w_n(Parambuff,Id)
        use pennes_model , ONLY: W_N
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(in)  :: Parambuff
        w_n = Parambuff
      end subroutine putparam_w_n
!----------------------------------------------------------------------
!
!    the following are various small routines that pass consitutive data
!      from hp3d to C++ for the parameter optimization of W_I
!      
!----------------------------------------------------------------------
      !get a parameter buffer from the w_i
      subroutine getparam_w_i(Parambuff,Id)
        use pennes_model , ONLY: W_I
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = w_i
      end subroutine getparam_w_i
      !get a parameter buffer from the w_i
      subroutine getparam_w_i_lb(Parambuff,Id)
        use pennes_model , ONLY: W_I_LB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = W_I_LB
      end subroutine getparam_w_i_lb
      !get a parameter buffer from the w_i
      subroutine getparam_w_i_ub(Parambuff,Id)
        use pennes_model , ONLY: W_I_UB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = W_I_UB
      end subroutine getparam_w_i_ub
      !put a parameter buffer into the w_i array
      subroutine putparam_w_i(Parambuff,Id)
        use pennes_model , ONLY: W_I
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(in)  :: Parambuff
        w_i = Parambuff
      end subroutine putparam_w_i
!----------------------------------------------------------------------
!
!    the following are various small routines that pass consitutive data
!      from hp3d to C++ for the parameter optimization of W_D
!      
!----------------------------------------------------------------------
      !get a parameter buffer from the w_d
      subroutine getparam_w_d(Parambuff,Id)
        use pennes_model , ONLY: W_D
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = w_d
      end subroutine getparam_w_d
      !get a w_d_lb
      subroutine getparam_w_d_lb(Parambuff,Id)
        use pennes_model , ONLY: W_D_LB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = W_D_LB
      end subroutine getparam_w_d_lb
      !get a w_d_ub
      subroutine getparam_w_d_ub(Parambuff,Id)
        use pennes_model , ONLY: W_D_UB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = W_D_UB
      end subroutine getparam_w_d_ub
      !put a parameter buffer into the w_d array
      subroutine putparam_w_d(Parambuff,Id)
        use pennes_model , ONLY: W_D
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(in)  :: Parambuff
        w_d = Parambuff
      end subroutine putparam_w_d
!----------------------------------------------------------------------
!
!    the following are various small routines that pass consitutive data
!      from hp3d to C++ for the parameter optimization of W_2
!      
!----------------------------------------------------------------------
      !get a parameter buffer from the w_2
      subroutine getparam_w_2(Parambuff,Id)
        use pennes_model , ONLY: W_2
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = w_2
      end subroutine getparam_w_2
      !get w_2_lb
      subroutine getparam_w_2_lb(Parambuff,Id)
        use pennes_model , ONLY: W_2_LB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = W_2_LB
      end subroutine getparam_w_2_lb
      !get w_2_ub
      subroutine getparam_w_2_ub(Parambuff,Id)
        use pennes_model , ONLY: W_2_UB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = W_2_UB
      end subroutine getparam_w_2_ub
      !put a parameter buffer into the w_2 array
      subroutine putparam_w_2(Parambuff,Id)
        use pennes_model , ONLY: W_2
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(in)  :: Parambuff
        w_2 = Parambuff
      end subroutine putparam_w_2
!----------------------------------------------------------------------
!
!    the following are various small routines that pass consitutive data
!      from hp3d to C++ for the parameter optimization of W_NI
!      
!----------------------------------------------------------------------
      !get a parameter buffer from the w_ni
      subroutine getparam_w_ni(Parambuff,Id)
        use pennes_model , ONLY: W_NI
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = w_ni
      end subroutine getparam_w_ni
      !get W_NID_LB
      subroutine getparam_w_ni_lb(Parambuff,Id)
        use pennes_model , ONLY: W_NID_LB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = W_NID_LB
      end subroutine getparam_w_ni_lb
      !get W_NID_UB
      subroutine getparam_w_ni_ub(Parambuff,Id)
        use pennes_model , ONLY: W_NID_MD
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = W_NID_MD
      end subroutine getparam_w_ni_ub
      !put a parameter buffer into the w_ni array
      subroutine putparam_w_ni(Parambuff,Id)
        use pennes_model , ONLY: W_NI
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(in)  :: Parambuff
        w_ni = Parambuff
      end subroutine putparam_w_ni
!----------------------------------------------------------------------
!
!    the following are various small routines that pass consitutive data
!      from hp3d to C++ for the parameter optimization of W_ID
!      
!----------------------------------------------------------------------
      !get a parameter buffer from the w_id
      subroutine getparam_w_id(Parambuff,Id)
        use pennes_model , ONLY: W_ID
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = w_id
      end subroutine getparam_w_id
      !get W_NID_MD
      subroutine getparam_w_id_lb(Parambuff,Id)
        use pennes_model , ONLY: W_NID_MD
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = W_NID_MD
      end subroutine getparam_w_id_lb
      !get W_NID_UB
      subroutine getparam_w_id_ub(Parambuff,Id)
        use pennes_model , ONLY: W_NID_UB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = W_NID_UB
      end subroutine getparam_w_id_ub
      !put a parameter buffer into the w_id array
      subroutine putparam_w_id(Parambuff,Id)
        use pennes_model , ONLY: W_ID
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(in)  :: Parambuff
        w_id = Parambuff
      end subroutine putparam_w_id
!----------------------------------------------------------------------
!
!    the following are various small routines that pass consitutive data
!      from hp3d to C++ for the parameter optimization of X_0
!      
!----------------------------------------------------------------------
      !get a parameter buffer from the x_0
      subroutine getparam_x_0(Parambuff,Id)
        use pennes_model , ONLY: X_0
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = x_0
      end subroutine getparam_x_0
      !get x_0_LB
      subroutine getparam_x_0_lb(Parambuff,Id)
        use pennes_model , ONLY: X_0_LB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = X_0_LB
      end subroutine getparam_x_0_lb
      !get x_0_UB
      subroutine getparam_x_0_ub(Parambuff,Id)
        use pennes_model , ONLY: X_0_UB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = X_0_UB
      end subroutine getparam_x_0_ub
      !put a parameter buffer into the x_0 array
      subroutine putparam_x_0(Parambuff,Id)
        use pennes_model , ONLY: X_0
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(in)  :: Parambuff
        x_0 = Parambuff
      end subroutine putparam_x_0
!----------------------------------------------------------------------
!
!    the following are various small routines that pass consitutive data
!      from hp3d to C++ for the parameter optimization of Y_0
!      
!----------------------------------------------------------------------
      !get a parameter buffer from the y_0
      subroutine getparam_y_0(Parambuff,Id)
        use pennes_model , ONLY: Y_0
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = y_0
      end subroutine getparam_y_0
      !get y_0_LB
      subroutine getparam_y_0_lb(Parambuff,Id)
        use pennes_model , ONLY: Y_0_LB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = Y_0_LB
      end subroutine getparam_y_0_lb
      !get y_0_UB
      subroutine getparam_y_0_ub(Parambuff,Id)
        use pennes_model , ONLY: Y_0_UB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = Y_0_UB
      end subroutine getparam_y_0_ub
      !put a parameter buffer into the y_0 array
      subroutine putparam_y_0(Parambuff,Id)
        use pennes_model , ONLY: Y_0
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(in)  :: Parambuff
        y_0 = Parambuff
      end subroutine putparam_y_0
!----------------------------------------------------------------------
!
!    the following are various small routines that pass consitutive data
!      from hp3d to C++ for the parameter optimization of Z_0
!      
!----------------------------------------------------------------------
      !get a parameter buffer from the z_0
      subroutine getparam_z_0(Parambuff,Id)
        use pennes_model , ONLY: Z_0
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = z_0
      end subroutine getparam_z_0
      !get z_0_LB
      subroutine getparam_z_0_lb(Parambuff,Id)
        use pennes_model , ONLY: Z_0_LB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = z_0_LB
      end subroutine getparam_z_0_lb
      !get z_0_UB
      subroutine getparam_z_0_ub(Parambuff,Id)
        use pennes_model , ONLY: Z_0_UB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = z_0_UB
      end subroutine getparam_z_0_ub
      !put a parameter buffer into the z_0 array
      subroutine putparam_z_0(Parambuff,Id)
        use pennes_model , ONLY: Z_0
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(in)  :: Parambuff
        z_0 = Parambuff
      end subroutine putparam_z_0
!----------------------------------------------------------------------
!
!    the following are various small routines that pass consitutive data
!      from hp3d to C++ for the parameter optimization of MU_A
!      
!----------------------------------------------------------------------
      !get a parameter buffer from the mu_a
      subroutine getparam_mu_a(Parambuff,Id)
        use pennes_model , ONLY: MU_A
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = mu_a
      end subroutine getparam_mu_a
      !get mu_a_LB
      subroutine getparam_mu_a_lb(Parambuff,Id)
        use pennes_model , ONLY: MU_A_LB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = mu_a_LB
      end subroutine getparam_mu_a_lb
      !get mu_a_UB
      subroutine getparam_mu_a_ub(Parambuff,Id)
        use pennes_model , ONLY: MU_A_UB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = mu_a_UB
      end subroutine getparam_mu_a_ub
      !put a parameter buffer into the mu_a array
      subroutine putparam_mu_a(Parambuff,Id)
        use pennes_model , ONLY: MU_A
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(in)  :: Parambuff
        mu_a = Parambuff
      end subroutine putparam_mu_a
!----------------------------------------------------------------------
!
!    the following are various small routines that pass consitutive data
!      from hp3d to C++ for the parameter optimization of MU_S
!      
!----------------------------------------------------------------------
      !get a parameter buffer from the mu_s
      subroutine getparam_mu_s(Parambuff,Id)
        use pennes_model , ONLY: MU_S
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = mu_s
      end subroutine getparam_mu_s
      !get mu_s_LB
      subroutine getparam_mu_s_lb(Parambuff,Id)
        use pennes_model , ONLY: MU_S_LB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = MU_S_LB
      end subroutine getparam_mu_s_lb
      !get mu_s_UB
      subroutine getparam_mu_s_ub(Parambuff,Id)
        use pennes_model , ONLY: MU_S_UB
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(out) :: Parambuff
        Parambuff = MU_S_UB
      end subroutine getparam_mu_s_ub
      !put a parameter buffer into the mu_s array
      subroutine putparam_mu_s(Parambuff,Id)
        use pennes_model , ONLY: MU_S
        implicit none
        PetscInt ,intent(in)  :: Id
        PetscScalar ,intent(in)  :: Parambuff
        mu_s = Parambuff
      end subroutine putparam_mu_s
!----------------------------------------------------------------------
!
!   subroutine  - setup_pennes_model  (latest revision: Aug 08)
!
!   purpose     - setup all pennes model parameters and
!                 setup the field parameters for optimization(if any)
!                 for and INITIAL UNREFINED MESH.
!----------------------------------------------------------------------
      subroutine setup_pennes_model(NfieldElem)
       use parse_ini
       use pennes_model
       implicit none
#include "finclude/petsc.h"
       ! total # of element amongst processors in hp3d
       PetscInt    , intent(in) :: NfieldElem         ! # of field elements
       character(len=MAXLEN) :: inifile   ! INI filename
       logical :: istat
       PetscInt :: irank,i
       PetscScalar :: w_0_const, k_0_const , w_0_tumor, k_0_tumor, locdamage
       PetscScalar :: getlocarrdam
       PetscErrorCode :: ierr

       ! initialize the blood perfusion and thermal conductivity 
       !  from the control file. All processors open files
       inifile = dfltINIFile 
       call openIniFile(inifile, istat) ! read the INI file
       !If unsuccessfully read,stop
       if(.not. istat) then
         write(*,*) 'could not open', inifile
         stop
       end if

       ! have exact solution? 0 = no
       call getIniInt('method', 'NEXACT', NEXACTVERIFNUMBER , 0,istat)

       ! get the constant parameter
       call getIniReal('perfusion','W_0',w_0_const,6.0d0,istat)
       ! get the constant parameter
       call getIniReal('thermal_conductivity','K_0',k_0_const,0.5d0,istat)
       call getIniPetscTruth('field','W_0_FIELD',W_0_FIELD)
       call getIniPetscTruth('field','K_0_FIELD',K_0_FIELD)
       call getIniReal('field','W_0_Tumor',w_0_tumor,w_0_const,istat)
       call getIniReal('field','K_0_Tumor',k_0_tumor,k_0_const,istat)

       ! set size of  thermal conductivities
       if(k_0_field .eqv. PETSC_TRUE ) then ! let K_0 vary spatially
           ! NfieldElem + 1 (for the constant field)
           allocate(K_0(0:NfieldElem),STAT=ierr)
           K_0(:) = k_0_tumor ! set the field portion
       else ! K_0 is spatially constant
           allocate(K_0(0:0),STAT=ierr)
       endif
       if(ierr.ne.0) then
         write(*,*) 'error allocating thermal conductivity'
         call abort
       end if
       !the first entry is always the constant portion
       K_0(0) = k_0_const 

       ! set size of  blood perfusion field
       if(w_0_field .eqv. PETSC_TRUE ) then ! let W_0 vary spatially
           ! NfieldElem + 1 (for the constant field)
           allocate(W_0(0:NfieldElem),STAT=ierr)
           W_0(:) = w_0_tumor ! set the field portion
       else ! W_0 is spatially constant
           allocate(W_0(0:0),STAT=ierr)
       endif
       if(ierr.ne.0) then
         write(*,*) 'error allocating blood perfusion'
         call abort
       end if
       !the first entry is always the constant portion
       W_0(0) = w_0_const 

       ! timestep section
       ! deltat for optimization problems
       call getIniReal('timestep','IDEAL_DT',IDEAL_DT, &
                                             dble(dfltideal_dt),istat)
       call getIniInt('timestep','ISTEPS_PER_IDEAL',ISTEPS_PER_IDEAL, &
                                             dfltistepsperideal,istat)

       ! anisotropy factor
       call getIniReal('optical','ANFACT',ANFACT,0.71d0, istat)

       ! Arrhenius Damage Model Data from Rylander  from
       ! @article{chin2001cdp,
       !   title={{Changes in dielectric properties of ex vivo bovine liver at 915 MHz  during heating}},
       !   author={Chin, L. and Sherar, M.},
       !   journal={Physics in Medicine and Biology},
       !   volume={46},
       !   number={1},
       !   pages={197--212},
       !   year={2001},
       !   publisher={IOP PUBLISHING LTD}
       ! }
       ! default universal gas const = Arr_R given in units of [kcal/mol/K] 
       call getIniReal('arrhenius','Arr_R' ,Arr_R ,1.98d-3,istat)
       ! default A = Arr_A given in units of [1/s]
       call getIniReal('arrhenius','Arr_A' ,Arr_A ,1.95d36,istat)
       ! default Ea = Arr_Ea given in units of [kcal/mol]
       call getIniReal('arrhenius','Arr_Ea',Arr_Ea,60.1d0,istat)

       ! Two State Damage Model Data
       ! default h = TS_h given in units of [K]
       call getIniReal('two_state','TS_h'    ,TS_h    ,7.0031d4,istat)
       ! scale TS_h for now default is too big and giving about zero for damage
       call getIniReal('two_state','TS_h'    ,TS_h    ,7.0031d2,istat)
       ! default alpha = TS_alpha given in units of [1/s]
       call getIniReal('two_state','TS_alpha',TS_alpha,4.93d-2 ,istat)
       ! default beta  = TS_beta given in units of [K]
       call getIniReal('two_state','TS_beta' ,TS_beta ,215.64d0,istat)

       ! initial temperature

       call getIniReal('INITIAL_CONDITION','U_INIT',U_INIT,310.0d0,istat)
       call getIniReal('INITIAL_CONDITION','U_PROBE',U_PROBE,294.0d0,istat)

       ! specific heat of blood

       call getIniReal('perfusion','C_BLOOD',C_BLOOD,0.0d0,istat)

       ! Arterial temperature

       call getIniReal('perfusion','U_ARTERY',U_A,330.0d0,istat)

       !read in tissue density and specific heat

       call getIniReal('material','RHO',RHO,0.0d0,istat)
       call getIniReal('material','SPECIFIC_HEAT',C_P,0.0d0,istat)

       !read in Cauchy boundary data

       call getIniReal('bc','NEWTON_COEFF', COEFF_COOL,100.0d0,istat)
       call getIniReal('bc','U_INFTY',U_INF,315.0d0,istat)

       !read in Neumann boundary data

       call getIniReal('bc','TEMP_FLUX',G_FLUX,1.0d2,istat)

       !read in electric conductivity
       call getIniReal('conductivity','S_0',S_0,0.0d0,istat)
       call getIniReal('conductivity','S_1',S_1,0.0d0,istat)

       ! close ini file
       call closeIniFile()

      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !@@@  echo model parameters
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      call mpi_comm_rank(PETSC_COMM_WORLD,irank,ierr)

      ! to get an idea of the expected damage to accumulate for on time interval
      locdamage =  getlocarrdam(U_INIT,IDEAL_DT)

      if(irank.eq.0) then 
        call printpetscint(   'setup_pennes_model: w_0_field    ='//char(0),w_0_field     ) 
        call printpetscint(   'setup_pennes_model: k_0_field    ='//char(0),k_0_field     )
        call printpetscscalar('setup_pennes_model: S_0          ='//char(0),S_0           )
        call printpetscscalar('setup_pennes_model: S_1          ='//char(0),S_1           )
        call printpetscscalar('setup_pennes_model: w_0_const    ='//char(0),w_0_const     )
        call printpetscscalar('setup_pennes_model: k_0_const    ='//char(0),k_0_const     )
        call printpetscscalar('setup_pennes_model: w_0_tumor    ='//char(0),w_0_tumor     )
        call printpetscscalar('setup_pennes_model: k_0_tumor    ='//char(0),k_0_tumor     )
        call printpetscscalar('setup_pennes_model: Arr_A        ='//char(0),Arr_A         )
        call printpetscscalar('setup_pennes_model: Arr_R        ='//char(0),Arr_R         )
        call printpetscscalar('setup_pennes_model: Arr_Ea       ='//char(0),Arr_Ea        )
        call printpetscscalar('setup_pennes_model: locdamage    ='//char(0),locdamage        )
        call printpetscscalar('setup_pennes_model: TS_h         ='//char(0),TS_h          )
        call printpetscscalar('setup_pennes_model: TS_alpha     ='//char(0),TS_alpha      )
        call printpetscscalar('setup_pennes_model: TS_beta      ='//char(0),TS_beta       )
        call printpetscscalar('setup_pennes_model: ANFACT       ='//char(0),ANFACT        )
        call printpetscscalar('setup_pennes_model: C_BLOOD      ='//char(0),C_BLOOD       )
        call printpetscscalar('setup_pennes_model: U_A          ='//char(0),U_A           )
        call printpetscscalar('setup_pennes_model: RHO          ='//char(0),RHO           )
        call printpetscscalar('setup_pennes_model: C_P          ='//char(0),C_P           )
        call printpetscscalar('setup_pennes_model: COEFF_COOL   ='//char(0),COEFF_COOL    )
        call printpetscscalar('setup_pennes_model: U_INF        ='//char(0),U_INF         )
        call printpetscscalar('setup_pennes_model: G_FLUX       ='//char(0),G_FLUX        )
        call printpetscscalar('setup_pennes_model: U_INIT       ='//char(0),U_INIT        )
        call printpetscscalar('setup_pennes_model: U_PROBE      ='//char(0),U_PROBE       )
      endif


      end subroutine setup_pennes_model
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine setup_power(Maxtime)
       use parse_ini
       use pennes_model
       implicit none
#include "finclude/petsc.h"
       PetscScalar , intent(in) :: Maxtime
       character(len=MAXLEN) :: inifile   ! INI filename
       character(len=MAXLEN) :: compfilelocation ! local location of data
       character(len=MAXLEN) :: powerfile ! text file containing power vs time 
       character(len=MAXLEN) :: buffer    ! scratch storage
       type(POWER) , allocatable , dimension(:) :: powfile
       logical :: istat
       PetscInt :: irank,numlines,numpow,i,j,jlo,jhi
       PetscScalar :: powerdata
       PetscErrorCode :: ierr

        ! initialize the blood perfusion and thermal conductivity 
        !  from the control file. All processors open files
        inifile = dfltINIFile 
        call openIniFile(inifile, istat) ! read the INI file
        !If unsuccessfully read,stop
        if(.not. istat) then
          write(*,*) 'could not open', inifile
          stop
        end if
        ! file location
        call getIniString('compexec','compfilelocation', compfilelocation,'files/',istat)
        call getIniReal('probe', 'power', powerdata, 0.0d0, istat)
        ! close ini file
        call closeIniFile()

      !  ! Try to open the file
      !  powerfile = trim(compfilelocation)//'/power.dat'
      !  open(unit=27, file=powerfile,iostat=ierr, form='formatted',  &
      !                access='sequential',status='old',action='read')
      !  ! Error checking
      !  if(ierr /= 0) then
      !      write(*,*) 'File I/O error!'
      !      write(*,*)' setup_pennes_model(): Unable to open ',    &
      !          trim(powerfile), ': ierr=', ierr
      !      write(*,*) ''
      !      call abort
      !  end if
      !  ! Read through the file to determine how much memory to allocate
      !  numlines = 0
      !  do
      !    read(unit=27, fmt='(A)', iostat=ierr) buffer
      !    if (ierr /= 0) exit   ! exit at EOF
      !    ! Increase the line counter
      !    numlines = numlines + 1
      !  end do
      !  if (numlines == 0) then
      !      write(*,*) 'File I/O error!'
      !      write(*,'(1X,3A)') ' setup_pennes_model: File ',trim(powerfile), ' is empty'
      !      close(unit=27) ; call abort
      !  end if
      !  ! Try to allocate memory to hold file
      !  allocate(powfile(0:numlines),stat=ierr)
      !  ! Error checking
      !  if(ierr /= 0) then
      !    write(*,*) 'Memory error!'
      !    write(*,*) ' setup_pennes_model: Unable to allocate memory'
      !    close(unit=27) ; call abort
      !  end if
      !  ! Now, read the file into memory
      !  rewind(unit=27)  ! rewind the file
      !  ! Read in the data, assume reading in data of the form:
      !  ! 
      !  ! 
      !  !               line 1     line 2 ... line n-1       line n
      !  !   |----P(t1)----|---P(t2)---|----------|-----P(tn)-----| 
      !  !  t=0           t1          t2   ...   tn-1             tn
      !  ! 
      !  ! 
      !  powfile(0)%time=0.0d0 ;  powfile(0)%power=0.0d0
      !  do i=1,numlines
      !    read(27,*) powfile(i)%time, powfile(i)%power
      !  end do
      !  ! All finished
      !  close(unit=27)

      !  call mpi_comm_rank(PETSC_COMM_WORLD,irank,ierr)
      !  if(irank.eq.0) then 
      !     write(buffer,'(a,a,I5,a)')' setup_pennes_model: [optical] ',  &
      !                            'size of power file = ', numlines , ' lines'
      !     call PetscPrintf(PETSC_COMM_SELF,trim(buffer)//char(10)//char(0),ierr)
      !  endif
      !  ! ensure that the final time goes all the way to the final
      !  !  time of interest to the fem calculations
      !  if(powfile(numlines)%time.lt.Maxtime)                            &
      !                        powfile(numlines)%time=Maxtime

        numpow=ceiling(Maxtime/IDEAL_DT)
        allocate(POW(0:numpow),stat=ierr)
        if (ierr /= 0) then
          write(*,*) 'Memory error!  setup_pennes_model.F' ; call abort
        end if
      !                      POW(0)%power SHOULD NEVER BE USED POW(0)%power 
      !                          is set to 1.0d0 for plotting purposes
        POW(0)%time=0.0d0 ;  POW(0)%power=1.0d0
        do i = 1, numpow
           POW(i)%time = IDEAL_DT*i
           POW(i)%power= powerdata
        enddo
      !  do i = 1, numlines-1
      !     jlo = int(powfile(i-1)%time/IDEAL_DT)+1
      !     jhi = int(powfile( i )%time/IDEAL_DT)
      !     do j =  jlo,jhi
      !        POW(j)%time =IDEAL_DT*j
      !        POW(j)%power=powfile( i )%power
      !     enddo
      !  enddo
      !     jlo = int(powfile(numlines-1)%time/IDEAL_DT)+1
      !     jhi = ceiling(powfile( numlines )%time/IDEAL_DT)
      !     do j =  jlo,jhi
      !        POW(j)%time =IDEAL_DT*j
      !        POW(j)%power=powfile( i )%power
      !     enddo
      !  deallocate(powfile,stat=ierr)
      end subroutine setup_power
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      ! set time window over which fem computations will be done 
      subroutine setfinaltime(FinalTime)
      use pennes_model, only : TAU
      implicit none
#include "finclude/petsc.h"
#include "cinout.blk"
        PetscScalar, intent(in ):: FinalTime ! final time instance for this 
                                             ! optimization step 
                                             ! (used in verification problems)
        ! final time (used in verification problems)
        TAU=FinalTime
      end subroutine setfinaltime
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      ! set time window over which fem computations will be done 
      subroutine setfemwindow(Ideallo,Idealhi,Nsteplo,Nstephi,FinalTime)
      use global_params, only : TAO_MONITOR
      use pennes_model, only : IDEAL_NZERO,IDEAL_NTIME,TAU
      implicit none
#include "finclude/petsc.h"
#include "cinout.blk"
        PetscInt, intent(in ):: Ideallo,Idealhi ! lower/upper time bounds on 
                                                ! ideal temp data
        PetscInt, intent(in ):: Nsteplo,Nstephi ! lower/upper time bounds on 
                                                ! fem computations 
        PetscInt :: ipetscrank,ierr
        PetscScalar, intent(in ):: FinalTime ! final time instance for this 
                                             ! optimization step 
                                             ! (used in verification problems)
        character(len=MAXLEN) :: buffer    ! scratch storage

        ! this stored for the routine manipcntrlgrad 
        IDEAL_NZERO = Ideallo   
        IDEAL_NTIME = Idealhi   

        ! final time (used in verification problems)
        TAU=FinalTime

        if(TAO_MONITOR .eqv. PETSC_TRUE)then
          call MPI_Comm_rank(PETSC_COMM_WORLD,ipetscrank,ierr);
          if(ipetscrank .eq.0 ) then
             write(buffer,'(a,2(a,i5),2(a,i5))')                              &
                      "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"//char(10),  &
                                 " IDEAL_NZERO= ", IDEAL_NZERO,               &
                                 " IDEAL_NTIME= ", IDEAL_NTIME,               &
                       char(10)//" Nsteplo= ",Nsteplo,                        &
                                 " Nstephi= ",Nstephi
             call PetscPrintf(PETSC_COMM_SELF,trim(buffer)//char(10)//char(0),ierr)
          endif
          write(NOUT,*)"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
          write(NOUT,'(2(a,i5))')" IDEAL_NZERO= ", IDEAL_NZERO, &
                                 " IDEAL_NTIME= ", IDEAL_NTIME
          write(NOUT,'(2(a,i5))')" Nsteplo= ",Nsteplo, &
                                 " Nstephi= ",Nstephi
        endif
      end subroutine setfemwindow
!----------------------------------------------------------------------
!     update the iteration count for each function/gradient evaluation
!----------------------------------------------------------------------
      subroutine update_iterno(NumOptVar,Global_Iter,Idopt,Ndof_control)
        use global_params , only: ITERNO, TAO_MONITOR
        use pennes_model
        implicit none
#include "cinout.blk"
#include "finclude/petsc.h"
        PetscInt,intent(in) :: NumOptVar,Global_Iter,Idopt,Ndof_control
        PetscInt :: jlo,jhi,jjj,ipetscrank
        character(len=MAXLEN) :: buffer    ! scratch storage
        PetscErrorCode :: ierr    
        ITERNO=ITERNO+1
        if(TAO_MONITOR .eqv. PETSC_TRUE )then
          call MPI_Comm_rank(PETSC_COMM_WORLD,ipetscrank,ierr);
          if(ipetscrank .eq.0 ) then
             write(buffer,'(a,5(a,i7))')                             &
               "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"//char(10),&
                      "optimization #", Idopt,                       &
                      " # of optimization parameters=",NumOptVar,    &
                      "  iter= ",Global_Iter,                        &
            char(10)//"# of control dofs=",Ndof_control,             &
                      "  func_eval= ",ITERNO
           call PetscPrintf(PETSC_COMM_SELF,trim(buffer)//char(10)//char(0),ierr)
          endif
          jlo=IDEAL_NZERO+1
          jhi=IDEAL_NTIME
          write(NOUT,*)"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          write(NOUT,'(3(a,i7))') "optimization #", Idopt,                 &
                   " # of optimization parameters=",NumOptVar,             &
                                         "  iter= ",Global_Iter
          write(NOUT,'(2(a,i7))') "# of control dofs=",Ndof_control,       &
                               "  func_eval= ",ITERNO
          ! 
          write(NOUT,'(2(a,E12.5))') "    RHO= ",RHO," C_P= ",C_P
          ! FIXME  the print out does not include updates from other procs
          ! write out parameter and check bounds 
          ! write out parameter and check bounds 
          !write(NOUT,'(2(a,E12.5))') "    K_0_LB= ",K_0_LB," K_0_UB= ",K_0_UB
          !do jjj = 0 , size(K_0) - 1 
          !  write(NOUT,'(a,i4,a,E12.5)') "    K_0(",jjj,")=",K_0(jjj)
          !  if(K_0(jjj) .lt. K_0_LB .or.  K_0(jjj) .gt. K_0_UB)            &
          !  write(NOUT,'(a,i7,a)') "***NOTE***K_0(",jjj,")out of bounds"
          !enddo
          ! write out parameter and check bounds 
          write(NOUT,'(3(a,E12.5))') "    K_1= ",K_1," K_1_LB= ",K_1_LB,   &
                                                  " K_1_UB= ",K_1_UB
          if(K_1 .lt. K_1_LB .or.  K_1 .gt. K_1_UB)                        &
             write(NOUT,'(a,i7,a)') "***NOTE***K_1 out of bounds"
          ! write out parameter and check bounds 
          write(NOUT,'(3(a,E12.5))') "    K_2= ",K_2," K_2_LB= ",K_2_LB,   &
                                                  " K_2_UB= ",K_2_UB
          if(K_2 .lt. K_2_LB .or.  K_2 .gt. K_2_UB)                        &
             write(NOUT,'(a,i7,a)') "***NOTE***K_2 out of bounds"
          write(NOUT,'(3(a,E12.5))') "    K_3= ",K_3," K_3_LB= ",K_3_LB,   &
                                                  " K_3_UB= ",K_3_UB
          ! write out parameter and check bounds 
          if(K_3 .lt. K_3_LB .or.  K_3 .gt. K_3_UB)                        &
             write(NOUT,'(a,i7,a)') "***NOTE***K_3 out of bounds"
          write(NOUT,'(2(a,E12.5))') "    W_0_LB= ",W_0_LB," W_0_UB= ",W_0_UB
          ! FIXME  the print out does not include updates from other procs
          ! write out parameter and check bounds 
          !do jjj = 0 , size(W_0) - 1 
          !  write(NOUT,'(a,i4,a,E12.5)') "    W_0(",jjj,")=",W_0(jjj)
          !  if(W_0(jjj) .lt. W_0_LB .or.  W_0(jjj) .gt. W_0_UB)            &
          !  write(NOUT,'(a,i7,a)') "***NOTE***W_0(",jjj,")out of bounds"
          !enddo
          ! write out parameter and check bounds 
          write(NOUT,'(3(a,E12.5))') "    W_N= ",W_N," W_N_LB= ",W_N_LB,   &
                                                  " W_N_UB= ",W_N_UB
          if(W_N .lt. W_N_LB .or.  W_N .gt. W_N_UB)                        &
             write(NOUT,'(a,i7,a)') "***NOTE***W_N out of bounds"
          ! write out parameter and check bounds 
          write(NOUT,'(3(a,E12.5))') "    W_I= ",W_I," W_I_LB= ",W_I_LB,   &
                                                  " W_I_UB= ",W_I_UB
          if(W_I .lt. W_I_LB .or.  W_I .gt. W_I_UB)                        &
             write(NOUT,'(a,i7,a)') "***NOTE***W_I out of bounds"
          ! write out parameter and check bounds 
          write(NOUT,'(3(a,E12.5))') "    W_D= ",W_D," W_D_LB= ",W_D_LB,   &
                                                  " W_D_UB= ",W_D_UB
          if(W_D .lt. W_D_LB .or.  W_D .gt. W_D_UB)                        &
             write(NOUT,'(a,i7,a)') "***NOTE***W_D out of bounds"
          ! write out parameter and check bounds 
          write(NOUT,'(3(a,E12.5))') "    W_2= ",W_2," W_2_LB= ",W_2_LB,   &
                                                  " W_2_UB= ",W_2_UB
          if(W_2 .lt. W_2_LB .or.  W_2 .gt. W_2_UB)                        &
             write(NOUT,'(a,i7,a)') "***NOTE***W_2 out of bounds"
          ! write out parameter and check bounds 
          write(NOUT,'(5(a,E12.5))') "    W_NI= ",W_NI, " W_ID= ",W_ID,    &
            " W_NID_LB= ",W_NID_LB, " W_NID_MD= ",                         &
              W_NID_MD ," W_NID_UB= ",W_NID_UB
          if(W_NI .lt. W_NID_LB .or.  W_NI .gt. W_NID_MD)                  & 
             write(NOUT,'(a,i7,a)') "***NOTE***W_NI out of bounds"
          if(W_ID .lt. W_NID_MD .or.  W_ID .gt. W_NID_UB)                  &
             write(NOUT,'(a,i7,a)') "***NOTE***W_ID out of bounds"
          ! write out parameter and check bounds 
          write(NOUT,'(3(a,E12.5))') "    X_0= ",X_0," X_0_LB= ",X_0_LB,   &
                                                  " X_0_UB= ",X_0_UB
          if(X_0 .lt. X_0_LB .or.  X_0 .gt. X_0_UB)                        &
             write(NOUT,'(a,i7,a)') "***NOTE***X_0 out of bounds"
          ! write out parameter and check bounds 
          write(NOUT,'(3(a,E12.5))') "    Y_0= ",Y_0," Y_0_LB= ",Y_0_LB,   &
                                                  " Y_0_UB= ",Y_0_UB
          if(Y_0 .lt. Y_0_LB .or.  Y_0 .gt. Y_0_UB)                        &
             write(NOUT,'(a,i7,a)') "***NOTE***Y_0 out of bounds"
          ! write out parameter and check bounds 
          write(NOUT,'(3(a,E12.5))') "    Z_0= ",Z_0," Z_0_LB= ",Z_0_LB,   &
                                                  " Z_0_UB= ",Z_0_UB
          if(Z_0 .lt. Z_0_LB .or.  Z_0 .gt. Z_0_UB)                        &
             write(NOUT,'(a,i7,a)') "***NOTE***Z_0 out of bounds"
          ! write out parameter and check bounds 
          write(NOUT,'(3(a,E12.5))') "    MU_A=",MU_A, " MU_A_LB= ",       &
                                       MU_A_LB," MU_A_UB= ",MU_A_UB
          if(MU_A .lt. MU_A_LB .or. MU_A .gt. MU_A_UB)                     &
             write(NOUT,'(a,i7,a)') "***NOTE***MU_A out of bounds"
          ! write out parameter and check bounds 
          write(NOUT,'(3(a,E12.5))') "    MU_S=",MU_S, " MU_S_LB= ",       &
                                       MU_S_LB," MU_S_UB= ",MU_S_UB
          if(MU_S .lt. MU_S_LB .or. MU_S .gt. MU_S_UB)                     &
             write(NOUT,'(a,i7,a)') "***NOTE***MU_S out of bounds"
          ! write out parameter and check bounds 
          write(NOUT,'(2(a,E12.5))') "    POW_LB= ",POW_LB,                &
                                     " POW_UB= ",POW_UB
          do jjj = jlo,jhi
            write(NOUT,'(a,i4,a,E12.5)') "    POW(",jjj,")=",              &
                                                POW(jjj)%power
            if(POW(jjj)%power.lt.POW_LB.or.POW(jjj)%power.gt.POW_UB)       &
            write(NOUT,'(a,i7,a)') "***NOTE***POW(",jjj,")out of bounds"
          enddo
          write(NOUT,*)"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          write(NOUT,*)" "
          call flush(NOUT)
        endif
      end subroutine update_iterno
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      ! used for forward rule integration in time of QOI
      subroutine adjloadtemp(Xpoint,Deltat,Idstep,Nstephi,Soln,Nsize, &  
                                         Adiff,Bconv,Creac,Dualsource)
        use pennes_model, only: RHO,C_P
        implicit none
        PetscInt,   intent(in) :: Idstep,Nstephi,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Deltat,Soln(0:Nsize-1)
        PetscScalar,intent(out):: Adiff,Bconv,Creac,Dualsource
        PetscScalar ::  Usln_prev,Usln,Usln_next,Qoiload,Qoiload_prev
        Usln         = Soln(0)
        Usln_prev    = Soln(1)
        Usln_next    = Soln(2)
        Qoiload      = Soln(3)
        Qoiload_prev = Soln(4)
        ! qoi = 1/2 int(rho * c_p *(t - t_exp)^2)
        Dualsource = RHO*C_P*(usln-qoiload)      
        ! inline
        call nonlinpennesjac(Deltat,Usln,Usln_next,Adiff,Bconv,Creac)
      end subroutine adjloadtemp
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      ! used for mid point rule integration in time of QOI
      subroutine adjloadtempmid(Xpoint,Deltat,Idstep,Nstephi,Soln,Nsize, &  
                                            Adiff,Bconv,Creac,Dualsource)
        use pennes_model, only: RHO,C_P
        implicit none
        PetscInt,   intent(in) :: Idstep,Nstephi,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Deltat,Soln(0:Nsize-1)
        PetscScalar,intent(out):: Adiff,Bconv,Creac,Dualsource
        PetscScalar :: Usln_prev   ,Usln   ,Usln_next,     &
                       Qoiload_prev,Qoiload,Qoiload_next 
        PetscScalar, parameter  :: onesixth=0.16666666666666d0, &
                                   onethird=0.33333333333333d0
        Usln         = Soln(0)
        Usln_prev    = Soln(1)
        Usln_next    = Soln(2)
        Qoiload      = Soln(3)
        Qoiload_prev = Soln(4)
        Qoiload_next = Soln(5)
        ! qoi = 1/2 int(rho * c_p *(t - t_exp)^2)
        if(Idstep .lt. Nstephi) then 
           Dualsource= RHO*C_P*onesixth*( (usln_prev-qoiload_prev)+   &
                                        4.0d0*(usln     -qoiload     )+   &
                                           (usln_next-qoiload_next))
        else ! slightly different for first time step backwards
           Dualsource= RHO*C_P*( onesixth*(usln_prev-qoiload_prev)      &
                                  +  onethird*(usln     -qoiload     ) )
        endif
        call nonlinpennesjac(Deltat,Usln,Usln_next,Adiff,Bconv,Creac)
      end subroutine adjloadtempmid
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine adjloaddam_arr(Xpoint,Deltat,Idstep,Nstephi,Soln,Nsize,&  
                                         Adiff,Bconv,Creac,Dualsource)
        use pennes_model , ONLY: Arr_A,Arr_Ea,Arr_R
        implicit none
        PetscInt,   intent(in) :: Idstep,Nstephi,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Deltat,Soln(0:Nsize-1)
        PetscScalar,intent(out):: Adiff,Bconv,Creac,Dualsource
        PetscScalar :: u_kmhalf , u_kphalf, dphidu_kmhalf, dphidu_kphalf, &
                       Usln_prev   ,Usln   ,Usln_next,        & 
                       Qoiload_zero,Qoiload_full 

        Usln         = Soln(0)
        Usln_prev    = Soln(1)
        Usln_next    = Soln(2)
        Qoiload_zero = Soln(6)
        Qoiload_full = Soln(7)
        ! compute solution at half time steps
        u_kmhalf = 0.5d0*(Usln + Usln_prev)
        u_kphalf = 0.5d0*(Usln + Usln_next)

        dphidu_kmhalf = & ! derivative of damage
               Arr_A*Arr_Ea/Arr_R/u_kmhalf**2/exp(Arr_Ea/Arr_R/u_kmhalf)

        if(Idstep .lt. Nstephi) then
           dphidu_kphalf =  & ! derivative of damage
               Arr_A*Arr_Ea/Arr_R/u_kphalf**2/exp(Arr_Ea/Arr_R/u_kphalf)
        else ! slightly different for first time step backwards
           dphidu_kphalf = 0.0d0
        endif
        Dualsource= (Qoiload_full-Qoiload_zero)*   &
                         0.5d0*(dphidu_kmhalf + dphidu_kphalf)
        ! inline
        call nonlinpennesjac(Deltat,Usln,Usln_next,Adiff,Bconv,Creac)
      end subroutine adjloaddam_arr
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine adjloaddam__ts(Xpoint,Deltat,Idstep,Nstephi,Soln,Nsize,&  
                                         Adiff,Bconv,Creac,Dualsource)
        use pennes_model, ONLY: TS_h,TS_alpha,TS_beta
        implicit none
        PetscInt,   intent(in) :: Idstep,Nstephi,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Deltat,Soln(0:Nsize-1)
        PetscScalar,intent(out):: Adiff,Bconv,Creac,Dualsource
        PetscScalar :: Usln_prev   ,Usln   ,Usln_next,    &
                       Qoiload_zero,Qoiload_full 
        PetscScalar :: t,u_kmhalf,u_kphalf,dphidu_kmhalf,dphidu_kphalf

        Usln         = Soln(0)
        Usln_prev    = Soln(1)
        Usln_next    = Soln(2)
        Qoiload_zero = Soln(6)
        Qoiload_full = Soln(7)

        ! compute solution at half time steps
        u_kmhalf = 0.5d0*(Usln + Usln_prev)
        u_kphalf = 0.5d0*(Usln + Usln_next)

        t=(Deltat*(Idstep) +Deltat*(Idstep-1))*0.5d0

        dphidu_kmhalf = &  ! derivative of damage
           TS_h*exp(-(TS_h/u_kmhalf+TS_alpha*t+TS_beta))/    &
             (1+exp(-(TS_h/u_kmhalf+TS_alpha*t+TS_beta)))**2/u_kmhalf**2

        if(Idstep .lt. Nstephi) then
           dphidu_kphalf =  & ! derivative of damage
           TS_h*exp(-(TS_h/u_kphalf+TS_alpha*t+TS_beta))/  &
             (1+exp(-(TS_h/u_kphalf+TS_alpha*t+TS_beta)))**2/u_kphalf**2
        else ! slightly different for first time step backwards
           dphidu_kphalf = 0.0d0
        endif
        Dualsource = (Qoiload_full-Qoiload_zero)*    &
                         0.5d0*(dphidu_kmhalf + dphidu_kphalf)
        ! inline
        call nonlinpennesjac(Deltat,Usln,Usln_next,Adiff,Bconv,Creac)
      end subroutine adjloaddam__ts
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine adjloadhsp(Xpoint,Deltat,Idstep,Nstephi,Soln,Nsize,&  
                                         Adiff,Bconv,Creac,Dualsource)
        use pennes_model, only: RHO,C_P
        implicit none
        PetscInt,   intent(in) :: Idstep,Nstephi,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Deltat,Soln(0:Nsize-1)
        PetscScalar,intent(out):: Adiff,Bconv,Creac,Dualsource
        PetscScalar :: Usln, Usln_prev,Usln_next
        Usln         = Soln(0)
        Usln_prev    = Soln(1)
        Usln_next    = Soln(2)
        ! inline
        call nonlinpennesjac(Deltat,Usln,Usln_next,Adiff,Bconv,Creac)
        Dualsource = 0.0d0
      end subroutine adjloadhsp
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine getdualsource(Xpoint,Deltat,Idstep,Nstephi,Soln,Nsize,&  
                                         Adiff,Bconv,Creac,Dualsource)
        use global_params, only: L,R
        use pennes_model, only: RHO,C_P,COEFF_COOL,K_0,POW, &
                                W_0,U_A,NEXACTVERIFNUMBER,TAU
        implicit none
        PetscInt,   intent(in) :: Idstep,Nstephi,Nsize
        PetscScalar,intent(in) :: Xpoint(3),Deltat,Soln(0:Nsize-1)
        PetscScalar,intent(out):: Adiff,Bconv,Creac,Dualsource
        PetscScalar :: t,xcoord,q0,Usln, Usln_prev,Usln_next
        Usln         = Soln(0)
        Usln_prev    = Soln(1)
        Usln_next    = Soln(2)

        ! inline
        call nonlinpennesjac(Deltat,Usln,Usln_next,Adiff,Bconv,Creac)

        xcoord = Xpoint(1)
        q0=POW(1)%power
            select case(NEXACTVERIFNUMBER) 
            case(1)
               Dualsource= -1.0d0/K_0(0)*exp(xcoord/K_0(0))
            case(2)
               t=(Deltat*(Idstep) +Deltat*(Idstep-1))/2.0d0
               Dualsource=(RHO*C_P*q0*(L-xcoord)**2-2*K_0(0)*q0*(TAU-t))
            case(8)
               t=Deltat*Idstep
               Dualsource= 0.0d0
!               Dualsource= -RHO*C_P*(3.0d0*(xpoint(1)-L)**3/
!     .              (3.0d0*K_0(0)+3.0d0*K_1*atan(K_2/exp(t)*L**2)
!     .              +COEFF_COOL*L)**2
!     .              *K_1*K_2/exp(t)/(1.0d0+K_2**2/exp(t)**2*L**4)*
!     .              (TAU-t)-((xpoint(1)-L)**3*COEFF_COOL/L**2/
!     .              (3.0d0*K_0(0)+3.0d0*K_1*atan(K_2/exp(t)*L**2)
!     .              +COEFF_COOL*L)+1.0d0)/COEFF_COOL)-(6.0d0*K_0(0)+
!     .              6.0d0*K_1*atan(K_2*(L-xpoint(1))**2*exp(-t)))
!     .              *(xpoint(1)-L)/L**2/(3.0d0*K_0(0)+
!     .              3.0d0*K_1*atan(K_2/exp(t)*L**2)+COEFF_COOL*L)*
!     .              (TAU-t)+ (W_0(1)+W_1*atan(W_2*((L-xpoint(1))**2* 
!     .              exp(-t)+K_3-W_3)))*((xpoint(1)-L)**3*COEFF_COOL
!     .              /L**2/ (3.0d0*K_0(0)+3.0d0*K_1*atan(K_2/exp(t)*L**2)
!     .              +COEFF_COOL*L)+1.0d0) /COEFF_COOL*(TAU-t)+
!     .              W_1*W_2/(1.0d0+W_2**2*((L-xpoint(1))**2
!     .              *exp(-t)+K_3-W_3)**2)*((L-xpoint(1))**2*exp(-t)
!     .              +K_3-U_A)*((xpoint(1)-L)**3*COEFF_COOL/L**2/
!     .              (3.0d0*K_0(0)+3.0d0*K_1*atan(K_2/exp(t)*L**2)+
!     .              COEFF_COOL*L)+1.0d0) /COEFF_COOL*(TAU-t)
            end select
      end subroutine getdualsource
!----------------------------------------------------------------------
!  
! The following are functions to add regularization terms to the QOI and 
! its gradient for W_0
!
!----------------------------------------------------------------------
      function regularization_w_0(Id,Scalefact)
         use pennes_model , only: W_0,W_0_LB,W_0_UB
         implicit none
         PetscScalar :: regularization_w_0 
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in) :: Scalefact
         regularization_w_0 =                                         &
                   +exp( Scalefact/(W_0_UB-W_0_LB)*(W_0(Id)-W_0_UB) ) &
                   +exp(-Scalefact/(W_0_UB-W_0_LB)*(W_0(Id)-W_0_LB) )
      end function regularization_w_0
!
      function regularizderiv_w_0(Id,Scalefact)
         use pennes_model , only: W_0,W_0_LB,W_0_UB
         implicit none
         PetscScalar :: regularizderiv_w_0
         PetscInt    , intent(in)  :: Id
         PetscScalar , intent(in)  :: Scalefact
         regularizderiv_w_0  =                                           &
                        Scalefact/(W_0_UB-W_0_LB) * (                    &
                   exp( Scalefact/(W_0_UB-W_0_LB)*(W_0(Id)-W_0_UB) )     &
                  -exp(-Scalefact/(W_0_UB-W_0_LB)*(W_0(Id)-W_0_LB) ) )
      end function regularizderiv_w_0
!----------------------------------------------------------------------
!  
! The following are functions to add regularization terms to the QOI and 
! its gradient for K_0
!
!----------------------------------------------------------------------
      function regularization_k_0(Id,Scalefact)
         use pennes_model , only: K_0,K_0_LB,K_0_UB
         implicit none
         PetscScalar :: regularization_k_0 
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in) :: Scalefact
         regularization_k_0 =                                            &
                   +exp( Scalefact/(K_0_UB-K_0_LB)*(K_0(Id)-K_0_UB) )    &
                   +exp(-Scalefact/(K_0_UB-K_0_LB)*(K_0(Id)-K_0_LB) )
      end function regularization_k_0
!
      function regularizderiv_k_0(Id, Scalefact)
         use pennes_model , only: K_0,K_0_LB,K_0_UB
         implicit none
         PetscScalar :: regularizderiv_k_0
         PetscInt    , intent(in)  :: Id
         PetscScalar , intent(in)  :: Scalefact
         regularizderiv_k_0 =                                             &
                         Scalefact/(K_0_UB-K_0_LB) * (                    &
                    exp( Scalefact/(K_0_UB-K_0_LB)*(K_0(Id)-K_0_UB) )     &
                   -exp(-Scalefact/(K_0_UB-K_0_LB)*(K_0(Id)-K_0_LB) ) )
      end function regularizderiv_k_0
!----------------------------------------------------------------------
!  
! The following are functions to add regularization terms to the QOI and 
! its gradient for the power term
!
!----------------------------------------------------------------------
     function regularization_pow(Id,Scalefact)
         use pennes_model , only: POW,POW_LB,POW_UB , IDEAL_NZERO
         implicit none
         PetscScalar :: regularization_pow 
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in)  :: Scalefact
         PetscInt :: jjj
         jjj = IDEAL_NZERO+1
         regularization_pow =                                              &
            +exp( Scalefact/(POW_UB-POW_LB)*(POW(jjj+Id)%power-POW_UB))    &
            +exp(-Scalefact/(POW_UB-POW_LB)*(POW(jjj+Id)%power-POW_LB))     
     end function regularization_pow                                        
!                                                                           
     function regularizderiv_pow(Id,Scalefact)                              
        use pennes_model , only: POW,POW_LB,POW_UB , IDEAL_NZERO                                  
        implicit none                                                       
        PetscScalar :: regularizderiv_pow                                   
        PetscInt    , intent(in) :: Id                                      
        PetscScalar , intent(in)  :: Scalefact                              
        PetscInt :: jjj                                                     
        jjj = IDEAL_NZERO+1                                                 
        regularizderiv_pow =                                               &
          +    Scalefact/(POW_UB-POW_LB) * (                               &
          exp( Scalefact/(POW_UB-POW_LB)*(POW(jjj+Id)%power-POW_UB) )      &
         -exp(-Scalefact/(POW_UB-POW_LB)*(POW(jjj+Id)%power-POW_LB) ) )
     end function regularizderiv_pow
!----------------------------------------------------------------------
!  
! The following are functions to add regularization to the OQI 
! and derivative of the regularization terms to the gradient for the W_N term
!
!----------------------------------------------------------------------
      function regularization_w_n(Id,Scalefact)
         use pennes_model , only: W_N,W_N_LB,W_N_UB
         implicit none
         PetscScalar :: regularization_w_n 
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in)  :: Scalefact
         regularization_w_n =                                   &
             +exp( Scalefact/(W_N_UB-W_N_LB)*(W_N-W_N_UB))      &
             +exp(-Scalefact/(W_N_UB-W_N_LB)*(W_N-W_N_LB))
      end function regularization_w_n
!
      function regularizderiv_w_n(Id,Scalefact)
         use pennes_model , only: W_N,W_N_LB,W_N_UB
         implicit none
         PetscScalar :: regularizderiv_w_n
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in)  :: Scalefact
         regularizderiv_w_n =                                      &
                 +    Scalefact/(W_N_UB-W_N_LB) * (                &
                 exp( Scalefact/(W_N_UB-W_N_LB)*(W_N-W_N_UB) )     &
                -exp(-Scalefact/(W_N_UB-W_N_LB)*(W_N-W_N_LB) ) )
      end function regularizderiv_w_n
!----------------------------------------------------------------------
!  
! The following are functions to add regularization to the OQI 
! and derivative of the regularization terms to the gradient for the W_I term
!
!----------------------------------------------------------------------
      function regularization_w_i(Id,Scalefact)
         use pennes_model , only: W_I,W_I_LB,W_I_UB
         implicit none
         PetscScalar :: regularization_w_i 
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in)  :: Scalefact
         regularization_w_i =                                     &
             +exp( Scalefact/(W_I_UB-W_I_LB)*(W_I-W_I_UB))        &
             +exp(-Scalefact/(W_I_UB-W_I_LB)*(W_I-W_I_LB))
      end function regularization_w_i
!
      function regularizderiv_w_i(Id,Scalefact)
         use pennes_model , only: W_I,W_I_LB,W_I_UB
         implicit none
         PetscScalar :: regularizderiv_w_i
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in)  :: Scalefact
         regularizderiv_w_i =                                       &
                 +    Scalefact/(W_I_UB-W_I_LB) * (                 &
                 exp( Scalefact/(W_I_UB-W_I_LB)*(W_I-W_I_UB) )      &
                -exp(-Scalefact/(W_I_UB-W_I_LB)*(W_I-W_I_LB) ) )
      end function regularizderiv_w_i
!----------------------------------------------------------------------
!  
! The following are functions to add regularization to the OQI 
! and derivative of the regularization terms to the gradient for the W_D term
!
!----------------------------------------------------------------------
      function regularization_w_d(Id,Scalefact)
         use pennes_model , only: W_D,W_D_LB,W_D_UB
         implicit none
         PetscScalar :: regularization_w_d 
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in)  :: Scalefact
         regularization_w_d =                                        &
             +exp( Scalefact/(W_D_UB-W_D_LB)*(W_D-W_D_UB))           &
             +exp(-Scalefact/(W_D_UB-W_D_LB)*(W_D-W_D_LB))
      end function regularization_w_d
!
      function regularizderiv_w_d(Id,Scalefact)
         use pennes_model , only: W_D,W_D_LB,W_D_UB
         implicit none
         PetscScalar :: regularizderiv_w_d
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in)  :: Scalefact
         regularizderiv_w_d =                                           &
                 +    Scalefact/(W_D_UB-W_D_LB) * (                     &
                 exp( Scalefact/(W_D_UB-W_D_LB)*(W_D-W_D_UB) )          &
                -exp(-Scalefact/(W_D_UB-W_D_LB)*(W_D-W_D_LB) ) )
      end function regularizderiv_w_d
!----------------------------------------------------------------------
!  
! The following are functions to add regularization to the OQI 
! and derivative of the regularization terms to the gradient for the W_2 term
!
!----------------------------------------------------------------------
      function regularization_w_2(Id,Scalefact)
         use pennes_model , only: W_2,W_2_LB,W_2_UB
         implicit none
         PetscScalar :: regularization_w_2 
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in)  :: Scalefact
         regularization_w_2 =                                      &
             +exp( Scalefact/(W_2_UB-W_2_LB)*(W_2-W_2_UB))         &
             +exp(-Scalefact/(W_2_UB-W_2_LB)*(W_2-W_2_LB))
      end function regularization_w_2
!
      function regularizderiv_w_2(Id,Scalefact)
         use pennes_model , only: W_2,W_2_LB,W_2_UB
         implicit none
         PetscScalar :: regularizderiv_w_2
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in)  :: Scalefact
         regularizderiv_w_2 =                                         &
                 +    Scalefact/(W_2_UB-W_2_LB) * (                   &
                 exp( Scalefact/(W_2_UB-W_2_LB)*(W_2-W_2_UB) )        &
                -exp(-Scalefact/(W_2_UB-W_2_LB)*(W_2-W_2_LB) ) )
      end function regularizderiv_w_2
!----------------------------------------------------------------------
!  
! The following are functions to add regularization to the OQI 
! and derivative of the regularization terms to the gradient for the W_NI term
!
!----------------------------------------------------------------------
      function regularization_w_ni(Id,Scalefact)
         use pennes_model , only: W_NI,W_NID_LB,W_NID_MD,W_NID_UB
         implicit none
         PetscScalar :: regularization_w_ni 
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in)  :: Scalefact
         regularization_w_ni =                                          &
             +exp( Scalefact/(W_NID_MD-W_NID_LB)*(W_NI-W_NID_MD))       &
             +exp(-Scalefact/(W_NID_MD-W_NID_LB)*(W_NI-W_NID_LB))
      end function regularization_w_ni
!
      function regularizderiv_w_ni(Id,Scalefact)
         use pennes_model , only: W_NI,W_NID_LB,W_NID_MD,W_NID_UB
         implicit none
         PetscScalar :: regularizderiv_w_ni
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in)  :: Scalefact
         regularizderiv_w_ni =                                               &
                 +    Scalefact/(W_NID_MD-W_NID_LB) * (                      &
                 exp( Scalefact/(W_NID_MD-W_NID_LB)*(W_NI-W_NID_MD) )        &
                -exp(-Scalefact/(W_NID_MD-W_NID_LB)*(W_NI-W_NID_LB) ) )
      end function regularizderiv_w_ni
!----------------------------------------------------------------------
!  
! The following are functions to add regularization to the OQI 
! and derivative of the regularization terms to the gradient for the W_ID term
!
!----------------------------------------------------------------------
      function regularization_w_id(Id,Scalefact)
         use pennes_model , only: W_ID,W_NID_LB,W_NID_MD,W_NID_UB
         implicit none
         PetscScalar :: regularization_w_id 
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in)  :: Scalefact
         regularization_w_id =                                           &
             +exp( Scalefact/(W_NID_UB-W_NID_MD)*(W_ID-W_NID_UB))        &
             +exp(-Scalefact/(W_NID_UB-W_NID_MD)*(W_ID-W_NID_MD))
      end function regularization_w_id
!
      function regularizderiv_w_id(Id,Scalefact)
         use pennes_model , only: W_ID,W_NID_LB,W_NID_MD,W_NID_UB
         implicit none
         PetscScalar :: regularizderiv_w_id
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in)  :: Scalefact
         regularizderiv_w_id =                                            &
                 +    Scalefact/(W_NID_UB-W_NID_MD) * (                   &
                 exp( Scalefact/(W_NID_UB-W_NID_MD)*(W_ID-W_NID_UB) )     &
                -exp(-Scalefact/(W_NID_UB-W_NID_MD)*(W_ID-W_NID_MD) ) )
      end function regularizderiv_w_id
!----------------------------------------------------------------------
!  
! The following are functions to add regularization to the OQI 
! and derivative of the regularization terms to the gradient for the K_1 term
!
!----------------------------------------------------------------------
      function regularization_k_1(Id,Scalefact)
         use pennes_model , only: K_1,K_1_LB,K_1_UB
         implicit none
         PetscScalar :: regularization_k_1 
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in)  :: Scalefact
         regularization_k_1 =                                             &
             +exp( Scalefact/(K_1_UB-K_1_LB)*(K_1-K_1_UB))                &
             +exp(-Scalefact/(K_1_UB-K_1_LB)*(K_1-K_1_LB))
      end function regularization_k_1
!
      function regularizderiv_k_1(Id,Scalefact)
         use pennes_model , only: K_1,K_1_LB,K_1_UB
         implicit none
         PetscScalar :: regularizderiv_k_1
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in)  :: Scalefact
         regularizderiv_k_1 =                                               &
                 +    Scalefact/(K_1_UB-K_1_LB) * (                         &
                 exp( Scalefact/(K_1_UB-K_1_LB)*(K_1-K_1_UB) )              &
                -exp(-Scalefact/(K_1_UB-K_1_LB)*(K_1-K_1_LB) ) )
      end function regularizderiv_k_1
!----------------------------------------------------------------------
!  
! The following are functions to add regularization to the OQI 
! and derivative of the regularization terms to the gradient for the K_2 term
!
!----------------------------------------------------------------------
      function regularization_k_2(Id,Scalefact)
         use pennes_model , only: K_2,K_2_LB,K_2_UB
         implicit none
         PetscScalar :: regularization_k_2 
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in)  :: Scalefact
         regularization_k_2 =                                               &
             +exp( Scalefact/(K_2_UB-K_2_LB)*(K_2-K_2_UB))                  &
             +exp(-Scalefact/(K_2_UB-K_2_LB)*(K_2-K_2_LB))
      end function regularization_k_2
!
      function regularizderiv_k_2(Id,Scalefact)
         use pennes_model , only: K_2,K_2_LB,K_2_UB
         implicit none
         PetscScalar :: regularizderiv_k_2
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in)  :: Scalefact
         regularizderiv_k_2 =                                                &
                 +    Scalefact/(K_2_UB-K_2_LB) * (                          &
                 exp( Scalefact/(K_2_UB-K_2_LB)*(K_2-K_2_UB) )               &
                -exp(-Scalefact/(K_2_UB-K_2_LB)*(K_2-K_2_LB) ) )
      end function regularizderiv_k_2
!----------------------------------------------------------------------
!  
! The following are functions to add regularization to the OQI 
! and derivative of the regularization terms to the gradient for the K_3 term
!
!----------------------------------------------------------------------
      function regularization_k_3(Id,Scalefact)
         use pennes_model , only: K_3,K_3_LB,K_3_UB
         implicit none
         PetscScalar :: regularization_k_3 
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in)  :: Scalefact
         regularization_k_3 =                                         &
             +exp( Scalefact/(K_3_UB-K_3_LB)*(K_3-K_3_UB))            &
             +exp(-Scalefact/(K_3_UB-K_3_LB)*(K_3-K_3_LB))
      end function regularization_k_3
!
      function regularizderiv_k_3(Id,Scalefact)
         use pennes_model , only: K_3,K_3_LB,K_3_UB
         implicit none
         PetscScalar :: regularizderiv_k_3
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in)  :: Scalefact
         regularizderiv_k_3 =                                           &
                 +    Scalefact/(K_3_UB-K_3_LB) * (                     &
                 exp( Scalefact/(K_3_UB-K_3_LB)*(K_3-K_3_UB) )          &
                -exp(-Scalefact/(K_3_UB-K_3_LB)*(K_3-K_3_LB) ) )
      end function regularizderiv_k_3
!----------------------------------------------------------------------
!  
! The following are functions to add regularization to the OQI 
! and derivative of the regularization terms to the gradient for the X_0 term
!
!----------------------------------------------------------------------
      function regularization_x_0(Id,Scalefact)
         use pennes_model , only: X_0,X_0_LB,X_0_UB
         implicit none
         PetscScalar :: regularization_x_0 
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in)  :: Scalefact
         regularization_x_0 =                                           &
             +exp( Scalefact/(X_0_UB-X_0_LB)*(X_0-X_0_UB))              &
             +exp(-Scalefact/(X_0_UB-X_0_LB)*(X_0-X_0_LB))
      end function regularization_x_0
!
      function regularizderiv_x_0(Id,Scalefact)
         use pennes_model , only: X_0,X_0_LB,X_0_UB
         implicit none
         PetscScalar :: regularizderiv_x_0
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in)  :: Scalefact
         regularizderiv_x_0 =                                              &
                 +    Scalefact/(X_0_UB-X_0_LB) * (                        &
                 exp( Scalefact/(X_0_UB-X_0_LB)*(X_0-X_0_UB) )             &
                -exp(-Scalefact/(X_0_UB-X_0_LB)*(X_0-X_0_LB) ) )
      end function regularizderiv_x_0
!----------------------------------------------------------------------
!  
! The following are functions to add regularization to the OQI 
! and derivative of the regularization terms to the gradient for the Y_0 term
!
!----------------------------------------------------------------------
      function regularization_y_0(Id,Scalefact)
         use pennes_model , only: Y_0,Y_0_LB,Y_0_UB
         implicit none
         PetscScalar :: regularization_y_0 
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in)  :: Scalefact
         regularization_y_0 =                                          &
             +exp( Scalefact/(Y_0_UB-Y_0_LB)*(Y_0-Y_0_UB))             &
             +exp(-Scalefact/(Y_0_UB-Y_0_LB)*(Y_0-Y_0_LB))
      end function regularization_y_0
!
      function regularizderiv_y_0(Id,Scalefact)
         use pennes_model , only: Y_0,Y_0_LB,Y_0_UB
         implicit none
         PetscScalar :: regularizderiv_y_0
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in)  :: Scalefact
         regularizderiv_y_0 =                                            &
                 +    Scalefact/(Y_0_UB-Y_0_LB) * (                      &
                 exp( Scalefact/(Y_0_UB-Y_0_LB)*(Y_0-Y_0_UB) )           &
                -exp(-Scalefact/(Y_0_UB-Y_0_LB)*(Y_0-Y_0_LB) ) )
      end function regularizderiv_y_0
!----------------------------------------------------------------------
!  
! The following are functions to add regularization to the OQI 
! and derivative of the regularization terms to the gradient for the Z_0 term
!
!----------------------------------------------------------------------
      function regularization_z_0(Id,Scalefact)
         use pennes_model , only: Z_0,Z_0_LB,Z_0_UB
         implicit none
         PetscScalar :: regularization_z_0 
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in)  :: Scalefact
         regularization_z_0 =                                   &
             +exp( Scalefact/(Z_0_UB-Z_0_LB)*(Z_0-Z_0_UB))      &
             +exp(-Scalefact/(Z_0_UB-Z_0_LB)*(Z_0-Z_0_LB))
      end function regularization_z_0
!
      function regularizderiv_z_0(Id,Scalefact)
         use pennes_model , only: Z_0,Z_0_LB,Z_0_UB
         implicit none
         PetscScalar :: regularizderiv_z_0
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in)  :: Scalefact
         regularizderiv_z_0 =                                         &
                 +    Scalefact/(Z_0_UB-Z_0_LB) * (                   &
                 exp( Scalefact/(Z_0_UB-Z_0_LB)*(Z_0-Z_0_UB) )        &
                -exp(-Scalefact/(Z_0_UB-Z_0_LB)*(Z_0-Z_0_LB) ) )
      end function regularizderiv_z_0
!----------------------------------------------------------------------
!  
! The following are functions to add regularization to the OQI 
! and derivative of the regularization terms to the gradient for the MU_A term
!
!----------------------------------------------------------------------
      function regularization_mu_a(Id,Scalefact)
         use pennes_model , only: MU_A,MU_A_LB,MU_A_UB
         implicit none
         PetscScalar :: regularization_mu_a 
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in)  :: Scalefact
         regularization_mu_a =                                        &
             +exp( Scalefact/(MU_A_UB-MU_A_LB)*(MU_A-MU_A_UB))        &
             +exp(-Scalefact/(MU_A_UB-MU_A_LB)*(MU_A-MU_A_LB))
      end function regularization_mu_a
!
      function regularizderiv_mu_a(Id,Scalefact)
         use pennes_model , only: MU_A,MU_A_LB,MU_A_UB
         implicit none
         PetscScalar :: regularizderiv_mu_a
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in)  :: Scalefact
         regularizderiv_mu_a =                                        &
                 +    Scalefact/(MU_A_UB-MU_A_LB) * (                 &
                 exp( Scalefact/(MU_A_UB-MU_A_LB)*(MU_A-MU_A_UB) )    &
                -exp(-Scalefact/(MU_A_UB-MU_A_LB)*(MU_A-MU_A_LB) ) )
      end function regularizderiv_mu_a
!----------------------------------------------------------------------
!  
! The following are functions to add regularization to the OQI 
! and derivative of the regularization terms to the gradient for the MU_S term
!
!----------------------------------------------------------------------
      function regularization_mu_s(Id,Scalefact)
         use pennes_model , only: MU_S,MU_S_LB,MU_S_UB
         implicit none
         PetscScalar :: regularization_mu_s 
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in)  :: Scalefact
         regularization_mu_s =                                        &
             +exp( Scalefact/(MU_S_UB-MU_S_LB)*(MU_S-MU_S_UB))        &
             +exp(-Scalefact/(MU_S_UB-MU_S_LB)*(MU_S-MU_S_LB))
      end function regularization_mu_s
!
      function regularizderiv_mu_s(Id,Scalefact)
         use pennes_model , only: MU_S,MU_S_LB,MU_S_UB
         implicit none
         PetscScalar :: regularizderiv_mu_s
         PetscInt    , intent(in) :: Id 
         PetscScalar , intent(in)  :: Scalefact
         regularizderiv_mu_s =                                           &
                 +    Scalefact/(MU_S_UB-MU_S_LB) * (                    &
                 exp( Scalefact/(MU_S_UB-MU_S_LB)*(MU_S-MU_S_UB) )       &
                -exp(-Scalefact/(MU_S_UB-MU_S_LB)*(MU_S-MU_S_LB) ) )
      end function regularizderiv_mu_s
!----------------------------------------------------------------------
!
!   function name   - getinittemp    (last revision: sept, 2008)
!
!   purpose         - get data at a point to initialize temperature field
!
!   arguments       - IN :  Xpoint - coordinate of the point
!                           Time   - time 
!                     OUT:  temperature value at the point
!----------------------------------------------------------------------
      function getinittemp(Xpoint,Time)
        use global_params  , only: L,R
        use pennes_model 
        implicit none
         PetscScalar,intent(in):: Xpoint(3),Time

         PetscScalar :: getinittemp
         PetscScalar :: radius,xcoord,t,q0

         getinittemp=U_INIT ! Default to U_INIT from input deck

         xcoord = Xpoint(1)
         radius=sqrt(Xpoint(1)**2+Xpoint(2)**2+Xpoint(3)**2)
         ! t= gettime(Idstep)
         t= Time
         q0=POW(1)%power
         select case(NEXACTVERIFNUMBER) 
         case(1)
            getinittemp= 0.0d0
!            getinittemp=(L-xcoord)**2/2.0d0+COEFF_COOL*xcoord+K_0(0)
         case(2)
            getinittemp= q0*(L-xcoord)**2*t+K_0(0)
         case(6)
            getinittemp= radius*t/K_2+K_3
         case(7)
            getinittemp= (L-Xpoint(1))**2*exp(-t)+K_3
         case(8)
            getinittemp= (L-Xpoint(1))**3*exp(-t)+K_3*Xpoint(1)
         case(9)
!            getinittemp= (L-xcoord)**2*exp(-t)/K_0(0)*q0*W_0(0)
!            getinittemp= (L-xcoord)**2*t/K_0(0)*q0*W_0(0)
         end select
      end function getinittemp
!----------------------------------------------------------------------
!   function name   - get_pennes_mass    (last revision: nov, 2009)
!   purpose         - get mass
!   arguments       - none
!----------------------------------------------------------------------
      function get_pennes_mass()
        use pennes_model , only: RHO,C_P
        implicit none
        PetscScalar :: get_pennes_mass 
        get_pennes_mass = RHO*C_P
      end function get_pennes_mass
!----------------------------------------------------------------------
!   function name   - get_probe_temp    (last revision: nov, 2009)
!   arguments       - none
!----------------------------------------------------------------------
      function get_probe_temp()
        use pennes_model , only: U_PROBE
        implicit none
        PetscScalar :: get_probe_temp 
        get_probe_temp = U_PROBE
      end function get_probe_temp
!----------------------------------------------------------------------
!   program name       - read_control_fortran
!----------------------------------------------------------------------
!   latest revision    - sept 08
!   purpose            - read in control parameters before each 
!                             optimizaiton problem solve
!----------------------------------------------------------------------
      subroutine read_control_fortran(ControlfileID,ImageSpacing)
!     
      use parse_ini
      use pennes_model
      use global_params,only: GroupID,VISUALASE_MAX_POW

      implicit none
#include "cinout.blk"
#include "finclude/petsc.h"
      PetscScalar, intent(in) :: ImageSpacing(3)
      PetscInt,    intent(in) :: ControlfileID
      character(len=MAXLEN) :: inifile   ! INI filename
      character(len=MAXLEN) :: tag   ! 
      ! file containing element wise perfusion/conductivity fields
      character(len=32) :: idcharfile,idgroup ! convert int to string
      logical :: istat,dfstat ! I/O status flag
      PetscInt :: irank,pixel_window
      PetscScalar :: dx_mrti, dy_mrti, dz_mrti
      PetscErrorCode :: ierr
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      ! all processors open files
      write(idcharfile ,*) ControlfileID
      write(idgroup,*) GroupID

      ! called by a compute group
      if(ControlfileID .ge. 0) then
        ! tag the output with groupID and ID of ini file
        tag      ="Group"//trim(adjustl(idgroup))// &
                  "control"//trim(adjustl(idcharfile))//".ini"
        inifile  ="files/control"//trim(adjustl(idcharfile))//".ini"
      else
        ! called by the server use the default file
        tag      ="control.ini"
        inifile = dfltINIFile 
      endif

      call openIniFile(inifile, istat)
      !If unsuccessfully read,stop
      if( .not. istat ) then
        write(*,*) 'could not open', inifile
        stop
      end if
      
      ! the lower bound on the power is set less than zero to allow
      ! a zero value of power other wise the penalty term will blow up
      ! the objective function
      call getIniReal('probe', 'POW_LB', POW_LB,-0.1d0, istat)
      call getIniReal('probe', 'POW_UB', POW_UB,VISUALASE_MAX_POW , istat)
     
      !  laser source term 
      call getIniReal('probe', 'X_0', X_0, 0.0d0, istat)
      call getIniReal('probe', 'Y_0', Y_0, 0.0d0, istat)
      call getIniReal('probe', 'Z_0', Z_0, 0.0d0, istat)
      call getIniReal('optical','MU_A', MU_A, 4.6d0, istat)
        if(MU_A.le.0.0d0) then
           write(*,*) "MU_A must be greater than ZERO!!!"
           call abort
        endif
      call getIniReal('optical', 'MU_S', MU_S, 1474.0d0, istat)
        if(MU_S.le.0.0d0) then
           write(*,*) "MU_S must be greater than ZERO!!!"
           call abort
        endif
      call getIniReal('optical','MU_A_LB', MU_A_LB,0.03d2,istat)
        if(MU_A_LB.le.0.0d0) then
           write(*,*) "MU_A_LB must be greater than ZERO!!!"
           call abort
        endif
      call getIniReal('optical','MU_A_UB',MU_A_UB, 16.0d2, istat)
        if(MU_A_UB.le.0.0d0) then
           write(*,*) "MU_A_UB must be greater than ZERO!!!"
           call abort
        endif
      call getIniReal('optical','MU_S_LB',MU_S_LB,2.0d2,istat)
        if(MU_S_LB.le.0.0d0) then
           write(*,*) "MU_S_LB must be greater than ZERO!!!"
           call abort
        endif
      call getIniReal('optical','MU_S_UB',MU_S_UB,2820.0d2,istat)
        if(MU_S_UB.le.0.0d0) then
           write(*,*) "MU_S_UB must be greater than ZERO!!!"
           call abort
        endif

      ! set the default bounds on the laser position to within 5% of the
      ! field of view 2.5% in the +/- direction
      call getIniInt( 'probe','pixel_window', pixel_window , 12 , istat)

      call getIniReal('probe','X_0_LB', X_0_LB, X_0 - pixel_window*ImageSpacing(1), istat)
      call getIniReal('probe','X_0_UB', X_0_UB, X_0 + pixel_window*ImageSpacing(1), istat)
      call getIniReal('probe','Y_0_LB', Y_0_LB, Y_0 - pixel_window*ImageSpacing(2), istat)
      call getIniReal('probe','Y_0_UB', Y_0_UB, Y_0 + pixel_window*ImageSpacing(2), istat)
      ! set the in-plane/out-plane default bounds on the laser position
      ! to within +/- one mrti slice thickness
      call getIniReal('probe','Z_0_LB', Z_0_LB, Z_0 - ImageSpacing(3), istat)
      call getIniReal('probe','Z_0_UB', Z_0_UB, Z_0 + ImageSpacing(3), istat)

      !read in blood perfusion
      call getIniReal('perfusion','W_N' ,W_N ,0.0d0,dfstat)
      call getIniReal('perfusion','W_I' ,W_I ,0.0d0,dfstat)
      call getIniReal('perfusion','W_D' ,W_D ,0.0d0,dfstat)
      call getIniReal('perfusion','W_2' ,W_2 ,0.0d0,dfstat)
      call getIniReal('perfusion','W_NI',W_NI,0.0d0,dfstat)
      call getIniReal('perfusion','W_ID',W_ID,0.0d0,dfstat)

      ! blood perfusion bounds
      call getIniReal('perfusion','W_0_LB', W_0_LB, 6.0d0 , istat)
      call getIniReal('perfusion','W_0_UB', W_0_UB, 50.0d0, istat)

      call getIniReal('perfusion','W_N_LB', W_N_LB, 4.0d0 , istat)
      call getIniReal('perfusion','W_N_UB', W_N_UB,20.0d0 , istat)

      call getIniReal('perfusion','W_I_LB', W_I_LB,20.0d0 , istat)
      call getIniReal('perfusion','W_I_UB', W_I_UB, 1.0d2 , istat)

      call getIniReal('perfusion','W_D_LB', W_D_LB, 0.0d0 , istat)
      call getIniReal('perfusion','W_D_UB', W_D_UB, 4.0d0 , istat)

      call getIniReal('perfusion','W_2_LB', W_2_LB, 1.0d-1, istat)
      call getIniReal('perfusion','W_2_UB', W_2_UB, 10.0d0, istat)

      call getIniReal('perfusion','W_NID_LB',W_NID_LB,300.0d0,istat)
      call getIniReal('perfusion','W_NID_MD',W_NID_MD,315.0d0,istat)
      call getIniReal('perfusion','W_NID_UB',W_NID_UB,335.0d0,istat)

      !error checking
      if( W_NID_LB .gt. W_NID_MD .or.  W_NID_MD .gt. W_NID_UB) then
          write(*,*) "blood perfusion range error!!!"
          write(*,*) "W_NID_LB = ", W_NID_LB
          write(*,*) "W_NID_MD = ", W_NID_MD
          write(*,*) "W_NID_UB = ", W_NID_UB
          call abort
      endif
      if( W_D_LB .lt. 0.0d0 .or.  W_D_UB .lt. W_N_LB .or.  W_N_UB .lt. W_I_LB) then
          write(*,*) "non physical blood perfusion detected!!!"
          write(*,*) "W_D_LB = ", W_D_LB
          write(*,*) "W_D_UB = ", W_D_UB
          write(*,*) "W_N_LB = ", W_N_LB
          write(*,*) "W_N_UB = ", W_N_UB
          write(*,*) "W_I_LB = ", W_I_LB
          call abort
      endif

      !read in thermal conductivity
      call getIniReal('thermal_conductivity','K_1',K_1,0.0d0,dfstat)
      call getIniReal('thermal_conductivity','K_2',K_2,0.0d0,dfstat)
      call getIniReal('thermal_conductivity','K_3',K_3,0.0d0,dfstat)

      !read in thermal conductivity bounds
      call getIniReal('thermal_conductivity','K_0_LB', K_0_LB, 0.40d0 , istat)
      call getIniReal('thermal_conductivity','K_1_LB', K_1_LB, 0.0d0  , istat)
      call getIniReal('thermal_conductivity','K_2_LB', K_2_LB, 1.0d-2 , istat)
      call getIniReal('thermal_conductivity','K_3_LB', K_3_LB, 305.0d0, istat)

      call getIniReal('thermal_conductivity','K_0_UB', K_0_UB, 0.72d0 , istat)
      call getIniReal('thermal_conductivity','K_1_UB', K_1_UB, 0.33d0 , istat)
      call getIniReal('thermal_conductivity','K_2_UB', K_2_UB, 10.0d0 , istat)
      call getIniReal('thermal_conductivity','K_3_UB', K_3_UB, 325.0d0, istat)

      !error checking
      if ( K_0_LB -  abs(K_1_UB) .lt. 0.0d0 .or. K_0_LB -  abs(K_1_LB) .lt. 0.0d0) then
          write(*,*) "non physical thermal conductivity detected!!!"
          write(*,*) "K_0_LB = ", K_0_LB
          write(*,*) "K_1_LB = ", K_1_LB
          write(*,*) "K_1_UB = ", K_1_UB
          call abort
      endif

      call closeIniFile()

      call mpi_comm_rank(PETSC_COMM_WORLD,irank,ierr)
      if(irank.eq.0) then 
        call printpetscscalar(trim(tag)//': POW_LB       ='//char(0),POW_LB        )
        call printpetscscalar(trim(tag)//': POW_UB       ='//char(0),POW_UB        )
        call printpetscint(   trim(tag)//': pixel_window ='//char(0),pixel_window  )
        call printpetscscalar(trim(tag)//': X_0          ='//char(0),X_0           )
        call printpetscscalar(trim(tag)//': Y_0          ='//char(0),Y_0           )
        call printpetscscalar(trim(tag)//': Z_0          ='//char(0),Z_0           )
        call printpetscscalar(trim(tag)//': X_0_LB       ='//char(0),X_0_LB        )
        call printpetscscalar(trim(tag)//': X_0_UB       ='//char(0),X_0_UB        )
        call printpetscscalar(trim(tag)//': Y_0_LB       ='//char(0),Y_0_LB        )
        call printpetscscalar(trim(tag)//': Y_0_UB       ='//char(0),Y_0_UB        )
        call printpetscscalar(trim(tag)//': Z_0_LB       ='//char(0),Z_0_LB        )
        call printpetscscalar(trim(tag)//': Z_0_UB       ='//char(0),Z_0_UB        )
        call printpetscscalar(trim(tag)//': MU_A         ='//char(0),MU_A          )
        call printpetscscalar(trim(tag)//': MU_A_LB      ='//char(0),MU_A_LB       )
        call printpetscscalar(trim(tag)//': MU_A_UB      ='//char(0),MU_A_UB       )
        call printpetscscalar(trim(tag)//': MU_S         ='//char(0),MU_S          )
        call printpetscscalar(trim(tag)//': MU_S_LB      ='//char(0),MU_S_LB       )
        call printpetscscalar(trim(tag)//': MU_S_UB      ='//char(0),MU_S_UB       )
        call printpetscscalar(trim(tag)//': W_N          ='//char(0),W_N           )
        call printpetscscalar(trim(tag)//': W_I          ='//char(0),W_I           )
        call printpetscscalar(trim(tag)//': W_D          ='//char(0),W_D           )
        call printpetscscalar(trim(tag)//': W_2          ='//char(0),W_2           )
        call printpetscscalar(trim(tag)//': W_NI         ='//char(0),W_NI          )
        call printpetscscalar(trim(tag)//': W_ID         ='//char(0),W_ID          )
        call printpetscscalar(trim(tag)//': W_0_LB       ='//char(0),W_0_LB        )
        call printpetscscalar(trim(tag)//': W_0_UB       ='//char(0),W_0_UB        )
        call printpetscscalar(trim(tag)//': W_N_LB       ='//char(0),W_N_LB        )
        call printpetscscalar(trim(tag)//': W_N_UB       ='//char(0),W_N_UB        )
        call printpetscscalar(trim(tag)//': W_I_LB       ='//char(0),W_I_LB        )
        call printpetscscalar(trim(tag)//': W_I_UB       ='//char(0),W_I_UB        )
        call printpetscscalar(trim(tag)//': W_D_LB       ='//char(0),W_D_LB        )
        call printpetscscalar(trim(tag)//': W_D_UB       ='//char(0),W_D_UB        )
        call printpetscscalar(trim(tag)//': W_2_LB       ='//char(0),W_2_LB        )
        call printpetscscalar(trim(tag)//': W_2_UB       ='//char(0),W_2_UB        )
        call printpetscscalar(trim(tag)//': W_NID_LB     ='//char(0),W_NID_LB      )
        call printpetscscalar(trim(tag)//': W_NID_MD     ='//char(0),W_NID_MD      )
        call printpetscscalar(trim(tag)//': W_NID_UB     ='//char(0),W_NID_UB      )
        call printpetscscalar(trim(tag)//': K_1          ='//char(0),K_1           )
        call printpetscscalar(trim(tag)//': K_2          ='//char(0),K_2           )
        call printpetscscalar(trim(tag)//': K_3          ='//char(0),K_3           )
        call printpetscscalar(trim(tag)//': K_0_LB       ='//char(0),K_0_LB        )
        call printpetscscalar(trim(tag)//': K_0_UB       ='//char(0),K_0_UB        )
        call printpetscscalar(trim(tag)//': K_1_LB       ='//char(0),K_1_LB        )
        call printpetscscalar(trim(tag)//': K_1_UB       ='//char(0),K_1_UB        )
        call printpetscscalar(trim(tag)//': K_2_LB       ='//char(0),K_2_LB        )
        call printpetscscalar(trim(tag)//': K_2_UB       ='//char(0),K_2_UB        )
        call printpetscscalar(trim(tag)//': K_3_LB       ='//char(0),K_3_LB        )
        call printpetscscalar(trim(tag)//': K_3_UB       ='//char(0),K_3_UB        )
      endif

      write(NOUT,*) 'read_control(',ControlfileID,'): X_0='   ,X_0
      write(NOUT,*) 'read_control(',ControlfileID,'): Y_0='   ,Y_0
      write(NOUT,*) 'read_control(',ControlfileID,'): Z_0='   ,Z_0
      write(NOUT,*) 'read_control(',ControlfileID,'): MU_A='  ,MU_A
      write(NOUT,*) 'read_control(',ControlfileID,'): MU_S='  ,MU_S
      write(NOUT,*) 'read_control(',ControlfileID,'): W_N='   ,W_N
      write(NOUT,*) 'read_control(',ControlfileID,'): W_I='   ,W_I
      write(NOUT,*) 'read_control(',ControlfileID,'): W_D='   ,W_D
      write(NOUT,*) 'read_control(',ControlfileID,'): W_2='   ,W_2
      write(NOUT,*) 'read_control(',ControlfileID,'): W_NI='  ,W_NI
      write(NOUT,*) 'read_control(',ControlfileID,'): W_ID='  ,W_ID
      write(NOUT,*) 'read_control(',ControlfileID,'): K_1='   ,K_1
      write(NOUT,*) 'read_control(',ControlfileID,'): K_2='   ,K_2
      write(NOUT,*) 'read_control(',ControlfileID,'): K_3='   ,K_3

      ! Error checking
      if(     K_0_LB  .ge.  K_0_UB  .or.  &
              K_1_LB  .ge.  K_1_UB  .or.  &
              K_2_LB  .ge.  K_2_UB  .or.  &
              K_3_LB  .ge.  K_3_UB  .or.  &
              X_0_LB  .ge.  X_0_UB  .or.  &
              Y_0_LB  .ge.  Y_0_UB  .or.  &
              Z_0_LB  .ge.  Z_0_UB  .or.  &
              POW_LB  .ge.  POW_UB  .or.  &
              MU_A_LB .ge.  MU_A_UB .or.  &
              MU_S_LB .ge.  MU_S_UB .or.  &
              W_0_LB  .ge.  W_0_UB  .or.  &
              W_N_LB  .ge.  W_N_UB  .or.  &
              W_I_LB  .ge.  W_I_UB  .or.  &
              W_D_LB  .ge.  W_D_UB  .or.  &
              W_2_LB  .ge.  W_2_UB) then
          write(*,*) "Invalid Parameter Bounds "
          write(*,*) "must have param_LB < param_UB "
          write(*,*) "param_LB==param_UB gives inf/nan in penalty term"
          call MPI_Barrier(PETSC_COMM_WORLD,ierr)
          call abort
      endif

      end subroutine read_control_fortran
!----------------------------------------------------------------------
!   program name       - init_compute_drone (latest revision - Sept 08)
!
!   purpose            - initialize the objective function specific 
!                        ideal values of compute drone 
!                        default to the arrhenius damage model for now
!----------------------------------------------------------------------
      subroutine init_compute_drone(Maxtime)
        use global_params, only: GroupID
!        use plot_params
        use pennes_model, only: Arr_R,Arr_A,Arr_Ea,U_INIT,X_0,Y_0,Z_0, &
                                CHAR_X_0  , CHAR_Y_0, CHAR_Z_0  , CHAR_RAD, &
                                IDEAL_X_0 ,IDEAL_Y_0 ,IDEAL_Z_0 ,IDEAL_RAD, &
                                IDEAL__IN ,IDEAL_OUT                       
        use parse_ini
      implicit none
#include "cinout.blk"
#include "finclude/petsc.h"
      PetscScalar , intent(in) :: Maxtime
      character(len=32)  :: idchargroup , qoiID
      character(len=MAXLEN),parameter::default_ini='files/control.ini'
      character(len=MAXLEN) :: objective
      PetscInt :: petscrank
      PetscErrorCode :: ierr
      PetscScalar :: dflt__in , dflt_out   

      logical :: istat ! I/O status flag

      call mpi_comm_rank(PETSC_COMM_WORLD,petscrank,ierr)
      ! all processors Try to read the INI file
      call openIniFile(default_ini, istat)
      !If unsuccessfully read,stop
      if( .not. istat ) then
        write(*,*) 'could not open', default_ini
        call abort
      endif

      write(idchargroup,*) GroupID
      qoiID="qoi_"//trim(adjustl(idchargroup))

      ! use laser source as default
      call getIniReal('probe','X_0',X_0, 0.0d0 ,istat)
      call getIniReal('probe','Y_0',Y_0, 0.0d0 ,istat)
      call getIniReal('probe','Z_0',Z_0, 0.0d0 ,istat)

      ! get the characteristic function which gives the domain to 
      ! integrate the qoi. The characteristic funtion is 1.0d0 over a ball of
      ! radius CHAR_RAD centered at (CHAR_X_0,CHAR_Y_0,CHAR_Z_0).
      ! 0.0d0 outside this ball. the DEFAULT IS setting the characteristic 
      ! function to 1.0d0 over THE ENTIRE DOMAIN
      call getIniReal(qoiID,'char_X_0',CHAR_X_0,  X_0  ,istat)
      call getIniReal(qoiID,'char_Y_0',CHAR_Y_0,  Y_0  ,istat)
      call getIniReal(qoiID,'char_Z_0',CHAR_Z_0,  Z_0  ,istat)
      call getIniReal(qoiID,'char_RAD',CHAR_RAD,1.0d+99,istat)
        
      ! get the function which gives the ideal temperature, HSP, or Damage
      ! field. the function equals IDEAL__IN on a ball of radius IDEAL_RAD
      ! centered at (IDEAL_X_0,IDEAL_Y_0,IDEAL_Z_0). outside this ball
      ! the function equals IDEAL__OUT
      call getIniReal(qoiID,'ideal_X_0',IDEAL_X_0,  X_0 ,istat)
      call getIniReal(qoiID,'ideal_Y_0',IDEAL_Y_0,  Y_0 ,istat)
      call getIniReal(qoiID,'ideal_Z_0',IDEAL_Z_0,  Z_0 ,istat)
      call getIniReal(qoiID,'ideal_RAD',IDEAL_RAD,1.0d0 ,istat)

      call getIniString(qoiID,'objective',objective,'calibration',istat)
      select case(trim(objective))
      case("hsp_control") ! hsp based optimization 
          dflt__in = 0.0d0
          dflt_out = 0.0d0
      case("dam_control_arr") !Arrhenius damage based optimization
          dflt__in = Maxtime * Arr_A * exp(-Arr_Ea/Arr_R/400.0d0)   ! kelvin
          dflt_out = Maxtime * Arr_A * exp(-Arr_Ea/Arr_R/U_INIT)    ! kelvin
      case("dam_control_two") !two state model damage based optimization
          dflt__in = 0.0d0
          dflt_out = 0.0d0
      case("calibration") !calibration to MRTI data 
          dflt__in = 0.0d0
          dflt_out = 0.0d0
      case("computeideal_mc") !compute ideal field from Monte Carlo data
          dflt__in = 0.0d0
          dflt_out = 0.0d0
      case("verifprob")
          dflt__in = 0.0d0
          dflt_out = 0.0d0
      case("temp_control")
          dflt__in = 350.d0   ! kelvin
          dflt_out = 290.d0   ! kelvin
      case default
        write(*,*) "unknown objective function" 
        call abort
      end select
      call getIniReal(qoiID,'ideal__IN',IDEAL__IN,dflt__in,istat)
      call getIniReal(qoiID,'ideal_OUT',IDEAL_OUT,dflt_out,istat)


      ! close ini file
      call closeIniFile()

      if(petscrank.eq.0) then 
        call printpetscscalar('init_compute_drone: [field]  CHAR_X_0='//char(0), CHAR_X_0)
        call printpetscscalar('init_compute_drone: [field]  CHAR_Y_0='//char(0), CHAR_Y_0)
        call printpetscscalar('init_compute_drone: [field]  CHAR_Z_0='//char(0), CHAR_Z_0)
        call printpetscscalar('init_compute_drone: [field]  CHAR_RAD='//char(0), CHAR_RAD)
        call printpetscscalar('init_compute_drone: [field] IDEAL_X_0='//char(0),IDEAL_X_0)
        call printpetscscalar('init_compute_drone: [field] IDEAL_Y_0='//char(0),IDEAL_Y_0)
        call printpetscscalar('init_compute_drone: [field] IDEAL_Z_0='//char(0),IDEAL_Z_0)
        call printpetscscalar('init_compute_drone: [field] IDEAL_RAD='//char(0),IDEAL_RAD)
        call printpetscscalar('init_compute_drone: [field] IDEAL__IN='//char(0),IDEAL__IN)
        call printpetscscalar('init_compute_drone: [field] IDEAL_OUT='//char(0),IDEAL_OUT)
      endif
      end subroutine init_compute_drone
!----------------------------------------------------------------------
!   purpose         - get the time steps contribution to the damage field
!                         arrhenius model calculation
!
!   arguments       - IN :  Xpoint - coordinate of the point
!                           Idstep - time step id 
!                     OUT:  ideal field value at the point
!----------------------------------------------------------------------
      function getlocarrdam(Usln,Deltat)
        use pennes_model , ONLY: Arr_A,Arr_Ea,Arr_R
        implicit none
         PetscScalar :: getlocarrdam
         PetscScalar,intent(in):: Usln,Deltat
         getlocarrdam = Arr_A * exp( - Arr_Ea / Arr_R / Usln) * deltat
      end function getlocarrdam
!----------------------------------------------------------------------
!
!   function name   - getidealfield    (last revision: jun, 2007)
!
!   purpose         - get ideal field data at a point
!
!   arguments       - IN :  Xpoint - coordinate of the point
!                     OUT:  ideal field value at the point
!----------------------------------------------------------------------
      function getidealfield(Xpoint)
        use pennes_model, ONLY: X_0, Y_0, Z_0, IDEAL_RAD, IDEAL__IN, IDEAL_OUT 
        implicit none
         PetscScalar :: getidealfield
         PetscScalar,intent(in):: Xpoint(3)
         PetscScalar :: radius

         radius = sqrt( (Xpoint(1)-X_0)**2 +  &
                        (Xpoint(2)-Y_0)**2 +  &
                        (Xpoint(3)-Z_0)**2 )

         getidealfield = IDEAL_OUT
         if(radius .lt. IDEAL_RAD) getidealfield = IDEAL__IN

      end function getidealfield
!----------------------------------------------------------------------
!
!   function name   - getidealfield    (last revision: jun, 2007)
!
!   purpose         - get ideal field data at a point
!
!   arguments       - IN :  Xpoint - coordinate of the point
!                     OUT:  ideal field value at the point
!----------------------------------------------------------------------
      subroutine setfieldid( IdField )
        use pennes_model , only:  ID_W_0,ID_K_0,K_0_FIELD,W_0_FIELD 
        implicit none
#include "finclude/petsc.h"
        PetscInt,intent(in):: IdField 
        if(K_0_FIELD .eqv. PETSC_TRUE) then
           ID_K_0 = IdField 
        else
           ID_K_0 =  0
        endif
        if(W_0_FIELD .eqv. PETSC_TRUE) then
           ID_W_0 = IdField 
        else
           ID_W_0 =  0
        endif
      end subroutine setfieldid
