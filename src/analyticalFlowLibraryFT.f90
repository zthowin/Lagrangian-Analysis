!!---------------------------------------------------------------------
!! Fortran 90/95 style module that contains an analytical flow library.
!!
!!
!! Author:  Zachariah Irwin
!!          University of Colorado Boulder
!! Version: March 2020
!!---------------------------------------------------------------------
MODULE ANALYTICAL_FLOW_LIBRARY_FT

  USE TOOLKIT_MATH

  IMPLICIT NONE

  CONTAINS

    !!--------------------------------------------------------------------------------------------------------
    !> @brief Returns velocity of a given point in steady Double Gyre flow. FLOW_MODEL == 1.
    !!
    !! @param[out]     V_LOC            double array, velocity of point
    !!                                  size: 3
    !! @param[in]      X|Y_IN           double, point coordinates
    !!--------------------------------------------------------------------------------------------------------
    SUBROUTINE DOUBLE_GYRE(X_IN, Y_IN, V_LOC)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN):: X_IN, Y_IN
  
      DOUBLE PRECISION, DIMENSION(3), INTENT(out):: V_LOC

      DOUBLE PRECISION, PARAMETER:: PI = 22.0D0/7.0D0
  
      V_LOC(1) = -PI*sin(PI*X_IN)*cos(PI*Y_IN)
      V_LOC(2) = PI*cos(PI*X_IN)*sin(PI*Y_IN)
      V_LOC(3) = 0.0D0

    END SUBROUTINE DOUBLE_GYRE

    !!---------------------------------------------------------------------------------------
    !> @brief Returns velocity of a given point in unsteady Double Gyre flow, FLOW_MODEL == 2
    !!
    !! @param[out]     V_LOC            double array, velocity of point
    !!                                  size: 3
    !! @param[in]      X|Y_IN           double, point coordinates
    !! @param[in]      T_IN             double, time associated with point
    !!---------------------------------------------------------------------------------------
    SUBROUTINE UNSTEADY_DOUBLE_GYRE(X_IN, Y_IN, T_IN, V_LOC)
      
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN):: X_IN, Y_IN, T_IN
  
      DOUBLE PRECISION, DIMENSION(2), INTENT(out):: V_LOC

      DOUBLE PRECISION, PARAMETER:: PI         = 22.0D0/7.0D0
      DOUBLE PRECISION, PARAMETER:: AMPLITUDE  = 0.10D0
      DOUBLE PRECISION, PARAMETER:: OMEGA      = 0.20D0*PI
      DOUBLE PRECISION, PARAMETER:: EPSILON    = 0.250D0
      
      DOUBLE PRECISION:: a, b, f, dfdx

      a     = EPSILON*sin(OMEGA*T_IN)
      b     = 1.0D0 - 2.0D0*EPSILON*sin(OMEGA*T_IN)
      f     = a*X_IN*X_IN + b*X_IN
      dfdx  = 2.0D0*a*X_IN + b

      V_LOC(1)  = -PI*AMPLITUDE*sin(PI*f)*cos(PI*Y_IN)
      V_LOC(2)  = PI*AMPLITUDE*cos(PI*f)*sin(PI*Y_IN)*dfdx

    END SUBROUTINE UNSTEADY_DOUBLE_GYRE

    !!--------------------------------------------------------------------------------------------------------
    !> @brief Returns velocity of a given point in the unsteady bickley jet, FLOW_MODEL == 3
    !!
    !! See: Haller et al. (2017) "A Critical Comparison of Lagrangian Methods for Coherent Structure Detection"
    !!
    !! @param[out]     V_LOC            double array, velocity of point
    !!                                  size: 3
    !! @param[in]      X|Y_IN           double, point coordinates
    !! @param[in]      T_IN             double, time associated with point
    !!--------------------------------------------------------------------------------------------------------
    SUBROUTINE BICKLEY_JET(X_IN, Y_IN, T_IN, V_LOC)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN):: X_IN, Y_IN, T_IN

      DOUBLE PRECISION, DIMENSION(3), INTENT(out):: V_LOC

      DOUBLE PRECISION, PARAMETER:: U      = 62.66D0
      DOUBLE PRECISION, PARAMETER:: L      = 1770000.0D0
      DOUBLE PRECISION, PARAMETER:: R      = 6371000.0D0
      DOUBLE PRECISION, PARAMETER:: K_N(3) = [2.0D0/R, 4.0D0/R, 6.0D0/R]
      DOUBLE PRECISION, PARAMETER:: E_N(3) = [0.075D0, 0.4D0, 0.3D0]
      DOUBLE PRECISION, PARAMETER:: C_N2   = 0.205D0*U
      DOUBLE PRECISION, PARAMETER:: C_N3   = 0.461D0*U
      DOUBLE PRECISION, PARAMETER:: C_N1   = C_N3+((SQRT(5.0D0)-1.0D0)/2.0D0)*(K_N(2)/K_N(1))*(C_N2-C_N3)

      V_LOC(1) = U*(1/(cosh(Y_IN/L)))**2&
      +2*U*tanh(Y_IN/L)*((1/cosh((Y_IN/L)))**2)*&
      ((E_N(1)*(cos(K_N(1)*C_N1*T_IN)*cos(K_N(1)*X_IN)+sin(K_N(1)*C_N1*T_IN)*sin(K_N(1)*X_IN)))&
      +(E_N(2)*(cos(K_N(2)*C_N2*T_IN)*cos(K_N(2)*X_IN)+sin(K_N(2)*C_N2*T_IN)*sin(K_N(2)*X_IN)))&
      +(E_N(3)*(cos(K_N(3)*C_N3*T_IN)*cos(K_N(3)*X_IN)+sin(K_N(3)*C_N3*T_IN)*sin(K_N(3)*X_IN))))

      V_LOC(2) = U*L*(1/cosh((Y_IN/L))**2)*&
      ((E_N(1)*(K_N(1)*cos(K_N(1)*X_IN)*sin(K_N(1)*C_N1*T_IN)-cos(K_N(1)*C_N1*T_IN)*K_N(1)*sin(K_N(1)*X_IN)))&
      +(E_N(2)*(K_N(2)*cos(K_N(2)*X_IN)*sin(K_N(2)*C_N2*T_IN)-cos(K_N(2)*C_N2*T_IN)*K_N(2)*sin(K_N(2)*X_IN)))&
      +(E_N(3)*(K_N(3)*cos(K_N(3)*X_IN)*sin(K_N(3)*C_N3*T_IN)-cos(K_N(3)*C_N3*T_IN)*K_N(3)*sin(K_N(3)*X_IN))))

      V_LOC(3) = 0.0D0

    END SUBROUTINE BICKLEY_JET

    !!------------------------------------------------------------------------------------
    !> @brief Returns velocity of a given point in the Lamb-Osseen vortex, FLOW_MODEL == 4
    !!
    !! Vortex core structure and layout configurable. Currently set to 3 vortices.
    !!
    !! @param[out]     V_LOC            double array, velocity of point
    !!                                  size: 3
    !! @param[in]      X|Y_IN           double, point coordinates
    !! @param[in]      T_IN             double, time associated with point
    !!-------------------------------------------------------------------------------------
    SUBROUTINE LAMB_OSSEEN(X_IN, Y_IN, T_IN, V_LOC)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN):: X_IN, Y_IN, T_IN

      DOUBLE PRECISION, DIMENSION(3), INTENT(out):: V_LOC

      DOUBLE PRECISION, PARAMETER:: NU     = 0.001D0
      DOUBLE PRECISION, PARAMETER:: ALPHA  = 1.25643D0
      DOUBLE PRECISION, PARAMETER:: U1     = 0.25D0
      DOUBLE PRECISION, PARAMETER:: R0     = 0.29D0/2.0D0
      DOUBLE PRECISION, PARAMETER:: XY1(2) = [R0, R0]
      DOUBLE PRECISION, PARAMETER:: XY2(2) = [-R0, R0]
      DOUBLE PRECISION, PARAMETER:: XY3(2) = [0.0D0, R0-SQRT((R0*2.0D0)**2-R0**2)]

      DOUBLE PRECISION:: RC, THETA_DOT_1, THETA_DOT_2, THETA_DOT_3
      DOUBLE PRECISION, DIMENSION(2):: RT_0, RT_1, RT_2, RT_3
      
      !!--------------------------------------
      !! Time dependent vs. fixed core radius.
      !!--------------------------------------
      !RC = SQRT(4.0D0*NU*abs(T_IN))
      RC = 0.005D0

      !!----------------------------------------------------------
      !! Get reference coordinates of point to all three vortices.
      !!----------------------------------------------------------
      CALL CARTESIAN_TO_POLAR(X_IN-XY1(1), Y_IN-XY1(2), RT_1)
      CALL CARTESIAN_TO_POLAR(X_IN-XY2(1), Y_IN-XY2(2), RT_2)
      CALL CARTESIAN_TO_POLAR(X_IN-XY3(1), Y_IN-XY3(2), RT_3)

      !!------------------------------
      !! Compute tangential velocities
      !!------------------------------
      THETA_DOT_1 = U1*(1.0D0 + &
      (0.50D0/ALPHA))*(RC/RT_1(1))*(1.0D0-exp(-ALPHA*(RT_1(1)**2/RC**2)))

      THETA_DOT_2 = U1*(1.0D0 + & 
      (0.50D0/ALPHA))*(RC/RT_2(1))*(1.0D0-exp(-ALPHA*(RT_2(1)**2/RC**2)))

      THETA_DOT_3 = U1*(1.0D0 + &
      (0.50D0/ALPHA))*(RC/RT_3(1))*(1.0D0-exp(-ALPHA*(RT_3(1)**2/RC**2)))

      !!-----------------------------
      !! Compute cartesian velocities
      !!-----------------------------
      V_LOC(1) = -(THETA_DOT_1*cos(RT_1(2)) + THETA_DOT_2*cos(RT_2(2)) + THETA_DOT_3*cos(RT_3(2)))
      V_LOC(2) = THETA_DOT_1*sin(RT_1(2)) + THETA_DOT_2*sin(RT_2(2)) + THETA_DOT_3*sin(RT_3(2))
      V_LOC(3) = 0.0D0
      
    END SUBROUTINE LAMB_OSSEEN

    !!-----------------------------------------------------------------------------
    !> @brief Returns velocity of a given point in steady ABC flow, FLOW_MODEL == 5
    !!
    !! Steady Arnold-Beltrami-Childress flow (direct solutions to Euler equations).
    !!
    !! @param[out]     V_LOC            double array, velocity of point
    !!                                  size: 3
    !! @param[in]      X|Y|Z_IN         double, point coordinates
    !!------------------------------------------------------------------------------
    SUBROUTINE ABC(X_IN, Y_IN, Z_IN, V_LOC)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN):: X_IN, Y_IN, Z_IN

      DOUBLE PRECISION, DIMENSION(3), INTENT(out):: V_LOC

      DOUBLE PRECISION, PARAMETER:: a_A = SQRT(3.0D0)
      DOUBLE PRECISION, PARAMETER:: a_B = SQRT(2.0D0)
      DOUBLE PRECISION, PARAMETER:: a_C = 1.0D0

      V_LOC(1)  = a_A*sin(Z_IN) + a_C*cos(Y_IN)
      V_LOC(2)  = a_B*sin(X_IN) + a_A*cos(Z_IN)
      V_LOC(3)  = a_C*sin(Y_IN) + a_B*cos(X_IN)

    END SUBROUTINE ABC

    !!-------------------------------------------------------------------------------
    !> @brief Returns velocity of a given point in unsteady ABC flow, FLOW_MODEL == 6
    !!
    !! Unteady Arnold-Beltrami-Childress flow (direct solutions to Euler equations).
    !!
    !! @param[out]     V_LOC            double array, velocity of point
    !!                                  size(3)
    !! @param[in]      X|Y|Z_IN         double, point coordinates
    !! @param[in]      T_IN             double, time associated with point
    !!--------------------------------------------------------------------------------
    SUBROUTINE UNSTEADY_ABC(X_IN, Y_IN, Z_IN, T_IN, V_LOC)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN):: X_IN, Y_IN, Z_IN, T_IN

      DOUBLE PRECISION, DIMENSION(3), INTENT(out):: V_LOC

      DOUBLE PRECISION, PARAMETER:: PI  = 22.0D0/7.0D0
      DOUBLE PRECISION, PARAMETER:: a_A = SQRT(3.0D0)
      DOUBLE PRECISION, PARAMETER:: a_B = SQRT(2.0D0)
      DOUBLE PRECISION, PARAMETER:: a_C = 1.0D0

      V_LOC(1) = (a_A + 0.50D0*T_IN*sin(PI*T_IN))*sin(Z_IN) + a_C*cos(Y_IN)
      V_LOC(2) = a_B*sin(X_IN) + (a_A + 0.50D0*T_IN*sin(PI*T_IN))*cos(Z_IN)
      V_LOC(3) = a_C*sin(Y_IN) + a_B*cos(X_IN)

    END SUBROUTINE UNSTEADY_ABC

END MODULE ANALYTICAL_FLOW_LIBRARY_FT
