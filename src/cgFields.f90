!!----------------------------------------------------------
!! Fortran 90/95 style module to perform Lagrangian analysis
!! on particles seeded on a Cartesian grid.
!!
!!
!! Author:   Zachariah Irwin
!!           University of Colorado, Boulder
!! Version:  March 2020
!!----------------------------------------------------------
MODULE CGFIELDS

  USE TOOLKIT_MATH

  IMPLICIT NONE

  CONTAINS

    !!--------------------------------------------------------------------------------------------------------
    !> @brief Calculates Lagrangian fields on cartesian grids
    !!
    !! LAGRANGIAN_MAIN handles all computational subroutines for every point in a cartesian grid by assessing
    !! whether or not the point has left the domain during the integration time step from T_N to T_N+1 (T_NP1)
    !! and, if it has, calculating the Right Cauchy-Green deformation tensor for that point, and its immediate
    !! neighbors, based on the grid at time T_N. If the point has not left the domain and the simulation has
    !! reached the end of the integration window, the Right Cauchy-Green deformation tensor is calculated based
    !! on the grid at T_N+1. From the input file, the user specifies what Lagrangian fields (currently FTLE, 
    !! elongational stretch, and longitudinal strain) are to be calculated.
    !!
    !! @param[in, out]   HAS_CG         integer array, flags for tracers that have a Cauchy Green value
    !!                                  size: NUM_P
    !! @param[in, out]   FTLE_T         double array, time-scaled FTLE
    !!                                  size: NUM_P
    !! @param[in, out]   FTLE_NO_T      double array, non-time-scaled FTLE
    !!                                  size: NUM_P
    !! @param[in, out]   STRETCH_1|2|3  double array, elongational stretch along unit vectors 1|2|3
    !!                                  size: NUM_P
    !! @param[in, out]   STRAIN_1|2|3   double array, Green-Lagrange strain along unit vectors 1|2|3
    !!                                  size: NUM_P
    !! @param[in, out]   CAUCHY_GREEN   double array, Right Cauchy-Green matrix
    !!                                  size: NUM_P x 3 x 3
    !! @param[in]        DIMN           integer, problem dimension (2D or 3D)
    !! @param[in]        NUM_P          integer, number of tracers
    !! @param[in]        REF_NUMS       integer, number of tracers in each dimension
    !!                                  size: 3
    !! @param[in]        FLAGS          integer, flags for which Lagrangian fields to compute
    !!                                  size: 3
    !! @param[in]        LEFT_GRID      integer array, flags for tracers which have left the flow domain
    !!                                  size: NUM_X x NUM_Y x NUM_Z
    !! @param[in]        REF_DIMNS      double array, grid spacing in referential configuration
    !!                                  size: 3
    !! @param[in]        TIMES          double array, FTLE window start/stop and simulation times
    !!                                  size: 4
    !! @param[in]        UNIT_VECTORS   double array, unit vectors along user-specified directions 1,2,3
    !!                                  size: 3x3
    !!--------------------------------------------------------------------------------------------------------
    SUBROUTINE LAGRANGIAN_MAIN(DIMN, NUM_P, REF_NUMS, REF_DIMNS, TIMES, PARTICLES_XYZ_TN, PARTICLES_XYZ_TNP1,&
                               FLAGS, HAS_CG, LEFT_GRID, UNIT_VECTORS, CAUCHY_GREEN, FTLE_T, FTLE_NO_T,& 
                               STRETCH_1, STRETCH_2, STRETCH_3, STRAIN_1, STRAIN_2, STRAIN_3)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: DIMN, NUM_P
      INTEGER, DIMENSION(3), INTENT(IN):: FLAGS, REF_NUMS
      INTEGER, DIMENSION(0:NUM_P-1), INTENT(IN):: LEFT_GRID
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: REF_DIMNS
      DOUBLE PRECISION, DIMENSION(4), INTENT(IN):: TIMES
      DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN):: UNIT_VECTORS
      DOUBLE PRECISION, DIMENSION(0:NUM_P-1, 3), INTENT(IN):: PARTICLES_XYZ_TN, PARTICLES_XYZ_TNP1
      
      INTEGER, DIMENSION(0:NUM_P-1), INTENT(INOUT):: HAS_CG
      DOUBLE PRECISION, DIMENSION(0:NUM_P-1), INTENT(INOUT):: FTLE_T, FTLE_NO_T, STRAIN_1, STRAIN_2, STRAIN_3,&
                                                              STRETCH_1, STRETCH_2, STRETCH_3
      DOUBLE PRECISION, DIMENSION(0:NUM_P-1, 3, 3), INTENT(INOUT):: CAUCHY_GREEN

      LOGICAL:: ONE, TWO, THREE
      INTEGER:: I, J, K, P, NUM_X, NUM_Y, NUM_Z
      DOUBLE PRECISION:: T_START, T_END, T_N, T_NP1
      DOUBLE PRECISION, DIMENSION(3):: UNIT_VEC_1, UNIT_VEC_2, UNIT_VEC_3

      DOUBLE PRECISION, PARAMETER:: EPS = 1.0e-12

      !f2py depend(NUM_P) PARTICLES_XYZ_TN, PARTICLES_XYZ_TNP1, HAS_CG, LEFT_GRID, CAUCHY_GREEN, FTLE_T, FTLE_NO_T, STRAIN_1, STRAIN_2, STRAIN_3, STRETCH_1, STRETCH_2, STRETCH_3

      !!--------------------------------
      !! Initialize reference dimensions
      !!--------------------------------
      NUM_X = REF_NUMS(1)
      NUM_Y = REF_NUMS(2)
      NUM_Z = REF_NUMS(3)

      !!--------------------------------------------
      !! Initialize FTLE window and simulation times
      !!--------------------------------------------
      T_START = TIMES(1)
      T_END   = TIMES(2)
      T_N     = TIMES(3)
      T_NP1   = TIMES(4)

      UNIT_VEC_1 = UNIT_VECTORS(1,:)
      UNIT_VEC_2 = UNIT_VECTORS(2,:)
      UNIT_VEC_3 = UNIT_VECTORS(3,:)

      !!----------------------------------------------------
      !! Flags to not compute along a unit vector direction
      !! if it is the 0 vector
      !!----------------------------------------------------
      IF (ALL(UNIT_VEC_1 .EQ. 0.0D0)) THEN
        ONE = .FALSE.
      ELSE
        ONE = .TRUE.
      END IF
      IF (ALL(UNIT_VEC_2 .EQ. 0.0D0)) THEN
        TWO = .FALSE.
      ELSE
        TWO = .TRUE.
      END IF
      IF (ALL(UNIT_VEC_3 .EQ. 0.0D0)) THEN
        THREE = .FALSE.
      ELSE
        THREE = .TRUE.
      END IF

      !!---------------------
      !! Loop over all points
      !!---------------------
      DO K = 0, NUM_Z-1
        DO J = 0, NUM_Y-1
          DO I = 0, NUM_X-1

            CALL LINEAR_INDEXING([I,J,K], DIMN, REF_NUMS, P)

            IF (LEFT_GRID(P) == 1) THEN

              !!---------------------------------------------------------------------------------------------
              !! Point has left the domain, calculate the CG matrix for the point and its immediate neighbors
              !!---------------------------------------------------------------------------------------------
              IF (HAS_CG(P) == -1) THEN

                CALL GET_CG_FOR_POINT(DIMN, NUM_P, I, J, K, REF_DIMNS, REF_NUMS, PARTICLES_XYZ_TN,&
                                      HAS_CG, CAUCHY_GREEN)

                IF (FLAGS(1) == 1) THEN
                  CALL GET_FTLE_FOR_POINT(NUM_P, P, CAUCHY_GREEN, T_START, T_NP1, FTLE_T, FTLE_NO_T)
                END IF

                IF (FLAGS(2) == 1) THEN
                  IF (ONE) THEN
                    CALL GET_STRETCH_FOR_POINT(NUM_P, P, UNIT_VEC_1, CAUCHY_GREEN, STRETCH_1)
                  END IF
                  IF (TWO) THEN
                    CALL GET_STRETCH_FOR_POINT(NUM_P, P, UNIT_VEC_2, CAUCHY_GREEN, STRETCH_2)
                  END IF
                  IF (THREE) THEN
                    CALL GET_STRETCH_FOR_POINT(NUM_P, P, UNIT_VEC_3, CAUCHY_GREEN, STRETCH_3)
                  END IF
                END IF

                IF (FLAGS(3) == 1 .AND. FLAGS(2) == 1) THEN
                  IF (ONE) THEN
                    STRAIN_1(P) = STRETCH_1(P) - 1.0D0
                  END IF
                  IF (TWO) THEN
                    STRAIN_2(P) = STRETCH_2(P) - 1.0D0
                  END IF
                  IF (THREE) THEN
                    STRAIN_3(P) = STRETCH_3(P) - 1.0D0
                  END IF
                ELSE IF (FLAGS(3) == 1 .AND. FLAGS(2) .NE. 1) THEN
                  IF (ONE) THEN
                    CALL GET_STRAIN_FOR_POINT(NUM_P, P, UNIT_VEC_1, CAUCHY_GREEN, STRAIN_1)
                  END IF
                  IF (TWO) THEN
                    CALL GET_STRAIN_FOR_POINT(NUM_P, P, UNIT_VEC_2, CAUCHY_GREEN, STRAIN_2)
                  END IF
                  IF (THREE) THEN
                    CALL GET_STRAIN_FOR_POINT(NUM_P, P, UNIT_VEC_3, CAUCHY_GREEN, STRAIN_3)
                  END IF
                END IF
              END IF

              !!-------------------
              !! Left neighbor (-x)
              !!-------------------
              IF (I > 0) THEN
                CALL LINEAR_INDEXING([I-1,J,K], DIMN, REF_NUMS, P)
                IF (HAS_CG(P) == -1) THEN
                  CALL GET_CG_FOR_POINT(DIMN, NUM_P, I-1, J, K, REF_DIMNS, REF_NUMS, PARTICLES_XYZ_TN,&
                                        HAS_CG, CAUCHY_GREEN)

                  IF (FLAGS(1) == 1) THEN
                    CALL GET_FTLE_FOR_POINT(NUM_P, P, CAUCHY_GREEN, T_START, T_NP1, FTLE_T, FTLE_NO_T)
                  END IF

                  IF (FLAGS(2) == 1) THEN
                    IF (ONE) THEN
                      CALL GET_STRETCH_FOR_POINT(NUM_P, P, UNIT_VEC_1, CAUCHY_GREEN, STRETCH_1)
                    END IF
                    IF (TWO) THEN
                      CALL GET_STRETCH_FOR_POINT(NUM_P, P, UNIT_VEC_2, CAUCHY_GREEN, STRETCH_2)
                    END IF
                    IF (THREE) THEN
                      CALL GET_STRETCH_FOR_POINT(NUM_P, P, UNIT_VEC_3, CAUCHY_GREEN, STRETCH_3)
                    END IF
                  END IF

                  IF (FLAGS(3) == 1 .AND. FLAGS(2) == 1) THEN
                    IF (ONE) THEN
                      STRAIN_1(P) = STRETCH_1(P) - 1.0D0
                    END IF
                    IF (TWO) THEN
                      STRAIN_2(P) = STRETCH_2(P) - 1.0D0
                    END IF
                    IF (THREE) THEN
                      STRAIN_3(P) = STRETCH_3(P) - 1.0D0
                    END IF
                  ELSE IF (FLAGS(3) == 1 .AND. FLAGS(2) .NE. 1) THEN
                    IF (ONE) THEN
                      CALL GET_STRAIN_FOR_POINT(NUM_P, P, UNIT_VEC_1, CAUCHY_GREEN, STRAIN_1)
                    END IF
                    IF (TWO) THEN
                      CALL GET_STRAIN_FOR_POINT(NUM_P, P, UNIT_VEC_2, CAUCHY_GREEN, STRAIN_2)
                    END IF
                    IF (THREE) THEN
                      CALL GET_STRAIN_FOR_POINT(NUM_P, P, UNIT_VEC_3, CAUCHY_GREEN, STRAIN_3)
                    END IF
                  END IF
                END IF
              END IF

              !!--------------------
              !! Right neighbor (+x)
              !!--------------------
              IF (I < NUM_X-1) THEN
                CALL LINEAR_INDEXING([I+1,J,K], DIMN, REF_NUMS, P)
                IF (HAS_CG(P) == -1) THEN
                  CALL GET_CG_FOR_POINT(DIMN, NUM_P, I+1, J, K, REF_DIMNS, REF_NUMS, PARTICLES_XYZ_TN,&
                                        HAS_CG, CAUCHY_GREEN)

                  IF (FLAGS(1) == 1) THEN
                    CALL GET_FTLE_FOR_POINT(NUM_P, P, CAUCHY_GREEN, T_START, T_NP1, FTLE_T, FTLE_NO_T)
                  END IF

                  IF (FLAGS(2) == 1) THEN
                    IF (ONE) THEN
                      CALL GET_STRETCH_FOR_POINT(NUM_P, P, UNIT_VEC_1, CAUCHY_GREEN, STRETCH_1)
                    END IF
                    IF (TWO) THEN
                      CALL GET_STRETCH_FOR_POINT(NUM_P, P, UNIT_VEC_2, CAUCHY_GREEN, STRETCH_2)
                    END IF
                    IF (THREE) THEN
                      CALL GET_STRETCH_FOR_POINT(NUM_P, P, UNIT_VEC_3, CAUCHY_GREEN, STRETCH_3)
                    END IF
                  END IF

                  IF (FLAGS(3) == 1 .AND. FLAGS(2) == 1) THEN
                    IF (ONE) THEN
                      STRAIN_1(P) = STRETCH_1(P) - 1.0D0
                    END IF
                    IF (TWO) THEN
                      STRAIN_2(P) = STRETCH_2(P) - 1.0D0
                    END IF
                    IF (THREE) THEN
                      STRAIN_3(P) = STRETCH_3(P) - 1.0D0
                    END IF
                  ELSE IF (FLAGS(3) == 1 .AND. FLAGS(2) .NE. 1) THEN
                    IF (ONE) THEN
                      CALL GET_STRAIN_FOR_POINT(NUM_P, P, UNIT_VEC_1, CAUCHY_GREEN, STRAIN_1)
                    END IF
                    IF (TWO) THEN
                      CALL GET_STRAIN_FOR_POINT(NUM_P, P, UNIT_VEC_2, CAUCHY_GREEN, STRAIN_2)
                    END IF
                    IF (THREE) THEN
                      CALL GET_STRAIN_FOR_POINT(NUM_P, P, UNIT_VEC_3, CAUCHY_GREEN, STRAIN_3)
                    END IF
                  END IF
                END IF
              END IF

              !!---------------------
              !! Bottom neighbor (-y)
              !!---------------------
              IF (J > 0) THEN
                CALL LINEAR_INDEXING([I,J-1,K], DIMN, REF_NUMS, P)
                IF (HAS_CG(P) == -1) THEN
                  CALL GET_CG_FOR_POINT(DIMN, NUM_P, I, J-1, K, REF_DIMNS, REF_NUMS, PARTICLES_XYZ_TN,&
                                        HAS_CG, CAUCHY_GREEN)

                  IF (FLAGS(1) == 1) THEN
                    CALL GET_FTLE_FOR_POINT(NUM_P, P, CAUCHY_GREEN, T_START, T_NP1, FTLE_T, FTLE_NO_T)
                  END IF

                  IF (FLAGS(2) == 1) THEN
                    IF (ONE) THEN
                      CALL GET_STRETCH_FOR_POINT(NUM_P, P, UNIT_VEC_1, CAUCHY_GREEN, STRETCH_1)
                    END IF
                    IF (TWO) THEN
                      CALL GET_STRETCH_FOR_POINT(NUM_P, P, UNIT_VEC_2, CAUCHY_GREEN, STRETCH_2)
                    END IF
                    IF (THREE) THEN
                      CALL GET_STRETCH_FOR_POINT(NUM_P, P, UNIT_VEC_3, CAUCHY_GREEN, STRETCH_3)
                    END IF
                  END IF

                  IF (FLAGS(3) == 1 .AND. FLAGS(2) == 1) THEN
                    IF (ONE) THEN
                      STRAIN_1(P) = STRETCH_1(P) - 1.0D0
                    END IF
                    IF (TWO) THEN
                      STRAIN_2(P) = STRETCH_2(P) - 1.0D0
                    END IF
                    IF (THREE) THEN
                      STRAIN_3(P) = STRETCH_3(P) - 1.0D0
                    END IF
                  ELSE IF (FLAGS(3) == 1 .AND. FLAGS(2) .NE. 1) THEN
                    IF (ONE) THEN
                      CALL GET_STRAIN_FOR_POINT(NUM_P, P, UNIT_VEC_1, CAUCHY_GREEN, STRAIN_1)
                    END IF
                    IF (TWO) THEN
                      CALL GET_STRAIN_FOR_POINT(NUM_P, P, UNIT_VEC_2, CAUCHY_GREEN, STRAIN_2)
                    END IF
                    IF (THREE) THEN
                      CALL GET_STRAIN_FOR_POINT(NUM_P, P, UNIT_VEC_3, CAUCHY_GREEN, STRAIN_3)
                    END IF
                  END IF
                END IF
              END IF

              !!------------------
              !! Top neighbor (+y)
              !!------------------
              IF (J < NUM_Y-1) THEN
                CALL LINEAR_INDEXING([I,J+1,K], DIMN, REF_NUMS, P)
                IF (HAS_CG(P) == -1) THEN
                  CALL GET_CG_FOR_POINT(DIMN, NUM_P, I, J+1, K, REF_DIMNS, REF_NUMS, PARTICLES_XYZ_TN,&
                                        HAS_CG, CAUCHY_GREEN)

                  IF (FLAGS(1) == 1) THEN
                    CALL GET_FTLE_FOR_POINT(NUM_P, P, CAUCHY_GREEN, T_START, T_NP1, FTLE_T, FTLE_NO_T)
                  END IF

                  IF (FLAGS(2) == 1) THEN
                    IF (ONE) THEN
                      CALL GET_STRETCH_FOR_POINT(NUM_P, P, UNIT_VEC_1, CAUCHY_GREEN, STRETCH_1)
                    END IF
                    IF (TWO) THEN
                      CALL GET_STRETCH_FOR_POINT(NUM_P, P, UNIT_VEC_2, CAUCHY_GREEN, STRETCH_2)
                    END IF
                    IF (THREE) THEN
                      CALL GET_STRETCH_FOR_POINT(NUM_P, P, UNIT_VEC_3, CAUCHY_GREEN, STRETCH_3)
                    END IF
                  END IF

                  IF (FLAGS(3) == 1 .AND. FLAGS(2) == 1) THEN
                    IF (ONE) THEN
                      STRAIN_1(P) = STRETCH_1(P) - 1.0D0
                    END IF
                    IF (TWO) THEN
                      STRAIN_2(P) = STRETCH_2(P) - 1.0D0
                    END IF
                    IF (THREE) THEN
                      STRAIN_3(P) = STRETCH_3(P) - 1.0D0
                    END IF
                  ELSE IF (FLAGS(3) == 1 .AND. FLAGS(2) .NE. 1) THEN
                    IF (ONE) THEN
                      CALL GET_STRAIN_FOR_POINT(NUM_P, P, UNIT_VEC_1, CAUCHY_GREEN, STRAIN_1)
                    END IF
                    IF (TWO) THEN
                      CALL GET_STRAIN_FOR_POINT(NUM_P, P, UNIT_VEC_2, CAUCHY_GREEN, STRAIN_2)
                    END IF
                    IF (THREE) THEN
                      CALL GET_STRAIN_FOR_POINT(NUM_P, P, UNIT_VEC_3, CAUCHY_GREEN, STRAIN_3)
                    END IF
                  END IF
                END IF
              END IF

              !!-------------------
              !! Back neighbor (-z)
              !!-------------------
              IF (K > 0) THEN
                CALL LINEAR_INDEXING([I,J,K-1], DIMN, REF_NUMS, P)
                IF (HAS_CG(P) == -1) THEN
                  CALL GET_CG_FOR_POINT(DIMN, NUM_P, I, J, K-1, REF_DIMNS, REF_NUMS,&
                                        PARTICLES_XYZ_TN, HAS_CG, CAUCHY_GREEN)

                  IF (FLAGS(1) == 1) THEN
                    CALL GET_FTLE_FOR_POINT(NUM_P, P, CAUCHY_GREEN, T_START, T_NP1, FTLE_T, FTLE_NO_T)
                  END IF

                  IF (FLAGS(2) == 1) THEN
                    IF (ONE) THEN
                      CALL GET_STRETCH_FOR_POINT(NUM_P, P, UNIT_VEC_1, CAUCHY_GREEN, STRETCH_1)
                    END IF
                    IF (TWO) THEN
                      CALL GET_STRETCH_FOR_POINT(NUM_P, P, UNIT_VEC_2, CAUCHY_GREEN, STRETCH_2)
                    END IF
                    IF (THREE) THEN
                      CALL GET_STRETCH_FOR_POINT(NUM_P, P, UNIT_VEC_3, CAUCHY_GREEN, STRETCH_3)
                    END IF
                  END IF

                  IF (FLAGS(3) == 1 .AND. FLAGS(2) == 1) THEN
                    IF (ONE) THEN
                      STRAIN_1(P) = STRETCH_1(P) - 1.0D0
                    END IF
                    IF (TWO) THEN
                      STRAIN_2(P) = STRETCH_2(P) - 1.0D0
                    END IF
                    IF (THREE) THEN
                      STRAIN_3(P) = STRETCH_3(P) - 1.0D0
                    END IF
                  ELSE IF (FLAGS(3) == 1 .AND. FLAGS(2) .NE. 1) THEN
                    IF (ONE) THEN
                      CALL GET_STRAIN_FOR_POINT(NUM_P, P, UNIT_VEC_1, CAUCHY_GREEN, STRAIN_1)
                    END IF
                    IF (TWO) THEN
                      CALL GET_STRAIN_FOR_POINT(NUM_P, P, UNIT_VEC_2, CAUCHY_GREEN, STRAIN_2)
                    END IF
                    IF (THREE) THEN
                      CALL GET_STRAIN_FOR_POINT(NUM_P, P, UNIT_VEC_3, CAUCHY_GREEN, STRAIN_3)
                    END IF
                  END IF
                END IF
              END IF

              !!--------------------
              !! Front neighbor (+z)
              !!--------------------
              IF (K < NUM_Z-1) THEN
                CALL LINEAR_INDEXING([I,J,K+1], DIMN, REF_NUMS, P)
                IF (HAS_CG(P) == -1) THEN
                  CALL GET_CG_FOR_POINT(DIMN, NUM_P, I, J, K+1, REF_DIMNS, REF_NUMS,&
                                        PARTICLES_XYZ_TN, HAS_CG, CAUCHY_GREEN)

                  IF (FLAGS(1) == 1) THEN
                    CALL GET_FTLE_FOR_POINT(NUM_P, P, CAUCHY_GREEN, T_START, T_NP1, FTLE_T, FTLE_NO_T)
                  END IF

                  IF (FLAGS(2) == 1) THEN
                    IF (ONE) THEN
                      CALL GET_STRETCH_FOR_POINT(NUM_P, P, UNIT_VEC_1, CAUCHY_GREEN, STRETCH_1)
                    END IF
                    IF (TWO) THEN
                      CALL GET_STRETCH_FOR_POINT(NUM_P, P, UNIT_VEC_2, CAUCHY_GREEN, STRETCH_2)
                    END IF
                    IF (THREE) THEN
                      CALL GET_STRETCH_FOR_POINT(NUM_P, P, UNIT_VEC_3, CAUCHY_GREEN, STRETCH_3)
                    END IF
                  END IF

                  IF (FLAGS(3) == 1 .AND. FLAGS(2) == 1) THEN
                    IF (ONE) THEN
                      STRAIN_1(P) = STRETCH_1(P) - 1.0D0
                    END IF
                    IF (TWO) THEN
                      STRAIN_2(P) = STRETCH_2(P) - 1.0D0
                    END IF
                    IF (THREE) THEN
                      STRAIN_3(P) = STRETCH_3(P) - 1.0D0
                    END IF
                  ELSE IF (FLAGS(3) == 1 .AND. FLAGS(2) .NE. 1) THEN
                    IF (ONE) THEN
                      CALL GET_STRAIN_FOR_POINT(NUM_P, P, UNIT_VEC_1, CAUCHY_GREEN, STRAIN_1)
                    END IF
                    IF (TWO) THEN
                      CALL GET_STRAIN_FOR_POINT(NUM_P, P, UNIT_VEC_2, CAUCHY_GREEN, STRAIN_2)
                    END IF
                    IF (THREE) THEN
                      CALL GET_STRAIN_FOR_POINT(NUM_P, P, UNIT_VEC_3, CAUCHY_GREEN, STRAIN_3)
                    END IF
                  END IF
                END IF
              END IF

            !!-----------------------------------------------------------------------------
            !! Point still in domain, only calculate Lagrangian fields if at the end of the
            !! Lagrangian window (big T integration window length)
            !!-----------------------------------------------------------------------------
            ELSE IF (LEFT_GRID(P) == -1 .AND. (ABS(T_NP1) >= ABS(T_END) - EPS)) THEN

              IF (HAS_CG(P) == -1) THEN
                CALL GET_CG_FOR_POINT(DIMN, NUM_P, I, J, K, REF_DIMNS, REF_NUMS, PARTICLES_XYZ_TN,&
                                      HAS_CG, CAUCHY_GREEN)

                IF (FLAGS(1) == 1) THEN
                  CALL GET_FTLE_FOR_POINT(NUM_P, P, CAUCHY_GREEN, T_START, T_END, FTLE_T, FTLE_NO_T)
                END IF

                IF (FLAGS(2) == 1) THEN
                    IF (ONE) THEN
                      CALL GET_STRETCH_FOR_POINT(NUM_P, P, UNIT_VEC_1, CAUCHY_GREEN, STRETCH_1)
                    END IF
                    IF (TWO) THEN
                      CALL GET_STRETCH_FOR_POINT(NUM_P, P, UNIT_VEC_2, CAUCHY_GREEN, STRETCH_2)
                    END IF
                    IF (THREE) THEN
                      CALL GET_STRETCH_FOR_POINT(NUM_P, P, UNIT_VEC_3, CAUCHY_GREEN, STRETCH_3)
                    END IF
                  END IF

                IF (FLAGS(3) == 1 .AND. FLAGS(2) == 1) THEN
                  IF (ONE) THEN
                    STRAIN_1(P) = STRETCH_1(P) - 1.0D0
                  END IF
                  IF (TWO) THEN
                    STRAIN_2(P) = STRETCH_2(P) - 1.0D0
                  END IF
                  IF (THREE) THEN
                    STRAIN_3(P) = STRETCH_3(P) - 1.0D0
                  END IF
                ELSE IF (FLAGS(3) == 1 .AND. FLAGS(2) .NE. 1) THEN
                  IF (ONE) THEN
                    CALL GET_STRAIN_FOR_POINT(NUM_P, P, UNIT_VEC_1, CAUCHY_GREEN, STRAIN_1)
                  END IF
                  IF (TWO) THEN
                    CALL GET_STRAIN_FOR_POINT(NUM_P, P, UNIT_VEC_2, CAUCHY_GREEN, STRAIN_2)
                  END IF
                  IF (THREE) THEN
                    CALL GET_STRAIN_FOR_POINT(NUM_P, P, UNIT_VEC_3, CAUCHY_GREEN, STRAIN_3)
                  END IF
                END IF
              END IF

            END IF

          END DO
        END DO
      END DO

    END SUBROUTINE LAGRANGIAN_MAIN

    !!--------------------------------------------------------------------------------------------------------
    !> @brief Calculates the Right Cauchy-Green tensor for a given point
    !!
    !! First, the deformation gradient F is calculated for a given point using numerical differentiation.
    !! Then F(Transpose)F is calculated for the given point to give the Cauchy-Green tensor at that point.
    !!
    !! @param[in, out]   HAS_CG         integer array, flags for points that have a Cauchy Green value
    !!                                  size: NUM_X x NUM_Y x NUM_Z
    !! @param[in, out]   CAUCHY_GREEN   double array, Right Cauchy-Green matrix
    !!                                  size: NUM_X x NUM_Y x NUM_Z x 3 x 3
    !! @param[in]        DIMN           integer, problem dimension (2D or 3D)
    !! @param[in]        I|J|K          integer, point's index in a given grid
    !! @param[in]        NUM_P          integer, number of total points in the simulation
    !!                                  (equivalent to REF_NUMS(1)*REF_NUMS(2)*REF_NUMS(3))
    !! @param[in]        REF_NUMS       integer array, number of points along each cartesian direction
    !!                                  size: 3
    !! @param[in]        REF_DIMNS      double array, grid spacing in referential configuration
    !!                                  size: 3
    !! @param[in]        PARTICLES_XYZ  double array, grid coordinates
    !!                                  size: NUM_P x 3
    !!--------------------------------------------------------------------------------------------------------
    SUBROUTINE GET_CG_FOR_POINT(DIMN, NUM_P, I, J, K, REF_DIMNS, REF_NUMS, PARTICLES_XYZ, HAS_CG,&
                                CAUCHY_GREEN)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: DIMN, NUM_P, I, J, K
      INTEGER, DIMENSION(3), INTENT(IN):: REF_NUMS
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: REF_DIMNS
      DOUBLE PRECISION, DIMENSION(0:NUM_P-1, 3), INTENT(IN):: PARTICLES_XYZ
      
      INTEGER, DIMENSION(0:NUM_P-1), INTENT(INOUT):: HAS_CG
      DOUBLE PRECISION, DIMENSION(0:NUM_P-1, 3, 3), INTENT(INOUT):: CAUCHY_GREEN

      INTEGER:: II, JJ, KK, P, P_I_MINUS, P_I_PLUS, P_J_MINUS, P_J_PLUS, P_K_MINUS, P_K_PLUS,&
                NUM_X, NUM_Y, NUM_Z
      DOUBLE PRECISION:: X_I_MINUS, Y_I_MINUS, Z_I_MINUS, X_I_PLUS, Y_I_PLUS, Z_I_PLUS,&
                         X_J_MINUS, Y_J_MINUS, Z_J_MINUS, X_J_PLUS, Y_J_PLUS, Z_J_PLUS,&
                         X_K_MINUS, Y_K_MINUS, Z_K_MINUS, X_K_PLUS, Y_K_PLUS, Z_K_PLUS,&
                         D_DX, D_DY, D_DZ, REF_DX, REF_DY, REF_DZ
      DOUBLE PRECISION, DIMENSION(3,3):: DEF_GRAD, PT_CG

      DOUBLE PRECISION, PARAMETER:: EPS = 1.0E-12

      !f2py depend(NUM_P) PARTICLES_XYZ, HAS_CG, CAUCHY_GREEN

      CALL LINEAR_INDEXING([I,J,K], DIMN, REF_NUMS, P)

      PT_CG = CAUCHY_GREEN(P, :, :)

      NUM_X = REF_NUMS(1)
      NUM_Y = REF_NUMS(2)
      NUM_Z = REF_NUMS(3)
      
      REF_DX = REF_DIMNS(1)
      REF_DY = REF_DIMNS(2)
      REF_DZ = REF_DIMNS(3)

      !!------------------------------
      !! Get IDs of neighboring points
      !!------------------------------
      IF (I < NUM_X-1) THEN
        CALL LINEAR_INDEXING([I+1,J,K], DIMN, REF_NUMS, P_I_PLUS)
      END IF
      IF (I > 0) THEN
        CALL LINEAR_INDEXING([I-1,J,K], DIMN, REF_NUMS, P_I_MINUS)
      END IF
      IF (J < NUM_Y-1) THEN
        CALL LINEAR_INDEXING([I,J+1,K], DIMN, REF_NUMS, P_J_PLUS)
      END IF
      IF (J > 0) THEN
        CALL LINEAR_INDEXING([I,J-1,K], DIMN, REF_NUMS, P_J_MINUS)
      END IF
      
      IF (DIMN == 3) THEN
        IF (K < NUM_Z-1) THEN
          CALL LINEAR_INDEXING([I,J,K+1], DIMN, REF_NUMS, P_K_PLUS)
        END IF
        IF (K > 0) THEN
          CALL LINEAR_INDEXING([I,J,K-1], DIMN, REF_NUMS, P_K_MINUS) 
        END IF
      END IF

      !!--------------------------------------
      !! Store entries of deformation gradient
      !!--------------------------------------
      IF (I == 0) THEN
        X_I_MINUS = PARTICLES_XYZ(P, 1)
        Y_I_MINUS = PARTICLES_XYZ(P, 2)
        Z_I_MINUS = PARTICLES_XYZ(P, 3)
        X_I_PLUS  = PARTICLES_XYZ(P_I_PLUS, 1)
        Y_I_PLUS  = PARTICLES_XYZ(P_I_PLUS, 2)
        Z_I_PLUS  = PARTICLES_XYZ(P_I_PLUS, 3)
        D_DX      = 1.0D0

      ELSE IF (I == NUM_X-1) THEN
        X_I_MINUS = PARTICLES_XYZ(P_I_MINUS, 1)
        Y_I_MINUS = PARTICLES_XYZ(P_I_MINUS, 2)
        Z_I_MINUS = PARTICLES_XYZ(P_I_MINUS, 3)
        X_I_PLUS  = PARTICLES_XYZ(P, 1)
        Y_I_PLUS  = PARTICLES_XYZ(P, 2)
        Z_I_PLUS  = PARTICLES_XYZ(P, 3)
        D_DX      = 1.0D0

      ELSE
        X_I_MINUS = PARTICLES_XYZ(P_I_MINUS, 1)
        Y_I_MINUS = PARTICLES_XYZ(P_I_MINUS, 2)
        Z_I_MINUS = PARTICLES_XYZ(P_I_MINUS, 3)
        X_I_PLUS  = PARTICLES_XYZ(P_I_PLUS, 1)
        Y_I_PLUS  = PARTICLES_XYZ(P_I_PLUS, 2)
        Z_I_PLUS  = PARTICLES_XYZ(P_I_PLUS, 3)
        D_DX      = 2.0D0

      END IF

      IF (J == 0) THEN
        X_J_MINUS = PARTICLES_XYZ(P, 1)
        Y_J_MINUS = PARTICLES_XYZ(P, 2)
        Z_J_MINUS = PARTICLES_XYZ(P, 3)
        X_J_PLUS  = PARTICLES_XYZ(P_J_PLUS, 1)
        Y_J_PLUS  = PARTICLES_XYZ(P_J_PLUS, 2)
        Z_J_PLUS  = PARTICLES_XYZ(P_J_PLUS, 3)
        D_DY      = 1.0D0

      ELSE IF (J == NUM_Y-1) THEN
        X_J_MINUS = PARTICLES_XYZ(P_J_MINUS, 1)
        Y_J_MINUS = PARTICLES_XYZ(P_J_MINUS, 2)
        Z_J_MINUS = PARTICLES_XYZ(P_J_MINUS, 3)
        X_J_PLUS  = PARTICLES_XYZ(P, 1)
        Y_J_PLUS  = PARTICLES_XYZ(P, 2)
        Z_J_PLUS  = PARTICLES_XYZ(P, 3)
        D_DY      = 1.0D0

      ELSE
        X_J_MINUS = PARTICLES_XYZ(P_J_MINUS, 1)
        Y_J_MINUS = PARTICLES_XYZ(P_J_MINUS, 2)
        Z_J_MINUS = PARTICLES_XYZ(P_J_MINUS, 3)
        X_J_PLUS  = PARTICLES_XYZ(P_J_PLUS, 1)
        Y_J_PLUS  = PARTICLES_XYZ(P_J_PLUS, 2)
        Z_J_PLUS  = PARTICLES_XYZ(P_J_PLUS, 3)
        D_DY      = 2.0D0

      END IF

      IF (DIMN == 3) THEN
        IF (K == 0) THEN
          X_K_MINUS = PARTICLES_XYZ(P, 1)
          Y_K_MINUS = PARTICLES_XYZ(P, 2)
          Z_K_MINUS = PARTICLES_XYZ(P, 3)
          X_K_PLUS  = PARTICLES_XYZ(P_K_PLUS, 1)
          Y_K_PLUS  = PARTICLES_XYZ(P_K_PLUS, 2)
          Z_K_PLUS  = PARTICLES_XYZ(P_K_PLUS, 3)
          D_DZ      = 1.0D0

        ELSE IF (K == NUM_Z-1) THEN
          X_K_MINUS = PARTICLES_XYZ(P_K_MINUS, 1)
          Y_K_MINUS = PARTICLES_XYZ(P_K_MINUS, 2)
          Z_K_MINUS = PARTICLES_XYZ(P_K_MINUS, 3)
          X_K_PLUS  = PARTICLES_XYZ(P, 1)
          Y_K_PLUS  = PARTICLES_XYZ(P, 2)
          Z_K_PLUS  = PARTICLES_XYZ(P, 3)
          D_DZ      = 1.0D0

        ELSE
          X_K_MINUS = PARTICLES_XYZ(P_K_MINUS, 1)
          Y_K_MINUS = PARTICLES_XYZ(P_K_MINUS, 2)
          Z_K_MINUS = PARTICLES_XYZ(P_K_MINUS, 3)
          X_K_PLUS  = PARTICLES_XYZ(P_K_PLUS, 1)
          Y_K_PLUS  = PARTICLES_XYZ(P_K_PLUS, 2)
          Z_K_PLUS  = PARTICLES_XYZ(P_K_PLUS, 3)
          D_DZ      = 2.0D0

        END IF

      !!-----------------------------
      !! Compute deformation gradient
      !!-----------------------------
        DEF_GRAD(1,3) = (X_K_PLUS - X_K_MINUS)/(D_DZ*REF_DZ)
        DEF_GRAD(2,3) = (Y_K_PLUS - Y_K_MINUS)/(D_DZ*REF_DZ)
        DEF_GRAD(3,1) = (Z_I_PLUS - Z_I_MINUS)/(D_DX*REF_DX)
        DEF_GRAD(3,2) = (Z_J_PLUS - Z_J_MINUS)/(D_DY*REF_DY)
        DEF_GRAD(3,3) = (Z_K_PLUS - Z_K_MINUS)/(D_DZ*REF_DZ)

      ELSE
        DEF_GRAD(1,3) = 0.0D0
        DEF_GRAD(2,3) = 0.0D0
        DEF_GRAD(3,1) = 0.0D0
        DEF_GRAD(3,2) = 0.0D0
        DEF_GRAD(3,3) = 0.0D0

      END IF

      DEF_GRAD(1,1) = (X_I_PLUS - X_I_MINUS)/(D_DX*REF_DX)
      DEF_GRAD(1,2) = (X_J_PLUS - X_J_MINUS)/(D_DY*REF_DY)
      DEF_GRAD(2,1) = (Y_I_PLUS - Y_I_MINUS)/(D_DX*REF_DX)
      DEF_GRAD(2,2) = (Y_J_PLUS - Y_J_MINUS)/(D_DY*REF_DY)

      !!--------------------------------------
      !! Build Cauchy Green deformation tensor
      !!--------------------------------------
      PT_CG = MATMUL(TRANSPOSE(DEF_GRAD), DEF_GRAD)

      CAUCHY_GREEN(P, :, :)  = PT_CG
      HAS_CG(P)              = 1

    END SUBROUTINE GET_CG_FOR_POINT

    !!--------------------------------------------------------------------------------------------------------
    !> @brief Calculates the FTLE for a given point.
    !!
    !! This subroutine makes usage of LAPACK's DGEEV to find the eigenvalues of the Cauchy-Green tensor.
    !! The Cauchy-Green tensor is passed into DGEEV which outputs two vectors containing the real and
    !! imaginary values of the Cauchy-Green tensor. The FTLE is then computed by taking the log of the
    !! square root of the maximum real eigenvalue. Documentation for DGEEV can be found at:
    !! http://www.netlib.org/lapack/explore-html/d9/d8e/group__double_g_eeigen_ga66e19253344358f5dee1e60502b9e96f.html
    !!
    !! @param[in, out]   CAUCHY_GREEN   double array, Right Cauchy-Green matrix
    !!                                  size: NUM_X x NUM_Y x NUM_Z x 3 x 3
    !! @param[in, out]   FTLE_T|NO_T    double array, FTLE time scaled|not time scaled
    !!                                  size: NUM_X x NUM_Y x NUM_Z
    !! @param[in]        NUM_P          integer, number of points in the simulation
    !! @param[in]        P              integer, point ID of interest
    !! @param[in]        T0             double, t_n    of integration
    !! @param[in]        T1             double, t_n+1  of integration
    !!--------------------------------------------------------------------------------------------------------
    SUBROUTINE GET_FTLE_FOR_POINT(NUM_P, P, CAUCHY_GREEN, T0, T1, FTLE_T, FTLE_NO_T)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: NUM_P, P
      DOUBLE PRECISION, INTENT(IN):: T0, T1
      
      DOUBLE PRECISION, DIMENSION(0:NUM_P-1, 3, 3), INTENT(IN):: CAUCHY_GREEN
      DOUBLE PRECISION, DIMENSION(0:NUM_P-1), INTENT(INOUT):: FTLE_T, FTLE_NO_T

      INTEGER:: INFO, LWORK, N, LDVL, LDVR
      DOUBLE PRECISION:: EIGENVALUE_MAX
      DOUBLE PRECISION, DIMENSION(3):: EIGENVALUE_VEC_REAL, EIGENVALUE_VEC_IMAG
      DOUBLE PRECISION, DIMENSION(1, 3):: VL, VR
      DOUBLE PRECISION, DIMENSION(MAX(1, 15)):: WORK
      DOUBLE PRECISION, DIMENSION(3, 3):: PT_CG

      DOUBLE PRECISION, PARAMETER:: EPS = 1.0e-12

      !f2py depend(NUM_P) CAUCHY_GREEN, FTLE_T, FTLE_NO_T

      !!----------------------------------------------------------------
      !! Make temporary Cauchy-Green tensor; DGEEV will modify the input
      !!----------------------------------------------------------------
      PT_CG = CAUCHY_GREEN(P, :, :)

      N     = 3
      LWORK = 5*3
      LDVL  = N
      LDVR  = N

      CALL DGEEV('N', 'N', N, PT_CG, N, EIGENVALUE_VEC_REAL, EIGENVALUE_VEC_IMAG, VL, LDVL,&
                 VR, LDVR, WORK, LWORK, INFO)

      IF (INFO .NE. 0) THEN
        print*, "DGEEV FAILED, cgFields.F90 LINE 791"
        print*, INFO
        STOP
      END IF

      EIGENVALUE_MAX = MAX(EIGENVALUE_VEC_REAL(1), EIGENVALUE_VEC_REAL(2), EIGENVALUE_VEC_REAL(3))

      ! CALL GET_MAX_EIGENVALUE(PT_CG, 3, EIGENVALUE_MAX)

      IF (EIGENVALUE_MAX > 1.0D0 + EPS) THEN
        FTLE_T(P)    = 0.50D0*LOG(EIGENVALUE_MAX)/ABS(T1 - T0)
        FTLE_NO_T(P) = 0.50D0*LOG(EIGENVALUE_MAX)

      ELSE
        FTLE_T(P)    = 0.0D0
        FTLE_NO_T(P) = 0.0D0

      END IF

    END SUBROUTINE GET_FTLE_FOR_POINT

    !!--------------------------------------------------------------------------------------------------------
    !> @brief Calculates the elongational stretch of a given point.
    !!
    !! The elongational stretch is the square root of the dot product of the unit vector of interest
    !! transposed with the Cauchy-Green tensor of a given point and dotted with the unit vector again.
    !!
    !! @param[in, out]   CAUCHY_GREEN   double array, Right Cauchy-Green matrix
    !!                                  size: NUM_X x NUM_Y x NUM_Z x 3 x 3
    !! @param[in, out]   STRETCH        double array, elongational stretch along a unit vector
    !!                                  size: NUM_P
    !! @param[in]        NUM_P          integer, number of points in the simulation
    !! @param[in]        P              integer, point ID of interest
    !! @param[in]        UNIT_VEC       double array, unit vector along user-specified direction
    !!                                  size: 3
    !!--------------------------------------------------------------------------------------------------------
    SUBROUTINE GET_STRETCH_FOR_POINT(NUM_P, P, UNIT_VEC, CAUCHY_GREEN, STRETCH)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: NUM_P, p
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: UNIT_VEC
      DOUBLE PRECISION, DIMENSION(0:NUM_P-1, 3, 3), INTENT(IN):: CAUCHY_GREEN
      DOUBLE PRECISION, DIMENSION(0:NUM_P-1), INTENT(INOUT):: STRETCH

      INTEGER:: II, JJ
      DOUBLE PRECISION:: PT_STRETCH
      DOUBLE PRECISION, DIMENSION(3,3):: PT_CG

      !f2py depend(NUM_P) CAUCHY_GREEN, STRETCH

      PT_STRETCH  = 0.0D0
      PT_CG       = CAUCHY_GREEN(P,:,:)

      DO II = 1, 3
        DO JJ = 1, 3
          PT_STRETCH = PT_STRETCH + SQRT(UNIT_VEC(JJ)*PT_CG(JJ, II)*UNIT_VEC(II))
        END DO
      END DO

      STRETCH(P) = PT_STRETCH

    END SUBROUTINE GET_STRETCH_FOR_POINT

    !!--------------------------------------------------------------------------------------------------------
    !> @brief Calculates the Green-Lagrange strain of a given point.
    !!
    !! The Green-Lagrange Strain is the square root of the dot product of the unit vector of interest
    !! transposed with the Cauchy-Green tensor of a given point and dotted with the unit vector again
    !! minus the identity.
    !!
    !! @param[in, out]   CAUCHY_GREEN   double array, Right Cauchy-Green matrix
    !!                                  size: NUM_X x NUM_Y x NUM_Z x 3 x 3
    !! @param[in, out]   STRAIN         double array, Green-Lagrange strain along a unit vector
    !!                                  size: NUM_X x NUM_Y x NUM_Z
    !! @param[in]        NUM_P          integer, number of points in the simulation
    !! @param[in]        P              integer, point ID of interest
    !! @param[in]        UNIT_VEC       double array, unit vector along user-specified direction
    !!                                  size: 3
    !!--------------------------------------------------------------------------------------------------------
    SUBROUTINE GET_STRAIN_FOR_POINT(NUM_P, P, UNIT_VEC, CAUCHY_GREEN, STRAIN)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: NUM_P, P
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: UNIT_VEC
      DOUBLE PRECISION, DIMENSION(0:NUM_P-1, 3, 3), INTENT(IN):: CAUCHY_GREEN
      DOUBLE PRECISION, DIMENSION(0:NUM_P-1), INTENT(INOUT):: STRAIN

      INTEGER:: II, JJ
      DOUBLE PRECISION:: PT_STRAIN
      DOUBLE PRECISION, DIMENSION(3,3):: PT_CG

      !f2py depend(NUM_P) CAUCHY_GREEN, STRAIN

      PT_STRAIN = 0.0D0
      PT_CG     = CAUCHY_GREEN(P,:,:)

      DO II = 1, 3
        DO JJ = 1, 3
          PT_STRAIN = PT_STRAIN + SQRT(UNIT_VEC(JJ)*PT_CG(JJ, II)*UNIT_VEC(II))
        END DO
      END DO

      STRAIN(P) = PT_STRAIN - 1.0D0

    END SUBROUTINE GET_STRAIN_FOR_POINT

END MODULE CGFIELDS
