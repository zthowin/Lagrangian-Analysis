!!-----------------------------------------------------------------------
!! Fortran 90/95 style module that contains a boundary condition library.
!!
!!
!! Author:  Zachariah Irwin
!!          University of Colorado Boulder
!! Version: November 2019
!!-----------------------------------------------------------------------
MODULE BOUNDARY_CONDITIONS

    IMPLICIT NONE

    CONTAINS

    !!--------------------------------------------------------------------------------------------------------
    !> @brief Decide which boundary condition subroutine to call.
    !!
    !! @param[out]        INT_STATUS    integer, flag corresponding to whether or not point left the domain
    !! @param[out]        X|Y|Z_OUT     double, point coordinates after domain check
    !! @param[in]         BC_TYPE       integer, flag corresponding to boundary condition
    !! @param[in]         X|Y|Z_IN      double, point coordinates at t_n
    !! @param[in]         X|Y|Z_TEMP    double, point coordinates at t_n+1
    !! @param[in]         BOUNDING_BOX  double array, flow domain bounding box
    !!                                  size: 6
    !!--------------------------------------------------------------------------------------------------------
    SUBROUTINE RESOLVE_BOUNDARY_CONDITIONS(BC_TYPE, BOUNDING_BOX, X_IN, Y_IN, Z_IN, X_TEMP, Y_TEMP, Z_TEMP,&
                                           X_OUT, Y_OUT, Z_OUT, INT_STATUS)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: BC_TYPE
      DOUBLE PRECISION, DIMENSION(6), INTENT(IN):: BOUNDING_BOX
      DOUBLE PRECISION, INTENT(IN):: X_IN, Y_IN, Z_IN, X_TEMP, Y_TEMP, Z_TEMP

      INTEGER, INTENT(OUT):: INT_STATUS
      DOUBLE PRECISION, INTENT(OUT):: X_OUT, Y_OUT, Z_OUT

      INT_STATUS  = -1
      
      IF (BC_TYPE == 0) THEN
        CALL RESET(BOUNDING_BOX, X_IN, Y_IN, Z_IN, X_TEMP, Y_TEMP, Z_TEMP, X_OUT, Y_OUT, Z_OUT, INT_STATUS)
      
      ELSE IF (BC_TYPE == 1) THEN
        CALL PUSHBACK(BOUNDING_BOX, X_IN, Y_IN, X_TEMP, Y_TEMP, X_OUT, Y_OUT)

      ELSE IF (BC_TYPE == 2) THEN
        CALL REFLECT(BOUNDING_BOX, X_TEMP, Y_TEMP, Z_TEMP, X_OUT, Y_OUT, Z_OUT)

      ELSE IF (BC_TYPE == 3) THEN
        CALL PERIODIC_X(BOUNDING_BOX, X_TEMP, Y_TEMP, X_OUT, Y_OUT)

      ELSE IF (BC_TYPE == 4) THEN
        CALL FIX(BOUNDING_BOX, X_TEMP, Y_TEMP, Z_TEMP, X_OUT, Y_OUT, Z_OUT, INT_STATUS)

      END IF

    END SUBROUTINE RESOLVE_BOUNDARY_CONDITIONS

    !!--------------------------------------------------------------------------------------------------------
    !> @brief If the tracer crosses the boundary, reset it to its previous position at t_n. BCType == 0.
    !!
    !! @param[in,out]     INT_STATUS    integer, flag corresponding to whether or not point left the domain
    !! @param[out]        X|Y|Z_OUT     double, point coordinates after domain check
    !! @param[in]         X|Y|Z_IN      double, point coordinates at t_n
    !! @param[in]         X|Y|Z_TEMP    double, point coordinates at t_n+1
    !! @param[in]         BOUNDING_BOX  double array, flow domain bounding box
    !!                                  size: 6
    !!--------------------------------------------------------------------------------------------------------
    SUBROUTINE RESET(BOUNDING_BOX, X_IN, Y_IN, Z_IN, X_TEMP, Y_TEMP, Z_TEMP, X_OUT, Y_OUT, Z_OUT, INT_STATUS)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN):: X_IN, Y_IN, Z_IN, X_TEMP, Y_TEMP, Z_TEMP
      DOUBLE PRECISION, DIMENSION(6), INTENT(IN):: BOUNDING_BOX

      INTEGER, INTENT(INOUT):: INT_STATUS
      DOUBLE PRECISION, INTENT(OUT):: X_OUT, Y_OUT, Z_OUT

      DOUBLE PRECISION:: EPS = 1.0E-10

      X_OUT    = X_TEMP
      Y_OUT    = Y_TEMP
      Z_OUT    = Z_TEMP

      IF ((X_TEMP > BOUNDING_BOX(2) + EPS) .OR. (X_TEMP < BOUNDING_BOX(1) - EPS) &
          .OR. (Y_TEMP > BOUNDING_BOX(4) + EPS) .OR. (Y_TEMP < BOUNDING_BOX(3) - EPS)&
          .OR. (Z_TEMP > BOUNDING_BOX(6) + EPS) .OR. (Z_TEMP < BOUNDING_BOX(5) - EPS)) THEN

        X_OUT       = X_IN
        Y_OUT       = Y_IN
        Z_OUT       = Z_IN
        INT_STATUS  = 1

      END IF

    END SUBROUTINE RESET

    !!--------------------------------------------------------------------------------------------------------
    !> @brief If the tracer crosses the boundary, push it back into the domain along a linear trajectory. 
    !! BCType == 1. **Currently only for 2D flows**
    !!
    !! @param[in,out]     INT_STATUS    integer, flag corresponding to whether or not point left the domain
    !! @param[out]        X|Y_OUT       double, point coordinates after domain check
    !! @param[in]         X|Y_IN        double, point coordinates at t_n
    !! @param[in]         X|Y_TEMP      double, point coordinates at t_n+1
    !! @param[in]         BOUNDING_BOX  double array, flow domain bounding box
    !!                                  size: 6
    !!--------------------------------------------------------------------------------------------------------
    SUBROUTINE PUSHBACK(BOUNDING_BOX, X_IN, Y_IN, X_TEMP, Y_TEMP, X_OUT, Y_OUT)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN):: X_IN, Y_IN, X_TEMP, Y_TEMP
      DOUBLE PRECISION, DIMENSION(6), INTENT(IN):: BOUNDING_BOX

      DOUBLE PRECISION, INTENT(OUT):: X_OUT, Y_OUT

      DOUBLE PRECISION, PARAMETER:: EPS = 1.0e-14

      DOUBLE PRECISION:: COUNTER

      X_OUT = X_TEMP
      Y_OUT = Y_TEMP
      
      !!----------
      !! X > X_max
      !!----------
      IF (X_TEMP > BOUNDING_BOX(2) .AND. (Y_TEMP < BOUNDING_BOX(4) .AND. Y_TEMP > BOUNDING_BOX(3))) THEN 
        X_OUT = BOUNDING_BOX(2) - EPS
        Y_OUT = ((X_OUT-X_IN)*(Y_TEMP-Y_IN))/(X_TEMP-X_IN) + Y_IN
        RETURN

      ELSE IF (X_TEMP > BOUNDING_BOX(2) .AND. Y_TEMP > BOUNDING_BOX(4)) THEN
        COUNTER = 1.0D0
        DO WHILE (Y_OUT > BOUNDING_BOX(4))
          X_OUT = BOUNDING_BOX(2) - EPS*COUNTER
          Y_OUT = (X_OUT-X_IN)*((Y_TEMP-Y_IN)/(X_TEMP-X_IN)) + Y_IN
          COUNTER = COUNTER + 1.0D0
        END DO
        RETURN

      ELSE IF (X_TEMP > BOUNDING_BOX(2) .AND. Y_TEMP < BOUNDING_BOX(3)) THEN
        COUNTER = 1.0D0
        DO WHILE (Y_OUT < BOUNDING_BOX(3))
          X_OUT = BOUNDING_BOX(2) - EPS*COUNTER
          Y_OUT = (X_OUT-X_IN)*((Y_TEMP-Y_IN)/(X_TEMP-X_IN)) + Y_IN
          COUNTER = COUNTER + 1.0D0
        END DO
        RETURN

      END IF

      !!----------
      !! X < X_min
      !!----------
      IF (X_TEMP < BOUNDING_BOX(1) .AND. (Y_TEMP < BOUNDING_BOX(4) .AND. Y_TEMP > BOUNDING_BOX(3))) THEN
        X_OUT = BOUNDING_BOX(1) + EPS
        Y_OUT = ((X_OUT-X_IN)*(Y_TEMP-Y_IN))/(X_TEMP-X_IN) + Y_IN
        RETURN

      ELSE IF (X_TEMP < BOUNDING_BOX(1) .AND. Y_TEMP > BOUNDING_BOX(4)) THEN
        COUNTER = 1.0D0
        DO WHILE (Y_OUT > BOUNDING_BOX(4))
          X_OUT = BOUNDING_BOX(1) + EPS*COUNTER
          Y_OUT = (X_OUT-X_IN)*((Y_TEMP-Y_IN)/(X_TEMP-X_IN)) + Y_IN
          COUNTER = COUNTER + 1.0D0
        END DO
        RETURN
          
      ELSE IF (X_TEMP < BOUNDING_BOX(1) .AND. Y_TEMP < BOUNDING_BOX(3)) THEN
        COUNTER = 1.0D0
        DO WHILE (Y_OUT < BOUNDING_BOX(3))
          X_OUT = BOUNDING_BOX(1) + EPS*COUNTER
          Y_OUT = (X_OUT-X_IN)*((Y_TEMP-Y_IN)/(X_TEMP-X_IN)) + Y_IN
          COUNTER = COUNTER + 1.0D0
        END DO
        RETURN

      END IF

      !!----------
      !! Y > Y_max
      !!----------
      IF (Y_TEMP > BOUNDING_BOX(4) .AND. (X_TEMP < BOUNDING_BOX(2) .AND. X_TEMP > BOUNDING_BOX(1))) THEN
        Y_OUT = BOUNDING_BOX(4) - EPS
        X_OUT = ((Y_OUT-Y_IN)*(X_TEMP-X_IN))/(Y_TEMP-Y_IN) + X_IN
        RETURN
      END IF

      !!----------
      !! Y < Y_min
      !!----------
      IF (Y_TEMP < BOUNDING_BOX(3) .AND. (X_TEMP < BOUNDING_BOX(2) .AND. X_TEMP > BOUNDING_BOX(1))) THEN
        Y_OUT = BOUNDING_BOX(3) + EPS
        X_OUT = ((Y_OUT-Y_IN)*(X_TEMP-X_IN))/(Y_TEMP-Y_IN) + X_IN
        RETURN
      END IF

    END SUBROUTINE PUSHBACK 

    !!--------------------------------------------------------------------------------------------------------
    !> @brief If the tracer crosses the boundary, impose an elastic collision. BCType == 2.
    !!
    !! @param[out]        X|Y|Z_OUT     double, point coordinates after domain check
    !! @param[in]         X|Y|Z_TEMP    double, point coordinates at t_n+1
    !! @param[in]         BOUNDING_BOX  double array, flow domain bounding box
    !!                                  size: 6
    !!--------------------------------------------------------------------------------------------------------
    SUBROUTINE REFLECT(BOUNDING_BOX, X_TEMP, Y_TEMP, Z_TEMP, X_OUT, Y_OUT, Z_OUT)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN):: X_TEMP, Y_TEMP, Z_TEMP
      DOUBLE PRECISION, DIMENSION(6), INTENT(IN):: BOUNDING_BOX

      DOUBLE PRECISION, INTENT(OUT):: X_OUT, Y_OUT, Z_OUT

      X_OUT = X_TEMP
      Y_OUT = Y_TEMP
      Z_OUT = Z_TEMP

      IF (X_TEMP > BOUNDING_BOX(2)) THEN
        X_OUT = X_TEMP - 2.0D0*abs(X_TEMP - BOUNDING_BOX(2))
      END IF

      IF (X_TEMP < BOUNDING_BOX(1)) THEN
        X_OUT = X_TEMP + 2.0D0*abs(BOUNDING_BOX(1) - X_TEMP)
      END IF

      IF (Y_TEMP > BOUNDING_BOX(4)) THEN
        Y_OUT = Y_TEMP - 2.0D0*abs(Y_TEMP - BOUNDING_BOX(4))
      END IF

      IF (Y_TEMP < BOUNDING_BOX(3)) THEN
        Y_OUT = Y_TEMP + 2.0D0*abs(BOUNDING_BOX(3) - Y_TEMP)
      END IF

      IF (Z_TEMP > BOUNDING_BOX(6)) THEN
        Z_OUT = Z_TEMP - 2.0D0*abs(Z_TEMP - BOUNDING_BOX(6))
      END IF

      IF (Z_TEMP < BOUNDING_BOX(5)) THEN
        Z_OUT = Z_TEMP + 2.0D0*abs(BOUNDING_BOX(5) - Z_TEMP)
      END IF

    END SUBROUTINE REFLECT

    !!--------------------------------------------------------------------------------------------------------
    !> @brief Periodic boundary condition for +x. BCType == 3. **Currently only for 2D flows**
    !!
    !! @param[out]        X|Y_OUT       double, point coordinates after domain check
    !! @param[in]         X|Y_OUT       double, point coordinates at t_n+1
    !! @param[in]         BOUNDING_BOX  double array, flow domain bounding box
    !!                                  size: 6
    !!--------------------------------------------------------------------------------------------------------
    SUBROUTINE PERIODIC_X(BOUNDING_BOX, X_TEMP, Y_TEMP, X_OUT, Y_OUT)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN):: X_TEMP, Y_TEMP
      DOUBLE PRECISION, DIMENSION(6), INTENT(IN):: BOUNDING_BOX

      DOUBLE PRECISION, INTENT(OUT):: X_OUT, Y_OUT

      DOUBLE PRECISION:: COUNTER

      X_OUT = X_TEMP
      Y_OUT = Y_TEMP

      IF (X_TEMP > BOUNDING_BOX(2)) THEN
        COUNTER = 1.0
        DO WHILE (X_OUT > BOUNDING_BOX(2))
          X_OUT = X_TEMP - COUNTER*BOUNDING_BOX(2)
          COUNTER = COUNTER + 1.0
        END DO

      ELSE IF (X_TEMP < BOUNDING_BOX(1)) THEN
        COUNTER = 1.0
        DO WHILE (X_OUT < BOUNDING_BOX(1))
          X_OUT = COUNTER*BOUNDING_BOX(2) + X_TEMP
          COUNTER = COUNTER + 1.0
        END DO

      END IF

      ! Reflective boundary condition for y-domain crossings
      IF (Y_TEMP > BOUNDING_BOX(4)) THEN
        Y_OUT = Y_TEMP - 2.0D0*abs(Y_TEMP - BOUNDING_BOX(4))

      ELSE IF (Y_TEMP < BOUNDING_BOX(3)) THEN
        Y_OUT = Y_TEMP + 2.0D0*abs(BOUNDING_BOX(3) - Y_TEMP)

      END IF

    END SUBROUTINE PERIODIC_X

    !!--------------------------------------------------------------------------------------------------------
    !> @brief If the tracer crosses the boundary, fix it on the boundary. BCType == 4.
    !!
    !! @param[in,out]     INT_STATUS    integer, flag corresponding to whether or not point left the domain
    !! @param[out]        X|Y|Z_OUT     double, point coordinates after domain check
    !! @param[in]         X|Y|Z_TEMP    double, point coordinates at t_n+1
    !! @param[in]         BOUNDING_BOX  double array, flow domain bounding box
    !!                                  size: 6
    !!--------------------------------------------------------------------------------------------------------
    SUBROUTINE FIX(BOUNDING_BOX, X_TEMP, Y_TEMP, Z_TEMP, X_OUT, Y_OUT, Z_OUT, INT_STATUS)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN):: X_TEMP, Y_TEMP, Z_TEMP
      DOUBLE PRECISION, DIMENSION(6), INTENT(IN):: BOUNDING_BOX

      INTEGER, INTENT(INOUT):: INT_STATUS
      DOUBLE PRECISION, INTENT(OUT):: X_OUT, Y_OUT, Z_OUT

      X_OUT    = X_TEMP
      Y_OUT    = Y_TEMP
      Z_OUT    = Z_TEMP

      IF (X_TEMP > BOUNDING_BOX(2)) THEN
        X_OUT = BOUNDING_BOX(2)
        INT_STATUS = 1

      ELSE IF (X_TEMP < BOUNDING_BOX(1)) THEN
        X_OUT = BOUNDING_BOX(1)
        INT_STATUS = 1

      END IF

      IF (Y_TEMP > BOUNDING_BOX(4)) THEN
        Y_OUT = BOUNDING_BOX(4)
        INT_STATUS = 1

      ELSE IF (Y_TEMP < BOUNDING_BOX(3)) THEN
        Y_OUT = BOUNDING_BOX(3)
        INT_STATUS = 1

      END IF

      IF (Z_TEMP > BOUNDING_BOX(6)) THEN
        Z_OUT = BOUNDING_BOX(6)
        INT_STATUS = 1

      ELSE IF (Z_TEMP < BOUNDING_BOX(5)) THEN
        Z_OUT = BOUNDING_BOX(5)
        INT_STATUS = 1

      END IF

    END SUBROUTINE FIX

END MODULE BOUNDARY_CONDITIONS
