!!--------------------------------------------------------------
!! Fortran 90/95 style module that contains a library of useful
!! mathematical operations.
!!
!! Certain subroutines translated from Shadden's FlowVC mymath.c
!!
!! Author:   Zachariah Irwin
!!           University of Colorado, Boulder
!! Version:  March 2020
!!--------------------------------------------------------------
MODULE TOOLKIT_MATH

  IMPLICIT NONE

  CONTAINS

    !!--------------------------------------------------------------------------------------------------------
    !> @brief Finds the maximum eigenvalue of a real, symmetric matrix A.
    !!
    !! GET_MAX_EIGENVALUE uses Householder reduction on a real, symmetric matrix A of size NxN to find the
    !! eigenvalues. It returns the maximum eigenvalue EIGVAL of A. A is modified on output.
    !! 
    !! @param[in, out]  A                  double array, matrix
    !!                                     size: N x N
    !! @param[out]      EIGVAL             double, maximum eigenvalue of A
    !! @param[in]       N                  integer, dimension of A
    !!--------------------------------------------------------------------------------------------------------
    SUBROUTINE GET_MAX_EIGENVALUE(A, N, EIGVAL)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: N
      DOUBLE PRECISION, DIMENSION(0:n-1,0:n-1), INTENT(INOUT):: A

      DOUBLE PRECISION, INTENT(OUT):: EIGVAL

      INTEGER:: l, k, j, i, m, iter, kk, kkk, kkkk
      DOUBLE PRECISION, DIMENSION(0:2):: e, d
      DOUBLE PRECISION:: scale, hh, h, g, f, s, r, p, dd, c, b

      DO i = 3-1, 1, -1

        l     = i - 1
        scale = 0.0D0
        h     = scale

        IF (l > 0) THEN

          DO k = 0, l, 1
            scale = scale + ABS(A(i,k))

            IF (scale == 0.0D0) THEN
              e(i) = a(i,l)

            ELSE
              A(i,k) = A(i,k) / scale
              h      = h + A(i,k)*A(i,k)

              f = A(i,l)

              IF (f > 0.0D0) THEN
                g = -SQRT(h)

              ELSE
                g = SQRT(h)

              END IF

              e(i)   = scale*g
              h      = h - f*g
              A(i,l) = f - g
              f      = 0.0D0

              DO j = 0, l, 1
                g = 0.0D0
                DO kk = 0, j, 1
                  g = g + A(j,kk)*A(i,kk)
                END DO
                DO kkk = j+1, l, 1
                  g = g + A(kkk,j)*A(i,kkk)
                END DO
                e(j) = g/h
                f    = f + e(j)*A(i,j)
              END DO

              hh = f/(h+h)

              DO j = 0, l, 1
                f    = a(i,j)
                e(j) = e(j) - hh*f
                g    = e(j)
                DO kkkk = 0, j
                  A(j,kkkk) = A(j,kkkk) - (f*e(kkkk) + g*A(i,kkkk))
                END DO
              END DO

            END IF

          END DO

        ELSE

          e(i) = A(i,l)
          d(i) = h

        END IF

      END DO

      e(0) = 0.0D0

      DO i = 0, 2, 1
        d(i) = A(i,i)
      END DO

      DO i = 1, 2, 1
        e(i-1) = e(i)
      END DO

      e(3-1) = 0.0D0

      DO l = 0, 2, 1
        iter = 0
        DO m = l, 3-1, 2
          dd = ABS(d(m)) + ABS(d(m+1))
            IF (ABS(e(m))+dd == dd) THEN
              exit
            END IF

          IF (m .NE. l) THEN
            IF (iter == 30) THEN
              print*, "Too many iterations in tqli"
              RETURN
            END IF
            g = (d(l+1)-d(l))/(2.0D0*e(l))
            CALL PYTHAG(g, 1.0D0, r)
            g = d(m) - d(l) + e(l)/(g+SIGN(r,g))
            s = 1.0D0
            c = s
            p = 0.0D0

            DO i = m-1, l, -1
              f = s*e(i)
              b = c*e(i)
              CALL PYTHAG(f, g, r)
              e(i+1) = r

              IF (r == 0.0D0) THEN
                d(i+1) = d(i+1) - p
                e(m) = 0.0D0
                exit
              END IF

              s      = f/r
              c      = g/r
              g      = d(i+1) - p
              r      = (d(i) - g)*s + 2.0D0*c*b
              p      = s*r
              d(i+1) = g + p
              g      = c*r - b
            END DO

            IF (r == 0.0D0 .AND. i >= l) THEN
              CYCLE

            ELSE
              d(l) = d(l) - p
              e(l) = g
              e(m) = 0.0D0

            END IF
          END IF

        END DO 
        iter = iter + 1

      END DO

      EIGVAL = MAX(d(0), d(1), d(2))

    END SUBROUTINE GET_MAX_EIGENVALUE

    !!--------------------------------------------------------------------------------------------------------
    !> @brief Finds the length of the third side of a triangle with side lengths A and B.
    !!
    !! 
    !! @param[out]  R                  double, length of third side of the triangle
    !! @param[in]   A                  double, length of one side of the triangle
    !! @param[in]   B                  double, length of one side of the triangle
    !!--------------------------------------------------------------------------------------------------------
    SUBROUTINE PYTHAG(A, B, R)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN):: A, B

      DOUBLE PRECISION, INTENT(OUT):: R

      DOUBLE PRECISION:: ABS_A, ABS_B

      ABS_A = ABS(A)
      ABS_B = ABS(B)

      IF (ABS_A > ABS_B) THEN
        R = ABS_A*SQRT(1.0D0 + (ABS_B/ABS_A)*(ABS_B/ABS_A))

      ELSE
        IF (ABS_B == 0.0D0) THEN
          R = 0.0D0

        ELSE
          R = ABS_B*SQRT(1.0D0 + (ABS_A/ABS_B)*(ABS_A*ABS_B))

        END IF

      END IF

    END SUBROUTINE PYTHAG

    !!--------------------------------------------------------------------------------------------------------
    !> @brief Converts polar coordinate data to a cartesian vector.
    !!
    !! 
    !! @param[out]  XY_OUT             double array, x and y coordinates
    !!                                 size: 2
    !! @param[in]   R_IN               double, r coordinate
    !! @param[in]   THeTA_IN           double, theta coordinate
    !!--------------------------------------------------------------------------------------------------------
    SUBROUTINE POLAR_TO_CARTESIAN(R_IN, THETA_IN, XY_OUT)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN):: R_IN
      DOUBLE PRECISION, INTENT(IN):: THETA_IN

      DOUBLE PRECISION, DIMENSION(2), INTENT(OUT):: XY_OUT

      XY_OUT(1) = R_IN*cos(THETA_IN)
      XY_OUT(2) = R_IN*sin(THETA_IN)

    END SUBROUTINE POLAR_TO_CARTESIAN

    !!--------------------------------------------------------------------------------------------------------
    !> @brief Converts cartesian coordinate data to a polar vector.
    !!
    !! 
    !! @param[out]  RT_OUT             double array, r and theta coordinates
    !!                                 size: 2
    !! @param[in]   X_IN               double, x coordinate
    !! @param[in]   Y_IN               double, y coordinate
    !!--------------------------------------------------------------------------------------------------------
    SUBROUTINE CARTESIAN_TO_POLAR(X_IN, Y_IN, RT_OUT)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN):: X_IN
      DOUBLE PRECISION, INTENT(IN):: Y_IN

      DOUBLE PRECISION, DIMENSION(2), INTENT(OUT):: RT_OUT

      RT_OUT(1) = SQRT(X_IN**2 + Y_IN**2)
      RT_OUT(2) = DATAN2(X_IN, Y_IN)

    END SUBROUTINE CARTESIAN_TO_POLAR

    !!-----------------------------------------------------------------------------------------------------
    !> @brief Utility function that converts Cartesian indices into a linear index to map into data arrays.
    !!
    !!
    !! @param[out]   INDEX_OUT     integer, linear index
    !! @param[in]    DIMN          integer, problem dimension (2D or 3D)
    !! @param[in]    INDICES       integer array, cartesian indices to be converted
    !!                             size: 3
    !! @param[in]    DIMNS         integer array, number of points along each dimension
    !!                             size: 3
    !!-----------------------------------------------------------------------------------------------------
    SUBROUTINE LINEAR_INDEXING(INDICES, DIMN, DIMNS, INDEX_OUT)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: DIMN
      INTEGER, DIMENSION(3), INTENT(IN):: INDICES, DIMNS
      
      INTEGER, INTENT(OUT):: INDEX_OUT

      IF (ANY(INDICES < 0)) THEN
        PRINT*, INDICES
        PRINT*, "INVALID INDICES"
        STOP
      END IF

      IF (DIMN .EQ. 2) THEN
        INDEX_OUT = INDICES(1) + DIMNS(1)*INDICES(2)
      ELSE
        INDEX_OUT = INDICES(1) + DIMNS(1)*(INDICES(2) + DIMNS(2)*INDICES(3))
      END IF

    END SUBROUTINE LINEAR_INDEXING

END MODULE TOOLKIT_MATH
