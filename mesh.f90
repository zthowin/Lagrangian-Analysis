!!------------------------------------------------------------------
!! Fortran 90/95 style module to perform mesh-based algorithms.
!!
!!
!! Author:   Zachariah Irwin
!!           University of Colorado, Boulder
!! Version:  March 2020
!!-----------------------------------------------------------------
MODULE MESH

  USE TOOLKIT_MATH

  IMPLICIT NONE

  CONTAINS

    !!--------------------------------------------------------------------------------------------------------
    !> @brief Finds the cell ID of a point inside a triangular/tetrahedral mesh
    !!
    !! PARTICLE_LOCATOR returns the cell ID that a given point lies within. The point is mapped back to an 
    !! initial cell's (CELL_ID) reference configuration. If the point has crossed a side/face of that cell, 
    !! the algorithm will search the adjacent cell(s) until the point is found.
    !!
    !! (see: Lohner,R. and Ambrosiano,J.(1990). "A Vectorized Particle-Tracer For Unstructured Grids")
    !!
    !! 
    !! @param[out]  LOC_ID             integer, updated cell ID associated with the point
    !! @param[out]  SHAPE_FUNCTIONS    double array, shape functions for point in corresponding element
    !!                                 size: 3
    !! @param[in]   DIMN               integer, problem dimension (2D or 3D)
    !! @param[in]   CELL_ID            integer, previous cell ID associated with the point 
    !! @param[in]   NUM_NODES          integer, number of nodes in the mesh
    !! @param[in]   NUM_CELLS          integer, number of cells in the mesh
    !! @param[in]   CONNECTIVITY       integer array, connectivity matrix
    !!                                 size: NUM_CELLS x 4
    !! @param[in]   ADJACENCY          integer array, adjacency matrix
    !!                                 size: NUM_CELLS x 4
    !! @param[in]   X                  double array, coordinates of point
    !!                                 size: 3
    !! @param[in]   COORDINATES        double array, node coordinates
    !!                                 size: NUM_NODES x 3
    !!--------------------------------------------------------------------------------------------------------
    SUBROUTINE PARTICLE_LOCATOR(DIMN, X, CELL_ID, NUM_NODES, NUM_CELLS, CONNECTIVITY, ADJACENCY, COORDINATES,&
                                LOC_ID, SHAPE_FUNCTIONS)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: DIMN, CELL_ID, NUM_NODES, NUM_CELLS
      INTEGER, DIMENSION(0:NUM_CELLS-1, 4), INTENT(IN):: CONNECTIVITY, ADJACENCY
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: X
      DOUBLE PRECISION, DIMENSION(0:NUM_NODES-1, 3), INTENT(IN):: COORDINATES

      INTEGER, INTENT(OUT):: LOC_ID
      DOUBLE PRECISION, DIMENSION(3), INTENT(OUT):: SHAPE_FUNCTIONS

      DOUBLE PRECISION, PARAMETER:: EPS = 1.0e-12

      LOGICAL:: IS_INSIDE, IS_CHECK
      INTEGER:: NOW_ID
      INTEGER, DIMENSION(4):: NODE_IDS
      DOUBLE PRECISION, DIMENSION(3):: NODE0, NODE1, NODE2, NODE3
      DOUBLE PRECISION:: X0, X1, X2, X3, Y0, Y1, Y2, Y3, Z0, Z1, Z2, Z3, A11, A12, A13, A21, A22, A23, A31, A32, A33,&
                         V, ZETA, ETA, XI, N1, N2, N3, N4, N_MIN

      !f2py depend(NUM_NODES) COORDINATES
      !f2py depend(NUM_CELLS) CONNECTIVITY, ADJACENCY

      IS_CHECK  = .TRUE.
      IS_INSIDE = .FALSE.
      NOW_ID    = CELL_ID

      !!------------------------------------
      !! Search for cell in which point lies
      !!------------------------------------
      DO WHILE (IS_CHECK)

        !!-----------------------------------------------------------------------
        !! A cell ID of -1 corresponds to a non-existent cell outside of the mesh
        !! Stop searching, return a cell ID of -1
        !!-----------------------------------------------------------------------
        IF (NOW_ID == -1) THEN
          IS_CHECK = .FALSE.
          LOC_ID   = NOW_ID
          RETURN
        END IF

        !!----------------------------
        !! Get current cell's node IDs
        !!----------------------------
        NODE_IDS = CONNECTIVITY(NOW_ID,:)
        NODE0    = COORDINATES(NODE_IDS(1),:)
        NODE1    = COORDINATES(NODE_IDS(2),:)
        NODE2    = COORDINATES(NODE_IDS(3),:)
        IF (DIMN .EQ. 3) THEN
          NODE3    = COORDINATES(NODE_IDS(4),:)
        END IF

        !!-------------------------------
        !! Store coordinates of each node
        !!-------------------------------
        X0 = NODE0(1)
        Y0 = NODE0(2)
        Z0 = NODE0(3)
        X1 = NODE1(1)
        Y1 = NODE1(2)
        Z1 = NODE1(3)
        X2 = NODE2(1)
        Y2 = NODE2(2)
        Z2 = NODE2(3)
        IF (DIMN == 3) THEN
          X3 = NODE3(1)
          Y3 = NODE3(2)
          Z3 = NODE3(3)
        END IF

        !!-----------------------------------------
        !! Evaluate shape functions, i.e. \phi_i(x)
        !!-----------------------------------------
        IF (DIMN == 2) THEN

          A11 = Y2 - Y0
          A12 = X0 - X2
          A21 = Y0 - Y1
          A22 = X1 - X0

          V     = (X1 - X0)*(Y2 - Y0) - (X2 - X0)*(Y1 - Y0)
          ETA   = (A11*(X(1) - X0) + A12*(X(2) - Y0)) / V
          ZETA  = (A21*(X(1) - X0) + A22*(X(2) - Y0)) / V;
          
          N1    = 1.0D0 - ZETA - ETA
          N2    = ETA
          N3    = ZETA

          IS_INSIDE = (N1>=0.0D0) .AND. (N1<=1.0D0) .AND. (N2>=0.0D0) .AND. (N2<=1.0D0) .AND.&
          (N3>=0.0D0) .AND. (N3<=1.0D0)

          !!-------------------------------------------------------
          !! If all 0 <= \phi_i(x) <= 1, point lies in current cell
          !!-------------------------------------------------------
          IF (IS_INSIDE) THEN
            IS_CHECK            = .FALSE.
            LOC_ID              = NOW_ID
            SHAPE_FUNCTIONS(1)  = ETA
            SHAPE_FUNCTIONS(2)  = ZETA
            SHAPE_FUNCTIONS(3)  = 0.0D0
            RETURN

          !!--------------------------------------------------
          !! Search in direction of the face the point crossed
          !!--------------------------------------------------
          ELSE
            N_MIN = MIN(N1, N2, N3)

            IF (N_MIN == N1) THEN
              NOW_ID = ADJACENCY(NOW_ID, 3)

            ELSE IF (N_MIN == N2) THEN
              NOW_ID = ADJACENCY(NOW_ID, 1)

            ELSE IF (N_MIN == N3) THEN
              NOW_ID = ADJACENCY(NOW_ID, 2)

            END IF

          END IF

        ELSE IF (DIMN == 3) THEN

          A11 = (Z3 - Z0)*(Y2 - Y3) - (Z2 - Z3)*(Y3 - Y0)
          A21 = (Z3 - Z0)*(Y0 - Y1) - (Z0 - Z1)*(Y3 - Y0)
          A31 = (Z1 - Z2)*(Y0 - Y1) - (Z0 - Z1)*(Y1 - Y2)
          A12 = (X3 - X0)*(Z2 - Z3) - (X2 - X3)*(Z3 - Z0)
          A22 = (X3 - X0)*(Z0 - Z1) - (X0 - X1)*(Z3 - Z0)
          A32 = (X1 - X2)*(Z0 - Z1) - (X0 - X1)*(Z1 - Z2)
          A13 = (Y3 - Y0)*(X2 - X3) - (Y2 - Y3)*(X3 - X0)
          A23 = (Y3 - Y0)*(X0 - X1) - (Y0 - Y1)*(X3 - X0)
          A33 = (Y1 - Y2)*(X0 - X1) - (Y0 - Y1)*(X1 - X2)

          V = (X1 - X0)*((Y2 - Y0)*(Z3 - Z0) - (Z2 - Z0)*(Y3 - Y0)) + &
          (X2 - X0)*((Y0 - Y1)*(Z3 - Z0) - (Z0 - Z1)*(Y3 - Y0)) + &
          (X3 - X0)*((Y1 - Y0)*(Z2 - Z0) - (Z1 - Z0)*(Y2 - Y0))

          ETA   = (A11*(X(1) - X0) + A12*(X(2) - Y0) + A13*(X(3) - Z0)) / V
          ZETA  = (A21*(X(1) - X0) + A22*(X(2) - Y0) + A23*(X(3) - Z0)) / V
          XI    = (A31*(X(1) - X0) + A32*(X(2) - Y0) + A33*(X(3) - Z0)) / V

          N1    = 1.0D0 - ETA - ZETA - XI
          N2    = ETA
          N3    = ZETA
          N4    = XI

          IS_INSIDE = (N1>=0.0D0) .AND. (N1<=1.0D0) .AND. (N2>=0.0D0) .AND. (N2<=1.0D0) .AND.&
          (N3>=0.0D0) .AND. (N3<=1.0D0) .AND. (N4>=0.0D0) .AND. (N4<=1.0D0)

          !!-------------------------------------------------------
          !! If all 0 <= \phi_i(x) <= 1, point lies in current cell
          !!-------------------------------------------------------
          IF (IS_INSIDE) THEN
            IS_CHECK            = .FALSE.
            LOC_ID              = NOW_ID
            SHAPE_FUNCTIONS(1)  = ETA
            SHAPE_FUNCTIONS(2)  = ZETA
            SHAPE_FUNCTIONS(3)  = XI
            RETURN

          !!--------------------------------------------------
          !! Search in direction of the face the point crossed
          !!--------------------------------------------------
          ELSE
            N_MIN = MIN(N1, N2, N3, N4)

            IF (N_MIN == N1) THEN
              NOW_ID = ADJACENCY(NOW_ID, 4)
            ELSE IF (N_MIN == N2) THEN
              NOW_ID = ADJACENCY(NOW_ID, 1)
            ELSE IF (N_MIN == N3) THEN
              NOW_ID = ADJACENCY(NOW_ID, 2)
            ELSE IF (N_MIN == N4) THEN
              NOW_ID = ADJACENCY(NOW_ID, 3)

            END IF

          END IF
        END IF

      END DO

    END SUBROUTINE PARTICLE_LOCATOR

    !!---------------------------------------------------------------------------------------
    !> @brief Generates initial seeding coordinates of Lagrangian tracers
    !!
    !! @param[out]  PARTICLE_COORDINATES   double array, coordinates of all points in domain
    !!                                     size: NUM_P x 3
    !! @param[in]   NUM_PARTICLES          integer, total number of points in domain
    !! @param[in]   REF_DIMNS              double array, spacing between points in each
    !!                                     direction
    !!                                     size: 3
    !! @param[in]   BOUNDING_BOX           double array, bounding box/boundaries of domain
    !!---------------------------------------------------------------------------------------
    SUBROUTINE SEED_POINTS_CARTESIAN(NUM_PARTICLES, REF_DIMNS, BOUNDING_BOX,&
                                     PARTICLE_COORDINATES)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: NUM_PARTICLES
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: REF_DIMNS
      DOUBLE PRECISION, DIMENSION(6), INTENT(IN):: BOUNDING_BOX

      DOUBLE PRECISION, DIMENSION(0:NUM_PARTICLES-1, 3), INTENT(OUT):: PARTICLE_COORDINATES

      INTEGER:: DIMN, P, I, J, K, NUM_X, NUM_Y, NUM_Z
      DOUBLE PRECISION:: DX, DY, DZ, X, Y, Z

      DX = REF_DIMNS(1)
      DY = REF_DIMNS(2)
      DZ = REF_DIMNS(3)

      IF (DZ .EQ. 0.0D0) THEN
        DIMN = 2
      ELSE
        DIMN = 3
      END IF

      NUM_X = INT((BOUNDING_BOX(2)-BOUNDING_BOX(1))/DX) + 1
      NUM_Y = INT((BOUNDING_BOX(4)-BOUNDING_BOX(3))/DY) + 1
      
      IF (DIMN .EQ. 3) THEN
        NUM_Z = INT((BOUNDING_BOX(6)-BOUNDING_BOX(5))/DZ) + 1
      ELSE
        NUM_Z = 1
      END IF

      DO K = 0, NUM_Z-1
        DO J = 0, NUM_Y-1
          DO I = 0, NUM_X-1

            CALL LINEAR_INDEXING([I,J,K], DIMN, [NUM_X, NUM_Y, NUM_Z], P)
            PARTICLE_COORDINATES(P, 1) = BOUNDING_BOX(1) + FLOAT(I)*DX
            PARTICLE_COORDINATES(P, 2) = BOUNDING_BOX(3) + FLOAT(J)*DY
            PARTICLE_COORDINATES(P, 3) = BOUNDING_BOX(5) + FLOAT(K)*DZ

          END DO
        END DO
      END DO
    END SUBROUTINE SEED_POINTS_CARTESIAN

END MODULE MESH
