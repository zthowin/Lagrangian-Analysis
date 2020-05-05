!!------------------------------------------------------------------
!! Fortran 90/95 style module to return an interpolation of a 
!! point's velocity.
!!
!!
!! Author:   Zachariah Irwin
!!           University of Colorado, Boulder
!! Version:  March 2020
!!-----------------------------------------------------------------
MODULE VELOCITY

  IMPLICIT NONE

  CONTAINS

    !!--------------------------------------------------------------------------------------------------------
    !> @brief Finds the velocity of a point inside a triangular/tetrahedral mesh
    !!
    !! VELOCITY_INTERPOLATOR returns the interpolated velocity of a point within a mesh given a set of shape
    !! functions evaluated at the point's coordinates and an element (cell) ID.
    !!
    !! @param[out]  V_LOC              double array, velocity of the point
    !!                                 size: 3
    !! @param[in]   DIMN               integer, problem dimension (2D or 3D)
    !! @param[in]   CELL_ID            integer, current cell ID associated with the point 
    !! @param[in]   NUM_NODES          integer, number of nodes in the mesh
    !! @param[in]   NUM_CELLS          integer, number of cells in the mesh
    !! @param[in]   T                  double, current time associated with point
    !! @param[in]   CONNECTIVITY       integer array, connectivity matrix
    !!                                 size: NUM_CELLS x 4
    !! @param[in]   ADJACENCY          integer array, adjacency matrix
    !!                                 size: NUM_CELLS x 4
    !! @param[in]   X                  double array, coordinates of point
    !!                                 size: 3
    !! @param[in]   WINDOW             double array, time data associated with loaded velocity data
    !!                                 size: 3
    !! @param[in]   SHAPE_FUNCTIONS    double array, shape functions for point in corresponding element
    !!                                 size: 3
    !! @param[in]   COORDINATES        double array, node coordinates
    !!                                 size: NUM_NODES x 3
    !! @param[in]   VEL_LOW|MID|UP     double array, velocity data
    !!                                 size: NUM_NODES x 3
    !!--------------------------------------------------------------------------------------------------------
    SUBROUTINE VELOCITY_INTERPOLATION(DIMN, X, CELL_ID, NUM_NODES, NUM_CELLS, CONNECTIVITY, SHAPE_FUNCTIONS,&
                                      VEL_LOW, VEL_MID, VEL_UP, T, WINDOW, V_LOC)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: DIMN, CELL_ID, NUM_NODES, NUM_CELLS
      INTEGER, DIMENSION(0:NUM_CELLS-1, 4), INTENT(IN):: CONNECTIVITY
      DOUBLE PRECISION, INTENT(IN):: T
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: X, WINDOW, SHAPE_FUNCTIONS
      DOUBLE PRECISION, DIMENSION(0:NUM_NODES-1, 3), INTENT(IN):: VEL_LOW, VEL_MID, VEL_UP

      DOUBLE PRECISION, DIMENSION(3), INTENT(OUT):: V_LOC

      DOUBLE PRECISION, PARAMETER:: EPS = 1.0e-10

      LOGICAL:: LOAD_UPPER, LOAD_LOWER
      INTEGER, DIMENSION(4):: NODE_IDS
      DOUBLE PRECISION:: ETA, ZETA, XI, N1, N2, N3, N4, T_LOC
      DOUBLE PRECISION, DIMENSION(3):: V_DOWN_LOC, V_UP_LOC

      !f2py depend(NUM_NODES) COORDINATES, VEL_LOW, VEL_MID, VEL_UP
      !f2py depend(NUM_CELLS) CONNECTIVITY

      !!------------------------------------------
      !! If point has left the mesh, hold it fixed
      !!------------------------------------------
      IF (CELL_ID == -1) THEN
        V_LOC = 0.00D0
        RETURN
      END IF

      LOAD_UPPER = .FALSE.
      LOAD_LOWER = .FALSE.

      !!--------------
      !! T <= Window 1
      !!--------------
      IF (T <= WINDOW(1)) THEN
        T_LOC       = 0.0D0
        LOAD_LOWER  = .TRUE.

      !!------------------------
      !! Window 2 < T < Window 3
      !!------------------------
      ELSE IF (T > WINDOW(2)) THEN
        T_LOC        = (T - WINDOW(2))/(WINDOW(3) - WINDOW(2))
        LOAD_UPPER  = .TRUE.

      !!------------------------
      !! Window 1 < T < Window 2
      !!------------------------
      ELSE
        T_LOC        = (T - WINDOW(1))/(WINDOW(2) - WINDOW(1))
        LOAD_LOWER  = .TRUE.

      END IF

      IF ((T_LOC > 1.0D0 + EPS) .OR. (T_LOC < 0.0D0 - EPS)) THEN
        PRINT*, "SIM TIME:", T
        PRINT*, "LOWEST WINDOW:", WINDOW(1)
        PRINT*, "MIDDLE WINDOW:", WINDOW(2)
        PRINT*, "HIGHEST WINDOW:", WINDOW(3)
        PRINT*, "T_LOC MUST BE BETWEEN 0 AND 1"
        STOP
      END IF

      !!----------------------------------
      !! Store cell node IDs and \phi_i(x)
      !!----------------------------------
      NODE_IDS = CONNECTIVITY(CELL_ID,:)
      ETA      = SHAPE_FUNCTIONS(1)
      ZETA     = SHAPE_FUNCTIONS(2)
      XI       = SHAPE_FUNCTIONS(3)
      
      IF (DIMN == 2) THEN
          
        N1    = 1.0D0 - ZETA - ETA
        N2    = ETA
        N3    = ZETA

        !!------------
        !! Interpolate
        !!------------
        IF (LOAD_UPPER) THEN
          V_DOWN_LOC(1) = VEL_MID(NODE_IDS(1),1)*N1 + VEL_MID(NODE_IDS(2),1)*N2 + VEL_MID(NODE_IDS(3),1)*N3
          V_DOWN_LOC(2) = VEL_MID(NODE_IDS(1),2)*N1 + VEL_MID(NODE_IDS(2),2)*N2 + VEL_MID(NODE_IDS(3),2)*N3
          V_DOWN_LOC(3) = 0.0D0

          V_UP_LOC(1)   = VEL_UP(NODE_IDS(1),1)*N1 + VEL_UP(NODE_IDS(2),1)*N2 + VEL_UP(NODE_IDS(3),1)*N3
          V_UP_LOC(2)   = VEL_UP(NODE_IDS(1),2)*N1 + VEL_UP(NODE_IDS(2),2)*N2 + VEL_UP(NODE_IDS(3),2)*N3
          V_UP_LOC(3)   = 0.0D0

        ELSE IF (LOAD_LOWER) THEN
          V_DOWN_LOC(1) = VEL_LOW(NODE_IDS(1),1)*N1 + VEL_LOW(NODE_IDS(2),1)*N2 + VEL_LOW(NODE_IDS(3),1)*N3
          V_DOWN_LOC(2) = VEL_LOW(NODE_IDS(1),2)*N1 + VEL_LOW(NODE_IDS(2),2)*N2 + VEL_LOW(NODE_IDS(3),2)*N3
          V_DOWN_LOC(3) = 0.0D0

          V_UP_LOC(1)   = VEL_MID(NODE_IDS(1),1)*N1 + VEL_MID(NODE_IDS(2),1)*N2 + VEL_MID(NODE_IDS(3),1)*N3
          V_UP_LOC(2)   = VEL_MID(NODE_IDS(1),2)*N1 + VEL_MID(NODE_IDS(2),2)*N2 + VEL_MID(NODE_IDS(3),2)*N3
          V_UP_LOC(3)   = 0.0D0

        END IF

        V_LOC           = (1.0D0 - T_LOC)*V_DOWN_LOC + T_LOC*V_UP_LOC

      ELSE IF (DIMN == 3) THEN    

        N1    = 1.0D0 - ETA - ZETA - XI
        N2    = ETA
        N3    = ZETA
        N4    = XI

        !!------------
        !! Interpolate
        !!------------
        IF (LOAD_UPPER) THEN

          V_DOWN_LOC(1) = VEL_MID(NODE_IDS(1),1)*N1 + VEL_MID(NODE_IDS(2),1)*N2 + VEL_MID(NODE_IDS(3),1)*N3 + &
                          VEL_MID(NODE_IDS(4),1)*N4
          V_DOWN_LOC(2) = VEL_MID(NODE_IDS(1),2)*N1 + VEL_MID(NODE_IDS(2),2)*N2 + VEL_MID(NODE_IDS(3),2)*N3 + &
                          VEL_MID(NODE_IDS(4),2)*N4
          V_DOWN_LOC(3) = VEL_MID(NODE_IDS(1),3)*N1 + VEL_MID(NODE_IDS(2),3)*N2 + VEL_MID(NODE_IDS(3),3)*N3 + &
                          VEL_MID(NODE_IDS(4),3)*N4

          V_UP_LOC(1)   = VEL_UP(NODE_IDS(1),1)*N1 + VEL_UP(NODE_IDS(2),1)*N2 + VEL_UP(NODE_IDS(3),1)*N3 + &
                          VEL_UP(NODE_IDS(4),1)*N4
          V_UP_LOC(2)   = VEL_UP(NODE_IDS(1),2)*N1 + VEL_UP(NODE_IDS(2),2)*N2 + VEL_UP(NODE_IDS(3),2)*N3 + &
                          VEL_UP(NODE_IDS(4),2)*N4
          V_UP_LOC(3)   = VEL_UP(NODE_IDS(1),3)*N1 + VEL_UP(NODE_IDS(2),3)*N2 + VEL_UP(NODE_IDS(3),3)*N3 + &
                          VEL_UP(NODE_IDS(4),3)*N4

        ELSE IF (LOAD_LOWER) THEN
          V_DOWN_LOC(1) = VEL_LOW(NODE_IDS(1),1)*N1 + VEL_LOW(NODE_IDS(2),1)*N2 + VEL_LOW(NODE_IDS(3),1)*N3 + &
                          VEL_LOW(NODE_IDS(4),1)*N4
          V_DOWN_LOC(2) = VEL_LOW(NODE_IDS(1),2)*N1 + VEL_LOW(NODE_IDS(2),2)*N2 + VEL_LOW(NODE_IDS(3),2)*N3 + &
                          VEL_LOW(NODE_IDS(4),2)*N4
          V_DOWN_LOC(3) = VEL_LOW(NODE_IDS(1),3)*N1 + VEL_LOW(NODE_IDS(2),3)*N2 + VEL_LOW(NODE_IDS(3),3)*N3 + &
                          VEL_LOW(NODE_IDS(4),3)*N4

          V_UP_LOC(1)   = VEL_MID(NODE_IDS(1),1)*N1 + VEL_MID(NODE_IDS(2),1)*N2 + VEL_MID(NODE_IDS(3),1)*N3 + &
                          VEL_MID(NODE_IDS(4),1)*N4
          V_UP_LOC(2)   = VEL_MID(NODE_IDS(1),2)*N1 + VEL_MID(NODE_IDS(2),2)*N2 + VEL_MID(NODE_IDS(3),2)*N3 + &
                          VEL_MID(NODE_IDS(4),2)*N4
          V_UP_LOC(3)   = VEL_MID(NODE_IDS(1),3)*N1 + VEL_MID(NODE_IDS(2),3)*N2 + VEL_MID(NODE_IDS(3),3)*N3 + &
                          VEL_MID(NODE_IDS(4),3)*N4
        END IF

        V_LOC           = (1.0D0 - T_LOC)*V_DOWN_LOC + T_LOC*V_UP_LOC

      END IF

    END SUBROUTINE VELOCITY_INTERPOLATION

END MODULE VELOCITY
