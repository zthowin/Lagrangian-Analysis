!!-----------------------------------------------------------------------------
!! Fortran 90/95/2003 style module that integrates a grid of cartesian tracers.
!!
!!
!! Author:  Zachariah Irwin
!!          University of Colorado, Boulder
!! Version: March 2020
!!-----------------------------------------------------------------------------
MODULE INTEGRATION

  USE ANALYTICAL_FLOW_LIBRARY_FT
  USE TOOLKIT_MATH
  USE BOUNDARY_CONDITIONS
  USE VELOCITY
  USE MESH
  
  IMPLICIT NONE

  CONTAINS

    !!--------------------------------------------------------------------------------------------------------
    !> @brief Prints the simulation integration parameters.
    !!
    !! @param[in]      FLOW_MODEL       integer, flag corresponding to the relevant flow model
    !! @param[in]      BC_TYPE          integer, flag corresponding to the relevant boundary condition
    !! @param[in]      INT_SCHEME       integer, flag corresponding to the relevant integration scheme
    !! @param[in]      INT_TIME         double, total integration time
    !! @param[in]      INT_TIME_STEP    double, integration time step
    !! @param[in]      BOUNDING_BOX     double array, bounding box/boundaries of the mesh domain
    !!                                  size: 6
    !!--------------------------------------------------------------------------------------------------------
    SUBROUTINE PRINT_INTEGRATION_PARAMETERS(FLOW_MODEL, BC_TYPE, INT_SCHEME, INT_TIME, INT_TIME_STEP, BOUNDING_BOX)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: FLOW_MODEL, BC_TYPE, INT_SCHEME
      DOUBLE PRECISION, INTENT(IN):: INT_TIME, INT_TIME_STEP
      DOUBLE PRECISION, DIMENSION(6), INTENT(IN):: BOUNDING_BOX

      IF (BOUNDING_BOX(5) == BOUNDING_BOX(6)) THEN

        PRINT*,"Flow Model Choice     :", FLOW_MODEL
        PRINT*,"Boundary Condition    :", BC_TYPE
        IF (INT_SCHEME == 1) THEN
          PRINT*, "Integration Scheme    :  RK4"
        ELSE IF (INT_SCHEME == 2) THEN
          PRINT*, "Integration Scheme    :  Euler"
        END IF
        PRINT*,"Integration Window    :", INT_TIME
        PRINT*,"Integration Time Step :", INT_TIME_STEP
        PRINT*,"Domain X Limits       :", BOUNDING_BOX(1), BOUNDING_BOX(2)
        PRINT*,"Domain Y Limits       :", BOUNDING_BOX(3), BOUNDING_BOX(4)

      ELSE IF (BOUNDING_BOX(5) < BOUNDING_BOX(6)) THEN

        PRINT*,"Flow Model Choice     :", FLOW_MODEL
        PRINT*,"Boundary Condition    :", BC_TYPE
        IF (INT_SCHEME == 1) THEN
          PRINT*, "Integration Scheme    :  RK4"
        ELSE IF (INT_SCHEME == 2) THEN
          PRINT*, "Integration Scheme    :  Euler"
        END IF
        PRINT*,"Integration Time      :", INT_TIME
        PRINT*,"Integration Time Step :", INT_TIME_STEP
        PRINT*,"Domain X Limits       :", BOUNDING_BOX(1), BOUNDING_BOX(2)
        PRINT*,"Domain Y Limits       :", BOUNDING_BOX(3), BOUNDING_BOX(4)
        PRINT*,"Domain Z Limits       :", BOUNDING_BOX(5), BOUNDING_BOX(6)

      END IF

    END SUBROUTINE PRINT_INTEGRATION_PARAMETERS
    
    !!--------------------------------------------------------------------------------------------------------
    !> @brief Makes data index and time windows synced to the number of time points in the simulation for
    !! appropriate velocity field loading
    !!
    !! @param[out]     DATA_INDEX_WINDOWS   integer array, data index windows synced to number of time points
    !!                                      in simulation
    !!                                      size: NUM_T+1 x 3
    !! @param[out]     DATA_TIME_WINDOWS    double array, data time windows synced to number of time points in
    !!                                      the simulation
    !!                                      size: NUM_T+1 x 3
    !! @param[in]      NUM_FILES            integer, number of velocity files in the data set
    !! @param[in]      DATA_INDEX_START     integer, index of first velocity file
    !! @param[in]      DATA_INDEX_END       integer, index of last velocity file
    !! @param[in]      DATA_INDEX_DELTA     integer, spacing between velocity file indices
    !! @param[in]      NUM_T                double, number of time points in the simulation
    !! @param[in]      DATA_START           double, timing of first velocity file
    !! @param[in]      DATA_END             double, timing of last velocity file
    !! @param[in]      DATA_DT              double, spacing between velocity file times
    !! @param[in]      SIM_START            double, simulation start time
    !! @param[in]      SIM_END              double, simulation stop time
    !! @param[in]      SIM_DT               double, simulation time step
    !!--------------------------------------------------------------------------------------------------------
    SUBROUTINE MAKE_DATA_WINDOWS(NUM_FILES, DATA_INDEX_START, DATA_INDEX_END, DATA_INDEX_DELTA,&
                                 NUM_T, DATA_START, DATA_END, DATA_DT, SIM_START, SIM_END, SIM_DT,&
                                 DATA_INDEX_WINDOWS, DATA_TIME_WINDOWS)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: NUM_FILES, DATA_INDEX_START, DATA_INDEX_END, DATA_INDEX_DELTA, NUM_T
      DOUBLE PRECISION, INTENT(IN):: DATA_START, DATA_END, DATA_DT, SIM_START, SIM_END, SIM_DT

      INTEGER, DIMENSION(0:NUM_T,3), INTENT(OUT):: DATA_INDEX_WINDOWS
      DOUBLE PRECISION, DIMENSION(0:NUM_T,3), INTENT(OUT):: DATA_TIME_WINDOWS

      INTEGER:: ID_LOW, ID_MID, ID_UP, I
      INTEGER, DIMENSION(0:NUM_FILES-1):: DATA_INDICES
      INTEGER, DIMENSION(0:NUM_T):: ID_LOW_ARR, ID_MID_ARR, ID_UP_ARR
      DOUBLE PRECISION:: T_LOW, T_MID, T_UP
      DOUBLE PRECISION, DIMENSION(0:NUM_FILES-1):: DATA_TIMING
      DOUBLE PRECISION, DIMENSION(0:NUM_T):: T_LOW_ARR, T_MID_ARR, T_UP_ARR
      DOUBLE PRECISION, DIMENSION(0:NUM_T):: INTEGRATION_TIMES

      !f2py depend(NUM_T) DATA_INDEX_WINDOWS, DATA_TIME_WINDOWS

      !!------------------------------------------------------------
      !! Build an array of file indexes and corresponding file times
      !!------------------------------------------------------------
      DO I = 0, NUM_FILES-1
        DATA_INDICES(I) = DATA_INDEX_START + I*DATA_INDEX_DELTA
        DATA_TIMING(I)  = DATA_START + I*DATA_DT
      END DO

      !!---------------------------------------------
      !! Build an array of integration time-instances
      !!---------------------------------------------
      DO I = 0, NUM_T
        INTEGRATION_TIMES(I) = SIM_START + I*SIM_DT
      END DO

      !!---------------------------------------------------------
      !! Build the lower, middle, and upper bound interval arrays
      !!---------------------------------------------------------
      DO I = 0, NUM_T
        T_LOW_ARR(I) = FLOOR(INTEGRATION_TIMES(I)/DATA_DT)*DATA_DT
        T_MID_ARR(I) = CEILING(INTEGRATION_TIMES(I)/DATA_DT)*DATA_DT
        T_UP_ARR(I)  = CEILING((INTEGRATION_TIMES(I) + DATA_DT)/DATA_DT)*DATA_DT

        ID_LOW_ARR(I) = T_LOW_ARR(I)/DATA_DT
        ID_MID_ARR(I) = T_MID_ARR(I)/DATA_DT
        ID_UP_ARR(I)  = T_UP_ARR(I)/DATA_DT

        !!------------------------------
        !! Adjust for periodic data sets
        !!------------------------------
        IF (ID_LOW_ARR(I) >= NUM_FILES) THEN
          ID_LOW_ARR(I) = MOD(ID_LOW_ARR(I), NUM_FILES)
        END IF

        IF (ID_MID_ARR(I) >= NUM_FILES) THEN
          ID_MID_ARR(I) = MOD(ID_MID_ARR(I), NUM_FILES)
        END IF

        IF (ID_UP_ARR(I) >= NUM_FILES) THEN
          ID_UP_ARR(I) = MOD(ID_UP_ARR(I), NUM_FILES)
        END IF

        !!---------------------------
        !! Set the data window array.
        !!---------------------------
        DATA_TIME_WINDOWS(I,1) = T_LOW_ARR(I)
        DATA_TIME_WINDOWS(I,2) = T_MID_ARR(I)
        DATA_TIME_WINDOWS(I,3) = T_UP_ARR(I)

        DATA_INDEX_WINDOWS(I,1) = ID_LOW_ARR(I)
        DATA_INDEX_WINDOWS(I,2) = ID_MID_ARR(I)
        DATA_INDEX_WINDOWS(I,3) = ID_UP_ARR(I)
      END DO

    END SUBROUTINE MAKE_DATA_WINDOWS

    !!---------------------------------------------------------------------
    !> @brief Loads a velocity file into an array
    !!
    !! @param[out]     VEL_OUT              double array, velocity field
    !!                                      size: NUM_NODES x 3
    !! @param[in]      FILE_INDEX           integer, index of velocity file
    !! @param[in]      NUM_NODES            integer, number of mesh nodes
    !! @param[in]      FILE_PREFIX          character, velocity file prefix
    !!                                      len: 1024
    !!---------------------------------------------------------------------
    SUBROUTINE LOAD_VELOCITY(FILE_INDEX, NUM_NODES, FILE_PREFIX, VEL_OUT)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: FILE_INDEX, NUM_NODES
      CHARACTER(LEN=1024), INTENT(IN):: FILE_PREFIX

      DOUBLE PRECISION, DIMENSION(0:NUM_NODES-1, 3), INTENT(INOUT):: VEL_OUT

      INTEGER:: UNIT, I
      CHARACTER(LEN=1024):: FILENAME, FILE_INDEX_STR

      !f2py depend(NUM_NODES) VEL_OUT

      UNIT = 10

      !!----------------------------------------
      !! Convert velocity file index into string
      !!----------------------------------------
      WRITE(FILE_INDEX_STR, '(I0)') FILE_INDEX

      !!-----------------------------------------------------
      !! Get velocity file name based on the prefix and index
      !!-----------------------------------------------------
      FILENAME = TRIM(FILE_PREFIX)//'.'//TRIM(FILE_INDEX_STR)//'.bin'

      !!-----------------------------------
      !! Prepare to read velocity file data
      !!-----------------------------------
      OPEN(NEWUNIT=UNIT, FILE=TRIM(FILENAME), STATUS='OLD', FORM='UNFORMATTED', ACCESS='STREAM')

      !!----------------------------
      !! Loop over all nodes in mesh
      !!----------------------------
      DO I = 1, NUM_NODES*3, 3
        
        !!-----------------------------------
        !! Do not skip bytes for node ID == 0
        !!-----------------------------------
        IF (I == 1) THEN
          READ(UNIT, POS=I) VEL_OUT((I/3), 1)
        ELSE
          READ(UNIT, POS=I*8 + 1) VEL_OUT((I/3), 1)
        END IF

        READ(UNIT, POS=(I+1)*8 + 1) VEL_OUT((I/3), 2)
        READ(UNIT, POS=(I+2)*8 + 1) VEL_OUT((I/3), 3)
      END DO
      CLOSE(UNIT)

    END SUBROUTINE LOAD_VELOCITY

    !!--------------------------------------------------------------------------------------------------------
    !> @brief Advects a grid of points from a given velocity field
    !!
    !! ADVECT_GRID_DATA is responsible for advecting an entire grid of cartesian point data. Each point is
    !! advected individually from T0 to T1 using either the Euler or RK4 integration scheme (fowards or
    !! backwards depends on sign of (T1 - T0). 
    !!
    !! The velocity field responsible for advection is obtained from simulation data. FLOW_MODEL == 0 for 
    !! simulation data.
    !!
    !! Currently, this subroutine and its dependencies can only interpolate triangular and tetrahedral meshes.
    !!
    !! @param[in, out]  PARTICLES_XYZ   double array, grid coordinates
    !!                                  size: NUM_P x 3
    !! @param[in, out]  LEFT_GRID       integer array, flags corresponding to integration status of a point
    !!                                  size: NUM_P   
    !! @param[in, out]  CELL_IDS        integer array, the cell IDs of every point in the cartesian grid
    !!                                  size: NUM_P
    !! @param[in]       NUM_P           integer, number of points in the simulation
    !! @param[in]       NUM_NODES       integer, number of nodes in the mesh
    !! @param[in]       NUM_CELLS       integer, number of cells in the mesh
    !! @param[in]       INT_FLAGS       integer array, flags corresponding to advection parameters
    !!                                  size: 3
    !! @param[in]       CONNECTIVITY    integer array, mesh connectivity matrix
    !!                                  size: NUM_CELLS x 4
    !! @param[in]       ADJACENCY       integer array, mesh adjacency matrix
    !!                                  size: NUM_CELLS x 4
    !! @param[in]       TIMES           double array, simulation time, simulation time step
    !!                                  size: 2
    !! @param[in]       WINDOW          double array, upper, middle, and lower time bounds of velocity data
    !!                                  size: 3
    !! @param[in]       BOUNDING_BOX    double array, bounding box/boundaries of the mesh domain
    !!                                  size: 6
    !! @param[in]       COORDINATES     double array, mesh node coordinates
    !!                                  size: NUM_NODES x 3
    !! @param[in]       VEL_LOW|MID|UP  double array, velocity data
    !!                                  size: NUM_NODES x 3
    !!--------------------------------------------------------------------------------------------------------
    SUBROUTINE ADVECT_GRID_DATA(INT_FLAGS, BOUNDING_BOX, NUM_P, PARTICLES_XYZ, CELL_IDS, LEFT_GRID, TIMES, &
                                WINDOW, NUM_NODES, NUM_CELLS, CONNECTIVITY, ADJACENCY, COORDINATES, &
                                VEL_LOW, VEL_MID, VEL_UP)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: NUM_P, NUM_NODES, NUM_CELLS
      INTEGER, DIMENSION(3), INTENT(IN):: INT_FLAGS
      INTEGER, DIMENSION(0:NUM_CELLS-1, 4), INTENT(IN):: CONNECTIVITY, ADJACENCY
      DOUBLE PRECISION, DIMENSION(2), INTENT(IN):: TIMES
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: WINDOW
      DOUBLE PRECISION, DIMENSION(6), INTENT(IN):: BOUNDING_BOX
      DOUBLE PRECISION, DIMENSION(0:NUM_NODES-1, 3), INTENT(IN):: COORDINATES, VEL_LOW, VEL_MID, VEL_UP

      INTEGER, DIMENSION(0:NUM_P-1), INTENT(INOUT):: CELL_IDS, LEFT_GRID
      DOUBLE PRECISION, DIMENSION(0:NUM_P-1, 3), INTENT(INOUT):: PARTICLES_XYZ

      DOUBLE PRECISION, PARAMETER:: EPS = 1.0e-12

      INTEGER:: P, DIMN, INT_SCHEME, BC_TYPE
      DOUBLE PRECISION:: T0, T1, DT
      DOUBLE PRECISION, DIMENSION(3):: XYZ

      !f2py depend(NUM_CELLS) CONNECTIVITY, ADJACENCY
      !f2py depend(NUM_NODES) COORDINATES, VEL_LOW, VEL_MID, VEL_UP
      !f2py depend(NUM_P) PARTICLES_XYZ, CELL_IDS, LEFT_GRID

      !!----------------------
      !! Initialize parameters
      !!----------------------
      DIMN        = INT_FLAGS(1)
      INT_SCHEME  = INT_FLAGS(2)
      BC_TYPE     = INT_FLAGS(3)

      T0 = TIMES(1)
      DT = TIMES(2)
      T1 = T0 + DT

      !!---------------------
      !! Loop over all points
      !!---------------------
      DO P = 0, NUM_P-1

        !!-------------------------------------------------------------
        !! CELL_ID == -1 corresponds to a point that is not in the mesh
        !! (and therefore not in the domain)
        !!-------------------------------------------------------------
        IF (CELL_IDS(P) == -1) THEN
          LEFT_GRID(P) = 1
        END IF

        !!-------------------------------------
        !! If point still in the domain, advect
        !!-------------------------------------
        IF (LEFT_GRID(P) == -1) THEN

          IF (INT_SCHEME == 1) THEN

            CALL RK4_DATA(DIMN, BC_TYPE, BOUNDING_BOX, PARTICLES_XYZ(P,:), CELL_IDS(P), LEFT_GRID(P),&
                          T0, T1, DT, WINDOW, NUM_NODES, NUM_CELLS, CONNECTIVITY, ADJACENCY, COORDINATES,&
                          VEL_LOW, VEL_MID, VEL_UP)

          ELSE IF (INT_SCHEME == 2) THEN

            CALL EULER_DATA(DIMN, BC_TYPE, BOUNDING_BOX, PARTICLES_XYZ(P,:), CELL_IDS(P), LEFT_GRID(P),&
                            T0, T1, DT, WINDOW, NUM_NODES, NUM_CELLS, CONNECTIVITY, ADJACENCY, COORDINATES,&
                            VEL_LOW, VEL_MID, VEL_UP)

          END IF

        END IF

      END DO

    END SUBROUTINE ADVECT_GRID_DATA

    !!--------------------------------------------------------------------------------------------------------
    !> @brief Advects a grid of points from a given velocity field
    !!
    !! ADVECT_GRID_ANALYTICAL is responsible for advecting an entire grid of cartesian point data. The velocity
    !! field responsible for advection is dependent on the FLOW_MODEL. See analyticalFlowLibraryFT.f90 for a
    !! complete list of analytical velocity fields.
    !!
    !! @param[in, out]  PARTICLES_XYZ   double array, grid coordinates
    !!                                  size: NUM_P x 3
    !! @param[in, out]  LEFT_GRID       integer array, flags corresponding to integration status of a point
    !!                                  size: NUM_P 
    !! @param[in]       FLOW_MODEL      integer, flag corresponding to the relevant analytical velocity field
    !! @param[in]       NUM_P           integer, number of points in the simulation
    !! @param[in]       INT_FLAGS       integer array, flags corresponding to advection parameters
    !!                                  size: 3
    !! @param[in]       TIMES           double array, simulation time, simulation time step
    !!                                  size: 2
    !! @param[in]       BOUNDING_BOX    double array, bounding box/boundaries of the mesh domain
    !!                                  size: 6
    !!--------------------------------------------------------------------------------------------------------
    SUBROUTINE ADVECT_GRID_ANALYTICAL(FLOW_MODEL, INT_FLAGS, BOUNDING_BOX, NUM_P, PARTICLES_XYZ, LEFT_GRID, TIMES)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: FLOW_MODEL, NUM_P
      INTEGER, DIMENSION(3), INTENT(IN):: INT_FLAGS
      DOUBLE PRECISION, DIMENSION(2), INTENT(IN):: TIMES
      DOUBLE PRECISION, DIMENSION(6), INTENT(IN):: BOUNDING_BOX
      
      INTEGER, DIMENSION(0:NUM_P-1), INTENT(INOUT):: LEFT_GRID
      DOUBLE PRECISION, DIMENSION(0:NUM_P-1, 3), INTENT(INOUT):: PARTICLES_XYZ
      
      INTEGER:: P, INT_SCHEME, BC_TYPE
      DOUBLE PRECISION:: T0, T1, DT

      !f2py depend(NUM_P) PARTICLES_XYZ, LEFT_GRID

      !!----------------------
      !! Initialize parameters
      !!----------------------
      INT_SCHEME = INT_FLAGS(2)
      BC_TYPE    = INT_FLAGS(3)

      T0 = TIMES(1)
      DT = TIMES(2)
      T1 = T0 + DT

      !!---------------------
      !! Loop over all points
      !!---------------------
      DO P = 0, NUM_P-1

        !!-------------------------------------
        !! If point still in the domain, advect
        !!-------------------------------------
        IF (LEFT_GRID(P) == -1) THEN
  
          IF (INT_SCHEME == 1) THEN  
            CALL RK4_ANALYTICAL(FLOW_MODEL, BC_TYPE, BOUNDING_BOX, PARTICLES_XYZ(P,:), LEFT_GRID(P),&
                                T0, T1, DT)
          
          ELSE IF (INT_SCHEME == 2) THEN
            CALL EULER_ANALYTICAL(FLOW_MODEL, BC_TYPE, BOUNDING_BOX, PARTICLES_XYZ(P,:), LEFT_GRID(P),&
                                  T0, T1, DT)

          END IF

        END IF

      END DO

    END SUBROUTINE ADVECT_GRID_ANALYTICAL

    !!--------------------------------------------------------------------------------------------------------
    !> @brief Advects a single point using forward Euler method
    !!
    !! EULER_DATA is responsible for advecting a single point from T0 to T1 given the coordinates X_IN, Y_IN, 
    !! Z_IN using forward Euler numerical integration method. The velocity field responsible for advection is 
    !! obtained from simulation data. 
    !! 
    !! Currently, this subroutine and its dependencies can only interpolate triangular and tetrahedral meshes.
    !!
    !! @param[in, out]  CELL_ID            integer, the cell ID of the point of interest
    !! @param[in, out]  INT_STATUS         integer, flag corresponding to whether or not point left domain
    !! @param[in, out]  XYZ                double, point coordinates
    !! @param[in]       DIMN               integer, problem dimension (2D or 3D)
    !! @param[in]       BC_TYPE            integer, flag corresponding to the relevant boundary condition
    !! @param[in]       NUM_NODES          integer, number of nodes in the mesh
    !! @param[in]       NUM_CELLS          integer, number of cells in the mesh
    !! @param[in]       T0                 double, t_n     of integration
    !! @param[in]       T1                 double, t_{n+1} of integration
    !! @param[in]       DT                 double, integration time step
    !! @param[in]       CONNECTIVITY       integer array, connectivity matrix
    !!                                     size: NUM_CELLS x 4
    !! @param[in]       ADJACENCY          integer array, adjacency matrix
    !!                                     size: NUM_CELLS x 4   
    !! @param[in]       WINDOW             double array, upper, middle, and lower time bounds of velocity data
    !!                                     size: 3
    !! @param[in]       BOUNDING_BOX       double array, bounding box/boundaries of the mesh domain
    !!                                     size: 6
    !! @param[in]       COORDINATES        double array, node coordinates
    !!                                     size: NUM_NODES x 3
    !! @param[in]       VEL_LOW|MID|UP     double array, velocity data
    !!                                     size: NUM_NODES x 3
    !!--------------------------------------------------------------------------------------------------------
    SUBROUTINE EULER_DATA(DIMN, BC_TYPE, BOUNDING_BOX, XYZ, CELL_ID, INT_STATUS, T0, T1, DT, WINDOW,&
                          NUM_NODES, NUM_CELLS, CONNECTIVITY, ADJACENCY, COORDINATES, VEL_LOW, VEL_MID, VEL_UP)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: DIMN, BC_TYPE, NUM_NODES, NUM_CELLS
      INTEGER, DIMENSION(0:NUM_CELLS-1, 4), INTENT(IN):: CONNECTIVITY, ADJACENCY
      DOUBLE PRECISION, INTENT(IN):: T0, T1, DT
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: WINDOW
      DOUBLE PRECISION, DIMENSION(6), INTENT(IN):: BOUNDING_BOX
      DOUBLE PRECISION, DIMENSION(0:NUM_NODES-1, 3), INTENT(IN):: COORDINATES, VEL_LOW, VEL_MID, VEL_UP

      INTEGER, INTENT(INOUT):: INT_STATUS, CELL_ID
      DOUBLE PRECISION, DIMENSION(3), INTENT(INOUT):: XYZ
    
      INTEGER:: LOC_ID, INT_STATUS_RESOLVED
      DOUBLE PRECISION:: X_TEMP, Y_TEMP, Z_TEMP, X_RESOLVED, Y_RESOLVED, Z_RESOLVED, T, H
      DOUBLE PRECISION, DIMENSION(3):: SHAPE_FUNCTIONS, VEL_LOC

      DOUBLE PRECISION, PARAMETER:: EPS = 1.0E-14

      !f2py depend(NUM_NODES) COORDINATES, VEL_LOW, VEL_MID, VEL_UP
      !f2py depend(NUM_CELLS) CONNECTIVITY, ADJACENCY

      !!------------------------------------------------
      !! Initialize modifiable variables from input data
      !!------------------------------------------------
      T           = T0
      H           = SIGN(DT, T1 - T0)
      X_TEMP      = XYZ(1)
      Y_TEMP      = XYZ(2)
      Z_TEMP      = XYZ(3)

      IF (INT_STATUS == 1) THEN
        PRINT*, "ATTEMPTING TO INTEGRATE POINT THAT HAS LEFT THE DOMAIN"
        STOP
      END IF

      DO WHILE ((T - T1) * (T1 - T0) < 0.0D0)
        
        IF ((T + H - T1) * (T + H - T0) > 0.0D0) THEN
          H = SIGN(T1 - T, T1 - T0)
        END IF

        !!------------------------------------------
        !! Check to see if point has left the domain
        !!------------------------------------------
        CALL RESOLVE_BOUNDARY_CONDITIONS(BC_TYPE, BOUNDING_BOX, X_TEMP, Y_TEMP, Z_TEMP, X_TEMP, Y_TEMP, Z_TEMP,&
                                         X_RESOLVED, Y_RESOLVED, Z_RESOLVED, INT_STATUS_RESOLVED)

        !!----------------------------------------------
        !! Stop integration if point has left the domain
        !!----------------------------------------------
        IF (INT_STATUS_RESOLVED == 1) THEN
          XYZ(1)     = X_RESOLVED
          XYZ(2)     = Y_RESOLVED
          XYZ(3)     = Z_RESOLVED
          INT_STATUS = 1
          RETURN
        END IF

        !!------------------------------
        !! Find cell in which point lies
        !!------------------------------
        CALL PARTICLE_LOCATOR(DIMN, [X_RESOLVED, Y_RESOLVED, Z_RESOLVED], CELL_ID, NUM_NODES, NUM_CELLS,&
                              CONNECTIVITY, ADJACENCY, COORDINATES, LOC_ID, SHAPE_FUNCTIONS)

        !!-------------------------------------------------------
        !! Interpolate point's velocity based on location in mesh
        !!-------------------------------------------------------
        CALL VELOCITY_INTERPOLATION(DIMN, [X_RESOLVED, Y_RESOLVED, Z_RESOLVED], LOC_ID, NUM_NODES, NUM_CELLS,&
                                    CONNECTIVITY, SHAPE_FUNCTIONS, VEL_LOW, VEL_MID, VEL_UP, T, WINDOW, VEL_LOC)

        !!-------------------------------------------------------
        !! Update point's coordinates using numerical integration
        !!-------------------------------------------------------
        CALL RESOLVE_BOUNDARY_CONDITIONS(BC_TYPE, BOUNDING_BOX, X_TEMP, Y_TEMP, Z_TEMP, &
                                         (X_TEMP + VEL_LOC(1)*H), (Y_TEMP + VEL_LOC(2)*H), &
                                         (Z_TEMP + VEL_LOC(3)*H), X_RESOLVED, Y_RESOLVED, Z_RESOLVED,&
                                         INT_STATUS_RESOLVED)

        !!------------------------
        !! Update points's cell ID
        !!------------------------
        CALL PARTICLE_LOCATOR(DIMN, [X_RESOLVED, Y_RESOLVED, Z_RESOLVED], CELL_ID, NUM_NODES, NUM_CELLS, &
                              CONNECTIVITY, ADJACENCY, COORDINATES, LOC_ID, SHAPE_FUNCTIONS)

        CELL_ID = LOC_ID

        !!----------------------------------------------
        !! Stop integration if point has left the domain
        !!----------------------------------------------
        IF (INT_STATUS_RESOLVED == 1) THEN
          XYZ(1)     = X_RESOLVED
          XYZ(2)     = Y_RESOLVED
          XYZ(3)     = Z_RESOLVED
          INT_STATUS = 1
          RETURN

        ELSE
          !!-----------------
          !! Periodic (+X) BC
          !!-----------------
          IF (BC_TYPE == 3) THEN
            X_TEMP = (X_TEMP+VEL_LOC(1)*H)
            Y_TEMP = Y_RESOLVED
            Z_TEMP = Z_RESOLVED

            ! Temporary workaround for single time step loops (T0 + DT = T1)
            XYZ(1) = X_TEMP
            XYZ(2) = Y_TEMP
            XYZ(3) = Z_TEMP

          ELSE
            X_TEMP  = X_RESOLVED
            Y_TEMP  = Y_RESOLVED
            Z_TEMP  = Z_RESOLVED

            ! Temporary workaround for single time step loops (T0 + DT = T1)
            XYZ(1) = X_TEMP
            XYZ(2) = Y_TEMP
            XYZ(3) = Z_TEMP

          END IF

        END IF

        T     = T + H
      
        END DO

    END SUBROUTINE EULER_DATA

    !!--------------------------------------------------------------------------------------------------------
    !> @brief Advects a single point using Runge-Kutta 4th order method
    !!
    !! RK4_DATA is responsible for advecting a single point from T0 to T1 given the coordinates XYZ
    !! using the 4th order Runge-Kutta numerical integration method. The velocity field responsible for 
    !! advection is obtained from simulation data. 
    !! 
    !! Currently, this subroutine and its dependencies can only interpolate triangular and tetrahedral meshes.
    !!
    !! @param[in, out]  CELL_ID            integer, the cell ID of the point
    !! @param[in, out]  INT_STATUS         integer, flag corresponding to whether or not point left domain
    !!                                     after integration
    !! @param[in, out]  XYZ                double, point coordinates
    !! @param[in]       DIMN               integer, problem dimension (2D or 3D)
    !! @param[in]       BC_TYPE            integer, flag corresponding to the relevant boundary condition
    !! @param[in]       NUM_NODES          integer, number of nodes in the mesh
    !! @param[in]       NUM_CELLS          integer, number of cells in the mesh
    !! @param[in]       T0                 double, t_n     of integration
    !! @param[in]       T1                 double, t_{n+1} of integration
    !! @param[in]       DT                 double, integration time step
    !! @param[in]       CONNECTIVITY       integer array, connectivity matrix
    !!                                     size: NUM_CELLS x 4
    !! @param[in]       ADJACENCY          integer array, adjacency matrix
    !!                                     size: NUM_CELLS x 4
    !! @param[in]       LEFT_GRID          integer array, flags corresponding to integration status of a point
    !!                                     size: NUM_X x NUM_Y x NUM_Z    
    !! @param[in]       WINDOW             double array, upper, middle, and lower time bounds of velocity data
    !!                                     size: 3
    !! @param[in]       BOUNDING_BOX       double array, bounding box/boundaries of the mesh domain
    !!                                     size: 6
    !! @param[in]       COORDINATES        double array, node coordinates
    !!                                     size: NUM_NODES x 3
    !! @param[in]       VEL_LOW|MID|UP     double array, velocity data
    !!                                     size: NUM_NODES x 3
    !!--------------------------------------------------------------------------------------------------------
    SUBROUTINE RK4_DATA(DIMN, BC_TYPE, BOUNDING_BOX, XYZ, CELL_ID, INT_STATUS, T0, T1, DT, WINDOW,&
                        NUM_NODES, NUM_CELLS, CONNECTIVITY, ADJACENCY, COORDINATES, VEL_LOW, VEL_MID, VEL_UP)
      
      IMPLICIT NONE

      INTEGER, INTENT(IN):: DIMN, BC_TYPE, NUM_NODES, NUM_CELLS
      INTEGER, DIMENSION(0:NUM_CELLS-1, 4), INTENT(IN):: CONNECTIVITY, ADJACENCY
      DOUBLE PRECISION, INTENT(IN):: T0, T1, DT
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: WINDOW
      DOUBLE PRECISION, DIMENSION(6), INTENT(IN):: BOUNDING_BOX
      DOUBLE PRECISION, DIMENSION(0:NUM_NODES-1, 3), INTENT(IN):: COORDINATES, VEL_LOW, VEL_MID, VEL_UP

      INTEGER, INTENT(INOUT):: INT_STATUS, CELL_ID
      DOUBLE PRECISION, DIMENSION(3), INTENT(INOUT):: XYZ
      
      INTEGER:: LOC_ID1, LOC_ID2, LOC_ID3, LOC_ID4, LOC_ID_OUT, INT_STATUS_RESOLVED
      DOUBLE PRECISION:: X_TEMP, Y_TEMP, Z_TEMP, X_RESOLVED, Y_RESOLVED, Z_RESOLVED, T, H
      DOUBLE PRECISION, DIMENSION(3):: K1, K2, K3, K4, SHAPE_FUNCTIONS, VEL_LOC

      DOUBLE PRECISION, PARAMETER:: EPS = 1.0e-12

      !f2py depend(NUM_NODES) COORDINATES, VEL_LOW, VEL_MID, VEL_UP
      !f2py depend(NUM_CELLS) CONNECTIVITY, ADJACENCY

      !!------------------------------------------------
      !! Initialize modifiable variables from input data
      !!------------------------------------------------
      T       = T0
      H       = SIGN(DT, T1 - T0)
      X_TEMP  = XYZ(1)
      Y_TEMP  = XYZ(2)
      Z_TEMP  = XYZ(3)

      IF (INT_STATUS == 1) THEN
        PRINT*, "ATTEMPTING TO INTEGRATE POINT THAT HAS LEFT THE DOMAIN"
        STOP
      END IF

      DO WHILE ((T - T1) * (T1 - T0) < 0.0D0)
        
        IF ((T + H - T1) * (T + H - T0) > 0.0D0) THEN
          H = SIGN(T1 - T, T1 - T0)
        END IF

        !!----------------------------------------------------
        ! Compute flow velocity at location of point (Stage 1)
        !!----------------------------------------------------
        !!------------------------------------------
        !! Check to see if point has left the domain
        !!------------------------------------------
        CALL RESOLVE_BOUNDARY_CONDITIONS(BC_TYPE, BOUNDING_BOX, X_TEMP, Y_TEMP, Z_TEMP, X_TEMP, Y_TEMP, Z_TEMP, &
                                         X_RESOLVED, Y_RESOLVED, Z_RESOLVED, INT_STATUS_RESOLVED)

        !!----------------------------------------------
        !! Stop integration if point has left the domain
        !!----------------------------------------------
        IF (INT_STATUS_RESOLVED == 1) THEN
          XYZ(1)     = X_RESOLVED
          XYZ(2)     = Y_RESOLVED
          XYZ(3)     = Z_RESOLVED
          INT_STATUS = 1
          RETURN
        END IF

        !!------------------------------
        !! Find cell in which point lies
        !!------------------------------
        CALL PARTICLE_LOCATOR(DIMN, [X_RESOLVED, Y_RESOLVED, Z_RESOLVED], CELL_ID, NUM_NODES, NUM_CELLS,&
                              CONNECTIVITY, ADJACENCY, COORDINATES, LOC_ID1, SHAPE_FUNCTIONS)

        !!-------------------------------------------------------
        !! Interpolate point's velocity based on location in mesh
        !!-------------------------------------------------------
        CALL VELOCITY_INTERPOLATION(DIMN, [X_RESOLVED, Y_RESOLVED, Z_RESOLVED], LOC_ID1, NUM_NODES, NUM_CELLS, &
                                    CONNECTIVITY, SHAPE_FUNCTIONS, VEL_LOW, VEL_MID, VEL_UP, T, WINDOW, VEL_LOC)
    
        K1  = H*VEL_LOC

        !!-----------------------------------------------------
        !! Compute flow velocity at location of point (Stage 2)
        !!-----------------------------------------------------
        !!------------------------------------------
        !! Check to see if point has left the domain
        !!------------------------------------------
        CALL RESOLVE_BOUNDARY_CONDITIONS(BC_TYPE, BOUNDING_BOX, X_TEMP, Y_TEMP, Z_TEMP, &
                                         (X_TEMP + 0.50D0*K1(1)), (Y_TEMP + 0.50D0*K1(2)), (Z_TEMP + 0.50D0*K1(3)),&
                                         X_RESOLVED, Y_RESOLVED, Z_RESOLVED, INT_STATUS_RESOLVED)
        
        !!----------------------------------------------
        !! Stop integration if point has left the domain
        !!----------------------------------------------
        IF (INT_STATUS_RESOLVED == 1) THEN
          XYZ(1)     = X_RESOLVED
          XYZ(2)     = Y_RESOLVED
          XYZ(3)     = Z_RESOLVED
          INT_STATUS = 1
          RETURN
        END IF

        !!------------------------------
        !! Find cell in which point lies
        !!------------------------------
        CALL PARTICLE_LOCATOR(DIMN, [X_RESOLVED, Y_RESOLVED, Z_RESOLVED], CELL_ID, NUM_NODES, NUM_CELLS, &
                              CONNECTIVITY, ADJACENCY, COORDINATES, LOC_ID2, SHAPE_FUNCTIONS)
        
        !!-------------------------------------------------------
        !! Interpolate point's velocity based on location in mesh
        !!-------------------------------------------------------
        CALL VELOCITY_INTERPOLATION(DIMN, [X_RESOLVED, Y_RESOLVED, Z_RESOLVED], LOC_ID2, NUM_NODES, NUM_CELLS, &
                                    CONNECTIVITY, SHAPE_FUNCTIONS, VEL_LOW, VEL_MID, VEL_UP, (T + 0.50D0*H), WINDOW, VEL_LOC)

        K2 = H*VEL_LOC

        !!-------------------------------------------------------
        !! Compute flow velocity at location of tracer (Stage 3).
        !!-------------------------------------------------------
        CALL RESOLVE_BOUNDARY_CONDITIONS(BC_TYPE, BOUNDING_BOX, X_TEMP, Y_TEMP, Z_TEMP, &
                                         (X_TEMP + 0.50D0*K2(1)), (Y_TEMP + 0.50D0*K2(2)), (Z_TEMP + 0.50D0*K2(3)), &
                                         X_RESOLVED, Y_RESOLVED, Z_RESOLVED, INT_STATUS_RESOLVED)

        !!----------------------------------------------
        !! Stop integration if point has left the domain
        !!----------------------------------------------
        IF (INT_STATUS_RESOLVED == 1) THEN
          XYZ(1)     = X_RESOLVED
          XYZ(2)     = Y_RESOLVED
          XYZ(3)     = Z_RESOLVED
          INT_STATUS = 1
          RETURN
        END IF

        !!------------------------------
        !! Find cell in which point lies
        !!------------------------------
        CALL PARTICLE_LOCATOR(DIMN, [X_RESOLVED, Y_RESOLVED, Z_RESOLVED], CELL_ID, NUM_NODES, NUM_CELLS, &
                              CONNECTIVITY, ADJACENCY, COORDINATES, LOC_ID3, SHAPE_FUNCTIONS)
        
        !!-------------------------------------------------------
        !! Interpolate point's velocity based on location in mesh
        !!-------------------------------------------------------
        CALL VELOCITY_INTERPOLATION(DIMN, [X_RESOLVED, Y_RESOLVED, Z_RESOLVED], LOC_ID3, NUM_NODES, NUM_CELLS, &
                                    CONNECTIVITY, SHAPE_FUNCTIONS, VEL_LOW, VEL_MID, VEL_UP, (T + 0.50D0*H), WINDOW, VEL_LOC)

        K3 = H*VEL_LOC

        !!-----------------------------------------------------
        !! Compute flow velocity at location of point (Stage 4)
        !!-----------------------------------------------------
        CALL RESOLVE_BOUNDARY_CONDITIONS(BC_TYPE, BOUNDING_BOX, X_TEMP, Y_TEMP, Z_TEMP, &
                                         (X_TEMP + K3(1)), (Y_TEMP + K3(2)), (Z_TEMP + K3(3)), &
                                         X_RESOLVED, Y_RESOLVED, Z_RESOLVED, INT_STATUS_RESOLVED)

        !!----------------------------------------------
        !! Stop integration if point has left the domain
        !!----------------------------------------------
        IF (INT_STATUS_RESOLVED == 1) THEN
          XYZ(1)     = X_RESOLVED
          XYZ(2)     = Y_RESOLVED
          XYZ(3)     = Z_RESOLVED
          INT_STATUS = 1
          RETURN
        END IF

        !!------------------------------
        !! Find cell in which point lies
        !!------------------------------
        CALL PARTICLE_LOCATOR(DIMN, [X_RESOLVED, Y_RESOLVED, Z_RESOLVED], CELL_ID, NUM_NODES, NUM_CELLS, &
                              CONNECTIVITY, ADJACENCY, COORDINATES, LOC_ID4, SHAPE_FUNCTIONS)
        
        !!-------------------------------------------------------
        !! Interpolate point's velocity based on location in mesh
        !!-------------------------------------------------------
        CALL VELOCITY_INTERPOLATION(DIMN, [X_RESOLVED, Y_RESOLVED, Z_RESOLVED], LOC_ID4, NUM_NODES, NUM_CELLS, &
                                    CONNECTIVITY, SHAPE_FUNCTIONS, VEL_LOW, VEL_MID, VEL_UP, (T + H), WINDOW, VEL_LOC)

        K4 = H*VEL_LOC

        !!-----------------------------------------------------
        !! Update point coordinates using numerical integration
        !!-----------------------------------------------------
        !!------------------------------------------
        !! Check to see if point has left the domain
        !!------------------------------------------
        CALL RESOLVE_BOUNDARY_CONDITIONS(BC_TYPE, BOUNDING_BOX, X_TEMP, Y_TEMP, Z_TEMP, &
                                         (X_TEMP + (1.0D0/6.0D0)*(K1(1) + 2.0D0*K2(1) + 2.0D0*K3(1) + K4(1))), &
                                         (Y_TEMP + (1.0D0/6.0D0)*(K1(2) + 2.0D0*K2(2) + 2.0D0*K3(2) + K4(2))), &
                                         (Z_TEMP + (1.0D0/6.0D0)*(K1(3) + 2.0D0*K2(3) + 2.0D0*K3(3) + K4(3))), &
                                         X_RESOLVED, Y_RESOLVED, Z_RESOLVED, INT_STATUS_RESOLVED)

        !!-----------------------
        !! Update point's cell ID
        !!-----------------------
        CALL PARTICLE_LOCATOR(DIMN, [X_RESOLVED, Y_RESOLVED, Z_RESOLVED], CELL_ID, NUM_NODES, NUM_CELLS, &
                              CONNECTIVITY, ADJACENCY, COORDINATES, LOC_ID_OUT, SHAPE_FUNCTIONS)

        CELL_ID = LOC_ID_OUT

        !!----------------------------------------------
        !! Stop integration if point has left the domain
        !!----------------------------------------------
        IF (INT_STATUS_RESOLVED == 1) THEN
          XYZ(1)     = X_RESOLVED
          XYZ(2)     = Y_RESOLVED
          XYZ(3)     = Z_RESOLVED
          INT_STATUS = 1
          RETURN

        ELSE
          !!-----------------
          !! Periodic BC (+x)
          !!-----------------
          IF (BC_TYPE == 3) THEN
            X_TEMP = X_TEMP+(1.0D0/6.0D0)*(K1(1) + 2.0D0*K2(1) + 2.0D0*K3(1) + K4(1))
            Y_TEMP = Y_RESOLVED
            Z_TEMP = Z_RESOLVED

            ! Temporary workaround for single time step loops (T0 + DT = T1).
            XYZ(1)  = X_TEMP
            XYZ(2)  = Y_TEMP
            XYZ(3)  = Z_TEMP

          ELSE
            X_TEMP  = X_RESOLVED
            Y_TEMP  = Y_RESOLVED
            Z_TEMP  = Z_RESOLVED

            ! Temporary workaround for single time step loops (T0 + DT = T1).
            XYZ(1)  = X_TEMP
            XYZ(2)  = Y_TEMP
            XYZ(3)  = Z_TEMP

          END IF

        END IF

        T = T + H
          
      END DO

    END SUBROUTINE RK4_DATA

    !!--------------------------------------------------------------------------------------------------------
    !> @brief Advects a single point using forward Euler method
    !!
    !! EULER_ANALYTICAL is responsible for advecting a single point from T0 to T1 given the coordinates
    !! XYZ using forward Euler numerical integration method. FLOW_MODEL is a flag that indicates 
    !! the analytical velocity field responsible for advecting the point.
    !! 
    !! @param[in, out]  INT_STATUS         integer, flag corresponding to whether or not point left domain after integration
    !! @param[in, out]  XYZ                double, point coordinates
    !! @param[in]       BC_TYPE            integer, flag corresponding to the relevant boundary condition
    !! @param[in]       FLOW_MODEL         integer, flag corresponding to the relevant analytical velocity field
    !! @param[in]       T0                 double, t_n     of integration
    !! @param[in]       T1                 double, t_{n+1} of integration
    !! @param[in]       DT                 double, integration time step
    !! @param[in]       BOUNDING_BOX       double array, bounding box/boundaries of the mesh domain
    !!                                     size: 6
    !!--------------------------------------------------------------------------------------------------------
    SUBROUTINE EULER_ANALYTICAL(FLOW_MODEL, BC_TYPE, BOUNDING_BOX, XYZ, INT_STATUS, T0, T1, DT)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: BC_TYPE, FLOW_MODEL
      DOUBLE PRECISION, INTENT(IN):: T0, T1, DT
      DOUBLE PRECISION, DIMENSION(6), INTENT(IN):: BOUNDING_BOX

      INTEGER, INTENT(INOUT):: INT_STATUS
      DOUBLE PRECISION, DIMENSION(3), INTENT(INOUT):: XYZ

      INTEGER:: INT_STATUS_RESOLVED
      DOUBLE PRECISION:: T, H, X_RESOLVED, Y_RESOLVED, Z_RESOLVED, X_TEMP, Y_TEMP, Z_TEMP
      DOUBLE PRECISION, DIMENSION(3):: VEL_LOC

      T      = T0
      H      = SIGN(DT, T1 - T0)
      X_TEMP = XYZ(1)
      Y_TEMP = XYZ(2)
      Z_TEMP = XYZ(3)
  
      IF (INT_STATUS == 1) THEN
        PRINT*, "ATTEMPTING TO INTEGRATE POINT THAT HAS LEFT THE DOMAIN"
        STOP
      END IF
      
      DO WHILE ((T - T1) * (T1 - T0) < 0.0D0)

        IF ((T + H - T1) * (T + H - T0) > 0.0D0) THEN
          H = SIGN(T1 - T, T1 - T0)
        END IF

        !!-------------------------------------------
        !! Compute flow velocity at location of point
        !!-------------------------------------------
        !!------------------------------------------
        !! Check to see if point has left the domain
        !!------------------------------------------
        CALL RESOLVE_BOUNDARY_CONDITIONS(BC_TYPE, BOUNDING_BOX, X_TEMP, Y_TEMP, Z_TEMP, X_TEMP, Y_TEMP, Z_TEMP, &
                                         X_RESOLVED, Y_RESOLVED, Z_RESOLVED, INT_STATUS_RESOLVED)
        
        !!----------------------------------------------
        !! Stop integration if point has left the domain
        !!----------------------------------------------
        IF (INT_STATUS_RESOLVED == 1) THEN
          XYZ(1)     = X_RESOLVED
          XYZ(2)     = Y_RESOLVED
          XYZ(3)     = Z_RESOLVED
          INT_STATUS = 1
        END IF

        IF (FLOW_MODEL == 1) THEN
          CALL DOUBLE_GYRE(X_RESOLVED, Y_RESOLVED, VEL_LOC)

        ELSE IF (FLOW_MODEL == 2) THEN
          CALL UNSTEADY_DOUBLE_GYRE(X_RESOLVED, Y_RESOLVED, T, VEL_LOC)

        ELSE IF (FLOW_MODEL == 3) THEN
          CALL BICKLEY_JET(X_RESOLVED, Y_RESOLVED, T, VEL_LOC)

        ELSE IF (FLOW_MODEL == 4) THEN
          CALL LAMB_OSSEEN(X_RESOLVED, Y_RESOLVED, T, VEL_LOC)

        ELSE IF (FLOW_MODEL == 5) THEN
          CALL ABC(X_TEMP, Y_TEMP, Z_TEMP, VEL_LOC)

        ELSE IF (FLOW_MODEL == 6) THEN
          CALL UNSTEADY_ABC(X_TEMP, Y_TEMP, Z_TEMP, T, VEL_LOC)

        END IF
        
        !!-----------------------------------------------------
        !! Update point coordinates using numerical integration
        !!-----------------------------------------------------
        !!------------------------------------------
        !! Check to see if point has left the domain
        !!------------------------------------------
        CALL RESOLVE_BOUNDARY_CONDITIONS(BC_TYPE, BOUNDING_BOX, X_TEMP, Y_TEMP, Z_TEMP, &
                                         (X_TEMP + VEL_LOC(1)*H), (Y_TEMP + VEL_LOC(2)*H), (Z_TEMP + VEL_LOC(3)*H), &
                                         X_RESOLVED, Y_RESOLVED, Z_RESOLVED, INT_STATUS_RESOLVED)
        
        !!----------------------------------------------
        !! Stop integration if point has left the domain
        !!----------------------------------------------
        IF (INT_STATUS_RESOLVED == 1) THEN
          XYZ(1)     = X_RESOLVED
          XYZ(2)     = Y_RESOLVED
          XYZ(3)     = Z_RESOLVED
          INT_STATUS = 1
          RETURN

        ELSE
          !!-----------------
          !! Periodic BC (+x)
          !!-----------------
          IF (BC_TYPE == 3) THEN
            X_TEMP = X_TEMP + (VEL_LOC(1)*H)
            Y_TEMP = Y_RESOLVED
            Z_TEMP = Z_RESOLVED

            ! Temporary workaround for single time step loops (T0 + DT = T1)
            XYZ(1) = X_TEMP + (VEL_LOC(1)*H)
            XYZ(2) = Y_RESOLVED
            XYZ(3) = Z_RESOLVED

          ELSE
            X_TEMP = X_RESOLVED
            Y_TEMP = Y_RESOLVED
            Z_TEMP = Z_RESOLVED
            
            ! Temporary workaround for single time step loops (T0 + DT = T1)
            XYZ(1) = X_RESOLVED
            XYZ(2) = Y_RESOLVED
            XYZ(3) = Z_RESOLVED
          
          END IF

        END IF

        T = T + H

      END DO

    END SUBROUTINE EULER_ANALYTICAL

    !!--------------------------------------------------------------------------------------------------------
    !> @brief Advects a single point using Runge-Kutta 4th order method.
    !!
    !! RK4_ANALYTICAL is responsible for advecting a single point from T0 to T1 given the coordinates
    !! X_IN, Y_IN, Z_IN using the 4th order Runge-Kutta numerical integration method. FLOW_MODEL is a 
    !! flag that indicates the analytical velocity field responsible for advecting the point.
    !! 
    !! @param[out]      INT_STATUS_RESOLVED     integer, flag corresponding to whether or not point left domain after integration
    !! @param[out]      X|Y|Z_OUT          double, point coordinates after integration
    !! @param[in]       INT_STATUS_IN      integer, flag corresponding to whether or not point left domain before integration
    !! @param[in]       BC_TYPE            integer, flag corresponding to the relevant boundary condition
    !! @param[in]       FLOW_MODEL         integer, flag corresponding to the relevant analytical velocity field
    !! @param[in]       T0                 double, t_n     of integration
    !! @param[in]       T1                 double, t_{n+1} of integration
    !! @param[in]       DT                 double, integration time step
    !! @param[in]       X|Y|Z_IN           double, point coordinates before integration
    !! @param[in]       BOUNDING_BOX       double array, bounding box/boundaries of the mesh domain
    !!                                     size: 6
    !!--------------------------------------------------------------------------------------------------------
    SUBROUTINE RK4_ANALYTICAL(FLOW_MODEL, BC_TYPE, BOUNDING_BOX, XYZ, INT_STATUS, T0, T1, DT)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: BC_TYPE, FLOW_MODEL
      DOUBLE PRECISION, INTENT(IN):: T0, T1, DT
      DOUBLE PRECISION, DIMENSION(6), INTENT(IN):: BOUNDING_BOX

      INTEGER, INTENT(INOUT):: INT_STATUS
      DOUBLE PRECISION, DIMENSION(3), INTENT(INOUT):: XYZ

      INTEGER:: INT_STATUS_RESOLVED
      DOUBLE PRECISION:: T, H, X_RESOLVED, Y_RESOLVED, Z_RESOLVED, X_TEMP, Y_TEMP, Z_TEMP
      DOUBLE PRECISION, DIMENSION(3):: VEL_LOC, K1, K2, K3, K4
      
      !!-------------------------------------------------
      !! Initialize modifiable variables from input data.
      !!-------------------------------------------------
      T      = T0
      H      = SIGN(DT, T1 - T0)
      X_TEMP = XYZ(1)
      Y_TEMP = XYZ(2)
      Z_TEMP = XYZ(3)

      IF (INT_STATUS == 1) THEN
        PRINT*, "ATTEMPTING TO INTEGRATE POINT THAT HAS LEFT THE DOMAIN"
        STOP
      END IF

      T = T0
      H = SIGN(DT, T1 - T0)
  
      DO WHILE ((T - T1) * (T1 - T0) < 0.0D0)

        IF ((T + H - T1) * (T + H - T0) > 0.0D0) THEN
          H = SIGN(T1 - T, T1 - T0)
        END IF

        !!-----------------------------------------------------
        !! Compute flow velocity at location of point (Stage 1)
        !!-----------------------------------------------------
        !!------------------------------------------
        !! Check to see if point has left the domain
        !!------------------------------------------
        CALL RESOLVE_BOUNDARY_CONDITIONS(BC_TYPE, BOUNDING_BOX, X_TEMP, Y_TEMP, Z_TEMP, X_TEMP, Y_TEMP, Z_TEMP, &
                                         X_RESOLVED, Y_RESOLVED, Z_RESOLVED, INT_STATUS_RESOLVED)
 
        !!----------------------------------------------
        !! Stop integration if point has left the domain
        !!----------------------------------------------
        IF (INT_STATUS_RESOLVED == 1) THEN
          XYZ(1)     = X_RESOLVED
          XYZ(2)     = Y_RESOLVED
          XYZ(3)     = Z_RESOLVED
          INT_STATUS = 1
        END IF

        IF (FLOW_MODEL == 1) THEN
          CALL DOUBLE_GYRE(X_RESOLVED, Y_RESOLVED, VEL_LOC)

        ELSE IF (FLOW_MODEL == 2 ) THEN
          CALL UNSTEADY_DOUBLE_GYRE(X_RESOLVED, Y_RESOLVED, T, VEL_LOC)

        ELSE IF (FLOW_MODEL == 3) THEN
          CALL BICKLEY_JET(X_RESOLVED, Y_RESOLVED, T, VEL_LOC)

        ELSE IF (FLOW_MODEL == 4) THEN
          CALL LAMB_OSSEEN(X_RESOLVED, Y_RESOLVED, T, VEL_LOC)

        ELSE IF (FLOW_MODEL == 5) THEN
          CALL ABC(X_TEMP, Y_TEMP, Z_TEMP, VEL_LOC)

        ELSE IF (FLOW_MODEL == 6) THEN
          CALL UNSTEADY_ABC(X_TEMP, Y_TEMP, Z_TEMP, T, VEL_LOC)

        END IF
    
        K1  = H*VEL_LOC

        !!-----------------------------------------------------
        !! Compute flow velocity at location of point (Stage 2)
        !!-----------------------------------------------------
        !!------------------------------------------
        !! Check to see if point has left the domain
        !!------------------------------------------
        CALL RESOLVE_BOUNDARY_CONDITIONS(BC_TYPE, BOUNDING_BOX, X_TEMP, Y_TEMP, Z_TEMP, &
                                         (X_TEMP + 0.50D0*K1(1)), (Y_TEMP + 0.50D0*K1(2)), (Z_TEMP + 0.50D0*K1(3)),&
                                         X_RESOLVED, Y_RESOLVED, Z_RESOLVED, INT_STATUS_RESOLVED)
        
        !!----------------------------------------------
        !! Stop integration if point has left the domain
        !!----------------------------------------------
        IF (INT_STATUS_RESOLVED == 1) THEN
          XYZ(1)     = X_RESOLVED
          XYZ(2)     = Y_RESOLVED
          XYZ(3)     = Z_RESOLVED
          INT_STATUS = 1
        END IF

        IF (FLOW_MODEL == 1) THEN
          CALL DOUBLE_GYRE(X_RESOLVED, Y_RESOLVED, VEL_LOC)

        ELSE IF (FLOW_MODEL == 2) THEN
          CALL UNSTEADY_DOUBLE_GYRE(X_RESOLVED, Y_RESOLVED, (T + 0.50D0*H), VEL_LOC)

        ELSE IF (FLOW_MODEL == 3) THEN
          CALL BICKLEY_JET(X_RESOLVED, Y_RESOLVED, (T + 0.50D0*H), VEL_LOC)

        ELSE IF (FLOW_MODEL == 4) THEN
          CALL LAMB_OSSEEN(X_RESOLVED, Y_RESOLVED, (T + 0.50D0*H), VEL_LOC)

        ELSE IF (FLOW_MODEL == 5) THEN
          CALL ABC(X_RESOLVED, Y_RESOLVED, Z_RESOLVED, VEL_LOC)

        ELSE IF (FLOW_MODEL == 6) THEN
          CALL UNSTEADY_ABC(X_RESOLVED, Y_RESOLVED, Z_RESOLVED, (T + 0.50D0*H), VEL_LOC)

        END IF

        K2  = H*VEL_LOC
        
        !!-----------------------------------------------------
        !! Compute flow velocity at location of point (Stage 3)
        !!-----------------------------------------------------
        !!------------------------------------------
        !! Check to see if point has left the domain
        !!------------------------------------------
        CALL RESOLVE_BOUNDARY_CONDITIONS(BC_TYPE, BOUNDING_BOX, X_TEMP, Y_TEMP, Z_TEMP, &
                                         (X_TEMP + 0.50D0*K2(1)), (Y_TEMP + 0.50D0*K2(2)), (Z_TEMP + 0.50D0*K2(3)),& 
                                         X_RESOLVED, Y_RESOLVED, Z_RESOLVED, INT_STATUS_RESOLVED)

        !!----------------------------------------------
        !! Stop integration if point has left the domain
        !!----------------------------------------------
        IF (INT_STATUS_RESOLVED == 1) THEN
          XYZ(1)     = X_RESOLVED
          XYZ(2)     = Y_RESOLVED
          XYZ(3)     = Z_RESOLVED
          INT_STATUS = 1
        END IF

        IF (FLOW_MODEL == 1) THEN
          CALL DOUBLE_GYRE(X_RESOLVED, Y_RESOLVED, VEL_LOC)

        ELSE IF (FLOW_MODEL == 2) THEN
          CALL UNSTEADY_DOUBLE_GYRE(X_RESOLVED, Y_RESOLVED, (T + 0.50D0*H), VEL_LOC)

        ELSE IF (FLOW_MODEL == 3) THEN
          CALL BICKLEY_JET(X_RESOLVED, Y_RESOLVED, (T + 0.50D0*H), VEL_LOC)

        ELSE IF (FLOW_MODEL == 4) THEN
          CALL LAMB_OSSEEN(X_RESOLVED, Y_RESOLVED, (T + 0.50D0*H), VEL_LOC)

        ELSE IF (FLOW_MODEL == 5) THEN
          CALL ABC(X_RESOLVED, Y_RESOLVED, Z_RESOLVED, VEL_LOC)

        ELSE IF (FLOW_MODEL == 6) THEN
          CALL UNSTEADY_ABC(X_RESOLVED, Y_RESOLVED, Z_RESOLVED, (T + 0.50D0*H), VEL_LOC)

        END IF

        K3  = H*VEL_LOC
      
        !!-----------------------------------------------------
        !! Compute flow velocity at location of point (Stage 4)
        !!-----------------------------------------------------
        !!------------------------------------------
        !! Check to see if point has left the domain
        !!------------------------------------------
        CALL RESOLVE_BOUNDARY_CONDITIONS(BC_TYPE, BOUNDING_BOX, X_TEMP, Y_TEMP, Z_TEMP, X_TEMP+K3(1), Y_TEMP+K3(2), Z_TEMP+K3(3),&
                                         X_RESOLVED, Y_RESOLVED, Z_RESOLVED, INT_STATUS_RESOLVED)

        !!----------------------------------------------
        !! Stop integration if point has left the domain
        !!----------------------------------------------
        IF (INT_STATUS_RESOLVED == 1) THEN
          XYZ(1)     = X_RESOLVED
          XYZ(2)     = Y_RESOLVED
          XYZ(3)     = Z_RESOLVED
          INT_STATUS = 1
        END IF

        IF (FLOW_MODEL == 1) THEN
          CALL DOUBLE_GYRE(X_RESOLVED, Y_RESOLVED, VEL_LOC)

        ELSE IF (FLOW_MODEL == 2) THEN
          CALL UNSTEADY_DOUBLE_GYRE(X_RESOLVED, Y_RESOLVED, (T + H), VEL_LOC)

        ELSE IF (FLOW_MODEL == 3) THEN
          CALL BICKLEY_JET(X_RESOLVED, Y_RESOLVED, (T + H), VEL_LOC)

        ELSE IF (FLOW_MODEL == 4) THEN
          CALL LAMB_OSSEEN(X_RESOLVED, Y_RESOLVED, (T + H), VEL_LOC)

        ELSE IF (FLOW_MODEL == 5) THEN
          CALL ABC(X_RESOLVED, Y_RESOLVED, Z_RESOLVED, VEL_LOC)

        ELSE IF (FLOW_MODEL == 6) THEN
          CALL UNSTEADY_ABC(X_RESOLVED, Y_RESOLVED, Z_RESOLVED, (T + H), VEL_LOC)

        END IF

        K4  = H*VEL_LOC

        !!-----------------------------------------------------
        !! Update point coordinates using numerical integration
        !!-----------------------------------------------------
        !!------------------------------------------
        !! Check to see if point has left the domain
        !!------------------------------------------
        CALL RESOLVE_BOUNDARY_CONDITIONS(BC_TYPE, BOUNDING_BOX, X_TEMP, Y_TEMP, Z_TEMP, &
                                         (X_TEMP + (1.0D0/6.0D0)*(K1(1) + 2.0D0*K2(1) + 2.0D0*K3(1) + K4(1))), &
                                         (Y_TEMP + (1.0D0/6.0D0)*(K1(2) + 2.0D0*K2(2) + 2.0D0*K3(2) + K4(2))), &
                                         (Z_TEMP + (1.0D0/6.0D0)*(K1(3) + 2.0D0*K2(3) + 2.0D0*K3(3) + K4(3))), &
                                         X_RESOLVED, Y_RESOLVED, Z_RESOLVED, INT_STATUS_RESOLVED)
        
        !!----------------------------------------------
        !! Stop integration if point has left the domain
        !!----------------------------------------------  
        IF (INT_STATUS_RESOLVED == 1) THEN
          XYZ(1)     = X_RESOLVED
          XYZ(2)     = Y_RESOLVED
          XYZ(3)     = Z_RESOLVED
          INT_STATUS = 1
          RETURN

        ELSE
          !!---------------s--
          !! Periodic BC (+x)
          !!-----------------
          IF (BC_TYPE == 3) THEN
            X_TEMP = X_TEMP + (1.0D0/6.0D0)*(K1(1) + 2.0D0*K2(1) + 2.0D0*K3(1) + K4(1))
            Y_TEMP = Y_RESOLVED
            Z_TEMP = Z_RESOLVED

            ! Temporary workaround for single time step loops (T0 + DT = T1)
            XYZ(1) = X_TEMP + (1.0D0/6.0D0)*(K1(1) + 2.0D0*K2(1) + 2.0D0*K3(1) + K4(1))
            XYZ(2) = Y_RESOLVED
            XYZ(3) = Z_RESOLVED

          ELSE
            X_TEMP = X_RESOLVED
            Y_TEMP = Y_RESOLVED
            Z_TEMP = Z_RESOLVED
            
            ! Temporary workaround for single time step loops (T0 + DT = T1)
            XYZ(1) = X_RESOLVED
            XYZ(2) = Y_RESOLVED
            XYZ(3) = Z_RESOLVED
          
          END IF

       END IF

       T = T + H

      END DO
    
    END SUBROUTINE RK4_ANALYTICAL

END MODULE INTEGRATION
