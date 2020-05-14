!!----------------------------------------------------------------
!! Fortran 90/95 style module to handle advection and lagrangian 
!! field computations directly from Fortran using NumPy generated 
!! arrays and pre-processing of VTK files from Python.
!!
!!
!! Author:   Zachariah Irwin
!!           University of Colorado, Boulder
!! Version:  March 2020
!!----------------------------------------------------------------
MODULE MAIN

  USE INTEGRATION
  USE CGFIELDS
  USE OUTPUTFORMAT

  IMPLICIT NONE

  CONTAINS

    !!------------------------------------------------------------------------------------------------------------
    !> @brief Performs advection and Lagrangian analysis for cartesian grids of tracers given a simulated velocity
    !!        field
    !!
    !! LAG_FIELDS_MESH handles all computational subroutines for every tracer in a cartesian grid. Velocity data
    !! is loaded from binary files into an array where each velocity field for a given point in time can be easily
    !! accessed from memory. Points are advected and Lagrangian analysis is performed at every time step. See
    !! integration.f90 and cgFields.f90 for further details.
    !!
    !! @param[in, out]   FTLE_T               double array, time-scaled FTLE
    !!                                        size: NUM_P
    !! @param[in, out]   FTLE_NO_T            double array, non-time-scaled FTLE
    !!                                        size: NUM_P
    !! @param[in, out]   STRETCH_1|2|3        double array, elongational stretch along unit vectors 1|2|3
    !!                                        size: NUM_P
    !! @param[in, out]   STRAIN_1|2|3         double array, Green-Lagrange strain along unit vectors 1|2|3
    !!                                        size: NUM_P
    !! @param[in, out]   CAUCHY_GREEN         double array, Right Cauchy-Green matrix
    !!                                        size: NUM_P x 3 x 3
    !! @param[in]        NUM_P                integer, number of tracers in the simulation
    !! @param[in]        NUM_NODES            integer, number of nodes in the mesh
    !! @param[in]        NUM_CELLS            integer, number of cells in the mesh  
    !! @param[in]        INT_FLAGS            integer array, flags corresponding to advection parameters
    !!                                        size: 3
    !! @param[in]        LAG_FLAGS            integer array, flags corresponding to which Lagrangian fields
    !!                                        are to be computed
    !!                                        size: 3
    !! @param[in]        FILE_ID_IN           integer array, start/stop/delta of velocity file indices
    !!                                        size: 3
    !! @param[in]        REF_NUMS             integer array, number of points along each direction
    !!                                        size: 3
    !! @param[in]        CELL_IDS_INIT        integer array, initial cell IDs for every point in the simulation
    !!                                        size: NUM_P
    !! @param[in]        CONNECTIVITY         integer array, mesh connectivity matrix
    !!                                        size: NUM_CELLS x 4
    !! @param[in]        ADJACENCY            integer array, mesh adjacency matrix
    !!                                        size: NUM_CELLS x 4
    !! @param[in]        REF_DIMNS            double array, grid spacing in referential configuration
    !!                                        size: 3
    !! @param[in]        SIM_ TIMES           double array, simulation times configured in input file
    !!                                        size: 3
    !! @param[in]        FILE_IN_TIMING       double array, start/stop/delta of velocity file times
    !!                                        size: 3
    !! @param[in]        FLOW_BBOX            double array, bounding box of flow domain
    !!                                        size: 6
    !! @param[in]        TRACER_BBOX          double array, bounding box for seeding tracers
    !!                                        size: 6
    !! @param[in]        UNIT_VECTORS         double array, unit vectors along user-specified directions 1,2,3
    !!                                        size: 3x3
    !! @param[in]        COORDINATES          double array, mesh node coordinates
    !!                                        size: NUM_NODES x 3
    !! @param[in]        FILE_IN_PREFIX       character, velocity file prefix, e.g. "data_vel"
    !!                                        len: 1024
    !!------------------------------------------------------------------------------------------------------------
    SUBROUTINE LAG_FIELDS_MESH(INT_FLAGS, LAG_FLAGS, NUM_P, REF_DIMNS, REF_NUMS, FLOW_BBOX, TRACER_BBOX,&
                               SIM_TIMES, FILE_ID_IN, FILE_IN_TIMING, FILE_IN_PREFIX, NUM_NODES, NUM_CELLS,&
                               CONNECTIVITY, ADJACENCY, COORDINATES, UNIT_VECTORS, FTLE_T, FTLE_NO_T,&
                               STRETCH_1, STRETCH_2, STRETCH_3, STRAIN_1, STRAIN_2, STRAIN_3, CAUCHY_GREEN,&
                               CELL_IDS_INIT)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: NUM_P, NUM_NODES, NUM_CELLS
      INTEGER, DIMENSION(3), INTENT(IN):: INT_FLAGS, LAG_FLAGS, FILE_ID_IN, REF_NUMS
      INTEGER, DIMENSION(0:NUM_CELLS-1, 4), INTENT(IN):: CONNECTIVITY, ADJACENCY
      INTEGER, DIMENSION(0:NUM_P-1), INTENT(IN):: CELL_IDS_INIT
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: REF_DIMNS, FILE_IN_TIMING, SIM_TIMES
      DOUBLE PRECISION, DIMENSION(6), INTENT(IN):: FLOW_BBOX, TRACER_BBOX
      DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN):: UNIT_VECTORS
      DOUBLE PRECISION, DIMENSION(0:NUM_NODES-1, 3), INTENT(IN):: COORDINATES
      CHARACTER(LEN=1024), INTENT(IN):: FILE_IN_PREFIX

      DOUBLE PRECISION, DIMENSION(0:NUM_P-1), INTENT(INOUT):: FTLE_T, FTLE_NO_T,&
                                                              STRAIN_1, STRAIN_2, STRAIN_3,&
                                                              STRETCH_1, STRETCH_2, STRETCH_3                                                   
      DOUBLE PRECISION, DIMENSION(0:NUM_P-1, 3, 3), INTENT(INOUT):: CAUCHY_GREEN

      LOGICAL:: IS_LOAD_FRAME, IS_STEADY_FLOW
      INTEGER:: I, DIMN, TIME_INDEX, FILE_ID_START, FILE_ID_STOP, FILE_ID_DELTA, NUM_FILES, NUM_T
      INTEGER, DIMENSION(0:NUM_P-1):: CELL_IDS, LEFT_GRID, HAS_CG
      INTEGER, DIMENSION(:,:), ALLOCATABLE:: DATA_INDEX_WINDOWS
      DOUBLE PRECISION:: T0, T1, DT, SIM_TIME, FILE_IN_START, FILE_IN_STOP, FILE_IN_DELTA
      DOUBLE PRECISION, DIMENSION(3):: WINDOW
      DOUBLE PRECISION, DIMENSION(0:NUM_NODES-1, 3):: VEL_LOW, VEL_MID, VEL_UP
      DOUBLE PRECISION, DIMENSION(0:NUM_P-1, 3):: PARTICLES_XYZ, PARTICLES_XYZ_OLD
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE:: DATA_TIME_WINDOWS
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE:: VEL_BIG

      DOUBLE PRECISION, PARAMETER:: EPS = 1.0e-12

      !f2py depend(NUM_CELLS) CONNECTIVITY, ADJACENCY
      !f2py depend(NUM_NODES) COORDINATES, VEL_LOW, VEL_MID, VEL_UP, VEL_BIG
      !f2py depend(NUM_P) PARTICLES_XYZ, PARTICLES_XYZ_OLD, CELL_IDS, CELL_IDS_INIT, LEFT_GRID, CAUCHY_GREEN, FTLE_T, FTLE_NO_T, STRETCH_1, STRETCH_2, STRETCH_3, STRAIN_1, STRAIN_2, STRAIN_3
      
      !!----------------------
      !! Initialize parameters
      !!----------------------
      DIMN           = INT_FLAGS(1)

      FILE_ID_START  = FILE_ID_IN(1)
      FILE_ID_STOP   = FILE_ID_IN(2)
      FILE_ID_DELTA  = FILE_ID_IN(3)

      FILE_IN_START  = FILE_IN_TIMING(1)
      FILE_IN_STOP   = FILE_IN_TIMING(2)
      FILE_IN_DELTA  = FILE_IN_TIMING(3)

      T0             = SIM_TIMES(1)
      T1             = SIM_TIMES(2)
      DT             = SIGN(SIM_TIMES(3), T1-T0)

      SIM_TIME       = T0
      TIME_INDEX     = 0

      ! Velocity files should be sequential (for now)
      IF (FILE_ID_DELTA /= 1 .AND. FILE_ID_DELTA /= 0) THEN
        print*, "ERROR: FILE DELTA != 1"
        STOP
      END IF

      !!----------------------------------------------------------
      !! If only one velocity file is provided, assume steady flow
      !!----------------------------------------------------------
      IS_STEADY_FLOW = (FILE_ID_DELTA == 0)
      
      IF (IS_STEADY_FLOW) THEN

        !!------------------------------
        !! Load velocity data into array
        !!------------------------------
        ALLOCATE(VEL_BIG(NUM_FILES, NUM_NODES, 3))

        CALL LOAD_VELOCITY(FILE_ID_START, NUM_NODES, FILE_IN_PREFIX, VEL_BIG(0,:,:))
        
        !!----------------------
        !! Fix data window times
        !!----------------------
        WINDOW(1) = 0.0
        WINDOW(2) = 1.0
        WINDOW(3) = 1.0

        !!---------------------------------------
        !! All velocities (across time) are fixed
        !!---------------------------------------
        VEL_LOW   = VEL_BIG(0,:,:)
        VEL_MID   = VEL_LOW
        VEL_UP    = VEL_LOW
      
      !!-----------------------------------------------------
      !! If more than one velocity file, assume unsteady flow
      !!-----------------------------------------------------
      ELSE 
        NUM_FILES = INT((FILE_ID_STOP - FILE_ID_START)/FILE_ID_DELTA)+1

        !!------------------------------------------------------------
        !! Load velocity data into array that can be indexed with time
        !!------------------------------------------------------------
        ALLOCATE(VEL_BIG(0:NUM_FILES-1, NUM_NODES, 3))

        DO I = 0, NUM_FILES-1
          CALL LOAD_VELOCITY(I, NUM_NODES, FILE_IN_PREFIX, VEL_BIG(I,:,:))
        END DO
      
      END IF

      NUM_T = INT((T1 - T0)/DT)

      !!--------------------------------------------------------
      !! Sync velocity file indices to time points in simulation
      !!--------------------------------------------------------
      ALLOCATE(DATA_INDEX_WINDOWS(0:NUM_T,3), DATA_TIME_WINDOWS(0:NUM_T,3))

      CALL MAKE_DATA_WINDOWS(NUM_FILES, FILE_ID_START, FILE_ID_STOP, FILE_ID_DELTA,&
                             NUM_T, FILE_IN_START, FILE_IN_STOP, FILE_IN_DELTA,&
                             T0, T1, DT, DATA_INDEX_WINDOWS, DATA_TIME_WINDOWS)

      !!--------------------------------
      !! Seed points in a cartesian grid
      !!--------------------------------
      CALL SEED_POINTS_CARTESIAN(NUM_P, REF_DIMNS, TRACER_BBOX, PARTICLES_XYZ)

      !!-----------------------------------------------------
      !! Assign point cell IDs to initial IDs
      !! Assume all points are in the domain
      !! Assume all points' CG tensor has not been calculated
      !!-----------------------------------------------------
      CELL_IDS   = CELL_IDS_INIT
      LEFT_GRID  = -1
      HAS_CG     = -1

      !!-----------------
      !! Begin simulation
      !!-----------------
      DO WHILE (ABS(T1-SIM_TIME) > EPS)

        !!-----------------------------------------------------------------------------------
        !! Initialize SIM_TIME tracer coordinates for Lagrangian module if points leave early
        !!-----------------------------------------------------------------------------------
        PARTICLES_XYZ_OLD = PARTICLES_XYZ

        !!--------------------------------------------------------------
        !! Determine if new velocity data needs to be loaded into memory
        !!--------------------------------------------------------------
        IS_LOAD_FRAME  = ((SIM_TIME == T0) .OR. &
                          ((DATA_TIME_WINDOWS(TIME_INDEX, 1) .NE. WINDOW(1)) .OR. &
                           (DATA_TIME_WINDOWS(TIME_INDEX, 2) .NE. WINDOW(2)) .OR. &
                           (DATA_TIME_WINDOWS(TIME_INDEX, 3) .NE. WINDOW(3)))&
                           .AND. (IS_STEADY_FLOW .EQV. .FALSE.))

        !!---------------------------------
        !! Load velocity fields if required
        !!---------------------------------
        IF (IS_LOAD_FRAME .AND. (IS_STEADY_FLOW .EQV. .FALSE.)) THEN

          ! Temporary fix for periodic data sets
          IF (SIM_TIME > (T1-T0)) THEN
            WINDOW(1) = DATA_TIME_WINDOWS(TIME_INDEX, 1) + (T1 - T0)
            WINDOW(2) = DATA_TIME_WINDOWS(TIME_INDEX, 2) + (T1 - T0)
            WINDOW(3) = DATA_TIME_WINDOWS(TIME_INDEX, 3) + (T1 - T0)
          ELSE
            WINDOW(1) = DATA_TIME_WINDOWS(TIME_INDEX, 1)
            WINDOW(2) = DATA_TIME_WINDOWS(TIME_INDEX, 2)
            WINDOW(3) = DATA_TIME_WINDOWS(TIME_INDEX, 3)
          END IF

          VEL_LOW = VEL_BIG(DATA_INDEX_WINDOWS(TIME_INDEX, 1), :, :)
          VEL_MID = VEL_BIG(DATA_INDEX_WINDOWS(TIME_INDEX, 2), :, :)
          VEL_UP  = VEL_BIG(DATA_INDEX_WINDOWS(TIME_INDEX, 3), :, :)
        END IF

        !!------------------------------------------------
        !! Advect tracers from SIM_TIME to (SIM_TIME + DT)
        !!------------------------------------------------
        CALL ADVECT_GRID_DATA(INT_FLAGS, FLOW_BBOX, NUM_P, PARTICLES_XYZ, CELL_IDS, LEFT_GRID,&
                              [SIM_TIME, DT], WINDOW, NUM_NODES, NUM_CELLS, CONNECTIVITY, ADJACENCY,&
                              COORDINATES, VEL_LOW, VEL_MID, VEL_UP)
        
        !!----------------------------
        !! Perform Lagrangian analysis
        !!----------------------------
        CALL LAGRANGIAN_MAIN(DIMN, NUM_P, REF_NUMS, REF_DIMNS, [T0, T1, SIM_TIME, SIM_TIME+DT],&
                             PARTICLES_XYZ_OLD, PARTICLES_XYZ, LAG_FLAGS, HAS_CG, LEFT_GRID,&
                             UNIT_VECTORS, CAUCHY_GREEN, FTLE_T, FTLE_NO_T,&
                             STRETCH_1, STRETCH_2, STRETCH_3, STRAIN_1, STRAIN_2, STRAIN_3)

        !!----------------------------------------------
        !! Update SIM_TIME, TIME_INDEX to next time step
        !!----------------------------------------------
        SIM_TIME   = SIM_TIME + DT
        TIME_INDEX = TIME_INDEX + 1

      END DO

      print*, "Deallocating arrays..."
      DEALLOCATE(DATA_INDEX_WINDOWS)
      DEALLOCATE(DATA_TIME_WINDOWS)
      DEALLOCATE(VEL_BIG)
      print*, "Finished Deallocating arrays."

    END SUBROUTINE LAG_FIELDS_MESH

    !!------------------------------------------------------------------------------------------------------------
    !> @brief Performs advection and Lagrangian analysis for cartesian grids given an analytical velocity field.
    !!
    !! LAG_FIELDS_ANLYT handles all computational subroutines for every point in a cartesian grid for analytical
    !! velocity fields ONLY. Points are advected and Lagrangian analysis is performed at every time step. 
    !! See integration.f90 and cgFields.f90 for further details.
    !!  
    !! @param[in, out]   FTLE_T               double array, time-scaled FTLE
    !!                                        size: NUM_P
    !! @param[in, out]   FTLE_NO_T            double array, non-time-scaled FTLE
    !!                                        size: NUM_P
    !! @param[in, out]   STRETCH_1|2|3        double array, elongational stretch along unit vectors 1|2|3
    !!                                        size: NUM_P
    !! @param[in, out]   STRAIN_1|2|3         double array, Green-Lagrange strain along unit vectors 1|2|3
    !!                                        size: NUM_P
    !! @param[in, out]   CAUCHY_GREEN         double array, Right Cauchy-Green matrix
    !!                                        size: NUM_P x 3 x 3
    !! @param[in]        FLOW_MODEL           integer, flag corresponding to relevant analytical velocity field
    !! @param[in]        NUM_P                integer, number of points in the simulation
    !! @param[in]        INT_FLAGS            integer array, flags corresponding to advection parameters
    !!                                        size: 3
    !! @param[in]        LAG_FLAGS            integer array, flags corresponding to which Lagrangian fields
    !!                                        are to be computed
    !!                                        size: 3
    !! @param[in]        REF_NUMS             integer array, number of points along each dimension
    !!                                        size: 3
    !! @param[in]        REF_DIMNS            double array, grid spacing in referential configuration
    !!                                        size: 3
    !! @param[in]        TIMES                double array, FTLE window start/stop and simulation times
    !!                                        size: 3
    !! @param[in]        FLOW_BBOX            double array, bounding box of flow domain
    !!                                        size: 6
    !! @param[in]        TRACER_BBOX          double array, bounding box for seeding tracers
    !!                                        size: 6
    !! @param[in]        UNIT_VECTORS         double array, unit vectors along user-specified directions 1,2,3
    !!                                        size: 3x3
    !!------------------------------------------------------------------------------------------------------------
    SUBROUTINE LAG_FIELDS_ANLYT(FLOW_MODEL, INT_FLAGS, LAG_FLAGS, NUM_P, REF_DIMNS, REF_NUMS, FLOW_BBOX, &
                                TRACER_BBOX, SIM_TIMES, UNIT_VECTORS, FTLE_T, FTLE_NO_T,  STRETCH_1, &
                                STRETCH_2, STRETCH_3, STRAIN_1, STRAIN_2, STRAIN_3, CAUCHY_GREEN)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: FLOW_MODEL, NUM_P
      INTEGER, DIMENSION(3), INTENT(IN):: INT_FLAGS, LAG_FLAGS, REF_NUMS
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: REF_DIMNS, SIM_TIMES
      DOUBLE PRECISION, DIMENSION(6), INTENT(IN):: FLOW_BBOX, TRACER_BBOX
      DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN):: UNIT_VECTORS

      DOUBLE PRECISION, DIMENSION(0:NUM_P-1), INTENT(INOUT):: FTLE_T, FTLE_NO_T,&
                                                              STRAIN_1, STRAIN_2, STRAIN_3,&
                                                              STRETCH_1, STRETCH_2, STRETCH_3                                                   
      DOUBLE PRECISION, DIMENSION(0:NUM_P-1, 3, 3), INTENT(INOUT):: CAUCHY_GREEN

      INTEGER:: I, DIMN
      INTEGER, DIMENSION(0:NUM_P-1):: LEFT_GRID, HAS_CG
      DOUBLE PRECISION:: T0, T1, DT, SIM_TIME
      DOUBLE PRECISION, DIMENSION(3):: WINDOW
      DOUBLE PRECISION, DIMENSION(0:NUM_P-1, 3):: PARTICLES_XYZ, PARTICLES_XYZ_OLD

      DOUBLE PRECISION, PARAMETER:: EPS = 1.0e-12

      !f2py depend(NUM_P) PARTICLES_XYZ, PARTICLES_XYZ_OLD, LEFT_GRID, CAUCHY_GREEN, FTLE_T, FTLE_NO_T, STRETCH_1, STRETCH_2, STRETCH_3, STRAIN_1, STRAIN_2, STRAIN_3
      
      !!----------------------
      !! Initialize parameters
      !!----------------------
      DIMN       = INT_FLAGS(1)

      T0         = SIM_TIMES(1)
      T1         = SIM_TIMES(2)
      DT         = SIGN(SIM_TIMES(3), T1-T0)

      SIM_TIME   = T0

      !!--------------------------------
      !! Seed points in a cartesian grid
      !!--------------------------------
      CALL SEED_POINTS_CARTESIAN(NUM_P, REF_DIMNS, TRACER_BBOX, PARTICLES_XYZ)

      !!-----------------------------------------------------
      !! Assume all points are in the domain
      !! Assume all points' CG tensor has not been calculated
      !!-----------------------------------------------------
      LEFT_GRID  = -1
      HAS_CG     = -1

      !!-----------------
      !! Begin simulation
      !!-----------------
      DO WHILE (ABS(T1 - SIM_TIME) > EPS)
        
        !!-----------------------------------------------------------------------------------
        !! Initialize SIM_TIME tracer coordinates for Lagrangian module if points leave early
        !!-----------------------------------------------------------------------------------
        PARTICLES_XYZ_OLD = PARTICLES_XYZ
        
        !!------------------------------------------------
        !! Advect tracers from SIM_TIME to (SIM_TIME + DT)
        !!------------------------------------------------
        CALL ADVECT_GRID_ANALYTICAL(FLOW_MODEL, INT_FLAGS, TRACER_BBOX, NUM_P, PARTICLES_XYZ,&
                                    LEFT_GRID, [SIM_TIME, DT])
        
        !!----------------------------
        !! Perform Lagrangian analysis
        !!----------------------------
        CALL LAGRANGIAN_MAIN(DIMN, NUM_P, REF_NUMS, REF_DIMNS, [T0, T1, SIM_TIME, SIM_TIME+DT],&
                             PARTICLES_XYZ_OLD, PARTICLES_XYZ, LAG_FLAGS, HAS_CG, LEFT_GRID,&
                             UNIT_VECTORS, CAUCHY_GREEN, FTLE_T, FTLE_NO_T, STRETCH_1, STRETCH_2,&
                             STRETCH_3, STRAIN_1, STRAIN_2, STRAIN_3)
        
        !!---------------------------------
        ! Update SIM_TIME to next time step
        !!---------------------------------
        SIM_TIME = SIM_TIME + DT

      END DO

    END SUBROUTINE LAG_FIELDS_ANLYT

    !!------------------------------------------------------------------------------------------------------------
    !> @brief Performs advection and with Fortran-based I/O for cartesian grids of tracers given a simulated 
    !! velocity field
    !!
    !!
    !! @param[in]        NUM_P                integer, number of tracers in the simulation
    !! @param[in]        NUM_NODES            integer, number of nodes in the mesh
    !! @param[in]        NUM_CELLS            integer, number of cells in the mesh  
    !! @param[in]        NUM_INJECTIONS       integer, number of tracer injections into the domain for the
    !!                                        length of the entire simulation
    !! @param[in]        INT_FLAGS            integer array, flags corresponding to advection parameters
    !!                                        size: 3
    !! @param[in]        FILE_ID_IN           integer array, start/stop/delta of velocity file indices
    !!                                        size: 3
    !! @param[in]        FILE_ID_OUT          integer array, start/stop/delta of vtp output files
    !!                                        size: 3
    !! @param[in]        REF_NUMS             integer array, number of points along each direction
    !!                                        size: 3
    !! @param[in]        CELL_IDS_INIT        integer array, initial cell IDs for every point in the simulation
    !!                                        size: NUM_P
    !! @param[in]        CONNECTIVITY         integer array, mesh connectivity matrix
    !!                                        size: NUM_CELLS x 4
    !! @param[in]        ADJACENCY            integer array, mesh adjacency matrix
    !!                                        size: NUM_CELLS x 4
    !! @param[in]        FILE_OUT_TIMING      double array, start/delta of vtp output files
    !!                                        size: 2
    !! @param[in]        REF_DIMNS            double array, grid spacing in referential configuration
    !!                                        size: 3
    !! @param[in]        TIMES                double array, simulation times and injection times
    !!                                        configured in input file
    !!                                        size: 5
    !! @param[in]        FILE_IN_TIMING       double array, start/stop/delta of velocity file times
    !!                                        size: 3
    !! @param[in]        FLOW_BBOX            double array, bounding box of flow domain
    !!                                        size: 6
    !! @param[in]        TRACER_BBOX          double array, bounding box for seeding tracers
    !!                                        size: 6
    !! @param[in]        COORDINATES          double array, mesh node coordinates
    !!                                        size: NUM_NODES x 3
    !! @param[in]        FILE_IN_PREFIX       character, velocity file prefix, e.g. "data_vel"
    !!                                        len: 1024
    !! @param[in]        FILE_OUT_PREFIX      character, vtp output file prefix
    !!                                        len: 1024
    !!------------------------------------------------------------------------------------------------------------
    SUBROUTINE ADVECT_SEED_CART(INT_FLAGS, FLOW_BBOX, TRACER_BBOX, NUM_P,&
                                NUM_INJECTIONS, REF_DIMNS, TIMES, FILE_ID_IN, FILE_IN_TIMING,&
                                FILE_IN_PREFIX, FILE_ID_OUT, FILE_OUT_TIMING, FILE_OUT_PREFIX,&
                                NUM_NODES, NUM_CELLS, CONNECTIVITY, ADJACENCY, COORDINATES,&
                                CELL_IDS_INIT)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: NUM_P, NUM_NODES, NUM_CELLS, NUM_INJECTIONS
      INTEGER, DIMENSION(3), INTENT(IN):: INT_FLAGS, FILE_ID_IN, FILE_ID_OUT
      INTEGER, DIMENSION(0:NUM_CELLS-1, 4), INTENT(IN):: CONNECTIVITY, ADJACENCY
      INTEGER, DIMENSION(0:NUM_P-1), INTENT(IN):: CELL_IDS_INIT
      DOUBLE PRECISION, DIMENSION(2), INTENT(IN):: FILE_OUT_TIMING
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: REF_DIMNS, FILE_IN_TIMING
      DOUBLE PRECISION, DIMENSION(5), INTENT(IN):: TIMES
      DOUBLE PRECISION, DIMENSION(6), INTENT(IN):: FLOW_BBOX, TRACER_BBOX
      DOUBLE PRECISION, DIMENSION(0:NUM_NODES-1, 3), INTENT(IN):: COORDINATES
      CHARACTER(LEN=1024), INTENT(IN):: FILE_IN_PREFIX, FILE_OUT_PREFIX

      LOGICAL:: IS_LOAD_FRAME, IS_STEADY_FLOW, STAT
      INTEGER:: I, J, DIMN, TIME_INDEX, NUM_FILES, FILE_ID_START, FILE_ID_STOP, FILE_ID_DELTA,&
                NUM_T, NUM_RELEASED, OUT_ID, DUMP_IDX, FILE_ID_OUT_START, FILE_ID_OUT_DELTA,&
                FILE_ID_OUT_STOP, STOP_SAVE, START_NEW, STOP_NEW, NUM_OUT_FILES, DLOC
      INTEGER, DIMENSION(:), ALLOCATABLE:: TRACERS_LEFT_GRID, TRACERS_LEFT_GRID_MVALLOC,&
                                           TRACERS_CELL_ID, TRACERS_CELL_ID_MVALLOC, OPF_PTS
      INTEGER, DIMENSION(:,:), ALLOCATABLE:: DATA_INDEX_WINDOWS
      DOUBLE PRECISION:: T_START, T_END, DT, SIM_TIME, INJECTION_DT, FILE_IN_START, FILE_IN_STOP,&
                         FILE_IN_DELTA, INJECTION_START, FILE_OUT_START, FILE_OUT_DELTA
      DOUBLE PRECISION, DIMENSION(3):: WINDOW
      DOUBLE PRECISION, DIMENSION(0:NUM_NODES-1, 3):: VEL_LOW, VEL_MID, VEL_UP
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: RELEASE_TIMES, DUMP_TIMES, RADII,&
                                                    RADII_MVALLOC
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: DATA_TIME_WINDOWS, TRACERS_XYZ,&
                                                      TRACERS_XYZ_MVALLOC, TRACERS_SEED
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE:: VEL_BIG
      CHARACTER(LEN=300):: VTKNAME, FILE_INDEX_STR, DTYPE
      CLASS(VTKmetaData), dimension(:), allocatable:: VTKMD

      DOUBLE PRECISION, PARAMETER:: EPS = 1.0e-12

      !f2py depend(NUM_CELLS) CONNECTIVITY, ADJACENCY
      !f2py depend(NUM_NODES) COORDINATES, VEL_LOW, VEL_MID, VEL_UP
      !f2py depend(NUM_P) CELL_IDS_INIT

      !!----------------------
      !! Initialize parameters
      !!----------------------
      DIMN              = INT_FLAGS(1)

      FILE_ID_START     = FILE_ID_IN(1)
      FILE_ID_STOP      = FILE_ID_IN(2)
      FILE_ID_DELTA     = FILE_ID_IN(3)

      FILE_IN_START     = FILE_IN_TIMING(1)
      FILE_IN_STOP      = FILE_IN_TIMING(2)
      FILE_IN_DELTA     = FILE_IN_TIMING(3)

      ! Velocity files should be seqeuntial, for now
      IF (FILE_ID_DELTA /= 1) THEN
        print*, "ERROR: FILE DELTA != 1"
        STOP
      END IF

      FILE_ID_OUT_START = FILE_ID_OUT(1)
      FILE_ID_OUT_STOP  = FILE_ID_OUT(2)
      FILE_ID_OUT_DELTA = FILE_ID_OUT(3)
      NUM_OUT_FILES     = INT((FILE_ID_OUT_STOP-FILE_ID_OUT_START)/FILE_ID_OUT_DELTA)+1
      OUT_ID            = FILE_ID_OUT_START
      DUMP_IDX          = 0

      FILE_OUT_START    = FILE_OUT_TIMING(1)
      FILE_OUT_DELTA    = FILE_OUT_TIMING(2)

      T_START           = TIMES(1)
      T_END             = TIMES(2)
      DT                = TIMES(3)
      INJECTION_START   = TIMES(4)
      INJECTION_DT      = TIMES(5)

      SIM_TIME          = T_START
      TIME_INDEX        = 0

      !!---------------------------------------
      !! See outputFormat.f90 for documentation
      !!---------------------------------------
      DTYPE             = 'DATASET POLYDATA'
      ALLOCATE(VTKMD(1), OPF_PTS(1))
      
      !!----------------------------------------------------------
      !! If only one velocity file is provided, assume steady flow
      !!----------------------------------------------------------
      IS_STEADY_FLOW = (FILE_ID_DELTA == 0)

      IF (IS_STEADY_FLOW) THEN

        NUM_FILES = 1
        
        !!------------------------------
        !! Load velocity data into array
        !!------------------------------
        ALLOCATE(VEL_BIG(NUM_FILES, NUM_NODES, 3))
        
        CALL LOAD_VELOCITY(FILE_ID_START, NUM_NODES, FILE_IN_PREFIX, VEL_BIG(FILE_ID_START,:,:))
        
        !!----------------------
        !! Fix data window times
        !!----------------------
        WINDOW(1) = 0.0
        WINDOW(2) = 1.0
        WINDOW(3) = 1.0

        !!---------------------------------------
        !! All velocities (across time) are fixed
        !!---------------------------------------
        VEL_LOW   = VEL_BIG(0,:,:)
        VEL_MID   = VEL_LOW
        VEL_UP    = VEL_LOW
      
      !!-----------------------------------------------------
      !! If more than one velocity file, assume unsteady flow
      !!-----------------------------------------------------
      ELSE
        NUM_FILES = INT((FILE_ID_STOP - FILE_ID_START)/FILE_ID_DELTA)+1
       
        !!------------------------------------------------------------
        !! Load velocity data into array that can be indexed with time
        !!------------------------------------------------------------
        ALLOCATE(VEL_BIG(0:NUM_FILES-1, NUM_NODES, 3))
        
        DO I = 0, NUM_FILES-1
          CALL LOAD_VELOCITY(I, NUM_NODES, FILE_IN_PREFIX, VEL_BIG(I,:,:))
        END DO
      
      END IF

      NUM_T = INT((T_END - T_START)/DT)

      !!--------------------------------------------------------
      !! Sync velocity file indices to time points in simulation
      !!--------------------------------------------------------
      ALLOCATE(DATA_INDEX_WINDOWS(0:NUM_T,3), DATA_TIME_WINDOWS(0:NUM_T,3))

      CALL MAKE_DATA_WINDOWS(NUM_FILES, FILE_ID_START, FILE_ID_STOP, FILE_ID_DELTA,&
                             NUM_T, FILE_IN_START, FILE_IN_STOP, FILE_IN_DELTA,&
                             T_START, T_END, DT, DATA_INDEX_WINDOWS, DATA_TIME_WINDOWS)
      
      !!-------------------------------------------
      !! Initialize array of tracer injection times
      !!-------------------------------------------
      ALLOCATE(RELEASE_TIMES(0:NUM_INJECTIONS-1))
      DO I = 0, NUM_INJECTIONS-1
        RELEASE_TIMES(I) = INJECTION_START + FLOAT(I)*INJECTION_DT
      END DO

      !!------------------------------------------------------------
      !! Initialize array of output times at which to write vtp data
      !!------------------------------------------------------------
      ALLOCATE(DUMP_TIMES(0:NUM_OUT_FILES-1))
      DO I = 0, NUM_OUT_FILES-1
        DUMP_TIMES(I) = FILE_OUT_START + FLOAT(I)*FILE_OUT_DELTA
      END DO


      !!--------------------------------
      !! Seed points in a cartesian grid
      !!--------------------------------
      ALLOCATE(TRACERS_SEED(0:NUM_P-1,3))
      CALL SEED_POINTS_CARTESIAN(NUM_P, REF_DIMNS, TRACER_BBOX, TRACERS_SEED)

      NUM_RELEASED  = 0

      !!-----------------
      !! Begin simulation
      !!-----------------
      DO WHILE (ABS(T_END-SIM_TIME) > EPS)

        !!--------------------------------------------------------------
        !! Determine if new velocity data needs to be loaded into memory
        !!--------------------------------------------------------------
        IS_LOAD_FRAME  = ((SIM_TIME == T_START) .OR. &
                          ((DATA_TIME_WINDOWS(TIME_INDEX, 1) .NE. WINDOW(1)) .OR. &
                           (DATA_TIME_WINDOWS(TIME_INDEX, 2) .NE. WINDOW(2)) .OR. &
                           (DATA_TIME_WINDOWS(TIME_INDEX, 3) .NE. WINDOW(3)))&
                           .AND. (IS_STEADY_FLOW .EQV. .FALSE.))

        !!---------------------------------
        !! Load velocity fields if required
        !!---------------------------------
        IF (IS_LOAD_FRAME .AND. (IS_STEADY_FLOW .EQV. .FALSE.)) THEN

          ! Temporary fix for periodic data sets
          IF (SIM_TIME > (T_END - T_START)) THEN
            WINDOW(1) = DATA_TIME_WINDOWS(TIME_INDEX, 1) + (T_END - T_START)
            WINDOW(2) = DATA_TIME_WINDOWS(TIME_INDEX, 2) + (T_END - T_START)
            WINDOW(3) = DATA_TIME_WINDOWS(TIME_INDEX, 3) + (T_END - T_START)
          ELSE
            WINDOW(1) = DATA_TIME_WINDOWS(TIME_INDEX, 1)
            WINDOW(2) = DATA_TIME_WINDOWS(TIME_INDEX, 2)
            WINDOW(3) = DATA_TIME_WINDOWS(TIME_INDEX, 3)
          END IF

          VEL_LOW = VEL_BIG(DATA_INDEX_WINDOWS(TIME_INDEX, 1), :, :)
          VEL_MID = VEL_BIG(DATA_INDEX_WINDOWS(TIME_INDEX, 2), :, :)
          VEL_UP  = VEL_BIG(DATA_INDEX_WINDOWS(TIME_INDEX, 3), :, :)
        END IF

        I = NUM_RELEASED

        !!------------------------------------------------
        !! Check to see if new tracers need to be released
        !!------------------------------------------------
        IF (I <= NUM_INJECTIONS - 1) THEN
          IF ((SIM_TIME >= (RELEASE_TIMES(I)-EPS)) .AND. (SIM_TIME <= (RELEASE_TIMES(I) + EPS))) THEN

            NUM_RELEASED = NUM_RELEASED + 1
            PRINT*, "Injecting", NUM_P, "particles at t = ", SIM_TIME

            ALLOCATE(TRACERS_XYZ_MVALLOC(0:NUM_P*NUM_RELEASED-1, 3))
            ALLOCATE(TRACERS_LEFT_GRID_MVALLOC(0:NUM_RELEASED*NUM_P-1))
            ALLOCATE(TRACERS_CELL_ID_MVALLOC(0:NUM_RELEASED*NUM_P-1))
            ALLOCATE(RADII_MVALLOC(NUM_P*NUM_RELEASED))

            !!------------------------------------------------------------------
            !! The following lines are for the first tracer injection wherein we
            !! do not need to save the data before overwriting with MOVE_ALLOC
            !!------------------------------------------------------------------
            IF (NUM_RELEASED == 1) THEN

              !!--------------------------------
              !! Copy initial tracer coordinates
              !!--------------------------------
              CALL MOVE_ALLOC(TRACERS_XYZ_MVALLOC, TRACERS_XYZ)
              TRACERS_XYZ = TRACERS_SEED

              !!-----------------------------
              !! Copy initial tracer cell IDs
              !!-----------------------------
              CALL MOVE_ALLOC(TRACERS_CELL_ID_MVALLOC, TRACERS_CELL_ID)
              TRACERS_CELL_ID = CELL_IDS_INIT

              !!-----------------------------------------
              !! Copy initial tracer integration statuses
              !!-----------------------------------------
              CALL MOVE_ALLOC(TRACERS_LEFT_GRID_MVALLOC, TRACERS_LEFT_GRID)
              TRACERS_LEFT_GRID = -1

              !!--------------------------------------------------------------
              !! Copy initial tracer radii
              !! NOTE: Finite-sized tracer integration with geometrical checks
              !!       has not yet been implemented. This line is present for
              !!       consistency with outputFormat.f90
              !!--------------------------------------------------------------
              CALL MOVE_ALLOC(RADII_MVALLOC, RADII)
              RADII = 0.0D0

            !!-------------------------------------------------------------
            !! The following lines are for the subsequent tracer injections
            !! wherein we do need to save the data before overwriting with 
            !! MOVE_ALLOC
            !!-------------------------------------------------------------
            ELSE

              !!------------------------------------------------------------
              !! Get IDs of last tracer released, first new tracer released,
              !! and last new tracer released, respectively.
              !!------------------------------------------------------------
              STOP_SAVE = NUM_P*(NUM_RELEASED-1)-1
              START_NEW = NUM_P*(NUM_RELEASED-1)
              STOP_NEW  = NUM_P*NUM_RELEASED-1

              !!-------------------------------------------------------------------
              !! Save previous coordinates, copy initial coordinates to new tracers
              !!-------------------------------------------------------------------
              TRACERS_XYZ_MVALLOC(0:STOP_SAVE, :) = TRACERS_XYZ
              CALL MOVE_ALLOC(TRACERS_XYZ_MVALLOC, TRACERS_XYZ)
              TRACERS_XYZ(START_NEW:STOP_NEW, :) = TRACERS_SEED

              !!-------------------------------------------------------------
              !! Save previous cell IDs, copy initial cell IDs to new tracers
              !!-------------------------------------------------------------
              TRACERS_CELL_ID_MVALLOC(0:STOP_SAVE) = TRACERS_CELL_ID
              CALL MOVE_ALLOC(TRACERS_CELL_ID_MVALLOC, TRACERS_CELL_ID)
              TRACERS_CELL_ID(START_NEW:STOP_NEW) = CELL_IDS_INIT

              !!--------------------------------------------------------------------
              !! Save previous integration statues, copy initial integration statues 
              !! to new tracers
              !!--------------------------------------------------------------------
              TRACERS_LEFT_GRID_MVALLOC(0:STOP_SAVE) = TRACERS_LEFT_GRID
              CALL MOVE_ALLOC(TRACERS_LEFT_GRID_MVALLOC, TRACERS_LEFT_GRID)
              TRACERS_LEFT_GRID(START_NEW:STOP_NEW) = -1

              !!---------------------------------
              !! No need to preserve tracer radii
              !!---------------------------------
              CALL MOVE_ALLOC(RADII_MVALLOC, RADII)
              RADII = 0.0D0

            END IF
          END IF
        END IF
        
        IF (NUM_RELEASED .GE. 1) THEN

          IF (DUMP_IDX > NUM_OUT_FILES - 1) THEN
            EXIT
          END IF

          !!---------------------------------------------
          !! If tracers have been released, the number
          !! of output files has not exceeded the desired
          !! number of output files, and the output file
          !! time matches the simulation time, write out
          !!---------------------------------------------
          IF ((SIM_TIME < DUMP_TIMES(DUMP_IDX) + EPS) .AND. &
              (SIM_TIME > DUMP_TIMES(DUMP_IDX) - EPS)) THEN

            WRITE(FILE_INDEX_STR, '(I0)') OUT_ID+DUMP_IDX
            VTKNAME = TRIM(FILE_OUT_PREFIX)//'.'//TRIM(FILE_INDEX_STR)//'.vtk'
            OPF_PTS(1) = NUM_P*NUM_RELEASED
            STAT = VTKMD(1)%createVTKMetaData(a_Npt=OPF_PTS, a_DType=DTYPE, a_Title = VTKNAME)
            
            CALL VTKParticleWriter(VTKMD(1), TRACERS_LEFT_GRID, TRACERS_XYZ, RADII, a_FnameWrite=VTKNAME)
            
            print*, "Wrote data to:", TRIM(VTKNAME)

            !!---------------------------------------------------------
            !! Update the number of output files that have been written
            !!---------------------------------------------------------
            DUMP_IDX = DUMP_IDX + 1

          END IF

          !!------------------------------------------------
          !! Advect tracers from SIM_TIME to (SIM_TIME + DT)
          !!------------------------------------------------
          CALL ADVECT_GRID_DATA(INT_FLAGS, FLOW_BBOX, NUM_P*NUM_RELEASED, TRACERS_XYZ,&
                                TRACERS_CELL_ID, TRACERS_LEFT_GRID, [SIM_TIME, DT], WINDOW,&
                                NUM_NODES, NUM_CELLS,&
                                CONNECTIVITY, ADJACENCY, COORDINATES, VEL_LOW,&
                                VEL_MID, VEL_UP)
        END IF

        !!----------------------------------------------
        !! Update SIM_TIME, TIME_INDEX to next time step
        !!----------------------------------------------
        SIM_TIME   = SIM_TIME + DT
        TIME_INDEX = TIME_INDEX + 1

      END DO
      print*, "Deallocating arrays..."
      DEALLOCATE(OPF_PTS)
      DEALLOCATE(VTKMD)
      DEALLOCATE(TRACERS_LEFT_GRID)
      DEALLOCATE(TRACERS_CELL_ID)
      DEALLOCATE(DATA_INDEX_WINDOWS)
      DEALLOCATE(RELEASE_TIMES)
      DEALLOCATE(DUMP_TIMES)
      DEALLOCATE(RADII)
      DEALLOCATE(DATA_TIME_WINDOWS)
      DEALLOCATE(TRACERS_SEED)
      DEALLOCATE(TRACERS_XYZ)
      DEALLOCATE(VEL_BIG)
      print*, "Finished Deallocating arrays."

    END SUBROUTINE ADVECT_SEED_CART

    !!------------------------------------------------------------------------------------------------------------
    !> @brief Performs advection and with Fortran-based I/O for VTK-specified positions of tracers given a 
    !! simulated velocity field
    !!
    !!
    !! @param[in]        NUM_P                integer, number of tracers in the simulation
    !! @param[in]        NUM_NODES            integer, number of nodes in the mesh
    !! @param[in]        NUM_CELLS            integer, number of cells in the mesh  
    !! @param[in]        NUM_INJECTIONS       integer, number of tracer injections into the domain for the
    !!                                        length of the entire simulation
    !! @param[in]        INT_FLAGS            integer array, flags corresponding to advection parameters
    !!                                        size: 3
    !! @param[in]        FILE_ID_IN           integer array, start/stop/delta of velocity file indices
    !!                                        size: 3
    !! @param[in]        FILE_ID_OUT          integer array, start/stop/delta of vtp output files
    !!                                        size: 3
    !! @param[in]        REF_NUMS             integer array, number of points along each direction
    !!                                        size: 3
    !! @param[in]        CELL_IDS_INIT        integer array, initial cell IDs for every point in the simulation
    !!                                        size: NUM_P
    !! @param[in]        CONNECTIVITY         integer array, mesh connectivity matrix
    !!                                        size: NUM_CELLS x 4
    !! @param[in]        ADJACENCY            integer array, mesh adjacency matrix
    !!                                        size: NUM_CELLS x 4
    !! @param[in]        FILE_OUT_TIMING      double array, start/delta of vtp output files
    !!                                        size: 2
    !! @param[in]        REF_DIMNS            double array, grid spacing in referential configuration
    !!                                        size: 3
    !! @param[in]        TIMES                double array, simulation times and injection times
    !!                                        configured in input file
    !!                                        size: 5
    !! @param[in]        FILE_IN_TIMING       double array, start/stop/delta of velocity file times
    !!                                        size: 3
    !! @param[in]        FLOW_BBOX            double array, bounding box of flow domain
    !!                                        size: 6
    !! @param[in]        COORDINATES          double array, mesh node coordinates
    !!                                        size: NUM_NODES x 3
    !! @param[in]        SEED_COORDINATES     double array, particle coordinates that were seeded along
    !!                                        a VTK polydata object in preProcess.py
    !!                                        size: NUM_P x 3
    !! @param[in]        FILE_IN_PREFIX       character, velocity file prefix, e.g. "data_vel"
    !!                                        len: 1024
    !! @param[in]        FILE_OUT_PREFIX      character, vtp output file prefix
    !!                                        len: 1024
    !!------------------------------------------------------------------------------------------------------------
    SUBROUTINE ADVECT_SEED_VTK(INT_FLAGS, FLOW_BBOX, SEED_COORDINATES, NUM_P,&
                                 NUM_INJECTIONS, TIMES, FILE_ID_IN, FILE_IN_TIMING,&
                                 FILE_IN_PREFIX, FILE_ID_OUT, FILE_OUT_TIMING, FILE_OUT_PREFIX,&
                                 NUM_NODES, NUM_CELLS, CONNECTIVITY, ADJACENCY, COORDINATES,&
                                 CELL_IDS_INIT)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: NUM_P, NUM_NODES, NUM_CELLS, NUM_INJECTIONS
      INTEGER, DIMENSION(3), INTENT(IN):: INT_FLAGS, FILE_ID_IN, FILE_ID_OUT
      INTEGER, DIMENSION(0:NUM_CELLS-1, 4), INTENT(IN):: CONNECTIVITY, ADJACENCY
      INTEGER, DIMENSION(0:NUM_P-1), INTENT(IN):: CELL_IDS_INIT
      DOUBLE PRECISION, DIMENSION(2), INTENT(IN):: FILE_OUT_TIMING
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: FILE_IN_TIMING
      DOUBLE PRECISION, DIMENSION(5), INTENT(IN):: TIMES
      DOUBLE PRECISION, DIMENSION(6), INTENT(IN):: FLOW_BBOX
      DOUBLE PRECISION, DIMENSION(0:NUM_P-1, 3), INTENT(IN):: SEED_COORDINATES
      DOUBLE PRECISION, DIMENSION(0:NUM_NODES-1, 3), INTENT(IN):: COORDINATES
      CHARACTER(LEN=1024), INTENT(IN):: FILE_IN_PREFIX, FILE_OUT_PREFIX

      LOGICAL:: IS_LOAD_FRAME, IS_STEADY_FLOW, STAT
      INTEGER:: I, J, DIMN, TIME_INDEX, NUM_FILES, FILE_ID_START, FILE_ID_STOP, FILE_ID_DELTA,&
                NUM_T, NUM_RELEASED, OUT_ID, DUMP_IDX, FILE_ID_OUT_START, FILE_ID_OUT_DELTA,&
                FILE_ID_OUT_STOP, STOP_SAVE, START_NEW, STOP_NEW, NUM_OUT_FILES, DLOC
      INTEGER, DIMENSION(:), ALLOCATABLE:: TRACERS_LEFT_GRID, TRACERS_LEFT_GRID_MVALLOC,&
                                           TRACERS_CELL_ID, TRACERS_CELL_ID_MVALLOC, OPF_PTS
      INTEGER, DIMENSION(:,:), ALLOCATABLE:: DATA_INDEX_WINDOWS
      DOUBLE PRECISION:: T_START, T_END, DT, SIM_TIME, INJECTION_DT, FILE_IN_START, FILE_IN_STOP,&
                         FILE_IN_DELTA, INJECTION_START, FILE_OUT_START, FILE_OUT_DELTA
      DOUBLE PRECISION, DIMENSION(3):: WINDOW
      DOUBLE PRECISION, DIMENSION(0:NUM_NODES-1, 3):: VEL_LOW, VEL_MID, VEL_UP
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: RELEASE_TIMES, DUMP_TIMES, RADII,&
                                                    RADII_MVALLOC
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: DATA_TIME_WINDOWS, TRACERS_XYZ,&
                                                      TRACERS_XYZ_MVALLOC, TRACERS_SEED
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE:: VEL_BIG
      CHARACTER(LEN=300):: VTKNAME, FILE_INDEX_STR, DTYPE
      CLASS(VTKmetaData), DIMENSION(:), ALLOCATABLE:: VTKMD

      DOUBLE PRECISION, PARAMETER:: EPS = 1.0e-12

      !f2py depend(NUM_CELLS) CONNECTIVITY, ADJACENCY
      !f2py depend(NUM_NODES) COORDINATES, VEL_LOW, VEL_MID, VEL_UP
      !f2py depend(NUM_P) CELL_IDS_INIT, SEED_COORDINATES

      !!----------------------
      !! Initialize parameters
      !!----------------------
      DIMN              = INT_FLAGS(1)

      FILE_ID_START     = FILE_ID_IN(1)
      FILE_ID_STOP      = FILE_ID_IN(2)
      FILE_ID_DELTA     = FILE_ID_IN(3)

      FILE_IN_START     = FILE_IN_TIMING(1)
      FILE_IN_STOP      = FILE_IN_TIMING(2)
      FILE_IN_DELTA     = FILE_IN_TIMING(3)

      ! Velocity files should be seqeuntial, for now
      IF (FILE_ID_DELTA /= 1) THEN
        print*, "ERROR: FILE DELTA != 1"
        STOP
      END IF

      FILE_ID_OUT_START = FILE_ID_OUT(1)
      FILE_ID_OUT_STOP  = FILE_ID_OUT(2)
      FILE_ID_OUT_DELTA = FILE_ID_OUT(3)
      NUM_OUT_FILES     = INT((FILE_ID_OUT_STOP-FILE_ID_OUT_START)/FILE_ID_OUT_DELTA)+1
      OUT_ID            = FILE_ID_OUT_START
      DUMP_IDX          = 0

      FILE_OUT_START    = FILE_OUT_TIMING(1)
      FILE_OUT_DELTA    = FILE_OUT_TIMING(2)

      T_START           = TIMES(1)
      T_END             = TIMES(2)
      DT                = TIMES(3)
      INJECTION_START   = TIMES(4)
      INJECTION_DT      = TIMES(5)

      SIM_TIME          = T_START
      TIME_INDEX        = 0

      !!---------------------------------------
      !! See outputFormat.f90 for documentation
      !!---------------------------------------
      DTYPE             = 'DATASET POLYDATA'
      ALLOCATE(VTKMD(1), OPF_PTS(1))
      
      !!----------------------------------------------------------
      !! If only one velocity file is provided, assume steady flow
      !!----------------------------------------------------------
      IS_STEADY_FLOW = (FILE_ID_DELTA == 0)

      IF (IS_STEADY_FLOW) THEN

        NUM_FILES = 1
        
        !!------------------------------
        !! Load velocity data into array
        !!------------------------------
        ALLOCATE(VEL_BIG(NUM_FILES, NUM_NODES, 3))
        
        CALL LOAD_VELOCITY(FILE_ID_START, NUM_NODES, FILE_IN_PREFIX, VEL_BIG(FILE_ID_START,:,:))
        
        !!----------------------
        !! Fix data window times
        !!----------------------
        WINDOW(1) = 0.0
        WINDOW(2) = 1.0
        WINDOW(3) = 1.0

        !!---------------------------------------
        !! All velocities (across time) are fixed
        !!---------------------------------------
        VEL_LOW   = VEL_BIG(0,:,:)
        VEL_MID   = VEL_LOW
        VEL_UP    = VEL_LOW
      
      !!-----------------------------------------------------
      !! If more than one velocity file, assume unsteady flow
      !!-----------------------------------------------------
      ELSE
        NUM_FILES = INT((FILE_ID_STOP - FILE_ID_START)/FILE_ID_DELTA)+1
       
        !!------------------------------------------------------------
        !! Load velocity data into array that can be indexed with time
        !!------------------------------------------------------------
        ALLOCATE(VEL_BIG(0:NUM_FILES-1, NUM_NODES, 3))
        
        DO I = 0, NUM_FILES-1
          CALL LOAD_VELOCITY(I, NUM_NODES, FILE_IN_PREFIX, VEL_BIG(I,:,:))
        END DO
      
      END IF

      NUM_T = INT((T_END - T_START)/DT)

      !!--------------------------------------------------------
      !! Sync velocity file indices to time points in simulation
      !!--------------------------------------------------------
      ALLOCATE(DATA_INDEX_WINDOWS(0:NUM_T,3), DATA_TIME_WINDOWS(0:NUM_T,3))

      CALL MAKE_DATA_WINDOWS(NUM_FILES, FILE_ID_START, FILE_ID_STOP, FILE_ID_DELTA,&
                             NUM_T, FILE_IN_START, FILE_IN_STOP, FILE_IN_DELTA,&
                             T_START, T_END, DT, DATA_INDEX_WINDOWS, DATA_TIME_WINDOWS)
      
      !!-------------------------------------------
      !! Initialize array of tracer injection times
      !!-------------------------------------------
      ALLOCATE(RELEASE_TIMES(0:NUM_INJECTIONS-1))
      DO I = 0, NUM_INJECTIONS-1
        RELEASE_TIMES(I) = INJECTION_START + FLOAT(I)*INJECTION_DT
      END DO

      !!------------------------------------------------------------
      !! Initialize array of output times at which to write vtp data
      !!------------------------------------------------------------
      ALLOCATE(DUMP_TIMES(0:NUM_OUT_FILES-1))
      DO I = 0, NUM_OUT_FILES-1
        DUMP_TIMES(I) = FILE_OUT_START + FLOAT(I)*FILE_OUT_DELTA
      END DO


      !!---------------------------------------------
      !! Seed points according to VTK polydata object
      !!---------------------------------------------
      ALLOCATE(TRACERS_SEED(0:NUM_P-1,3))
      TRACERS_SEED = SEED_COORDINATES

      NUM_RELEASED  = 0

      !!-----------------
      !! Begin simulation
      !!-----------------
      DO WHILE (ABS(T_END-SIM_TIME) > EPS)

        !!--------------------------------------------------------------
        !! Determine if new velocity data needs to be loaded into memory
        !!--------------------------------------------------------------
        IS_LOAD_FRAME  = ((SIM_TIME == T_START) .OR. &
                          ((DATA_TIME_WINDOWS(TIME_INDEX, 1) .NE. WINDOW(1)) .OR. &
                           (DATA_TIME_WINDOWS(TIME_INDEX, 2) .NE. WINDOW(2)) .OR. &
                           (DATA_TIME_WINDOWS(TIME_INDEX, 3) .NE. WINDOW(3)))&
                           .AND. (IS_STEADY_FLOW .EQV. .FALSE.))

        !!---------------------------------
        !! Load velocity fields if required
        !!---------------------------------
        IF (IS_LOAD_FRAME .AND. (IS_STEADY_FLOW .EQV. .FALSE.)) THEN

          ! Temporary fix for periodic data sets
          IF (SIM_TIME > (T_END - T_START)) THEN
            WINDOW(1) = DATA_TIME_WINDOWS(TIME_INDEX, 1) + (T_END - T_START)
            WINDOW(2) = DATA_TIME_WINDOWS(TIME_INDEX, 2) + (T_END - T_START)
            WINDOW(3) = DATA_TIME_WINDOWS(TIME_INDEX, 3) + (T_END - T_START)
          ELSE
            WINDOW(1) = DATA_TIME_WINDOWS(TIME_INDEX, 1)
            WINDOW(2) = DATA_TIME_WINDOWS(TIME_INDEX, 2)
            WINDOW(3) = DATA_TIME_WINDOWS(TIME_INDEX, 3)
          END IF

          VEL_LOW = VEL_BIG(DATA_INDEX_WINDOWS(TIME_INDEX, 1), :, :)
          VEL_MID = VEL_BIG(DATA_INDEX_WINDOWS(TIME_INDEX, 2), :, :)
          VEL_UP  = VEL_BIG(DATA_INDEX_WINDOWS(TIME_INDEX, 3), :, :)
        END IF

        I = NUM_RELEASED

        !!------------------------------------------------
        !! Check to see if new tracers need to be released
        !!------------------------------------------------
        IF (I <= NUM_INJECTIONS - 1) THEN
          IF ((SIM_TIME >= (RELEASE_TIMES(I)-EPS)) .AND. (SIM_TIME <= (RELEASE_TIMES(I) + EPS))) THEN

            NUM_RELEASED = NUM_RELEASED + 1
            PRINT*, "Injecting", NUM_P, "particles at t = ", SIM_TIME

            ALLOCATE(TRACERS_XYZ_MVALLOC(0:NUM_P*NUM_RELEASED-1, 3))
            ALLOCATE(TRACERS_LEFT_GRID_MVALLOC(0:NUM_RELEASED*NUM_P-1))
            ALLOCATE(TRACERS_CELL_ID_MVALLOC(0:NUM_RELEASED*NUM_P-1))
            ALLOCATE(RADII_MVALLOC(NUM_P*NUM_RELEASED))

            !!------------------------------------------------------------------
            !! The following lines are for the first tracer injection wherein we
            !! do not need to save the data before overwriting with MOVE_ALLOC
            !!------------------------------------------------------------------
            IF (NUM_RELEASED == 1) THEN

              !!--------------------------------
              !! Copy initial tracer coordinates
              !!--------------------------------
              CALL MOVE_ALLOC(TRACERS_XYZ_MVALLOC, TRACERS_XYZ)
              TRACERS_XYZ = TRACERS_SEED

              !!-----------------------------
              !! Copy initial tracer cell IDs
              !!-----------------------------
              CALL MOVE_ALLOC(TRACERS_CELL_ID_MVALLOC, TRACERS_CELL_ID)
              TRACERS_CELL_ID = CELL_IDS_INIT

              !!-----------------------------------------
              !! Copy initial tracer integration statuses
              !!-----------------------------------------
              CALL MOVE_ALLOC(TRACERS_LEFT_GRID_MVALLOC, TRACERS_LEFT_GRID)
              TRACERS_LEFT_GRID = -1

              !!--------------------------------------------------------------
              !! Copy initial tracer radii
              !! NOTE: Finite-sized tracer integration with geometrical checks
              !!       has not yet been implemented. This line is present for
              !!       consistency with outputFormat.f90
              !!--------------------------------------------------------------
              CALL MOVE_ALLOC(RADII_MVALLOC, RADII)
              RADII = 0.0D0

            !!-------------------------------------------------------------
            !! The following lines are for the subsequent tracer injections
            !! wherein we do need to save the data before overwriting with 
            !! MOVE_ALLOC
            !!-------------------------------------------------------------
            ELSE

              !!------------------------------------------------------------
              !! Get IDs of last tracer released, first new tracer released,
              !! and last new tracer released, respectively.
              !!------------------------------------------------------------
              STOP_SAVE = NUM_P*(NUM_RELEASED-1)-1
              START_NEW = NUM_P*(NUM_RELEASED-1)
              STOP_NEW  = NUM_P*NUM_RELEASED-1

              !!-------------------------------------------------------------------
              !! Save previous coordinates, copy initial coordinates to new tracers
              !!-------------------------------------------------------------------
              TRACERS_XYZ_MVALLOC(0:STOP_SAVE, :) = TRACERS_XYZ
              CALL MOVE_ALLOC(TRACERS_XYZ_MVALLOC, TRACERS_XYZ)
              TRACERS_XYZ(START_NEW:STOP_NEW, :) = TRACERS_SEED

              !!-------------------------------------------------------------
              !! Save previous cell IDs, copy initial cell IDs to new tracers
              !!-------------------------------------------------------------
              TRACERS_CELL_ID_MVALLOC(0:STOP_SAVE) = TRACERS_CELL_ID
              CALL MOVE_ALLOC(TRACERS_CELL_ID_MVALLOC, TRACERS_CELL_ID)
              TRACERS_CELL_ID(START_NEW:STOP_NEW) = CELL_IDS_INIT

              !!--------------------------------------------------------------------
              !! Save previous integration statues, copy initial integration statues 
              !! to new tracers
              !!--------------------------------------------------------------------
              TRACERS_LEFT_GRID_MVALLOC(0:STOP_SAVE) = TRACERS_LEFT_GRID
              CALL MOVE_ALLOC(TRACERS_LEFT_GRID_MVALLOC, TRACERS_LEFT_GRID)
              TRACERS_LEFT_GRID(START_NEW:STOP_NEW) = -1

              !!---------------------------------
              !! No need to preserve tracer radii
              !!---------------------------------
              CALL MOVE_ALLOC(RADII_MVALLOC, RADII)
              RADII = 0.0D0

            END IF
          END IF
        END IF
        
        IF (NUM_RELEASED .GE. 1) THEN

          IF (DUMP_IDX > NUM_OUT_FILES - 1) THEN
            EXIT
          END IF

          !!---------------------------------------------
          !! If tracers have been released, the number
          !! of output files has not exceeded the desired
          !! number of output files, and the output file
          !! time matches the simulation time, write out
          !!---------------------------------------------
          IF ((SIM_TIME < DUMP_TIMES(DUMP_IDX) + EPS) .AND. &
              (SIM_TIME > DUMP_TIMES(DUMP_IDX) - EPS)) THEN

            WRITE(FILE_INDEX_STR, '(I0)') OUT_ID+DUMP_IDX
            VTKNAME = TRIM(FILE_OUT_PREFIX)//'.'//TRIM(FILE_INDEX_STR)//'.vtk'
            OPF_PTS(1) = NUM_P*NUM_RELEASED
            STAT = VTKMD(1)%createVTKMetaData(a_Npt=OPF_PTS, a_DType=DTYPE, a_Title = VTKNAME)
            
            CALL VTKParticleWriter(VTKMD(1), TRACERS_LEFT_GRID, TRACERS_XYZ, RADII, a_FnameWrite=VTKNAME)
            
            print*, "Wrote data to:", TRIM(VTKNAME)

            !!---------------------------------------------------------
            !! Update the number of output files that have been written
            !!---------------------------------------------------------
            DUMP_IDX = DUMP_IDX + 1

          END IF

          !!------------------------------------------------
          !! Advect tracers from SIM_TIME to (SIM_TIME + DT)
          !!------------------------------------------------
          CALL ADVECT_GRID_DATA(INT_FLAGS, FLOW_BBOX, NUM_P*NUM_RELEASED, TRACERS_XYZ,&
                                TRACERS_CELL_ID, TRACERS_LEFT_GRID, [SIM_TIME, DT], WINDOW,&
                                NUM_NODES, NUM_CELLS,&
                                CONNECTIVITY, ADJACENCY, COORDINATES, VEL_LOW,&
                                VEL_MID, VEL_UP)
        END IF

        !!----------------------------------------------
        !! Update SIM_TIME, TIME_INDEX to next time step
        !!----------------------------------------------
        SIM_TIME   = SIM_TIME + DT
        TIME_INDEX = TIME_INDEX + 1

      END DO
      print*, "Deallocating arrays..."
      DEALLOCATE(OPF_PTS)
      DEALLOCATE(VTKMD)
      DEALLOCATE(TRACERS_LEFT_GRID)
      DEALLOCATE(TRACERS_CELL_ID)
      DEALLOCATE(DATA_INDEX_WINDOWS)
      DEALLOCATE(RELEASE_TIMES)
      DEALLOCATE(DUMP_TIMES)
      DEALLOCATE(RADII)
      DEALLOCATE(DATA_TIME_WINDOWS)
      DEALLOCATE(TRACERS_SEED)
      DEALLOCATE(TRACERS_XYZ)
      DEALLOCATE(VEL_BIG)
      print*, "Finished Deallocating arrays."

    END SUBROUTINE ADVECT_SEED_VTK

    !!------------------------------------------------------------------------------------------------------------
    !> @brief Performs advection and with Fortran-based I/O for cartesian grids of tracers given a simulated 
    !! velocity field
    !!
    !!
    !! @param[in]        FLOW_MODEL           integer, flag corresponding to analytical flow model
    !! @param[in]        NUM_P                integer, number of tracers in the simulation
    !! @param[in]        NUM_INJECTIONS       integer, number of tracer injections into the domain for the
    !!                                        length of the entire simulation
    !! @param[in]        INT_FLAGS            integer array, flags corresponding to advection parameters
    !!                                        size: 3
    !! @param[in]        FILE_ID_OUT          integer array, start/stop/delta of vtp output files
    !!                                        size: 3
    !! @param[in]        REF_NUMS             integer array, number of points along each direction
    !!                                        size: 3
    !! @param[in]        FILE_OUT_TIMING      double array, start/delta of vtp output files
    !!                                        size: 2
    !! @param[in]        REF_DIMNS            double array, grid spacing in referential configuration
    !!                                        size: 3
    !! @param[in]        TIMES                double array, simulation times and injection times
    !!                                        configured in input file
    !!                                        size: 5
    !! @param[in]        FLOW_BBOX            double array, bounding box of flow domain
    !!                                        size: 6
    !! @param[in]        TRACER_BBOX          double array, bounding box for seeding tracers
    !!                                        size: 6
    !! @param[in]        FILE_OUT_PREFIX      character, vtp output file prefix
    !!                                        len: 1024
    !!------------------------------------------------------------------------------------------------------------
    SUBROUTINE ADVECT_SEED_ANLYT(FLOW_MODEL, INT_FLAGS, FLOW_BBOX, TRACER_BBOX, NUM_P,&
                                 NUM_INJECTIONS, REF_DIMNS, TIMES, FILE_ID_OUT, FILE_OUT_TIMING,&
                                 FILE_OUT_PREFIX)

      IMPLICIT NONE

      INTEGER, INTENT(IN):: FLOW_MODEL, NUM_P, NUM_INJECTIONS
      INTEGER, DIMENSION(3), INTENT(IN):: INT_FLAGS, FILE_ID_OUT
      DOUBLE PRECISION, DIMENSION(2), INTENT(IN):: FILE_OUT_TIMING
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: REF_DIMNS
      DOUBLE PRECISION, DIMENSION(5), INTENT(IN):: TIMES
      DOUBLE PRECISION, DIMENSION(6), INTENT(IN):: FLOW_BBOX, TRACER_BBOX
      CHARACTER(LEN=1024), INTENT(IN):: FILE_OUT_PREFIX

      LOGICAL:: STAT
      INTEGER:: I, J, DIMN, NUM_RELEASED, OUT_ID, DUMP_IDX, FILE_ID_OUT_START, FILE_ID_OUT_DELTA,&
                FILE_ID_OUT_STOP, STOP_SAVE, START_NEW, STOP_NEW, NUM_OUT_FILES
      INTEGER, DIMENSION(:), ALLOCATABLE:: TRACERS_LEFT_GRID, TRACERS_LEFT_GRID_MVALLOC, OPF_PTS
      DOUBLE PRECISION:: T_START, T_END, DT, SIM_TIME, INJECTION_DT, INJECTION_START,&
                         FILE_OUT_START, FILE_OUT_DELTA
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: RELEASE_TIMES, DUMP_TIMES, RADII, RADII_MVALLOC
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: TRACERS_XYZ, TRACERS_XYZ_MVALLOC, TRACERS_SEED
      CHARACTER(LEN=300):: VTKNAME, FILE_INDEX_STR, DTYPE
      CLASS(VTKmetaData), dimension(:), allocatable:: VTKMD

      DOUBLE PRECISION, PARAMETER:: EPS = 1.0e-12

      !!----------------------
      !! Initialize parameters
      !!----------------------
      DIMN              = INT_FLAGS(1)

      FILE_ID_OUT_START = FILE_ID_OUT(1)
      FILE_ID_OUT_STOP  = FILE_ID_OUT(2)
      FILE_ID_OUT_DELTA = FILE_ID_OUT(3)
      NUM_OUT_FILES     = INT((FILE_ID_OUT_STOP-FILE_ID_OUT_START)/FILE_ID_OUT_DELTA)+1
      OUT_ID            = FILE_ID_OUT_START
      DUMP_IDX          = 0

      FILE_OUT_START    = FILE_OUT_TIMING(1)
      FILE_OUT_DELTA    = FILE_OUT_TIMING(2)

      T_START           = TIMES(1)
      T_END             = TIMES(2)
      DT                = TIMES(3)
      INJECTION_START   = TIMES(4)
      INJECTION_DT      = TIMES(5)

      SIM_TIME          = T_START

      !!---------------------------------------
      !! See outputFormat.f90 for documentation
      !!---------------------------------------
      DTYPE             = 'DATASET POLYDATA'
      ALLOCATE(VTKMD(1), OPF_PTS(1))
      
      !!-------------------------------------------
      !! Initialize array of tracer injection times
      !!-------------------------------------------
      ALLOCATE(RELEASE_TIMES(0:NUM_INJECTIONS-1))
      DO I = 0, NUM_INJECTIONS-1
        RELEASE_TIMES(I) = INJECTION_START + FLOAT(I)*INJECTION_DT
      END DO

      !!------------------------------------------------------------
      !! Initialize array of output times at which to write vtp data
      !!------------------------------------------------------------
      ALLOCATE(DUMP_TIMES(0:NUM_OUT_FILES-1))
      DO I = 0, NUM_OUT_FILES-1
        DUMP_TIMES(I) = FILE_OUT_START + FLOAT(I)*FILE_OUT_DELTA
      END DO

      !!--------------------------------
      !! Seed points in a cartesian grid
      !!--------------------------------
      ALLOCATE(TRACERS_SEED(0:NUM_P-1,3))
      CALL SEED_POINTS_CARTESIAN(NUM_P, REF_DIMNS, TRACER_BBOX, TRACERS_SEED)

      NUM_RELEASED  = 0

      !!-----------------
      !! Begin simulation
      !!-----------------
      DO WHILE (ABS(T_END-SIM_TIME) > EPS)

        I = NUM_RELEASED

        !!------------------------------------------------
        !! Check to see if new tracers need to be released
        !!------------------------------------------------
        IF (I <= NUM_INJECTIONS - 1) THEN
          IF ((SIM_TIME >= (RELEASE_TIMES(I)-EPS)) .AND. (SIM_TIME <= (RELEASE_TIMES(I) + EPS))) THEN

            NUM_RELEASED = NUM_RELEASED + 1
            PRINT*, "Injecting", NUM_P, "particles at t = ", SIM_TIME

            ALLOCATE(TRACERS_XYZ_MVALLOC(0:NUM_P*NUM_RELEASED-1, 3))
            ALLOCATE(TRACERS_LEFT_GRID_MVALLOC(0:NUM_RELEASED*NUM_P-1))
            ALLOCATE(RADII_MVALLOC(NUM_P*NUM_RELEASED))

            !!------------------------------------------------------------------
            !! The following lines are for the first tracer injection wherein we
            !! do not need to save the data before overwriting with MOVE_ALLOC
            !!------------------------------------------------------------------
            IF (NUM_RELEASED == 1) THEN

              !!--------------------------------
              !! Copy initial tracer coordinates
              !!--------------------------------
              CALL MOVE_ALLOC(TRACERS_XYZ_MVALLOC, TRACERS_XYZ)
              TRACERS_XYZ = TRACERS_SEED

              !!-----------------------------------------
              !! Copy initial tracer integration statuses
              !!-----------------------------------------
              CALL MOVE_ALLOC(TRACERS_LEFT_GRID_MVALLOC, TRACERS_LEFT_GRID)
              TRACERS_LEFT_GRID = -1

              !!--------------------------------------------------------------
              !! Copy initial tracer radii
              !! NOTE: Finite-sized tracer integration with geometrical checks
              !!       has not yet been implemented. This line is present for
              !!       consistency with outputFormat.f90
              !!--------------------------------------------------------------
              CALL MOVE_ALLOC(RADII_MVALLOC, RADII)
              RADII = 0.0D0

            !!-------------------------------------------------------------
            !! The following lines are for the subsequent tracer injections
            !! wherein we do need to save the data before overwriting with 
            !! MOVE_ALLOC
            !!-------------------------------------------------------------
            ELSE

              !!------------------------------------------------------------
              !! Get IDs of last tracer released, first new tracer released,
              !! and last new tracer released, respectively.
              !!------------------------------------------------------------
              STOP_SAVE = NUM_P*(NUM_RELEASED-1)-1
              START_NEW = NUM_P*(NUM_RELEASED-1)
              STOP_NEW  = NUM_P*NUM_RELEASED-1

              !!-------------------------------------------------------------------
              !! Save previous coordinates, copy initial coordinates to new tracers
              !!-------------------------------------------------------------------
              TRACERS_XYZ_MVALLOC(0:STOP_SAVE, :) = TRACERS_XYZ
              CALL MOVE_ALLOC(TRACERS_XYZ_MVALLOC, TRACERS_XYZ)
              TRACERS_XYZ(START_NEW:STOP_NEW, :) = TRACERS_SEED

              !!--------------------------------------------------------------------
              !! Save previous integration statues, copy initial integration statues 
              !! to new tracers
              !!--------------------------------------------------------------------
              TRACERS_LEFT_GRID_MVALLOC(0:STOP_SAVE) = TRACERS_LEFT_GRID
              CALL MOVE_ALLOC(TRACERS_LEFT_GRID_MVALLOC, TRACERS_LEFT_GRID)
              TRACERS_LEFT_GRID(START_NEW:STOP_NEW) = -1

              !!---------------------------------
              !! No need to preserve tracer radii
              !!---------------------------------
              CALL MOVE_ALLOC(RADII_MVALLOC, RADII)
              RADII = 0.0D0

            END IF
          END IF
        END IF
        
        IF (NUM_RELEASED .GE. 1) THEN

          IF (DUMP_IDX > NUM_OUT_FILES - 1) THEN
            EXIT
          END IF

          !!---------------------------------------------
          !! If tracers have been released, the number
          !! of output files has not exceeded the desired
          !! number of output files, and the output file
          !! time matches the simulation time, write out
          !!---------------------------------------------
          IF ((SIM_TIME < DUMP_TIMES(DUMP_IDX) + EPS) .AND. &
              (SIM_TIME > DUMP_TIMES(DUMP_IDX) - EPS)) THEN

            WRITE(FILE_INDEX_STR, '(I0)') OUT_ID+DUMP_IDX
            VTKNAME = TRIM(FILE_OUT_PREFIX)//'.'//TRIM(FILE_INDEX_STR)//'.vtk'
            OPF_PTS(1) = NUM_P*NUM_RELEASED
            STAT = VTKMD(1)%createVTKMetaData(a_Npt=OPF_PTS, a_DType=DTYPE, a_Title = VTKNAME)
            
            CALL VTKParticleWriter(VTKMD(1), TRACERS_LEFT_GRID, TRACERS_XYZ, RADII, a_FnameWrite=VTKNAME)
            
            print*, "Wrote data to:", TRIM(VTKNAME)

            !!---------------------------------------------------------
            !! Update the number of output files that have been written
            !!---------------------------------------------------------
            DUMP_IDX = DUMP_IDX + 1

          END IF

          !!------------------------------------------------
          !! Advect tracers from SIM_TIME to (SIM_TIME + DT)
          !!------------------------------------------------
          CALL ADVECT_GRID_ANALYTICAL(FLOW_MODEL, INT_FLAGS, FLOW_BBOX, NUM_P*NUM_RELEASED,&
                                      TRACERS_XYZ, TRACERS_LEFT_GRID, [SIM_TIME, DT])

        END IF

        !!----------------------------------------------
        !! Update SIM_TIME to next time step
        !!----------------------------------------------
        SIM_TIME   = SIM_TIME + DT

      END DO
      print*, "Deallocating arrays..."
      DEALLOCATE(OPF_PTS)
      DEALLOCATE(VTKMD)
      DEALLOCATE(TRACERS_LEFT_GRID)
      DEALLOCATE(RELEASE_TIMES)
      DEALLOCATE(DUMP_TIMES)
      DEALLOCATE(RADII)
      DEALLOCATE(TRACERS_SEED)
      DEALLOCATE(TRACERS_XYZ)
      print*, "Finished Deallocating arrays."

    END SUBROUTINE ADVECT_SEED_ANLYT

END MODULE MAIN
