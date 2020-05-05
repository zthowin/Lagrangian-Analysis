!---------------------------------------------------------------------------------
!> @file File containing source code for formatted data read/write functionalities
!!
!! @brief Module containing all functions and subroutines that comprise an API to 
!! transport data from a raw .dat format to a formatted file-type readable by
!! visualization libraries like Paraview, Visit, or TecPlot
!!
!! @author Debanjan Mukherjee
!!
!! @todo Add subroutine to read VTK data from a file for interpolating continuous
!! field data (or perhaps use a Python script to do this and call it from Fortran)
!<--------------------------------------------------------------------------------
module outputFormat

  implicit none

  !----------------------------------------------
  !> The user defined data type to define the 
  !! data structure for holding the vtk metadata
  !<---------------------------------------------
  type VTKmetaData
    private
    character(len=100):: m_DType        
    integer, dimension(:),allocatable:: m_numPts
    integer:: m_Scalars
    integer:: m_Vectors
    integer:: m_Tensors
    character(len=100):: m_Title
    character(len=100), dimension(:), allocatable:: m_Variables
    contains
    procedure, public:: createVTKmetaData
    ! procedure, public:: F_TO_VTK
    procedure, public:: VTKParticleWriter
    procedure, public:: VTKContinuumWriter
  end type VTKmetaData

  contains

  !-----------------------------------------------------------------
  !> @brief Function for constructing an object of type VTK metadata
  !!
  !! @param[in] a_NPt The number of points
  !! @param[in] a_DType The VTK data type classifier
  !! @param[in] (optional) a_S The number of scalar data fields
  !! @param[in] (optional) a_V The number of vector data fields
  !! @param[in] (optional) a_T The number of tensor data fields
  !! @param[in] a_Title The file title for the header
  !! @param[in] a_Variables (optional) The names of all variables
  !! each of length 100 characters (scalars first, then vectors, 
  !! then tensors)
  !<-----------------------------------------------------------------
  function createVTKmetaData(this, a_NPt, a_DType, a_S, a_V, a_T, a_Title, a_Variables)

    implicit none
    class(VTKmetaData):: this
    character(len=100), intent(in):: a_DType
    integer, dimension(:), allocatable, intent(in):: a_NPt
    integer, optional, intent(in):: a_S
    integer, optional, intent(in):: a_V
    integer, optional, intent(in):: a_T
    character(len=100), intent(in):: a_Title
    character(len=100), dimension(:), allocatable, optional, intent(in):: a_Variables
    logical:: createVTKmetaData
    
    integer:: k1, countVar, numData
    character(len=10):: tempStr

    countVar = 1

    this%m_DType    = a_DType

    if ( allocated(this%m_numPts) ) then
      deallocate(this%m_numPts)
      allocate(this%m_numPts(size(a_NPt,1)))
      this%m_numPts   = a_NPt
    else
      allocate(this%m_numPts(size(a_NPt,1)))
      this%m_numPts   = a_NPt
    end if
    
    if ( present(a_S) ) then
      this%m_Scalars  = a_S
    else
      this%m_Scalars  = 0
    end if

    if ( present(a_V) ) then
      this%m_Vectors  = a_V
    else
      this%m_Vectors  = 0
    end if

    if ( present(a_T) ) then
      this%m_Tensors  = a_T
    else
      this%m_Tensors = 0
    end if
    
    if ( (this%m_Scalars == 0) .and. (this%m_Vectors == 0) .and. (this%m_Tensors == 0) ) then
      if (trim(a_DType) /= 'DATASET POLYDATA') then
        print*,"No Data To Write Into The VTK File. Seems Like Wrong Inputs Are Used"
        stop
      end if
    end if

    this%m_Title    = a_Title

    if ( present(a_Variables) ) then
      if ( .not.allocated(this%m_Variables) ) then
        allocate(this%m_Variables(size(a_Variables,1)))
      else 
        deallocate(this%m_Variables)
        allocate(this%m_Variables(size(a_Variables,1)))
      end if        
      this%m_Variables  = a_Variables
    else
      numData = this%m_Scalars + &
                this%m_Vectors + &
                this%m_Tensors 

      if ( .not.allocated(this%m_Variables) ) then
        allocate(this%m_Variables(numData))
      else
        deallocate(this%m_Variables)
        allocate(this%m_Variables(numData))
      end if
      
      do k1 = 1, this%m_Scalars
        write(tempStr,*) k1
        this%m_Variables(countVar) = 'scalar'//trim(adjustl(tempStr))
        countVar = countVar + 1 
      end do
      do k1 = 1, this%m_Vectors
        write(tempStr,*) k1
        this%m_Variables(countVar) = 'vector'//trim(adjustl(tempStr))
        countVar = countVar + 1 
      end do
      do k1 = 1, this%m_Tensors
        write(tempStr,*) k1
        this%m_Variables(countVar) = 'tensor'//trim(adjustl(tempStr))
        countVar = countVar + 1 
      end do
    end if
  
    createVTKmetaData = .true.

  end function createVTKmetaData

  !--------------------------------------------------------------
  !> @brief Subroutine to write a Paraview readable .vtk
  !! file from a .dat file using a 'STRUCTURED GRID' format.
  !!
  !! Dat file data is read in form of a matrix, where every 
  !! column denotes a field. The data format can be changed 
  !! by changing the headers. 
  !! For convenient usage, try formatting .dat file as:
  !! x, y, z, scalar1, vector11, vector12, vector13, tensor11 ...
  !<--------------------------------------------------------------
  ! subroutine F_TO_VTK(this, a_FnameRead, a_FunitRead, a_FnameWrite, a_FunitWrite)
    
  !   use commonOperations, only: int2str, digit_to_ch
  
  !   implicit none
  !   class(VTKmetaData):: this
  !   character(len=40), intent(in):: a_FnameRead 
  !   character(len=40), intent(in):: a_FnameWrite
  !   integer, intent(in):: a_FunitRead
  !   integer, intent(in):: a_FunitWrite

  !   character(len=20):: title
  !   character(len=20):: scalardata
  !   character(len=20):: vectordata 
  !   character(len=20):: tensordata

  !   character(len=40), parameter:: dataset = 'DATASET STRUCTURED_GRID'
  !   character(len=40), parameter:: vtkGenInfo = '# vtk DataFile Version 2.0'
  !   ! In case more fields are in the data file, define more data labels 
  !   ! by changing the number in the label: e.g. scalardata_2 etc.

  !   character(len=10):: n_x_str, n_y_str, n_z_str, n_pt_str, num_lines_str

  !   character(len=60):: header_DIM, header_POINT, header_PDATA, header_SCAL
    
  !   integer:: n_cols, n_x, n_y, n_z, n_pt
  !   double precision, dimension(:,:), allocatable:: data_matrix
  !   integer:: num_lines, n_params
  !   integer:: k1, k2
  !   integer:: IO, err
    
  !   if ( size(this%m_numPts,1) /= 3 ) then
  !     print*,"Wrong dimension of m_numPts"
  !     stop
  !   end if

  !   n_x     = this%m_numPts(1)    
  !   n_y     = this%m_numPts(2)
  !   n_z     = this%m_numPts(3)
  !   n_pt    = n_x*n_y*n_z

  !   n_params = 3 + this%m_Scalars + 3*this%m_Vectors + 9*this%m_Tensors
    
  !   ! With the knowledge that each line has (n_params) number of columns:
  !   !--------------------------------------------------------------------
  !   allocate(data_matrix(n_pt, n_params), stat = err)
  !   if (err /= 0) then
  !     print *,'Could not allocate the data matrix'
  !   end if

  !   ! Now read and store data into the allocated matrix:
  !   !---------------------------------------------------
  !   open(unit = a_FunitRead, file = a_FnameRead)
  
  !   do k1 = 1,n_pt
  !     !do k2 = 1, n_params
  !       read(unit=a_FunitRead,fmt='(3E20.10)', iostat = IO)(data_matrix(k1,k2), k2=1, n_params)
  !     !end do
  !   end do

  !   close(unit = a_FunitRead)

  !   ! Now create headers for the VTK file:
  !   !-------------------------------------

  !   ! 0. Convert numbers to strings to add into headers:
  !   !---------------------------------------------------
  !   call int2str(n_x, n_x_str)
  !   call int2str(n_y, n_y_str)
  !   call int2str(n_z, n_z_str)
  !   call int2str(n_pt, n_pt_str)
    
  !   ! 1. Create the main headers:
  !   !-------------------------------------------------------------------
  !   title         = trim(this%m_Title)
  !   header_DIM    = 'DIMENSIONS ' // n_x_str // ' ' // n_y_str // ' ' // n_z_str
  !   header_POINT  = 'POINTS ' // n_pt_str //' float'
  !   header_PDATA  = 'POINT_DATA '// n_pt_str

  !   ! 2. Write the headers and the point coordinates:
  !   !------------------------------------------------
  !   open(unit = a_FunitWrite, file = a_FnameWrite)

  !   write(unit = a_FunitWrite, fmt = '(A)') vtkGenInfo
  !   write(unit = a_FunitWrite, fmt = '(A)') title
  !   write(unit = a_FunitWrite, fmt = '(A)') 'ASCII'
  !   write(unit = a_FunitWrite, fmt = '(A)') dataset
  !   write(unit = a_FunitWrite, fmt = '(A)') header_DIM
  !   write(unit = a_FunitWrite, fmt = '(A)') header_POINT

  !   ! Write the point coordinates:
  !   !-----------------------------
  !   do k1 = 1, n_pt
  !     write(unit = a_FunitWrite, fmt = '(3E20.10)') data_matrix(k1,1), data_matrix(k1,2), data_matrix(k1,3)
  !   end do

  !   ! Mandatory blank line between two sections of the file:
  !   !-------------------------------------------------------
  !   write(unit = a_FunitWrite, fmt = '(a)') ''

  !   ! Start writing data:
  !   !--------------------
  !   write(unit=a_FunitWrite, fmt='(A)') header_PDATA    

  !   do k1 = 1, this%m_Scalars 
  !     ! 2. Create appropriate header for data 
  !     !---------------------------------------
  !     header_SCAL   = 'SCALARS '//trim(this%m_Variables(k1))//' float 1'
  !     write(unit = a_FunitWrite, fmt = '(A)') trim(header_SCAL)
  !     write(unit = a_FunitWrite, fmt = '(A)') 'LOOKUP_TABLE default'
  !     do k2 = 1, n_pt
  !       write(unit = a_FunitWrite, fmt = '(E20.10)') data_matrix(k2,3+k1)
  !     end do
  !   end do

  !   ! Write all vector data here (to be completed later)
  !   !------------------------------------------------------

  !   ! Write all tensor data here (to be completed later)
  !   !------------------------------------------------------

  !   close(unit=a_FunitWrite)
        
  ! end subroutine F_TO_VTK

  !------------------------------------------------------------
  !> @brief Subroutine to write particle data into a VTK file
  !! 
  !! @note The program that is used for testing the 
  !! capabilities of the subroutine is in the source file
  !! testVTKparticle.f90
  !!
  !! @note The number of points is equivalent to the number of 
  !! particles in the ensemble being written into the VTK file
  !!
  !! @param[in] a_PCoords The vector of particle coordinates
  !! @param[in] a_PRads The vector of particle radius values
  !! @param[in] a_PScalars (optional) The vector of scalar 
  !! field data based on particles
  !! @param[in] a_PVectors (optional) The vector of vector 
  !! field data based on particles
  !! @param[in] a_FUnitWrite (optional) The integer file unit
  !! specified for the output VTK file
  !! @param[in] a_FnameWrite (optional) The name of the file
  !! specified for the output VTK file. If this is not specified 
  !! then the data is saved as 'particle-data.vtk'
  !<------------------------------------------------------------ 
  subroutine VTKParticleWriter(this, a_PLeftGrid, a_PCoords, a_PRads, a_PScalars, &
              a_PVectors, a_FUnitWrite, a_FnameWrite)

    implicit none 
    class(VTKmetaData):: this
    integer, dimension(:), allocatable, intent(in):: a_PLeftGrid
    double precision, dimension(:,:), allocatable, intent(in):: a_PCoords
    double precision, dimension(:), allocatable, intent(in):: a_PRads
    double precision, dimension(:,:), allocatable, optional, intent(in):: a_PScalars
    double precision, dimension(:,:,:), allocatable, optional, intent(in):: a_PVectors
    integer, optional, intent(in):: a_FUnitWrite
    character(len=256), optional, intent(in):: a_FnameWrite
    
    character(len=256):: fName
    character(len=20):: tempStr
    character(len=40), parameter:: vtkGenInfo = '# vtk DataFile Version 2.0'
    character(len=40), parameter:: vtkLookUp  = 'LOOKUP_TABLE default'
    character(len=60):: header_DIM, header_POINT, header_PDATA, header_SCAL, header_VECT

    integer:: fUnit, k1, countVar
    
    !---------------------------------------------------------------------
    !> Perform some basic checks on the consistency of the input arguments
    !---------------------------------------------------------------------
    if ( this%m_DType /= 'DATASET POLYDATA' ) then
      print*,"Particle Data Type Is Supported Currently Using VTK-POLYDATA"
      stop
    end if

    if ( size(this%m_numPts) /= 1 ) then
      print*,"Polydata should have a single specified number of points"
      stop
    end if
      
    !-------------------------------------------------------
    !> Resolve the naming of the output file, and assignment
    !! of an appropriate output file unit for I/O
    !<------------------------------------------------------
    if ( present(a_FUnitWrite) ) then
      fUnit = a_FUnitWrite
    else
      fUnit = 10
    end if
  
    if ( present(a_FnameWrite) ) then
      fName = a_FnameWrite
    else
      fname = 'particle-data.vtk'
    end if

    !----------------------------------------
    !> Open the file for writing data into it
    !<---------------------------------------
    open(unit=fUnit, file=trim(fName), action='WRITE')
    
    !-------------------------------------------------
    !> Write the generic file headers for the vtk file
    !<------------------------------------------------
    write(unit=fUnit,fmt='(A)') vtkGenInfo
    write(unit=fUnit,fmt='(A)') trim(this%m_Title)
    write(unit=fUnit,fmt='(A)') 'ASCII'
    write(unit=fUnit,fmt='(A)') trim(this%m_DType)
    
    !----------------------------------------------------------
    !> Write the point coordinates with all appropriate headers
    !<---------------------------------------------------------
    write(tempStr,*) this%m_NumPts(1)
    header_POINT = 'POINTS '//trim(adjustl(tempStr))//' float'
    
    write(unit=fUnit,fmt='(A)') trim(header_POINT)
    do k1 = 0, this%m_numPts(1)-1
      !if (a_PLeftGrid(k1) /= 1) then
        write(unit=fUnit,fmt='(3F20.10)') a_PCoords(k1,1), a_PCoords(k1,2), a_PCoords(k1,3)
      ! Temporary workaround for points that have left domain that we don't want to display
      !else
       ! write(unit=fUnit,fmt='(3F20.10)') -1.0d0, 0.0d0, 0.0d0
      !end if

    end do    

    !-------------------------------------
    !> Write the header for the point data
    !-------------------------------------
    header_PDATA = 'POINT_DATA '//trim(adjustl(tempStr))
    write(unit=fUnit,fmt='(A)') header_PDATA
    
    !----------------------------------------------
    !> The first scalar should be the radius values
    !! with all appropriate data headers
    !<---------------------------------------------
    header_SCAL = 'SCALARS radius float 1'
    write(unit=fUnit, fmt='(A)') header_SCAL
    write(unit=fUnit, fmt='(A)') vtkLookUp
    do k1 = 1, this%m_numPts(1)
      write(unit=fUnit,fmt='(E20.10)') a_PRads(k1)
    end do

    !---------------------------------------------------
    !> Then write all the other scalar data presented to 
    !! the subroutine for writing
    !<--------------------------------------------------
    if ( present(a_PScalars) ) then
      do countVar = 1, this%m_Scalars 
        header_SCAL = 'SCALARS '//trim(this%m_Variables(countVar))//' float 1'
        write(unit=fUnit,fmt='(A)') header_SCAL
        write(unit=fUnit,fmt='(A)') vtkLookUp
        do k1 = 1, this%m_numPts(1)
          write(unit=fUnit, fmt='(E20.10)') a_PScalars(k1,countVar)
        end do
      end do
    end if

    !---------------------------------------------------
    !> Then write all the other vector data presented to 
    !! the subroutine for writing
    !<--------------------------------------------------
    if ( present(a_PVectors) ) then
    
      do countVar = 1, this%m_Vectors
        header_VECT = 'VECTORS '//trim(this%m_Variables(this%m_Scalars+countVar))//' float'
        write(unit=fUnit,fmt='(A)') header_VECT
        !write(unit=fUnit,fmt='(A)') vtkLookUp  ! NO NEED FOR LOOKUP_TABLE TO BE THERE
        do k1 = 1, this%m_numPts(1)
          write(unit=fUnit, fmt='(3E20.10)') a_PVectors(k1,countVar,1:3)        
        end do
      end do
      
    end if

    close(unit=fUnit)

  end subroutine VTKParticleWriter

  !-----------------------------------------------------------------------------
  !> @brief Subroutine to wrote a Paraview readable VTK file for post-processing
  !! continnum field data using a STRUCTURED GRID formatted data file
  !!
  !! @param[in] a_CCoords The coordinates of the nodal points for the continuum
  !! grid over which data is to be visualized
  !! @param[in] a_CScalars (optional) The scalar data for visualization over the 
  !! continuum grid
  !! @param[in] a_CVectors (optional) The vector data for visualization over the
  !! continuum grid
  !! @param[in] a_CTensors (optional) The tensor data for visualization over the
  !! continuum grid
  !<----------------------------------------------------------------------------
  subroutine VTKContinuumWriter(this, a_CCoords, a_CScalars, a_CVectors, a_CTensors, &
              a_FUnitWrite, a_FnameWrite)

    implicit none
  
    class(VTKmetaData):: this
    double precision, dimension(:,:), allocatable, intent(in):: a_CCoords
    double precision, dimension(:,:), allocatable, optional, intent(in):: a_CScalars
    double precision, dimension(:,:,:), allocatable, optional, intent(in):: a_CVectors
    double precision, dimension(:,:,:,:), allocatable, optional, intent(in):: a_CTensors
    integer, optional, intent(in):: a_FUnitWrite
    character(len=100), optional, intent(in):: a_FnameWrite

    character(len=100):: fName
    character(len=40), parameter:: vtkGenInfo = '# vtk DataFile Version 2.0'
    character(len=40), parameter:: vtkLookUp  = 'LOOKUP_TABLE default'
    character(len=100):: header_DIM, header_POINT, header_PDATA, &
                          header_SCAL, header_VEC, header_TEN

    integer:: countVar, k1, totalPts, fUnit
    character(len=20):: numStr1, numStr2, numStr3, numStr
  
    !---------------------------------------------------------------------
    !> Perform some basic checks on the consistency of the input arguments
    !---------------------------------------------------------------------
    if ( this%m_DType /= 'DATASET STRUCTURED_GRID' ) then
      print*,"Continnum Data Is Supported Currently Using VTK-STRUCTURED-GRID"
      stop
    end if

    !-------------------------------------------------------
    !> Resolve the naming of the output file, and assignment
    !! of an appropriate output file unit for I/O
    !<------------------------------------------------------
    if ( present(a_FUnitWrite) ) then
      fUnit = a_FUnitWrite
    else
      fUnit = 10
    end if
  
    if ( present(a_FnameWrite) ) then
      fName = a_FnameWrite
    else
      fname = 'continuum-data.vtk'
    end if

    !----------------------------------------
    !> Open the file for writing data into it
    !<---------------------------------------
    open(unit=fUnit, file=trim(fName), action='WRITE')
    
    !-------------------------------------------------
    !> Write the generic file headers for the vtk file
    !<------------------------------------------------
    write(unit=fUnit,fmt='(A)') vtkGenInfo
    write(unit=fUnit,fmt='(A)') trim(this%m_Title)
    write(unit=fUnit,fmt='(A)') 'ASCII'
    write(unit=fUnit,fmt='(A)') trim(this%m_DType)

    !----------------------------------------------------------
    !> Write the point coordinates with all appropriate headers
    !<---------------------------------------------------------
    totalPts = (this%m_NumPts(1))*(this%m_NumPts(2))*(this%m_NumPts(3))
    write(numStr1,*) this%m_NumPts(1)
    write(numStr2,*) this%m_NumPts(2)
    write(numStr3,*) this%m_NumPts(3)
    write(numStr,*) totalPts

    header_DIM = 'DIMENSIONS '//trim(adjustl(numStr1))//' '//&
                  trim(adjustl(numStr2))//' '//trim(adjustl(numStr3))
  
    header_POINT = 'POINTS '//trim(adjustl(numStr))//' float'
    
    write(unit=fUnit,fmt='(A)') trim(header_DIM)
    write(unit=fUnit,fmt='(A)') trim(header_POINT)
    do k1 = 1, totalPts
      write(unit=fUnit,fmt='(3E20.10)') a_CCoords(k1,1), a_CCoords(k1,2), a_CCoords(k1,3)
    end do    
    
    header_PDATA = 'POINT_DATA '//trim(adjustl(numStr))
    write(unit=fUnit,fmt='(A)') header_PDATA

    !---------------------------------------------------
    !> Then write all the other scalar data presented to 
    !! the subroutine for writing
    !<--------------------------------------------------
    if ( present(a_CScalars) ) then
      do countVar = 1, this%m_Scalars
        header_SCAL = 'SCALARS '//trim(this%m_Variables(countVar))//' float 1'
        write(unit=fUnit, fmt='(A)') header_SCAL
        write(unit=fUnit, fmt='(A)') vtkLookUp
        do k1 = 1, totalPts
          write(unit=fUnit, fmt='(E20.10)') a_CScalars(countVar,k1)
        end do
      end do
    end if

    !---------------------------------------------------
    !> Then write all the other vector data presented to 
    !! the subroutine for writing
    !<--------------------------------------------------
    if ( present(a_CVectors) ) then
      do countVar = 1, this%m_Vectors
        header_VEC = 'VECTORS '//trim(this%m_Variables(this%m_Scalars+countVar))//' float'
        write(unit=fUnit, fmt='(A)') header_VEC
        write(unit=fUnit, fmt='(A)') vtkLookUp
        do k1 = 1, totalPts
          write(unit=fUnit, fmt='(3E20.10)') a_CVectors(countVar,k1,1:3)
        end do
      end do
    end if

    !---------------------------------------------------
    !> Then write all the other tensor data presented to 
    !! the subroutine for writing
    !<--------------------------------------------------
    if ( present(a_CTensors) ) then
      do countVar = 1, this%m_Tensors
        header_TEN = 'TENSORS '//trim(this%m_Variables(this%m_Scalars+this%m_Vectors+countVar))//' float'
        write(unit=fUnit, fmt='(A)') header_TEN
        write(unit=fUnit, fmt='(A)') vtkLookUp
        do k1 = 1, totalPts
          write(unit=fUnit,fmt='(3E20.10)') a_CTensors(countVar,k1,1,1:3)
          write(unit=fUnit,fmt='(3E20.10)') a_CTensors(countVar,k1,2,1:3)
          write(unit=fUnit,fmt='(3E20.10)') a_CTensors(countVar,k1,3,1:3)
        end do
      end do
    end if

  end subroutine VTKContinuumWriter

  !--------------------------------------------------------
  !> @brief Subroutine to convert a data file containing
  !! particle positions and radius into a VTK Polydata file
  !! for visualization
  !!
  !! @note This is particularly meant for quick generation 
  !! of a mesh of particles at verious specified locations
  !! and then convert all of them into a VTK file
  !!
  !! @param[in] a_DatUnit Integer valued file unit for the 
  !! the dat file to be read
  !! @param[in] a_DatName The name of the '.dat' file to be 
  !! read (include the extension)
  !! @param[in] a_VTKUnit Integer values file unit for the
  !! VTK file to be written
  !! @param[in] a_VTKName The name of the '.vtk' file to be
  !! created (include the extension)
  !<-------------------------------------------------------
  subroutine dataToParticleVTK(a_DatUnit, a_DatName, a_VTKUnit, a_VTKName)

    implicit none
    integer, intent(in):: a_DatUnit
    character(len=*), intent(in):: a_DatName
    integer, intent(in):: a_VTKUnit
    character(len=*), intent(in):: a_VTKName

    integer:: k1, npt, check
    double precision:: X, Y, Z, R
    double precision, dimension(:,:), allocatable:: temp, temptemp

    character(len=40), parameter:: vtkGenInfo = '# vtk DataFile Version 2.0'
    character(len=40), parameter:: vtkLookUp  = 'LOOKUP_TABLE default'
    character(len=40), parameter:: vtkDType   = 'DATASET POLYDATA'

    character(len=200):: vtkTitle, vtkPData, vtkPoint, vtkScal, tempStr

    open(unit=a_DatUnit, file=trim(a_DatName), action='READ')

    npt = 0

    do
      read(unit=a_DatUnit, fmt='(4E20.10)',iostat=check) X,Y,Z,R
      if ( check > 0 ) then
        print*,"File not read correctly"
        stop
      elseif ( check < 0 ) then
        exit
      else
        npt = npt + 1
        if ( .not.allocated(temp) ) then
          allocate(temp(1,4)) 
          temp(1,1:4) = (/X,Y,Z,R/)
        else
          call move_alloc(temp, temptemp)
          allocate(temp(size(temptemp,1)+1,4))
          temp(1:size(temptemp,1),1:4) = temptemp
          temp(size(temptemp,1)+1,1:4) = (/X,Y,Z,R/)
          deallocate(temptemp)
        end if
      end if
    end do

    close(unit=a_DatUnit)

    vtkTitle = trim('Converted Data From '//trim(a_DatName))

    write(tempStr,*) npt
    vtkPoint  = 'POINTS '//trim(adjustl(tempStr))//' float'
    vtkPdata  = 'POINT_DATA '//trim(adjustl(tempStr))
    vtkScal   = 'SCALARS radius float 1'

    open(unit=a_VTKUnit, file=trim(a_VTKName), action='WRITE')
    
    write(unit=a_VTKUnit,fmt='(A)') trim(vtkGenInfo)
    write(unit=a_VTKUnit,fmt='(A)') trim(vtkTitle)
    write(unit=a_VTKUnit,fmt='(A)') 'ASCII'
    write(unit=a_VTKUnit,fmt='(A)') trim(vtkDType)

    write(unit=a_VTKUnit,fmt='(A)') trim(vtkPoint)

    do k1 = 1, npt
      write(unit=a_VTKUnit,fmt='(3E20.10)') temp(k1,1), temp(k1,2), temp(k1,3)
    end do

    write(unit=a_VTKUnit,fmt='(A)') trim(vtkPdata)
    write(unit=a_VTKUnit,fmt='(A)') trim(vtkScal)
    write(unit=a_VTKUnit,fmt='(A)') trim(vtkLookUp)

    do k1 = 1, npt
      write(unit=a_VTKUnit,fmt='(E20.10)') temp(k1,4)
    end do

    close(unit=a_VTKUnit)
    
  end subroutine dataToParticleVTK

end module outputFormat
