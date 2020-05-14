!!--------------------------------------------------------------
!! Example code from https://modelingguru.nasa.gov/docs/DOC-2343
!!--------------------------------------------------------------
module matrixMult

  use matMult_loop1
  use matMult_loop2
  use matMult_func

  implicit none

  contains

      subroutine matrixMult_routine(A, B, C, n1, n2, n3, iCase)

        implicit none

        integer, intent(in) :: n1, n2, n3
        integer, intent(in) :: iCase
        double precision, dimension(n1, n2), intent(in):: A
        double precision, dimension(n2, n3), intent(in):: B

        double precision, dimension(n1, n3), intent(out):: C

        select case (iCase)
           case (1)
              call matMult_loop1_routine(A, B, C, n1, n2, n3)
           case (2)
              call matMult_loop2_routine(A, B, C, n1, n2, n3)
           case (3)
              call matMult_func_routine(A, B, C, n1, n2, n3)
        end select

      end subroutine matrixMult_routine

end module matrixMult
