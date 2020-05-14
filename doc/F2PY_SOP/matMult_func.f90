!!--------------------------------------------------------------
!! Example code from https://modelingguru.nasa.gov/docs/DOC-2343
!!--------------------------------------------------------------
module matMult_func

  implicit none

	contains

      subroutine matMult_func_routine(A, B, C, n1, n2, n3)

      	implicit none

        integer, intent(in):: n1, n2, n3
        double precision, dimension(n1, n2), intent(in):: A
        double precision, dimension(n2, n3), intent(in):: B

        double precision, dimension(n1, n3), intent(out):: C

        C = matmul(A,B)

      end subroutine matMult_func_routine

end module matMult_func