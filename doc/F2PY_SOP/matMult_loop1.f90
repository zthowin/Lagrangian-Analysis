!!--------------------------------------------------------------
!! Example code from https://modelingguru.nasa.gov/docs/DOC-2343
!!--------------------------------------------------------------
module matMult_loop1

  implicit none

  contains

    subroutine matMult_loop1_routine(A, B, C, n1, n2, n3)

      implicit none

      integer, intent(in):: n1, n2, n3
      double precision, dimension(n1, n2), intent(in):: A
      double precision, dimension(n2, n3), intent(in):: B

      double precision, dimension(n1, n3), intent(out):: C

      integer:: i, j, k

      do j = 1, n3
         do i = 1, n1
            C(i, j) = 0.0d0
            do k = 1, n2
               C(i, j) = C(i, j) + A(i, k)*B(k, j)
            end do
         end do
      end do

    end subroutine matMult_loop1_routine

end module matMult_loop1
