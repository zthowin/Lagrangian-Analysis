!file: callback_full.f90
module callback_full

    implicit none

    contains

        subroutine calculate(func,x,n)
            implicit none

            external func
            double precision:: func
            integer:: n,i
            double precision:: x(n)
            do i=1,n
                 x(i) = func(x(i))
            end do
        end subroutine calculate

end module callback_full