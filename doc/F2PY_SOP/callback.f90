!file: callback.f90
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