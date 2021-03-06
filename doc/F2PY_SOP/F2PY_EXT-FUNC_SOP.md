This is a detailed SOP on how to write Fortran subroutines that call external Python functions 
using NumPy's f2py module. Example code is included. MCVE codes used for the examples are
callback.f90, callback_sig.pyf, callback_full.f90, callback_full_sig.pyf.


Additional resources:

- https://docs.scipy.org/doc/numpy/f2py/signature-file.html
- https://github.com/thehackerwithin/PyTrieste/wiki/F2Py
- https://docs.scipy.org/doc/numpy/f2py/python-usage.html 
  (note: there are many inconsistencies in this example code)
-----------

*Author:*       &nbsp; Zachariah Irwin

*Institution:*  &nbsp; University of Colorado, Boulder

*Last Edit:*    &nbsp; January 2020

-----------

## General Procedure for .f90 as a Subroutine

 1. Write the .f90 subroutine that will call the Python function. The arguments should contain the Python function (exact name not required) and the variable that is being passed between the two languages.

    ***EXAMPLE:***

    ```
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
     ```


 2. Use the command shell to generate the signature file with: 

    `f2py -h [signature_file_name].pyf [.f90_file_name].f90`
     
    ***EXAMPLE:*** 

    `f2py -h callback_sig.pyf callback.f90`

    ***OUTPUT FILE:***

    ```
    !    -*- f90 -*-
    ! Note: the context of this file is case sensitive.

    python module calculate__user__routines 
      interface calculate_user_interface 
        function func(x_i_) ! in callback.f90:calculate:unknown_interface
          real :: x_i_
          double precision :: func
        end function func
      end interface calculate_user_interface
    end python module calculate__user__routines
    subroutine calculate(func,x,n) ! in callback.f90
      use calculate__user__routines
        external func
        double precision dimension(n) :: x
        integer, optional,check(len(x)>=n),depend(x) :: n=len(x)
    end subroutine calculate

    ! This file was auto-generated with f2py (version:2).
    ! See http://cens.ioc.ee/projects/f2py2e/
    ```

    The first module describes the parameters of the external Python function. The second module describes the parameters of the .f90 subroutine. Make the following modifications with respect to the .pyf:

     LINE 6: `function func(x) result (x)! in callback.f90:calculate:unknown_interface`

     LINE 7: `double precision :: x`

     LINE 8: REMOVE: `double precision :: func`

     LINE 12: INSERT: `python module callback`

     LINE 12: INSERT: `interface`

     LINE 15: `double precision dimension(n), intent(in, out) :: x`

     LINE 16: `integer, intent(in):: n`

     LINE 18: INSERT: `end interface`

     LINE 18: INSERT: `end python module callback`


 3. Recompile with: `f2py -c [signature_file_name].pyf [.f90_file_name].f90`\

    ***TO TEST FUNCTIONALITY:***

      1. Open up a Python terminal

      2. Import the compiled module (whose name was given by the user on LINE 12 of the .pyf)

      3. Define the external Python function

      4. Call the Fortran subroutine as follows: `[module_name].[subroutine_name]([python function], args)`

    ***EXAMPLE:***

      ```
      >>>import numpy as np
      >>>import callback
      >>>def foo(x): return x**2
      >>>x = np.asarray(range(5), order='F', dtype=np.float64)
      >>>print(callback.calculate(foo, x, 5))
      [ 0. 1. 4. 9. 16.]
      ```
<div style="page-break-after: always;"></div>
<div style="page-break-after: always;"></div>

## General Procedure for .f90 as a Module
 1. Write the .f90 subroutine that will call the Python function. The arguments should contain the Python function (exact name not required) and the variable that is being passed between the two languages.

    ***EXAMPLE:***

      ```
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
      ```
 2. Generate the signature file with:

      ```f2py -h [signature_file_name].pyf [.f90_file_name].f90```

      ***EXAMPLE:***

        `f2py -h callback_full_sig.pyf callback_full.f90`

      ***OUTPUT FILE:***

      ```
      !    -*- f90 -*-
      ! Note: the context of this file is case sensitive.

      python module calculate__user__routines 
        interface calculate_user_interface 
          function func(x_i_) ! in callback_full.f90:callback_full:calculate:unknown_interface
            real :: x_i_
            double precision :: func
          end function func
        end interface calculate_user_interface
      end python module calculate__user__routines
      module callback_full ! in callback_full.f90
        subroutine calculate(func,x,n) ! in callback_full.f90:callback_full
          use calculate__user__routines
          external func
          double precision dimension(n) :: x
          integer, optional,check(len(x)>=n),depend(x) :: n=len(x)
        end subroutine calculate
      end module callback_full

      ! This file was auto-generated with f2py (version:2).
      ! See http://cens.ioc.ee/projects/f2py2e/
      ```

      The first module describes the parameters of the external Python function. The second module describes the parameters of the .f90 subroutine(s). Make the following modifications with respect to the .pyf:

       LINE 6: `function func(x) result (x)! in callback.f90:calculate:unknown_interface`

       LINE 7: `double precision :: x`

       LINE 8: REMOVE: `double precision :: func`

       LINE 12: INSERT: `python module callback_full`

       LINE 12: INSERT: `interface`

       LINE 16: `double precision dimension(n), intent(in, out) :: x`

       LINE 17: `integer, intent(in):: n`

       LINE 19: INSERT: `end interface`

       LINE 19: INSERT: `end python module callback_full`


 3. Recompile with:

    `f2py -c [signature_file_name].pyf [.f90_file_name].f90`

    ***TO TEST FUNCTIONALITY:***

      1. Open up a Python terminal
      2. Import the compiled module (whose name was given by the user on LINE 124 of the example code)
      3. Define the external Python function
      4. Call the Fortran subroutine as follows: 

          `[python_module_name].[fortran_module_name].[subroutine_name]([python function], args)`

      ***EXAMPLE:***

      ```
      >>>import numpy as np
      >>>import callback_full
      >>>def foo(x): return x**2
      >>>x = np.asarray(range(5), order='F', dtype=np.float64)
      >>>print(callback_full.callback_full.calculate(foo, x, 5))
      [ 0. 1. 4. 9. 16.]
      ```

