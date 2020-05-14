This is a detailed SOP on how to link multiple Fortran modules together, each containing their own subroutines, in a shared library using gfortran and NumPy's f2py module. Example code is included. MCVE codes used for the example are matMult_loop1.f90, matMult_loop2.f90, matMult_func.f90, matrixMult.f90.


*Author:* &nbsp;       Zachariah Irwin

*Institution:* &nbsp;  University of Colorado, Boulder

*Last Edit:* &nbsp;    January 2020

----------------------------------------------------------------------------------------------

# General Steps

*Note that if you are using Python 3, all `f2py` statements should be changed to `f2py3` -- otherwise the `.so` file will not be read by Python 3.*

1. Include all `use` statements below the module name and before `implicit none` statement.

2. Compile the individual modules with `gfortran -c -fPIC module.f90`. This should produce
a `.mod` and a `.o` file for *each* `.f90`.

3. Link modules into a shared library `libNAME.a` with `ar crs libNAME.a module1.o module2.o`.

4. Create the signature file for the main module that uses the submodules with `f2py -m main_module_name -h sigFileName.pyf main_module.f90`.

5. Link the library to the main module with `f2py -c --fcompiler=gfortran sigFileName.pyf main_module.f90 -L{Full directory to shared library} -lNAME`. The -l flag will automatically wrap `NAME` with the prefix `lib` and the `.a` extension..

## Instructions for Example Code
------------------------------

1. Place all .f90 files in the same directory.

2. Use the command shell to compile `matMult_loop1.f90`,`matMult_loop2.f90`, `matMult_func.f90`:

     ```
     gfortran -c -fPIC matMult_loop1.f90
     gfortran -c -fPIC matMult_loop2.f90
     gfortran -c -fPIC matMult_func.f90
     ```
3. Use the command shell to create the shared static libary:

     `ar crs libMatMultiplication.a matMult_loop1.o matMult_loop2.o matMult_func.o`
    
4. Use the command shell to create the signature file for `matrixMult.f90`:

     `f2py -m multMatrices -h sgnFile.pyf matrixMult.f90`

5. Use the command shell to create and link the Python module `multMatrices` to `libMatMultiplication.a`:

     `f2py -c --fcompiler=gnu95 sgnFile.pyf matrixMult.f90 -L{Full directory to shared library} -lMatMultiplication`

6. Run the Python script `f2py_matrixMult.py`:

    `python3 f2py_matrixMult.py 3 3 3`

   Output:

   ```
   Loop1: time for AB(3, 3) = A(3, 3) B(3, 3) is 1.71661376953e-05 s
   Loop2: time for AB(3, 3) = A(3, 3) B(3, 3) is 4.05311584473e-06 s
   matmul function: time for AB(3, 3) = A(3, 3) B(3, 3) is 1.69277191162e-05 s
   ```

