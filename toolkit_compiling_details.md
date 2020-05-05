## General steps for compiling
As of March 2020 there are 9 Fortran modules which handle the bulk of the Lagrangian toolkit computations. Most are linked together using static libraries. There are three modules which are independent of one another (mymath, boundaryConditions, velocity) which can be compiled using the following syntax: `gfortran -c -fPIC module.f90`. The `-fPIC` flag is to ensure the generated `.o` file will allow for static linking (discussed later). During compilation, gfortran will generate  `.mod` Fortran module files which can be read by other Fortran modules with `USE MODULE_NAME` where `MODULE_NAME` is the name of compiled module `MODULE.mod`.

The remaining five modules are dependent on one another (analyticalFlowLibraryFT, mesh, integration, cgFields, main). The dependencies are listed as follows:

1. mesh: mymath
2. analyticalFlowLibraryFT: mymath
3. integration: mymath, analyticalFlowLibraryFT, boundaryConditions, velocity, mesh
4. cgFields: mymath
5. main: mymath, analyticalFlowLibraryFT, boundaryConditions, velocity, mesh, integration, cgFields, outputFormat

These five modules require an extra step in the gfortran and f2py compilations to be able to utilize one another's subroutines iff we want to be able to call them from Python. The extra step is the creation of the static library. The general syntax is as follows.

`ar crs libLIBRARY_NAME.a MODULE_1.o MODULE_2.o MODULE_3.o`

Once the static library has been created, we are ready to start compiling the Fortran modules we wish to import into Python. This is a two-step process.

1. Create the signature files which are automatically generated interfaces from f2py for the "main" module in question that we wish to import into and call from Python.

    `f2py -m MODULE_NAME -h SIGNATURE_FILE_NAME.pyf MODULE.f90 --overwrite signature`

    The `--overwrite signature` flag is added in case `SIGNATURE_FILE_NAME.pyf` already exists in the local directory.
3. Then in one step, we will compile `MODULE.f90` by linking it to the static library and the signature file.
`f2py -c --fcompiler=gnu95 SIGNATURE_FILE_NAME.pyf MODULE.f90 -L/path/to/static/library/ -lLIBRARY_NAME`

    The `--fcompiler` flag specifies what compiler to use (we have not tested anything other than GNU as of March 2020). The `-L` flag specifies that a library is to be linked to the compiled module. The `-l` will automatically wrap `LIBRARY_NAME` with `lib` and the `.a` extension. Alternatively you can write out the name of the static library file including its extension.

We have provided a simple bash script `compile.sh` which you can run in the folder where the source files are held to automatically generate the Python-callable Fortran modules. Please note that we have not included instructions on how to compile the static libraries from LAPACK whose subroutines are called from `cgFields.f90`. Instead, we have provided the `make.inc` which was used to `make` the LAPACK libraries.

Let us break down `compile.sh` line by line.
- Lines 1-4: Remove all previous compiled filetypes. This is simply to ensure that any edits made to the Fortran files are being updated.
- Lines 5-9: Compile Fortran modules which we will not import into fortran. We have added the static library linking flag `-fPIC` as well as the `-O3` optimization flag.
- Line 10: Create static library for `mesh.f90` which uses `math.f90` for some mesh operations.
- Line 11: Compile `mesh.f90`
- Line 12: Create static library for `integration.f90` which uses a multitude of modules in its subroutines.
- Lines 13-14: Compile Fortran modules.
- Line 15: Create static library for `main.f90`. This module has ties to every other sub-module and therefore needs to link all of them together.
- Lines 16-18: Create signature files/interface blocks for the modules we wish to call from Python.
- Lines 19-21: Compile modules we wish to call from Python. Note that this requires a recompilation of `mesh.f90` and `integration.f90` which does not affect the previously created static libraries. We have added flags to supress f2py verbosity and some warnings from gfortran.

## Calling Fortran module subroutines from Python
After compiling the Fortran modules, f2py will produce `.so` files which correspond to Python-importable Fortran modules. The name of the `.so` was specified after the `-m` flag when producing the signature file/interface block. This is the name of the module that needs to be imported into Python. The syntax for calling a Fortran module from Python is

`MODULE_NAME_FROM_-M.MODULE_NAME.SUBROUTINE_NAME(**args)`

As an example, suppose we have a Fortran file `my_file.f90` that contains a module `My_Module`. We will create the signature file/interface block and link that to the python module name `my_f2py_module`:

`f2py -m my_f2py_module -h my_f2py_module_SignatureFile.pyf`

Then we will wrap `my_file.f90` into a f2py module that can be called from Python in the following compilation step:

`f2py -c --fcompiler=gnu95 my_f2py_module_SignatureFile.pyf my_file.f90`

Then we will import the wrapped module into Python:

`import my_f2py_module`

We can call a specific subroutine `unique_subroutine` as follows:

`my_f2py_module.my_module.unique_subroutine(**args)`

Note that f2py is case sensitive and all subroutines/modules must be called in lower-case. Also note that all NumPy arrays which are passed into these subroutines need to be Fortran-contigious which you can supply with the argument `order='F'` when creating said arrays.