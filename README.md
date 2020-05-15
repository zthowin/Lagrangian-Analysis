## 1. Software Requirements

This package requires the following software dependencies:

* GNU Fortran compiler (7.5.0+)
* Python3 installation (3.6+)
* Scientifc Python stack: Numpy, Scipy (1.18+, 1.3.1+)
* VTK Library with Python bindings (8.1.2+)
* LAPACK Library (instructions for installation in Section 2) (3.8.0)

For Mac users, it is recommended to install all Python dependencies through HomeBrew.

For Linux users, use the Advanced Packaging Tool (e.g. sudo apt-get...)

This package has not yet been tested on Windows OS nor Mac OS. It is advised that Windows users run this package on an Ubuntu virtual machine as Windows does not have GNU compilers by default.

## 2. LAPACK Installation Instructions

First, download `make-LAPACK.inc` from this repository and rename it `make.inc`

Then, download the LAPACK .tar.gz for version 3.8.0 from https://www.netlib.org/lapack/ and extract it into a separate folder.

Copy `make.inc` into the LAPACK root directory (there should be a make.inc.example in this folder which is a general template -- the `make.inc` provided in this repository has been customized with optimization and static library linking flags, the latter of which must be included in the options and loader flags).

Navigate to the root directory and run the makefile to install LAPACK.

```bash
cd lapack-3.8.0/
make
```

## 3. Installation Instructions
**Note:** We will use `Lagrangian-Toolkit` as the placeholder for the name of the directory including the full package where the package was installed. You will need to edit this accordingly.

Download or clone the library. Copy the .a files from the LAPACK folder into the `src` directory for the toolkit library.

Navigate to the `src` directory and run the `makefile`:

```bash
cd /Lagrangian-Toolkit/src/
make
```

## 4. Usage
It is recommended that a common project directory structure is maintained as follows:

* root
  * flow-data-Directory
  * output-data-Directory

## 5. Examples

We have included several examples demonstrating the toolkit's capabilities, located in the `tests and demos` folder. There are five examples provided. The folder `double_gyre_demo` contains an example computation of FTLE fields for a canonical double gyre flow velocity field. The folder `clot_demo` contains an example computation of FTLE fields as well as an example computation of tracer particles trajectories for a two-dimensional artery with a fictitious clot embedded within the channel.

### Double Gyre FTLE Example

**Note:** We will use `Lagrangian-Toolkit` as the placeholder for the name of the directory including the full package where the package was installed. You will need to edit this accordingly.

* Once installation is complete, navigate to the example directory `tests-and-demos/double_gyre_demo`
* Create the folder `double_gyre_demo/output`
* Navigate to `double_gyre_demo` and open the provided file `input-DG-DEMO.dat`
* Replace the correct pathname in `Lagrangian-Toolkit` for entry `Project directory`:

```bash
Project directory: /Lagrangian-Toolkit/tests-and-demos/double_gyre_demo
```

* Save this edited configuration file.
* Navigate to `Lagrangian-Toolkit` and run `codeLagrangian_v2.py` with the argument `/Lagrangian-Toolkit/tests-and-demos/double_gyre_demo/input-DG.dat`

```bash
cd /Lagrangian-Toolkit/
python3 /src/codeLagrangian_v2.py /Lagrangian-Toolkit/test-and-demos/double_gyre_demo/input-DG.dat
```

This input file has computed the forward FTLE for steady double-gyre flow for 1.0s of advection time. There are two more input files in `double_gyre_demo` called `input-USDG.dat` and `input-BWDUSDG.dat` which will compute the forward and backward FTLE, respectively, for the unsteady double-gyre flow for 10.0s of advection time. Replace the project directory pathnames for both input files and run them as provided in the above syntax with their respective input file names.

### Clot FTLE Example

**Note:** We will use `Lagrangian-Toolkit` as the placeholder for the name of the directory including the full package where the package was installed. You will need to edit this accordingly.

* Navigate to the example directory `tests-and-demos/clot_demo`
* Create the directories `clot_demo/velocity/` and `clot_demo/output`
* Download `clot-velocity.zip` from [Google Drive](https://drive.google.com/file/d/1QELOELqdnk_DjXjSWpK6zE2UzJx_h4wp/view?usp=sharing)
* Unzip `clot-velocity.zip` and move the .vtu's to `clot_demo/velocity`
* Navigate to `clot_demo` and open the provided file `input-clot-FTLE.dat`
* Replace the correct pathname in `Lagrangian-Toolkit` for entry `Project directory`:

```bash
Project directory: /Lagrangian-Toolkit/tests-and-demos/clot_demo/
```

* Save this edited configuration file.
* Navigate to `Lagrangian-Toolkit` and run `codeLagrangian_v2.py` with the argument `/Lagrangian-Toolkit/test-and-demos/clot_demo/input-clot-FTLE.dat`

```bash
cd /Lagrangian-Toolkit/
python3 src/codeLagrangian_v2.py /Lagrangian-Toolkit/tests-and-demos/clot_demo/input-clot-FTLE.dat
```

### Clot Tracer Particle Example

**Note:** We will use `Lagrangian-Toolkit` as the placeholder for the name of the directory including the full package where the package was installed. You will need to edit this accordingly.

* Navigate to the example directory `tests-and-demos/clot_demo`
* Create the directories `clot_demo/velocity/` and `clot_demo/output`
* Download `clot-velocity.zip` from [Google Drive](https://drive.google.com/file/d/1QELOELqdnk_DjXjSWpK6zE2UzJx_h4wp/view?usp=sharing)
* Unzip `clot-velocity.zip` and move the .vtu's to `clot_demo/velocity`
* Navigate to `clot_demo` and open the provided file `input-clot-TRACER.dat`
* Replace the correct pathname in `Lagrangian-Toolkit` for entry `Project directory`:

```bash
Project directory: /Lagrangian-Toolkit/tests-and-demos/clot_demo/
```

* Save this edited configuration file.
* Navigate to `Lagrangian-Toolkit` and run `codeTracersAdvection_v1.py` with the argument `/Lagrangian-Toolkit/tests-and-demos/clot_demo/input-clot-TRACER.dat`

```bash
cd /Lagrangian-Toolkit/
python3 /src/codeTracersAdvection_v1.py /Lagrangian-Toolkit/tests-and-demos/clot_demo/input-clot-TRACER.dat
```

### Other Demos
There are two other demos--the calculations of FTLE fields for three interacting Lamb-Osseen vortices and a Bickley jet perturbed by a Rosby wave--included in `tests-and-demos` which the user may run in a similar fashion to the above examples. However, the user is warned that both demos are much more computationally expensive.
