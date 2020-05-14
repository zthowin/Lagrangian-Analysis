#----------------------------------------------
# Makefile for the Lagrangian analysis toolkit.
#
# Author:       Zachariah Irwin
# Institution:  University of Colorado Boulder
# Last Edit:    May 2020
#----------------------------------------------
ifeq ($(mode),debug)
	FCFLAGS = -fPIC -O0 -fbacktrace -fcheck=all
	PCFLAGS = --verbose --opt='-O0'
	PSFLAGS = --verbose
else
	mode = release
	FCFLAGS = -fPIC -O3 -w
	PCFLAGS = --quiet --opt='-O3' --f90flags=-w
	PSFLAGS = --quiet
endif

#----------
# Compilers
#----------
FC = gfortran
PC = f2py3

#-------------
# Source files
#-------------
SRC = boundaryConditions.f90 mymath.f90 analyticalFlowLibraryFT.f90 velocity.f90 outputFormat.f90 cgFields.f90 mesh.f90 integration.f90
LIBMESHO = mymath.o
LIBINTO = velocity.o boundaryConditions.o analyticalFlowLibraryFT.o mymath.o mesh.o
LIBMAINO = integration.o cgFields.o mymath.o analyticalFlowLibraryFT.o boundaryConditions.o velocity.o mesh.o outputFormat.o
LIBLAPACKA = -llapack -lrefblas -ltmglib
PREFIXES = mesh integration main


.PHONY: all
all: info preliminaries libraries f2py

info: 
ifneq ($(mode),release)
ifneq ($(mode),debug)
		@echo "Invalid build mode."
		@echo "Please use 'make mode=release' or 'make mode=debug'"
		@exit 1
endif
endif
	@echo "Building in mode "$(mode)" mode..."

#------------------------
# Compile Fortran modules
#------------------------
preliminaries:
		@echo "Compiling Fortran modules..."
		@$(FC) -c $(FCFLAGS) $(SRC)
		@echo "Finished compiling Fortran modules."
libraries:
		@echo "Linking static libraries..."
		@ar crs libmesh.a $(LIBMESHO)
		@ar crs libintegration.a $(LIBINTO)
		@ar crs libmain.a $(LIBMAINO)
		@echo "Finished linking static libraries."

#-------------------
# Build f2py modules
#-------------------
f2py:
		@echo "Writing signature files (interface blocks)..."
		@for prefix in $(PREFIXES); do \
			$(PC) $(PSFLAGS) -m $${prefix} -h $${prefix}SgnFile.pyf $${prefix}.f90; \
		done
		@echo "Finished writing signature files (interface blocks)."

		@echo "Creating f2py modules..."
		@for prefix in $(PREFIXES); do \
			if [ $$prefix = "main" ] ; then \
				$(PC) $(PCFLAGS) -c $${prefix}SgnFile.pyf $${prefix}.f90 -L. -l$${prefix} $(LIBLAPACKA) ; \
			else \
				$(PC) $(PCFLAGS) -c $${prefix}SgnFile.pyf $${prefix}.f90 -L. -l$$prefix ; \
			fi \
		done
		@echo "Finished creating f2py modules."
		@echo "Cleaning up..."
		@rm -vf *.o


.PHONY: clean
clean:
	@rm -vf *.so *.o *.mod *.pyf
	@for prefix in $(PREFIXES) ; do \
		rm -vf lib$${prefix}.a ; \
	done
