BPP = bpp

CC = gcc
CXX = g++
FC = gfortran

CCFLAGS= -fPIC -g 

FFLAGS = -O3 -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none

FORTRAN_LIBRARIES =-lgfortran

LIBS = $(FORTRAN_LIBRARIES)


# rule to compile .C files 
%.o : %.C
	$(CXX) $(CCFLAGS) -o $@ -c $<

# rule to compile f90 files
# NOTE: convert single to double precision when compiling

%.o : %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

# run bpp to generate bppExamle.C from bppExample.bC
bppExample.C: bppExample.bC; $(BPP) -clean -quiet bppExample.bC

bppExampleFiles = bppExample.o 
bppExample: $(bppExampleFiles)
	$(CXX) $(CCFLAGS) -o $@ $(bppExampleFiles) $(LIBS)


test1f90:
	bpp -splitFortranLines=0 test1.bf90
	gfortran -fPIC -fdefault-real-8 -fdefault-double-8  -c test1.f90
	gfortran -o test1 test1.o
