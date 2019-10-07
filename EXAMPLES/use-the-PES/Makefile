#===================================================================================
#*******************************  A U T O S U R F  *********************************
#===================================================================================
#
#    AUTOSURF Program Suite
#
#    Copyright (c) 2018 Missouri University of Science and Technology
#
#    This file is part of AUTOSURF.
#
#-----------------------------------------------------------------------------------
#===================================================================================


#===================================================================================
# "evaluatePES.f90"  -->  minimalistic example of how to use the constructed PES...
#-----------------------------------------------------------------------------------
# This Makefile will create the executable: ejec.x
#-----------------------------------------------------------------------------------

# Comment one of the following lines:
F90 = ifort    # compile with ifort
#F90 = gfortran # compile with gfortran

CMPLFLG = -c -O3
OBJS = PES_rigid4D.o evaluatePES.o 
.SUFFIXES: .f .f90 .o

all : $(OBJS)
	$(F90) $(OBJS) -o ejec.x
	rm *.o 
	rm *.mod 

$(OBJS) :
.f90.o:
	$(F90) $(CMPLFLG) $<
.f.o:
	$(F90) $(CMPLFLG) $<
