.SUFFIXES:
.SUFFIXES: .f .f90 .o


#sgi
#FL = f90
#CMPLFLG = -c -O3
#LNKFLG = -o $(TARGET)  

#alpha
#FL = fort
#CMPLFLG = -c -fast 
#LNKFLG = -o $(TARGET) -lcxml
#-warn argument_checking -C

#P4
#FL = /opt/intel/fc/9.1.043/bin/ifort
FL = mpif90
CMPLFLG = -c -O3 -lmpi
#CMPLFLG = -C  
#LNKFLG = -o $(TARGET) -L/opt/intel/Compiler/11.1/064/mkl/lib/em64t -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -static-intel
#LNKFLG = -o $(TARGET) -L/opt/intel/mkl/10.1.2.024/lib/em64t -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core
LNKFLG = -o $(TARGET) -L/opt/intel/compilerpro-12.0.2.137/mkl/lib/intel64/ -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
TARGET = IMLS.x

OBJS = MODULES_rigid4D.o AUTOSURF-PES_rigid4D.o FUNCTIONS_rigid4D.o diag.o SUBROUTINES_rigid4D.o 
$(TARGET) : $(OBJS)
	$(FL) $(OBJS) $(LNKFLG)

# this is a freeform

$(OBJS) :

#equivalent to
#%.f : %.o
.f90.o:
	$(FL) $(CMPLFLG) $<

.f.o:
	$(FL) $(CMPLFLG) $<
