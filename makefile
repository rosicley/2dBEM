MANSEC           = KSP
CLEANFILES       = rhs.vtk solution.vtk
NP               = 1
CSOURCES 				 = $(wildcard *.cpp src/*.cpp src/mesh/*.cpp)
FCOMPILER        = gfortran -O2
CXXFLAGS        += -w
CURRENT_DIR      = $(shell pwd)

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test

s:$(GCH_FILES:.h=.gch) $(CSOURCES:.cpp=.o)
	@-${CLINKER} -o $@ $^ ${PETSC_KSP_LIB} -std=c++0x 

debug: $(CSOURCES:.cpp=.o)
	@-${CLINKER} -o $@ $^ ${PETSC_KSP_LIB} -g
	@gdb debug

clear:
	@$ rm *.o *~ s *.vtu mirror* src/mesh/*.o src/*.o $(CURRENT_DIR)/src/*.gch

run1:
	@$ mpirun -np 1 ./s

run2:
	@$ mpirun -np 2 ./s

run3:
	@$ mpirun -np 3 ./s

run4:
	@$ mpirun -np 4 ./s

run5:
	@$ mpirun -np 5 ./s

run1t:
	@$ mpirun -np 1 ./s -ksp_view -ksp_monitor_singular_value

run6:
	@$ mpirun -np 6 ./s

run7:
	@$ mpirun -np 7 ./s

run8:
	@$ mpirun -np 8 ./s

run12:
	@$ mpirun -np 12 ./s

run10:
	@$ mpirun -np 10 ./s