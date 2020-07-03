DEST      = .
EXTHDRS         =
HDRS            =
MAKEFILE  = Makefile


FC        = ifort
FFLAGS    =  -nologo -O2 -debug none -nocheck -heap-arrays -mkl
DFLAGS    =  -check bounds -traceback
CFLAGS    =  -c
PFLAGS    =  -qopenmp 
MODS      = -mod .


OUT       = -o $(basename $(@F)).exe


.SUFFIXES: .exe .mod .obj .f90
.f90.obj:
	$(FC) $(FFLAGS) $(DFLAGS) $(PFLAGS) $(CFLAGS) $<

OBJ0 = LIBRARY.f90 time.f90 wavenumbers.f90 HamiltonianVII_SMO.f90 TshoppingII.f90 TshoppingIV_SMO.f90 readfile.f90 
OBJ1 = plotparameters.f90 HMdata.f90 HamiltonianVII_SMODATA.f90 

all:  TMDC 

TMDC:  module1.f90 module2.f90 paixumodule.f90 main.f90  $(OBJ0) $(OBJ1)  
	$(FC) $(FFLAGS) $(DFLAGS) $(PFLAGS)  $(LDFLAGS) $(OUT) $^
	
plot:module1.f90  plotparameters.f90   
	$(FC) $(FFLAGS) $(DFLAGS) $(PFLAGS) $(LDFLAGS) $(OUT) $^
module1: module1.f90
	$(FC) $(FFLAGS) $(DFLAGS) $(PFLAGS) $(CFLAGS) $^
module2: module2.f90
	$(FC) $(FFLAGS) $(DFLAGS) $(PFLAGS) $(CFLAGS) $^
paixumodule: paixumodule.f90
	$(FC) $(FFLAGS) $(DFLAGS) $(PFLAGS) $(CFLAGS) $^

	  
cleanup:
	@rm -f fort.*
	@rm -f *.log
	@rm -f *.out
	@rm -f *.eps
	@rm -f *.obj
	@rm -f *.mod
	@rm -f *.plt
	@rm -f *.py
	@rm -f *.dat    
	@rm -f *.pdf 

cleanall:
	@rm -f fort.*
	@rm -f *.log
	@rm -f *.out
	@rm -f *.eps
	@rm -f *.obj
	@rm -f *.mod
	@rm -f *.exe
	@rm -f *.plt
	@rm -f *.py
	@rm -f *.dat    
	@rm -f *.pdf    
